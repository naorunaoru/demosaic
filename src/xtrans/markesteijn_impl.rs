use alloc::vec;
use crate::cfa::{CfaPattern, Channel};
use crate::lab::rgb_to_lab;

const TS: usize = 64; // Tile size (matches LibRaw's LIBRAW_AHD_TILE)
const OVERLAP: usize = 16; // Tile overlap margin
const NDIR: usize = 4; // Number of directions for 1-pass

/// Markesteijn 1-pass demosaicing for X-Trans sensors.
///
/// Ported from LibRaw's `xtrans_interpolate(1)`. Processes the image in
/// 64x64 tiles with 16px overlap for cache efficiency. For each tile:
///
/// 1. Interpolate green in 4 directions using hex-neighbor weights
/// 2. Interpolate red/blue using directional variance analysis
/// 3. Convert all 4 directions to CIELab
/// 4. Compute gradient magnitude per direction
/// 5. Vote for smoothest directions in 3x3 neighborhoods
/// 6. Average selected directions → final RGB
pub fn demosaic(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    output: &mut [f32],
) {
    let npix = width * height;

    // Precompute flat CFA lookup (avoids modulo in hot loops)
    let cfa_lut = build_cfa_lut(cfa);

    // Precompute R/B neighbor offsets per CFA position
    let rb_offsets = build_rb_offsets(cfa);

    // Build hex neighbor offset LUT
    let allhex = build_hex_lut(cfa);

    // Pass 1: compute green bounds (min/max from hex neighbors)
    let mut green_bounds = vec![[0.0f32; 2]; npix];
    compute_green_bounds(input, width, height, &cfa_lut, &allhex, &mut green_bounds);

    // Border: simple nearest-color interpolation for 8px border
    border_interpolate(input, width, height, &cfa_lut, output, 8);

    // Pre-allocate tile buffers (reused across all tiles)
    let tpix = TS * TS;
    let mut rgb = vec![[0.0f32; 3]; NDIR * tpix];
    let mut lab = vec![[0.0f32; 3]; tpix];
    let mut drv = vec![0.0f32; NDIR * tpix];
    let mut homo = vec![0u8; NDIR * tpix];

    // Process tiles
    let stride = TS - OVERLAP;
    for top in (0..height).step_by(stride) {
        for left in (0..width).step_by(stride) {
            let tile_h = TS.min(height - top);
            let tile_w = TS.min(width - left);
            if tile_h < 16 || tile_w < 16 {
                continue;
            }

            // Clear buffers for this tile via memset
            let n = NDIR * tile_h * tile_w;
            // SAFETY: [f32;3], f32, u8 are all valid when zero-initialized
            unsafe {
                core::ptr::write_bytes(rgb.as_mut_ptr(), 0, n);
                core::ptr::write_bytes(drv.as_mut_ptr(), 0, n);
                core::ptr::write_bytes(homo.as_mut_ptr(), 0, n);
            }

            process_tile(
                input, width, height, &cfa_lut, &allhex, &rb_offsets, &green_bounds,
                top, left, tile_h, tile_w, output,
                &mut rgb, &mut lab, &mut drv, &mut homo,
            );
        }
    }
}

/// Flat CFA lookup: cfa_lut[row % 6][col % 6] = channel index (0/1/2)
type CfaLut = [[u8; 6]; 6];

fn build_cfa_lut(cfa: &CfaPattern) -> CfaLut {
    let mut lut = [[0u8; 6]; 6];
    for y in 0..6 {
        for x in 0..6 {
            lut[y][x] = cfa.color_at(y, x) as u8;
        }
    }
    lut
}

#[inline(always)]
fn cfa_color(lut: &CfaLut, y: usize, x: usize) -> u8 {
    lut[y % 6][x % 6]
}

/// Precomputed R/B neighbor offsets for each CFA position.
/// For each (row%6, col%6), stores offsets to nearest R neighbors and B neighbors.
/// Max 4 neighbors each within ±2 distance.
struct RbOffsets {
    /// offsets[row%6][col%6][target_ch (0=R, 1=B)][idx] = (dy, dx), count
    offsets: [[[[i8; 2]; 6]; 2]; 36], // [pos][target_ch][neighbor_idx] = (dy,dx)
    counts: [[u8; 2]; 36],             // [pos][target_ch] = count
}

fn build_rb_offsets(cfa: &CfaPattern) -> RbOffsets {
    let mut rb = RbOffsets {
        offsets: [[[[0i8; 2]; 6]; 2]; 36],
        counts: [[0u8; 2]; 36],
    };

    for y in 0..6 {
        for x in 0..6 {
            let pos = y * 6 + x;
            // Find R neighbors (target=0) and B neighbors (target=1)
            for target in 0..2u8 {
                let target_ch = if target == 0 { Channel::Red } else { Channel::Blue };
                let mut k = 0usize;
                for dy in -2..=2i32 {
                    for dx in -2..=2i32 {
                        if k >= 6 { break; }
                        let ny = ((y as i32 + dy + 600) % 6) as usize;
                        let nx = ((x as i32 + dx + 600) % 6) as usize;
                        if cfa.color_at(ny, nx) == target_ch {
                            rb.offsets[pos][target as usize][k] = [dy as i8, dx as i8];
                            k += 1;
                        }
                    }
                }
                rb.counts[pos][target as usize] = k as u8;
            }
        }
    }
    rb
}

/// Hex neighborhood offset lookup table.
type HexLut = [[[[(i32, i32); 8]; 2]; 3]; 3];

fn build_hex_lut(cfa: &CfaPattern) -> HexLut {
    let mut lut = [[[[(0i32, 0i32); 8]; 2]; 3]; 3];
    for row in 0..3i32 {
        for col in 0..3i32 {
            let color = cfa.color_at(row as usize, col as usize);
            let mut k = 0;

            for dy in -2..=2i32 {
                for dx in -2..=2i32 {
                    if dy == 0 && dx == 0 { continue; }
                    let nc = cfa.color_at((row + dy + 600) as usize, (col + dx + 600) as usize);
                    if color != Channel::Green {
                        if nc == Channel::Green && (dy.abs() + dx.abs() <= 3) && k < 8 {
                            lut[row as usize][col as usize][0][k] = (dy, dx);
                            k += 1;
                        }
                    } else if nc != Channel::Green && (dy.abs() + dx.abs() <= 2) && k < 8 {
                        lut[row as usize][col as usize][0][k] = (dy, dx);
                        k += 1;
                    }
                }
            }

            let n = k;
            lut[row as usize][col as usize][0][..n]
                .sort_by(|a, b| (a.0*a.0 + a.1*a.1).cmp(&(b.0*b.0 + b.1*b.1)));

            for i in 0..4.min(n) {
                lut[row as usize][col as usize][1][i] =
                    lut[row as usize][col as usize][0][i];
            }
        }
    }
    lut
}

fn compute_green_bounds(
    input: &[f32],
    width: usize,
    height: usize,
    cfa_lut: &CfaLut,
    allhex: &HexLut,
    bounds: &mut [[f32; 2]],
) {
    for y in 2..height - 2 {
        for x in 2..width - 2 {
            let idx = y * width + x;
            if cfa_color(cfa_lut, y, x) == 1 { // Green
                bounds[idx] = [input[idx], input[idx]];
                continue;
            }

            let hex = &allhex[y % 3][x % 3][1];
            let mut lo = f32::MAX;
            let mut hi = f32::MIN;

            for &(dy, dx) in hex.iter() {
                if dy == 0 && dx == 0 { break; }
                let ny = (y as i32 + dy) as usize;
                let nx = (x as i32 + dx) as usize;
                if ny < height && nx < width {
                    let v = input[ny * width + nx];
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }

            bounds[idx] = [lo, hi];
        }
    }
}

fn border_interpolate(
    input: &[f32],
    width: usize,
    height: usize,
    cfa_lut: &CfaLut,
    output: &mut [f32],
    border: usize,
) {
    let npix = width * height;
    let hb = height.min(border);
    let wb = width.min(border);

    // Top and bottom strips
    for y in (0..hb).chain(height.saturating_sub(border)..height) {
        for x in 0..width {
            interpolate_border_pixel(input, width, height, cfa_lut, output, npix, y, x);
        }
    }
    // Left and right strips (excluding corners already done)
    for y in hb..height.saturating_sub(border) {
        for x in (0..wb).chain(width.saturating_sub(border)..width) {
            interpolate_border_pixel(input, width, height, cfa_lut, output, npix, y, x);
        }
    }
}

#[inline]
fn interpolate_border_pixel(
    input: &[f32], width: usize, height: usize,
    cfa_lut: &CfaLut, output: &mut [f32], npix: usize,
    y: usize, x: usize,
) {
    let idx = y * width + x;
    let mut rgb = [0.0f32; 3];
    let mut count = [0u32; 3];

    let y_lo = y.saturating_sub(1);
    let y_hi = (y + 1).min(height - 1);
    let x_lo = x.saturating_sub(1);
    let x_hi = (x + 1).min(width - 1);

    for ny in y_lo..=y_hi {
        for nx in x_lo..=x_hi {
            let ch = cfa_color(cfa_lut, ny, nx) as usize;
            rgb[ch] += input[ny * width + nx];
            count[ch] += 1;
        }
    }

    for c in 0..3 {
        output[c * npix + idx] = if count[c] > 0 {
            rgb[c] / count[c] as f32
        } else {
            0.0
        };
    }
}

#[allow(clippy::too_many_arguments)]
fn process_tile(
    input: &[f32],
    width: usize,
    height: usize,
    cfa_lut: &CfaLut,
    allhex: &HexLut,
    rb_offsets: &RbOffsets,
    green_bounds: &[[f32; 2]],
    top: usize,
    left: usize,
    tile_h: usize,
    tile_w: usize,
    output: &mut [f32],
    rgb: &mut [[f32; 3]],
    lab: &mut [[f32; 3]],
    drv: &mut [f32],
    homo: &mut [u8],
) {
    let npix = width * height;
    let tpix = tile_h * tile_w;

    // Step 1: Green interpolation in 4 directions
    green_interpolation(
        input, width, height, cfa_lut, allhex, green_bounds,
        top, left, tile_h, tile_w, rgb,
    );

    // Step 2: Red/Blue interpolation
    rb_interpolation(
        input, width, height, cfa_lut, rb_offsets,
        top, left, tile_h, tile_w, rgb,
    );

    // Step 3: Homogeneity-directed selection
    for d in 0..NDIR {
        let base = d * tpix;

        // CIELab conversion
        for i in 0..tpix {
            let [r, g, b] = rgb[base + i];
            lab[i] = rgb_to_lab(r, g, b);
        }

        // Gradient magnitude (unrolled cardinal directions)
        for ty in 1..tile_h.saturating_sub(1) {
            let row_off = ty * tile_w;
            for tx in 1..tile_w.saturating_sub(1) {
                let ti = row_off + tx;
                let lc = lab[ti];

                let right = lab[ti + 1];
                let left_n = lab[ti - 1];
                let down = lab[ti + tile_w];
                let up = lab[ti - tile_w];

                let grad = (lc[0] - right[0]).abs() + (lc[1] - right[1]).abs() + (lc[2] - right[2]).abs()
                    + (lc[0] - left_n[0]).abs() + (lc[1] - left_n[1]).abs() + (lc[2] - left_n[2]).abs()
                    + (lc[0] - down[0]).abs() + (lc[1] - down[1]).abs() + (lc[2] - down[2]).abs()
                    + (lc[0] - up[0]).abs() + (lc[1] - up[1]).abs() + (lc[2] - up[2]).abs();

                drv[base + ti] = grad;
            }
        }
    }

    // Homogeneity voting
    for d in 0..NDIR {
        let base = d * tpix;
        for ty in 2..tile_h.saturating_sub(2) {
            let row_off = ty * tile_w;
            for tx in 2..tile_w.saturating_sub(2) {
                let ti = row_off + tx;
                let mut votes = 0u8;

                for &noff in &[
                    -(tile_w as i32) - 1, -(tile_w as i32), -(tile_w as i32) + 1,
                    -1, 0, 1,
                    tile_w as i32 - 1, tile_w as i32, tile_w as i32 + 1,
                ] {
                    let ni = (ti as i32 + noff) as usize;
                    // Find max gradient across all directions at this neighbor
                    let d0 = drv[ni];
                    let d1 = drv[tpix + ni];
                    let d2 = drv[2 * tpix + ni];
                    let d3 = drv[3 * tpix + ni];
                    let max_grad = d0.max(d1).max(d2).max(d3);
                    let threshold = max_grad * 0.875; // max - max/8
                    if drv[base + ni] <= threshold {
                        votes += 1;
                    }
                }

                homo[base + ti] = votes;
            }
        }
    }

    // Final averaging: select directions with best homogeneity
    let margin = OVERLAP / 2;
    for ty in margin..tile_h.saturating_sub(margin) {
        let iy = top + ty;
        if iy >= height { break; }

        for tx in margin..tile_w.saturating_sub(margin) {
            let ix = left + tx;
            if ix >= width { break; }

            let ti = ty * tile_w + tx;
            let img_idx = iy * width + ix;

            let h0 = homo[ti];
            let h1 = homo[tpix + ti];
            let h2 = homo[2 * tpix + ti];
            let h3 = homo[3 * tpix + ti];
            let max_homo = h0.max(h1).max(h2).max(h3);
            let threshold = max_homo.saturating_sub(max_homo >> 3);

            let mut sum = [0.0f32; 3];
            let mut cnt = 0u32;

            if h0 >= threshold { let v = rgb[ti]; sum[0] += v[0]; sum[1] += v[1]; sum[2] += v[2]; cnt += 1; }
            if h1 >= threshold { let v = rgb[tpix + ti]; sum[0] += v[0]; sum[1] += v[1]; sum[2] += v[2]; cnt += 1; }
            if h2 >= threshold { let v = rgb[2*tpix + ti]; sum[0] += v[0]; sum[1] += v[1]; sum[2] += v[2]; cnt += 1; }
            if h3 >= threshold { let v = rgb[3*tpix + ti]; sum[0] += v[0]; sum[1] += v[1]; sum[2] += v[2]; cnt += 1; }

            if cnt > 0 {
                let inv = 1.0 / cnt as f32;
                output[img_idx] = sum[0] * inv;
                output[npix + img_idx] = sum[1] * inv;
                output[2 * npix + img_idx] = sum[2] * inv;
            }
        }
    }
}

fn green_interpolation(
    input: &[f32],
    width: usize,
    height: usize,
    cfa_lut: &CfaLut,
    allhex: &HexLut,
    green_bounds: &[[f32; 2]],
    top: usize,
    left: usize,
    tile_h: usize,
    tile_w: usize,
    rgb: &mut [[f32; 3]],
) {
    let tpix = tile_h * tile_w;

    for ty in 0..tile_h {
        let iy = top + ty;
        if iy >= height { break; }

        for tx in 0..tile_w {
            let ix = left + tx;
            if ix >= width { break; }

            let ti = ty * tile_w + tx;
            let img_idx = iy * width + ix;
            let color = cfa_color(cfa_lut, iy, ix);

            if color == 1 { // Green
                let v = input[img_idx];
                for d in 0..NDIR {
                    rgb[d * tpix + ti][1] = v;
                }
            } else if iy >= 3 && iy + 3 < height && ix >= 3 && ix + 3 < width {
                let hex = &allhex[iy % 3][ix % 3][0];
                let [lo, hi] = green_bounds[img_idx];

                for d in 0..NDIR {
                    let h_idx = d * 2;
                    if h_idx + 1 >= 8 {
                        rgb[d * tpix + ti][1] = input[img_idx];
                        continue;
                    }

                    let (dy0, dx0) = hex[h_idx];
                    let (dy1, dx1) = hex[h_idx + 1];

                    let ny0 = (iy as i32 + dy0) as usize;
                    let nx0 = (ix as i32 + dx0) as usize;
                    let ny1 = (iy as i32 + dy1) as usize;
                    let nx1 = (ix as i32 + dx1) as usize;

                    if ny0 < height && nx0 < width && ny1 < height && nx1 < width {
                        let interp = (input[ny0 * width + nx0] + input[ny1 * width + nx1]) * 0.5;
                        rgb[d * tpix + ti][1] = interp.max(lo).min(hi);
                    } else {
                        rgb[d * tpix + ti][1] = input[img_idx];
                    }
                }
            } else {
                let v = input[img_idx];
                for d in 0..NDIR {
                    rgb[d * tpix + ti][1] = v;
                }
            }

            // Set known color channel for all directions
            let ch = color as usize;
            let v = input[img_idx];
            for d in 0..NDIR {
                rgb[d * tpix + ti][ch] = v;
            }
        }
    }
}

fn rb_interpolation(
    input: &[f32],
    width: usize,
    height: usize,
    cfa_lut: &CfaLut,
    rb_offsets: &RbOffsets,
    top: usize,
    left: usize,
    tile_h: usize,
    tile_w: usize,
    rgb: &mut [[f32; 3]],
) {
    let tpix = tile_h * tile_w;

    for d in 0..NDIR {
        let base = d * tpix;
        for ty in 2..tile_h.saturating_sub(2) {
            let iy = top + ty;
            if iy + 2 >= height { break; }

            for tx in 2..tile_w.saturating_sub(2) {
                let ix = left + tx;
                if ix + 2 >= width { break; }

                let ti = ty * tile_w + tx;
                let color = cfa_color(cfa_lut, iy, ix);

                // Determine which channels need interpolation
                // color 0=R: need B (target_idx=1)
                // color 2=B: need R (target_idx=0)
                // color 1=G: need both R and B
                let targets: &[u8] = match color {
                    0 => &[1],      // R pixel: need B
                    2 => &[0],      // B pixel: need R
                    _ => &[0, 1],   // G pixel: need R and B
                };

                let pos = (iy % 6) * 6 + (ix % 6);

                for &target in targets {
                    let count = rb_offsets.counts[pos][target as usize] as usize;
                    let ch = if target == 0 { 0 } else { 2 }; // R=0, B=2

                    let mut sum = 0.0f32;
                    let mut cnt = 0u32;

                    for k in 0..count {
                        let [dy, dx] = rb_offsets.offsets[pos][target as usize][k];
                        let nty = (ty as i32 + dy as i32) as usize;
                        let ntx = (tx as i32 + dx as i32) as usize;
                        if nty < tile_h && ntx < tile_w {
                            let ni = nty * tile_w + ntx;
                            let ny = (iy as i32 + dy as i32) as usize;
                            let nx = (ix as i32 + dx as i32) as usize;
                            let diff = input[ny * width + nx] - rgb[base + ni][1]; // value - green
                            sum += diff;
                            cnt += 1;
                        }
                    }

                    let green = rgb[base + ti][1];
                    rgb[base + ti][ch] = if cnt > 0 {
                        (green + sum / cnt as f32).max(0.0)
                    } else {
                        green
                    };
                }
            }
        }
    }
}

/// Markesteijn 3-pass demosaicing: initial 1-pass + two median refinement iterations.
///
/// The refinement passes median-filter color differences in a 5×5 window,
/// which smooths zipper artifacts and color fringing from the initial
/// direction-selection step.
pub fn demosaic3(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    output: &mut [f32],
) {
    // Run the 1-pass algorithm first
    demosaic(input, width, height, cfa, output);

    // Two additional median-based refinement passes
    let cfa_lut = build_cfa_lut(cfa);
    let npix = width * height;
    let mut snap = vec![0.0f32; 3 * npix];

    for _ in 0..2 {
        snap.copy_from_slice(&output[..3 * npix]);
        refine_pass(&snap, output, width, height, npix, &cfa_lut);
    }
}

/// Single refinement pass: reads from `snap`, writes updated values to `output`.
///
/// For each non-known channel at each pixel:
/// - R/B refinement: median of (target - green) from neighbors where target is CFA-known,
///   then refined = green_here + median.
/// - Green refinement at R/B sites: median of (known_ch - green) from neighbors where
///   the same channel is CFA-known, then refined_green = known_here - median.
fn refine_pass(
    snap: &[f32],
    output: &mut [f32],
    width: usize,
    height: usize,
    npix: usize,
    cfa_lut: &CfaLut,
) {
    for y in 3..height.saturating_sub(3) {
        for x in 3..width.saturating_sub(3) {
            let idx = y * width + x;
            let known = cfa_color(cfa_lut, y, x) as usize;

            for c in 0..3usize {
                if c == known {
                    continue;
                }

                let mut diffs = [0.0f32; 25];
                let mut n = 0usize;

                if c != 1 {
                    // Refine R or B: median of (c - green) from neighbors where c is CFA-known
                    for ky in -2isize..=2 {
                        for kx in -2isize..=2 {
                            let ny = (y as isize + ky) as usize;
                            let nx = (x as isize + kx) as usize;
                            if cfa_color(cfa_lut, ny, nx) as usize == c {
                                let nidx = ny * width + nx;
                                diffs[n] = snap[c * npix + nidx] - snap[npix + nidx];
                                n += 1;
                            }
                        }
                    }
                    if n >= 3 {
                        diffs[..n].sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
                        output[c * npix + idx] = snap[npix + idx] + diffs[n / 2];
                    }
                } else {
                    // Refine green at non-green site: median of (known_ch - green) from
                    // neighbors that share the same CFA channel as this pixel
                    for ky in -2isize..=2 {
                        for kx in -2isize..=2 {
                            let ny = (y as isize + ky) as usize;
                            let nx = (x as isize + kx) as usize;
                            if cfa_color(cfa_lut, ny, nx) as usize == known {
                                let nidx = ny * width + nx;
                                diffs[n] = snap[known * npix + nidx] - snap[npix + nidx];
                                n += 1;
                            }
                        }
                    }
                    if n >= 3 {
                        diffs[..n].sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
                        // G = known_value - median(known - G)
                        output[npix + idx] = snap[known * npix + idx] - diffs[n / 2];
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_lut_non_zero() {
        let cfa = CfaPattern::xtrans_default();
        let lut = build_hex_lut(&cfa);

        for r in 0..3 {
            for c in 0..3 {
                let has_nonzero = lut[r][c][0].iter().any(|&(dy, dx)| dy != 0 || dx != 0);
                assert!(has_nonzero, "hex LUT [{r}][{c}] is all zeros");
            }
        }
    }

    #[test]
    fn demosaic_small_image() {
        let cfa = CfaPattern::xtrans_default();
        let w = 72;
        let h = 72;
        let mut input = vec![0.0f32; w * h];

        for y in 0..h {
            for x in 0..w {
                input[y * w + x] = match cfa.color_at(y, x) {
                    Channel::Red => 0.7,
                    Channel::Green => 0.5,
                    Channel::Blue => 0.3,
                };
            }
        }

        let mut output = vec![0.0f32; 3 * w * h];
        demosaic(&input, w, h, &cfa, &mut output);

        for y in 16..h - 16 {
            for x in 16..w - 16 {
                let idx = y * w + x;
                let r = output[idx];
                let g = output[w * h + idx];
                let b = output[2 * w * h + idx];
                assert!(r > 0.0, "R=0 at ({y},{x})");
                assert!(g > 0.0, "G=0 at ({y},{x})");
                assert!(b > 0.0, "B=0 at ({y},{x})");
                assert!(r <= 1.0 && g <= 1.0 && b <= 1.0,
                    "value > 1.0 at ({y},{x}): R={r} G={g} B={b}");
            }
        }
    }

    #[test]
    fn demosaic3_small_image() {
        let cfa = CfaPattern::xtrans_default();
        let w = 72;
        let h = 72;
        let mut input = vec![0.0f32; w * h];

        // Spatially varying input — gradient with per-channel offsets
        // This ensures the 1-pass result has non-trivial color differences
        // that the median refinement can actually smooth.
        for y in 0..h {
            for x in 0..w {
                let base = (y as f32 / h as f32) * 0.5 + (x as f32 / w as f32) * 0.3;
                input[y * w + x] = match cfa.color_at(y, x) {
                    Channel::Red => (base + 0.3).min(1.0),
                    Channel::Green => (base + 0.1).min(1.0),
                    Channel::Blue => base.min(1.0),
                };
            }
        }

        let mut output1 = vec![0.0f32; 3 * w * h];
        let mut output3 = vec![0.0f32; 3 * w * h];
        demosaic(&input, w, h, &cfa, &mut output1);
        demosaic3(&input, w, h, &cfa, &mut output3);

        // 3-pass should produce valid output
        for y in 16..h - 16 {
            for x in 16..w - 16 {
                let idx = y * w + x;
                let r = output3[idx];
                let g = output3[w * h + idx];
                let b = output3[2 * w * h + idx];
                assert!(r > 0.0, "3-pass R=0 at ({y},{x})");
                assert!(g > 0.0, "3-pass G=0 at ({y},{x})");
                assert!(b > 0.0, "3-pass B=0 at ({y},{x})");
            }
        }

        // 3-pass should differ from 1-pass (refinement changed something)
        let differs = (0..3 * w * h).any(|i| (output1[i] - output3[i]).abs() > 1e-6);
        assert!(differs, "3-pass output identical to 1-pass — refinement had no effect");
    }
}
