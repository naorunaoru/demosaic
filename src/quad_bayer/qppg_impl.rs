use alloc::vec;
use crate::cfa::{CfaPattern, Channel};

/// Quad-PPG: gradient-based direct demosaicing for Quad Bayer CFA patterns.
///
/// Two-pass algorithm adapted from PPG for the 4x4 Quad Bayer pattern:
/// 1. Directional green interpolation with asymmetric distance-weighted
///    averaging and block-averaged Laplacian correction.
/// 2. R/B interpolation via color differences from nearby same-color neighbors.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    // Intermediate green channel at full resolution.
    let mut green = vec![0.0f32; npix];

    // Pass 1: green interpolation.
    interpolate_green(input, width, height, cfa, &mut green);

    // Pass 2: R/B via color differences, writing all three channels to output.
    interpolate_rb(input, &green, width, height, cfa, output, npix);
}

/// Corner type within a 2x2 same-color block.
#[derive(Clone, Copy)]
enum BlockCorner {
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

/// Classify a non-green pixel's position within its 2x2 same-color block.
fn block_corner(row: usize, col: usize, cfa: &CfaPattern) -> BlockCorner {
    let r4 = row % 4;
    let c4 = col % 4;
    // R block occupies (0..2, 0..2), B block occupies (2..4, 2..4) for RGGB.
    // Generalize: within whichever 2x2 block this pixel is in, find the local offset.
    let ch = cfa.color_at(row, col);
    // Find which corner of the 4x4 tile starts the 2x2 block containing this pixel.
    // For the known color, the 2x2 block is at (r4 & !1, c4 & !1) if that also has
    // the same color, but we need a simpler approach.
    let br = r4 & 1; // 0 = top row of block, 1 = bottom row
    let bc = c4 & 1; // 0 = left col of block, 1 = right col
    // But we need to check: is (r4 ^ 1, c4) the same color? If not, the block
    // boundary is at a different position. For quad Bayer, every 2x2 block of
    // the same color is aligned to even rows/cols within its quadrant.

    // Simpler: for any non-green pixel, its 2x2 block starts at even coordinates
    // within the 4x4 tile. R block starts at (0,0), B block starts at (2,2).
    // So the local offset is (r4 % 2, c4 % 2).
    let _ = ch; // used conceptually above
    match (br, bc) {
        (0, 0) => BlockCorner::TopLeft,
        (0, 1) => BlockCorner::TopRight,
        (1, 0) => BlockCorner::BottomLeft,
        (1, 1) => BlockCorner::BottomRight,
        _ => unreachable!(),
    }
}

/// Pass 1: interpolate green at non-green sites using gradient-adaptive
/// direction selection with Laplacian correction.
fn interpolate_green(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    green: &mut [f32],
) {
    let w = width;
    let h = height;

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;

            if cfa.color_at(y, x) == Channel::Green {
                green[idx] = input[idx];
                continue;
            }

            let center = input[idx];
            let corner = block_corner(y, x, cfa);

            let get = |dy: i32, dx: i32| -> f32 {
                let ny = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
                let nx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
                input[ny * w + nx]
            };

            // Directional offsets depend on corner type.
            // (g_near, g_far): offsets to nearest green along each axis
            // (mate): offset to within-block same-color neighbor
            // (prev1, prev2, next1, next2): offsets to adjacent same-color block pair
            let (g_near_h, g_far_h, mate_h, prev_h1, prev_h2, next_h1, next_h2) =
                match corner {
                    BlockCorner::TopLeft | BlockCorner::BottomLeft =>
                        (-1i32, 2i32, 1i32, -4i32, -3i32, 4i32, 5i32),
                    BlockCorner::TopRight | BlockCorner::BottomRight =>
                        (1, -2, -1, 3, 4, -5, -4),
                };
            let (g_near_v, g_far_v, mate_v, prev_v1, prev_v2, next_v1, next_v2) =
                match corner {
                    BlockCorner::TopLeft | BlockCorner::TopRight =>
                        (-1i32, 2i32, 1i32, -4i32, -3i32, 4i32, 5i32),
                    BlockCorner::BottomLeft | BlockCorner::BottomRight =>
                        (1, -2, -1, 3, 4, -5, -4),
                };

            // --- Horizontal green estimate ---
            let gl = get(0, g_near_h);
            let gr = get(0, g_far_h);
            let dl = g_near_h.unsigned_abs() as f32;
            let dr = g_far_h.unsigned_abs() as f32;
            let g_h = (dr * gl + dl * gr) / (dl + dr);

            // Block-averaged same-color Laplacian correction (horizontal)
            let c_here_h = (center + get(0, mate_h)) * 0.5;
            let c_prev_h = (get(0, prev_h1) + get(0, prev_h2)) * 0.5;
            let c_next_h = (get(0, next_h1) + get(0, next_h2)) * 0.5;
            let corr_h = (2.0 * c_here_h - c_prev_h - c_next_h) * 0.25;

            let h_grad = (gl - gr).abs()
                + (2.0 * c_here_h - c_prev_h - c_next_h).abs();

            // --- Vertical green estimate ---
            let gu = get(g_near_v, 0);
            let gd = get(g_far_v, 0);
            let du = g_near_v.unsigned_abs() as f32;
            let dd = g_far_v.unsigned_abs() as f32;
            let g_v = (dd * gu + du * gd) / (du + dd);

            let c_here_v = (center + get(mate_v, 0)) * 0.5;
            let c_prev_v = (get(prev_v1, 0) + get(prev_v2, 0)) * 0.5;
            let c_next_v = (get(next_v1, 0) + get(next_v2, 0)) * 0.5;
            let corr_v = (2.0 * c_here_v - c_prev_v - c_next_v) * 0.25;

            let v_grad = (gu - gd).abs()
                + (2.0 * c_here_v - c_prev_v - c_next_v).abs();

            // Direction selection
            green[idx] = if h_grad < v_grad {
                g_h + corr_h
            } else if v_grad < h_grad {
                g_v + corr_v
            } else {
                (g_h + corr_h + g_v + corr_v) * 0.5
            };
        }
    }
}

/// Pass 2: interpolate R and B channels using color differences.
fn interpolate_rb(
    input: &[f32],
    green: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    output: &mut [f32],
    npix: usize,
) {
    let w = width;
    let h = height;

    let r_off = 0;
    let g_off = npix;
    let b_off = 2 * npix;

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let ch = cfa.color_at(y, x);

            // Write green.
            output[g_off + idx] = green[idx];

            // Write the CFA-known channel.
            output[ch as usize * npix + idx] = input[idx];

            match ch {
                Channel::Green => {
                    // Determine which green block we are in.
                    let r4 = y % 4;
                    let c4 = x % 4;
                    // Top-right green block: r4 < 2 && c4 >= 2
                    // Bottom-left green block: r4 >= 2 && c4 < 2
                    let in_top_right = r4 < 2 && c4 >= 2;

                    if in_top_right {
                        // R from horizontal (left), B from vertical (below)
                        let cd_r = color_diff_axis(
                            input, green, w, h, cfa, y, x, Channel::Red, true,
                        );
                        output[r_off + idx] = green[idx] + cd_r;
                        let cd_b = color_diff_axis(
                            input, green, w, h, cfa, y, x, Channel::Blue, false,
                        );
                        output[b_off + idx] = green[idx] + cd_b;
                    } else {
                        // R from vertical (above), B from horizontal (right)
                        let cd_r = color_diff_axis(
                            input, green, w, h, cfa, y, x, Channel::Red, false,
                        );
                        output[r_off + idx] = green[idx] + cd_r;
                        let cd_b = color_diff_axis(
                            input, green, w, h, cfa, y, x, Channel::Blue, true,
                        );
                        output[b_off + idx] = green[idx] + cd_b;
                    }
                }
                Channel::Red => {
                    // Need B from diagonal neighbors.
                    let cd = color_diff_diagonal(
                        input, green, w, h, cfa, y, x, Channel::Blue,
                    );
                    output[b_off + idx] = green[idx] + cd;
                }
                Channel::Blue => {
                    // Need R from diagonal neighbors.
                    let cd = color_diff_diagonal(
                        input, green, w, h, cfa, y, x, Channel::Red,
                    );
                    output[r_off + idx] = green[idx] + cd;
                }
            }
        }
    }
}

/// Color difference from axis-aligned neighbors of the target channel.
#[allow(clippy::too_many_arguments)]
#[inline]
fn color_diff_axis(
    input: &[f32],
    green: &[f32],
    w: usize,
    h: usize,
    cfa: &CfaPattern,
    y: usize,
    x: usize,
    target: Channel,
    horizontal: bool,
) -> f32 {
    let mut sum = 0.0f32;
    let mut wt = 0.0f32;
    for k in -4i32..=4 {
        if k == 0 { continue; }
        let (ny, nx) = if horizontal {
            (y, (x as i32 + k).clamp(0, w as i32 - 1) as usize)
        } else {
            ((y as i32 + k).clamp(0, h as i32 - 1) as usize, x)
        };
        if cfa.color_at(ny, nx) == target {
            let nidx = ny * w + nx;
            let d = k.unsigned_abs() as f32;
            let weight = 1.0 / (d * d);
            sum += (input[nidx] - green[nidx]) * weight;
            wt += weight;
        }
    }
    if wt > 0.0 { sum / wt } else { 0.0 }
}

/// Color difference from diagonal neighbors of the target channel.
#[allow(clippy::too_many_arguments)]
#[inline]
fn color_diff_diagonal(
    input: &[f32],
    green: &[f32],
    w: usize,
    h: usize,
    cfa: &CfaPattern,
    y: usize,
    x: usize,
    target: Channel,
) -> f32 {
    let mut sum = 0.0f32;
    let mut wt = 0.0f32;
    for dy in -3i32..=3 {
        for dx in -3i32..=3 {
            if dy == 0 || dx == 0 { continue; }
            let ny = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
            let nx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
            if cfa.color_at(ny, nx) == target {
                let nidx = ny * w + nx;
                let d2 = (dy * dy + dx * dx) as f32;
                let weight = 1.0 / d2;
                sum += (input[nidx] - green[nidx]) * weight;
                wt += weight;
            }
        }
    }
    if wt > 0.0 { sum / wt } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use super::*;

    #[test]
    fn solid_color_reconstruction() {
        for cfa in &[
            CfaPattern::quad_bayer_rggb(),
            CfaPattern::quad_bayer_bggr(),
            CfaPattern::quad_bayer_grbg(),
            CfaPattern::quad_bayer_gbrg(),
        ] {
            let w = 32;
            let h = 32;
            let input = vec![0.5f32; w * h];
            let mut output = vec![0.0f32; 3 * w * h];

            demosaic(&input, w, h, cfa, &mut output);

            for y in 6..h - 6 {
                for x in 6..w - 6 {
                    let idx = y * w + x;
                    for c in 0..3 {
                        let v = output[c * w * h + idx];
                        assert!(
                            (v - 0.5).abs() < 1e-4,
                            "ch {c} at ({y},{x}) = {v}, expected 0.5"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn known_channel_preserved() {
        for cfa in &[
            CfaPattern::quad_bayer_rggb(),
            CfaPattern::quad_bayer_bggr(),
            CfaPattern::quad_bayer_grbg(),
            CfaPattern::quad_bayer_gbrg(),
        ] {
            let w = 32;
            let h = 32;
            let mut input = vec![0.0f32; w * h];
            for y in 0..h {
                for x in 0..w {
                    input[y * w + x] = match cfa.color_at(y, x) {
                        Channel::Red => 0.8,
                        Channel::Green => 0.5,
                        Channel::Blue => 0.3,
                    };
                }
            }
            let mut output = vec![0.0f32; 3 * w * h];
            demosaic(&input, w, h, cfa, &mut output);

            for y in 0..h {
                for x in 0..w {
                    let idx = y * w + x;
                    let ch = cfa.color_at(y, x) as usize;
                    assert_eq!(
                        output[ch * w * h + idx], input[idx],
                        "known channel mismatch at ({y},{x}) ch={ch}"
                    );
                }
            }
        }
    }
}
