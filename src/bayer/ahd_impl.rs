use alloc::vec;
use crate::cfa::{CfaPattern, Channel};
use crate::lab::rgb_to_lab;

/// Adaptive Homogeneity-Directed demosaicing for Bayer CFA patterns.
///
/// Three-pass algorithm:
/// 1. Interpolate green using horizontal-only neighbors (Laplacian-corrected).
/// 2. Interpolate green using vertical-only neighbors (Laplacian-corrected).
/// 3. Complete R/B for both directions via color differences, convert to
///    CIELab, measure 3×3 homogeneity, and pick the more homogeneous
///    direction per pixel.
///
/// Reference: K. Hirakawa, T.W. Parks, "Adaptive Homogeneity-Directed
/// Demosaicing Algorithm", IEEE TIP, 2005.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;
    let g_off = npix;
    let b_off = 2 * npix;

    // Two full RGB buffers for directional estimates (planar CHW).
    let mut rgb_h = vec![0.0f32; 3 * npix];
    let mut rgb_v = vec![0.0f32; 3 * npix];

    // Pass 1+2: directional green interpolation.
    interpolate_green(input, width, height, cfa, &mut rgb_h, &mut rgb_v, npix);

    // Pass 3a: complete R/B for both directions.
    interpolate_rb(input, width, height, cfa, &mut rgb_h, npix);
    interpolate_rb(input, width, height, cfa, &mut rgb_v, npix);

    // Pass 3b: homogeneity voting.
    let mut homo_h = vec![0u8; npix];
    let mut homo_v = vec![0u8; npix];
    compute_homogeneity(&rgb_h, &rgb_v, width, height, npix, &mut homo_h, &mut homo_v);

    // Pass 3c: direction selection.
    for i in 0..npix {
        let hh = homo_h[i];
        let hv = homo_v[i];

        let (r, g, b) = if hh > hv {
            (rgb_h[i], rgb_h[g_off + i], rgb_h[b_off + i])
        } else if hv > hh {
            (rgb_v[i], rgb_v[g_off + i], rgb_v[b_off + i])
        } else {
            (
                (rgb_h[i] + rgb_v[i]) * 0.5,
                (rgb_h[g_off + i] + rgb_v[g_off + i]) * 0.5,
                (rgb_h[b_off + i] + rgb_v[b_off + i]) * 0.5,
            )
        };
        output[i] = r;
        output[g_off + i] = g;
        output[b_off + i] = b;
    }
}

/// Interpolate green channel in horizontal and vertical directions independently.
fn interpolate_green(
    input: &[f32],
    w: usize,
    h: usize,
    cfa: &CfaPattern,
    rgb_h: &mut [f32],
    rgb_v: &mut [f32],
    npix: usize,
) {
    let g_off = npix;

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let ch = cfa.color_at(y, x);
            let c_off = ch as usize * npix;

            // Copy the CFA-known channel to both buffers.
            rgb_h[c_off + idx] = input[idx];
            rgb_v[c_off + idx] = input[idx];

            if ch == Channel::Green {
                rgb_h[g_off + idx] = input[idx];
                rgb_v[g_off + idx] = input[idx];
            } else {
                let center = input[idx];

                let get = |dy: i32, dx: i32| -> f32 {
                    let ny = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
                    let nx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
                    input[ny * w + nx]
                };

                // Horizontal green estimate.
                let g_left = get(0, -1);
                let g_right = get(0, 1);
                let c_left2 = get(0, -2);
                let c_right2 = get(0, 2);
                rgb_h[g_off + idx] = (g_left + g_right) * 0.5
                    + (2.0 * center - c_left2 - c_right2) * 0.25;

                // Vertical green estimate.
                let g_top = get(-1, 0);
                let g_bot = get(1, 0);
                let c_top2 = get(-2, 0);
                let c_bot2 = get(2, 0);
                rgb_v[g_off + idx] = (g_top + g_bot) * 0.5
                    + (2.0 * center - c_top2 - c_bot2) * 0.25;
            }
        }
    }
}

/// Interpolate R and B channels using color differences from the green channel
/// already present in `rgb`.
fn interpolate_rb(
    input: &[f32],
    w: usize,
    h: usize,
    cfa: &CfaPattern,
    rgb: &mut [f32],
    npix: usize,
) {
    let g_off = npix;
    let b_off = 2 * npix;

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let ch = cfa.color_at(y, x);

            match ch {
                Channel::Green => {
                    let adj_ch = cfa.color_at(y, x ^ 1);

                    // Horizontal neighbors share the row-color.
                    let h_off = adj_ch as usize * npix;
                    let lx = x.saturating_sub(1);
                    let rx = (x + 1).min(w - 1);
                    let li = y * w + lx;
                    let ri = y * w + rx;
                    let cd_h = (input[li] - rgb[g_off + li]
                        + input[ri] - rgb[g_off + ri]) * 0.5;
                    rgb[h_off + idx] = rgb[g_off + idx] + cd_h;

                    // Vertical neighbors share the column-color.
                    let v_off = if adj_ch == Channel::Red { b_off } else { 0 };
                    let ty = y.saturating_sub(1);
                    let by = (y + 1).min(h - 1);
                    let ti = ty * w + x;
                    let bi = by * w + x;
                    let cd_v = (input[ti] - rgb[g_off + ti]
                        + input[bi] - rgb[g_off + bi]) * 0.5;
                    rgb[v_off + idx] = rgb[g_off + idx] + cd_v;
                }
                Channel::Red => {
                    let cd = diagonal_color_diff(input, rgb, w, h, cfa, y, x, Channel::Blue, g_off);
                    rgb[b_off + idx] = rgb[g_off + idx] + cd;
                }
                Channel::Blue => {
                    let cd = diagonal_color_diff(input, rgb, w, h, cfa, y, x, Channel::Red, g_off);
                    rgb[idx] = rgb[g_off + idx] + cd;
                }
            }
        }
    }
}

/// Average color difference from diagonal neighbors matching the target channel.
#[allow(clippy::too_many_arguments)]
#[inline]
fn diagonal_color_diff(
    input: &[f32],
    rgb: &[f32],
    w: usize,
    h: usize,
    cfa: &CfaPattern,
    y: usize,
    x: usize,
    target: Channel,
    g_off: usize,
) -> f32 {
    let mut sum = 0.0f32;
    let mut count = 0u32;
    for &(dy, dx) in &[(-1i32, -1i32), (-1, 1), (1, -1), (1, 1)] {
        let ny = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
        let nx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
        if cfa.color_at(ny, nx) == target {
            let nidx = ny * w + nx;
            sum += input[nidx] - rgb[g_off + nidx];
            count += 1;
        }
    }
    if count > 0 { sum / count as f32 } else { 0.0 }
}

/// Convert both directional RGB buffers to CIELab and compute per-pixel
/// homogeneity scores in a 3×3 neighborhood.
fn compute_homogeneity(
    rgb_h: &[f32],
    rgb_v: &[f32],
    w: usize,
    h: usize,
    npix: usize,
    homo_h: &mut [u8],
    homo_v: &mut [u8],
) {
    let g_off = npix;
    let b_off = 2 * npix;

    // Pre-compute Lab for both directions.
    let mut lab_h = vec![[0.0f32; 3]; npix];
    let mut lab_v = vec![[0.0f32; 3]; npix];

    for i in 0..npix {
        lab_h[i] = rgb_to_lab(rgb_h[i], rgb_h[g_off + i], rgb_h[b_off + i]);
        lab_v[i] = rgb_to_lab(rgb_v[i], rgb_v[g_off + i], rgb_v[b_off + i]);
    }

    let l_thresh = 2.0f32;
    let c_thresh_sq = 2.0f32 * 2.0f32;

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let center_h = lab_h[idx];
            let center_v = lab_v[idx];

            let mut count_h = 0u8;
            let mut count_v = 0u8;

            for ky in -1i32..=1 {
                let ny = (y as i32 + ky).clamp(0, h as i32 - 1) as usize;
                for kx in -1i32..=1 {
                    if ky == 0 && kx == 0 {
                        continue;
                    }
                    let nx = (x as i32 + kx).clamp(0, w as i32 - 1) as usize;
                    let nidx = ny * w + nx;

                    // Horizontal direction homogeneity.
                    let nh = lab_h[nidx];
                    let dl = (center_h[0] - nh[0]).abs();
                    let da = center_h[1] - nh[1];
                    let db = center_h[2] - nh[2];
                    if dl < l_thresh && da * da + db * db < c_thresh_sq {
                        count_h += 1;
                    }

                    // Vertical direction homogeneity.
                    let nv = lab_v[nidx];
                    let dl = (center_v[0] - nv[0]).abs();
                    let da = center_v[1] - nv[1];
                    let db = center_v[2] - nv[2];
                    if dl < l_thresh && da * da + db * db < c_thresh_sq {
                        count_v += 1;
                    }
                }
            }

            homo_h[idx] = count_h;
            homo_v[idx] = count_v;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_channel_preserved() {
        for cfa in &[
            CfaPattern::bayer_rggb(),
            CfaPattern::bayer_bggr(),
            CfaPattern::bayer_grbg(),
            CfaPattern::bayer_gbrg(),
        ] {
            let w = 16;
            let h = 16;
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

            // AHD blends two directions, so the known channel may not be
            // exactly preserved at pixels where both directions tie. Check
            // that the known channel is within a very tight tolerance.
            for y in 0..h {
                for x in 0..w {
                    let idx = y * w + x;
                    let ch = cfa.color_at(y, x) as usize;
                    let expected = input[idx];
                    let got = output[ch * w * h + idx];
                    assert!(
                        (got - expected).abs() < 1e-5,
                        "known channel mismatch at ({y},{x}) ch={ch}: expected {expected}, got {got}"
                    );
                }
            }
        }
    }

    #[test]
    fn uniform_color_reconstructed() {
        for cfa in &[
            CfaPattern::bayer_rggb(),
            CfaPattern::bayer_bggr(),
            CfaPattern::bayer_grbg(),
            CfaPattern::bayer_gbrg(),
        ] {
            let w = 16;
            let h = 16;
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

            for y in 3..h - 3 {
                for x in 3..w - 3 {
                    let idx = y * w + x;
                    let r = output[idx];
                    let g = output[w * h + idx];
                    let b = output[2 * w * h + idx];
                    assert!((r - 0.8).abs() < 0.05, "R at ({y},{x}) = {r}");
                    assert!((g - 0.5).abs() < 0.05, "G at ({y},{x}) = {g}");
                    assert!((b - 0.3).abs() < 0.05, "B at ({y},{x}) = {b}");
                }
            }
        }
    }
}
