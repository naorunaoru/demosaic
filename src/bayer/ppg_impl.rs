use alloc::vec;
use crate::cfa::{CfaPattern, Channel};

/// Patterned Pixel Grouping demosaicing for Bayer CFA patterns.
///
/// Two-pass algorithm:
/// 1. Directional green interpolation with Laplacian correction — computes
///    horizontal and vertical gradients, picks the direction with lower
///    gradient (or averages when tied).
/// 2. R/B interpolation via color differences (R−G, B−G) from neighbors
///    where the target channel is CFA-known.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    // Intermediate green channel at full resolution.
    let mut green = vec![0.0f32; npix];

    // Pass 1: green interpolation.
    interpolate_green(input, width, height, cfa, &mut green);

    // Pass 2: R/B via color differences, writing all three channels to output.
    interpolate_rb(input, &green, width, height, cfa, output, npix);
}

/// Pass 1: interpolate green at non-green sites using gradient-adaptive
/// Laplacian-corrected directional interpolation.
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

            let get = |dy: i32, dx: i32| -> f32 {
                let ny = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
                let nx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
                input[ny * w + nx]
            };

            // Cardinal green neighbors (distance 1).
            let g_left = get(0, -1);
            let g_right = get(0, 1);
            let g_top = get(-1, 0);
            let g_bot = get(1, 0);

            // Same-color neighbors at distance 2 (for Laplacian correction).
            let c_left2 = get(0, -2);
            let c_right2 = get(0, 2);
            let c_top2 = get(-2, 0);
            let c_bot2 = get(2, 0);

            // Horizontal gradient.
            let h_grad = (g_left - g_right).abs()
                + (2.0 * center - c_left2 - c_right2).abs();

            // Vertical gradient.
            let v_grad = (g_top - g_bot).abs()
                + (2.0 * center - c_top2 - c_bot2).abs();

            // Green estimates with Laplacian correction.
            let g_h = (g_left + g_right) * 0.5
                + (2.0 * center - c_left2 - c_right2) * 0.25;
            let g_v = (g_top + g_bot) * 0.5
                + (2.0 * center - c_top2 - c_bot2) * 0.25;

            green[idx] = if h_grad < v_grad {
                g_h
            } else if v_grad < h_grad {
                g_v
            } else {
                (g_h + g_v) * 0.5
            };
        }
    }
}

/// Pass 2: interpolate R and B channels using color differences (value − green).
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
                    // Determine which color is in the same row vs same column.
                    let adj_ch = cfa.color_at(y, x ^ 1);

                    // Horizontal neighbors share the row-color.
                    let h_color = adj_ch as usize;
                    let h_cd = color_diff_cardinal(input, green, w, h, y, x, true);
                    output[h_color * npix + idx] = green[idx] + h_cd;

                    // Vertical neighbors share the column-color.
                    let v_color = if adj_ch == Channel::Red { b_off } else { r_off };
                    let v_cd = color_diff_cardinal(input, green, w, h, y, x, false);
                    output[v_color + idx] = green[idx] + v_cd;
                }
                Channel::Red => {
                    // Need B — diagonal neighbors.
                    let cd = color_diff_diagonal(input, green, w, h, cfa, y, x, Channel::Blue);
                    output[b_off + idx] = green[idx] + cd;
                }
                Channel::Blue => {
                    // Need R — diagonal neighbors.
                    let cd = color_diff_diagonal(input, green, w, h, cfa, y, x, Channel::Red);
                    output[r_off + idx] = green[idx] + cd;
                }
            }
        }
    }
}

/// Average color difference (input − green) from the two cardinal neighbors
/// in a given direction.
#[inline]
fn color_diff_cardinal(
    input: &[f32],
    green: &[f32],
    width: usize,
    height: usize,
    y: usize,
    x: usize,
    horizontal: bool,
) -> f32 {
    let (i0, i1) = if horizontal {
        let lx = x.saturating_sub(1);
        let rx = (x + 1).min(width - 1);
        (y * width + lx, y * width + rx)
    } else {
        let ty = y.saturating_sub(1);
        let by = (y + 1).min(height - 1);
        (ty * width + x, by * width + x)
    };
    let cd0 = input[i0] - green[i0];
    let cd1 = input[i1] - green[i1];
    (cd0 + cd1) * 0.5
}

/// Average color difference from diagonal neighbors that match the target channel.
#[allow(clippy::too_many_arguments)]
#[inline]
fn color_diff_diagonal(
    input: &[f32],
    green: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    y: usize,
    x: usize,
    target: Channel,
) -> f32 {
    let mut sum = 0.0f32;
    let mut count = 0u32;
    for &(dy, dx) in &[(-1i32, -1i32), (-1, 1), (1, -1), (1, 1)] {
        let ny = (y as i32 + dy).clamp(0, height as i32 - 1) as usize;
        let nx = (x as i32 + dx).clamp(0, width as i32 - 1) as usize;
        if cfa.color_at(ny, nx) == target {
            let nidx = ny * width + nx;
            sum += input[nidx] - green[nidx];
            count += 1;
        }
    }
    if count > 0 { sum / count as f32 } else { 0.0 }
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
