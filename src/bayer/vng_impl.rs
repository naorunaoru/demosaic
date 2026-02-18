use alloc::vec;
use crate::cfa::{CfaPattern, Channel};

/// Variable Number of Gradients (VNG) demosaicing for Bayer CFA patterns.
///
/// Two-pass algorithm:
/// 1. Green interpolation at non-green sites using 8-directional gradient
///    analysis — computes gradients in N, NE, E, SE, S, SW, W, NW directions,
///    selects directions with gradient below an adaptive threshold, and averages
///    the Laplacian-corrected green estimates from the selected directions.
/// 2. R/B interpolation via color differences (R−G, B−G) from neighbors
///    where the target channel is CFA-known.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    let mut green = vec![0.0f32; npix];

    // Pass 1: green interpolation via VNG direction selection.
    interpolate_green(input, width, height, cfa, &mut green);

    // Pass 2: R/B via color differences.
    interpolate_rb(input, &green, width, height, cfa, output, npix);
}

// Direction indices.
const N: usize = 0;
const NE: usize = 1;
const E: usize = 2;
const SE: usize = 3;
const S: usize = 4;
const SW: usize = 5;
const W: usize = 6;
const NW: usize = 7;

/// Pass 1: interpolate green at non-green sites using 8-directional gradient
/// analysis with adaptive threshold selection (Variable Number of Gradients).
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

            // Neighborhood values (5x5 window).
            let n1 = get(-1, 0);
            let n2 = get(-2, 0);
            let s1 = get(1, 0);
            let s2 = get(2, 0);
            let e1 = get(0, 1);
            let e2 = get(0, 2);
            let w1 = get(0, -1);
            let w2 = get(0, -2);
            let ne1 = get(-1, 1);
            let ne2 = get(-2, 2);
            let nw1 = get(-1, -1);
            let nw2 = get(-2, -2);
            let se1 = get(1, 1);
            let se2 = get(2, 2);
            let sw1 = get(1, -1);
            let sw2 = get(2, -2);

            // Compute 8 directional gradients.
            //
            // Cardinal gradients use 3 terms spanning the 5x5 window:
            //   distance-2 vs center, distance-1 vs opposite-distance-1,
            //   and a lateral term for robustness.
            //
            // Diagonal gradients use 2 terms:
            //   distance-2 vs center, distance-1 vs opposite-distance-1.
            let mut grad = [0.0f32; 8];
            grad[N] = (n2 - center).abs() + (n1 - s1).abs() + (nw1 - sw1).abs();
            grad[S] = (s2 - center).abs() + (s1 - n1).abs() + (se1 - ne1).abs();
            grad[E] = (e2 - center).abs() + (e1 - w1).abs() + (ne1 - nw1).abs();
            grad[W] = (w2 - center).abs() + (w1 - e1).abs() + (sw1 - se1).abs();
            grad[NE] = (ne2 - center).abs() + (ne1 - sw1).abs();
            grad[NW] = (nw2 - center).abs() + (nw1 - se1).abs();
            grad[SE] = (se2 - center).abs() + (se1 - nw1).abs();
            grad[SW] = (sw2 - center).abs() + (sw1 - ne1).abs();

            // Adaptive threshold: select directions with gradient <= 1.5 * minimum.
            let min_grad = grad.iter().copied().fold(f32::INFINITY, f32::min);
            let threshold = min_grad * 1.5;

            // Green estimates per direction with Laplacian correction.
            //
            // Cardinal: nearest green neighbor + (center - same_color_at_distance_2) * 0.5
            // Diagonal: average of two nearest greens + Laplacian correction from
            //           two same-color pixels at distance 2.
            let g_est = [
                /* N  */ n1 + (center - n2) * 0.5,
                /* NE */ (n1 + e1) * 0.5 + (2.0 * center - n2 - e2) * 0.25,
                /* E  */ e1 + (center - e2) * 0.5,
                /* SE */ (s1 + e1) * 0.5 + (2.0 * center - s2 - e2) * 0.25,
                /* S  */ s1 + (center - s2) * 0.5,
                /* SW */ (s1 + w1) * 0.5 + (2.0 * center - s2 - w2) * 0.25,
                /* W  */ w1 + (center - w2) * 0.5,
                /* NW */ (n1 + w1) * 0.5 + (2.0 * center - n2 - w2) * 0.25,
            ];

            // Average green estimates from qualifying directions.
            let mut sum = 0.0f32;
            let mut count = 0u32;
            for d in 0..8 {
                if grad[d] <= threshold {
                    sum += g_est[d];
                    count += 1;
                }
            }

            green[idx] = if count > 0 {
                sum / count as f32
            } else {
                // Fallback: average all estimates (shouldn't happen with 1.5x threshold).
                g_est.iter().sum::<f32>() * 0.125
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
                    let adj_ch = cfa.color_at(y, x ^ 1);

                    let h_color = adj_ch as usize;
                    let h_cd = color_diff_cardinal(input, green, w, h, y, x, true);
                    output[h_color * npix + idx] = green[idx] + h_cd;

                    let v_color = if adj_ch == Channel::Red { b_off } else { r_off };
                    let v_cd = color_diff_cardinal(input, green, w, h, y, x, false);
                    output[v_color + idx] = green[idx] + v_cd;
                }
                Channel::Red => {
                    let cd = color_diff_diagonal(input, green, w, h, cfa, y, x, Channel::Blue);
                    output[b_off + idx] = green[idx] + cd;
                }
                Channel::Blue => {
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
