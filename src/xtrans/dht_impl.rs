use alloc::vec;
use crate::cfa::{CfaPattern, Channel};

/// Directional Homogeneity Test (DHT) demosaicing for X-Trans CFA patterns.
///
/// 2-direction algorithm (H/V) — lighter than Markesteijn's 4 directions but
/// higher quality than bilinear. Well-suited for CPU fallback when the WebGPU
/// DHT shader is unavailable.
///
/// Pipeline:
///   1. Interpolate green in horizontal and vertical directions independently
///      using inverse-distance-weighted averaging of green neighbors within ±3.
///   2. Select/blend H vs V green via local homogeneity (3×3 gradient test).
///   3. Interpolate R/B via color-difference in a 5×5 window.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;
    let w = width as isize;
    let h = height as isize;

    // Pass 1: directional green interpolation
    let mut green_h = vec![0.0f32; npix];
    let mut green_v = vec![0.0f32; npix];

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let val = input[idx];

            if cfa.color_at(y, x) == Channel::Green {
                green_h[idx] = val;
                green_v[idx] = val;
                continue;
            }

            // Horizontal: inverse-distance-weighted green neighbors in same row, ±3
            let mut sum_h = 0.0f32;
            let mut wt_h = 0.0f32;
            for kx in -3isize..=3 {
                let nx = x as isize + kx;
                if nx >= 0 && nx < w {
                    let nx = nx as usize;
                    if cfa.color_at(y, nx) == Channel::Green {
                        let d = kx.unsigned_abs() as f32;
                        let weight = 1.0 / (1.0 + d * d);
                        sum_h += input[y * width + nx] * weight;
                        wt_h += weight;
                    }
                }
            }
            green_h[idx] = if wt_h > 0.0 { sum_h / wt_h } else { val };

            // Vertical: inverse-distance-weighted green neighbors in same column, ±3
            let mut sum_v = 0.0f32;
            let mut wt_v = 0.0f32;
            for ky in -3isize..=3 {
                let ny = y as isize + ky;
                if ny >= 0 && ny < h {
                    let ny = ny as usize;
                    if cfa.color_at(ny, x) == Channel::Green {
                        let d = ky.unsigned_abs() as f32;
                        let weight = 1.0 / (1.0 + d * d);
                        sum_v += input[ny * width + x] * weight;
                        wt_v += weight;
                    }
                }
            }
            green_v[idx] = if wt_v > 0.0 { sum_v / wt_v } else { val };
        }
    }

    // Pass 2: homogeneity selection + R/B color-difference interpolation
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let ch = cfa.color_at(y, x);

            let gh = green_h[idx];
            let gv = green_v[idx];

            // Homogeneity: sum of absolute green differences in 3×3 window
            let mut hom_h = 0.0f32;
            let mut hom_v = 0.0f32;
            for ky in -1isize..=1 {
                for kx in -1isize..=1 {
                    if ky == 0 && kx == 0 {
                        continue;
                    }
                    let ny = (y as isize + ky).clamp(0, h - 1) as usize;
                    let nx = (x as isize + kx).clamp(0, w - 1) as usize;
                    let nidx = ny * width + nx;
                    hom_h += (green_h[nidx] - gh).abs();
                    hom_v += (green_v[nidx] - gv).abs();
                }
            }

            // Smooth blending: direction with lower variation gets more weight
            let alpha = hom_v / (hom_h + hom_v + 1e-6);
            let green = alpha * gh + (1.0 - alpha) * gv;

            let mut rgb = [0.0f32; 3];

            if ch == Channel::Green {
                rgb[1] = input[idx];
                rgb[0] = interp_color_diff(input, &green_h, &green_v, width, height, cfa, y, x, Channel::Red);
                rgb[2] = interp_color_diff(input, &green_h, &green_v, width, height, cfa, y, x, Channel::Blue);
            } else {
                rgb[ch as usize] = input[idx];
                rgb[1] = green;
                let other = if ch == Channel::Red { Channel::Blue } else { Channel::Red };
                rgb[other as usize] = interp_color_diff(input, &green_h, &green_v, width, height, cfa, y, x, other);
            }

            output[idx] = rgb[0];
            output[npix + idx] = rgb[1];
            output[2 * npix + idx] = rgb[2];
        }
    }
}

/// Interpolate a missing channel via color-difference in a 5×5 window.
///
/// At each neighbor that has `target_ch` in the CFA, compute (cfa_value - green_avg)
/// and average with inverse-Manhattan-distance weighting. Add the result to the
/// local green estimate.
#[inline]
fn interp_color_diff(
    input: &[f32],
    green_h: &[f32],
    green_v: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    y: usize,
    x: usize,
    target_ch: Channel,
) -> f32 {
    let w = width as isize;
    let h = height as isize;
    let here = y * width + x;
    let green_here = (green_h[here] + green_v[here]) * 0.5;

    let mut cd_sum = 0.0f32;
    let mut cd_wt = 0.0f32;

    for ky in -2isize..=2 {
        let ny = y as isize + ky;
        if ny < 0 || ny >= h {
            continue;
        }
        let ny = ny as usize;
        for kx in -2isize..=2 {
            let nx = x as isize + kx;
            if nx < 0 || nx >= w {
                continue;
            }
            let nx = nx as usize;
            if cfa.color_at(ny, nx) == target_ch {
                let nidx = ny * width + nx;
                let green_neighbor = (green_h[nidx] + green_v[nidx]) * 0.5;
                let cd = input[nidx] - green_neighbor;
                let d = (ky.unsigned_abs() + kx.unsigned_abs()) as f32;
                let weight = 1.0 / (1.0 + d);
                cd_sum += cd * weight;
                cd_wt += weight;
            }
        }
    }

    if cd_wt > 0.0 {
        green_here + cd_sum / cd_wt
    } else {
        green_here
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_color_preserved() {
        let cfa = CfaPattern::xtrans_default();
        let w = 24;
        let h = 24;
        let mut input = vec![0.0f32; w * h];

        // Each CFA site gets its "true" channel value
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
        demosaic(&input, w, h, &cfa, &mut output);

        // Interior pixels should reconstruct close to true values
        for y in 3..h - 3 {
            for x in 3..w - 3 {
                let idx = y * w + x;
                let r = output[idx];
                let g = output[w * h + idx];
                let b = output[2 * w * h + idx];
                assert!((r - 0.8).abs() < 0.1, "R at ({y},{x}) = {r}");
                assert!((g - 0.5).abs() < 0.1, "G at ({y},{x}) = {g}");
                assert!((b - 0.3).abs() < 0.1, "B at ({y},{x}) = {b}");
            }
        }
    }

    #[test]
    fn known_channel_exact() {
        // Verify that the known CFA channel is preserved exactly
        let cfa = CfaPattern::xtrans_default();
        let w = 18;
        let h = 18;
        let mut input = vec![0.0f32; w * h];
        for i in 0..w * h {
            input[i] = (i as f32) / (w * h) as f32;
        }

        let mut output = vec![0.0f32; 3 * w * h];
        demosaic(&input, w, h, &cfa, &mut output);

        for y in 0..h {
            for x in 0..w {
                let idx = y * w + x;
                let ch = cfa.color_at(y, x) as usize;
                let out_val = output[ch * w * h + idx];
                let in_val = input[idx];
                assert_eq!(out_val, in_val, "known channel at ({y},{x}) ch={ch}");
            }
        }
    }
}
