use crate::cfa::{CfaPattern, Channel};

/// Bilinear remosaic: convert Quad Bayer to Standard Bayer.
///
/// For pixels where quad and standard Bayer agree on the color, the value
/// passes through unchanged. For mismatched pixels, the target color is
/// interpolated from surrounding same-color Quad Bayer pixels using
/// inverse-distance-squared weighting in a 7x7 window (radius 3).
pub fn remosaic(
    input: &[f32],
    width: usize,
    height: usize,
    quad_cfa: &CfaPattern,
    bayer_cfa: &CfaPattern,
    output: &mut [f32],
) {
    let w = width;
    let h = height;

    for r in 0..h {
        for c in 0..w {
            let idx = r * w + c;
            let target = bayer_cfa.color_at(r, c);
            let source = quad_cfa.color_at(r, c);

            if source == target {
                output[idx] = input[idx];
            } else {
                output[idx] = interpolate(input, w, h, quad_cfa, r, c, target);
            }
        }
    }
}

/// Interpolate a target color at position (r, c) from surrounding same-color
/// Quad Bayer pixels using inverse-distance-squared weighting.
#[inline]
fn interpolate(
    input: &[f32],
    w: usize,
    h: usize,
    quad_cfa: &CfaPattern,
    r: usize,
    c: usize,
    target: Channel,
) -> f32 {
    const RADIUS: i32 = 3;

    let mut weighted_sum = 0.0f32;
    let mut weight_sum = 0.0f32;

    for dr in -RADIUS..=RADIUS {
        let nr = (r as i32 + dr).clamp(0, h as i32 - 1) as usize;
        for dc in -RADIUS..=RADIUS {
            if dr == 0 && dc == 0 {
                continue;
            }
            let nc = (c as i32 + dc).clamp(0, w as i32 - 1) as usize;
            if quad_cfa.color_at(nr, nc) == target {
                let dist_sq = (dr * dr + dc * dc) as f32;
                let wt = 1.0 / dist_sq;
                weighted_sum += wt * input[nr * w + nc];
                weight_sum += wt;
            }
        }
    }

    if weight_sum > 0.0 {
        weighted_sum / weight_sum
    } else {
        input[r * w + c]
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use super::*;

    #[test]
    fn passthrough_pixels_preserved() {
        let quad_cfa = CfaPattern::quad_bayer_rggb();
        let bayer_cfa = quad_cfa.quad_to_bayer().unwrap();
        let (w, h) = (16, 16);
        let mut input = vec![0.0f32; w * h];

        // Fill with distinct values per pixel
        for i in 0..w * h {
            input[i] = i as f32;
        }

        let mut output = vec![0.0f32; w * h];
        remosaic(&input, w, h, &quad_cfa, &bayer_cfa, &mut output);

        // Check pass-through pixels are exact
        for r in 0..h {
            for c in 0..w {
                let idx = r * w + c;
                if quad_cfa.color_at(r, c) == bayer_cfa.color_at(r, c) {
                    assert_eq!(output[idx], input[idx],
                        "pass-through pixel ({r},{c}) changed");
                }
            }
        }
    }

    #[test]
    fn uniform_color_roundtrip() {
        // Uniform per-channel values should remosaic perfectly
        let quad_cfa = CfaPattern::quad_bayer_rggb();
        let bayer_cfa = quad_cfa.quad_to_bayer().unwrap();
        let (w, h) = (32, 32);
        let mut input = vec![0.0f32; w * h];
        for r in 0..h {
            for c in 0..w {
                input[r * w + c] = match quad_cfa.color_at(r, c) {
                    Channel::Red => 0.8,
                    Channel::Green => 0.5,
                    Channel::Blue => 0.3,
                };
            }
        }

        let mut output = vec![0.0f32; w * h];
        remosaic(&input, w, h, &quad_cfa, &bayer_cfa, &mut output);

        // Interior pixels should have the correct target color value
        for r in 4..h - 4 {
            for c in 4..w - 4 {
                let idx = r * w + c;
                let expected = match bayer_cfa.color_at(r, c) {
                    Channel::Red => 0.8,
                    Channel::Green => 0.5,
                    Channel::Blue => 0.3,
                };
                assert!(
                    (output[idx] - expected).abs() < 0.05,
                    "pixel ({r},{c}) = {}, expected {expected}",
                    output[idx]
                );
            }
        }
    }
}
