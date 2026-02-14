use crate::cfa::CfaPattern;

/// Bilinear demosaicing for Bayer CFA patterns.
///
/// For each pixel, interpolates missing channels by averaging available
/// neighbors in a 3x3 window, weighted by their CFA position.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let mut rgb = [0.0f32; 3];
            let mut count = [0u32; 3];

            // Gather from 3x3 neighborhood
            let y_lo = y.saturating_sub(1);
            let y_hi = (y + 1).min(height - 1);
            let x_lo = x.saturating_sub(1);
            let x_hi = (x + 1).min(width - 1);

            for ny in y_lo..=y_hi {
                for nx in x_lo..=x_hi {
                    let ch = cfa.color_at(ny, nx) as usize;
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
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use super::*;
    use crate::cfa::Channel;

    #[test]
    fn known_pixels_preserved() {
        // 4x4 RGGB pattern with uniform value per channel
        let cfa = CfaPattern::bayer_rggb();
        let w = 4;
        let h = 4;
        let mut input = vec![0.0f32; w * h];

        // Fill: R=0.8, G=0.5, B=0.3
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

        // Center pixel (1,1) is Green in RGGB — its green value should be
        // the average of greens in 3x3 neighborhood
        let idx = 1 * w + 1;
        let g = output[w * h + idx]; // G plane
        assert!((g - 0.5).abs() < 0.01, "G at (1,1) = {g}");

        // All pixels should have reasonable interpolated values
        for i in 0..w * h {
            let r = output[i];
            let g = output[w * h + i];
            let b = output[2 * w * h + i];
            assert!(r >= 0.0 && r <= 1.0, "R out of range at {i}: {r}");
            assert!(g >= 0.0 && g <= 1.0, "G out of range at {i}: {g}");
            assert!(b >= 0.0 && b <= 1.0, "B out of range at {i}: {b}");
        }
    }
}
