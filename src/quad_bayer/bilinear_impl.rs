use crate::cfa::CfaPattern;

/// Bilinear demosaicing for Quad Bayer CFA patterns.
///
/// For each pixel, averages all same-color neighbors within a 5x5 window.
/// The Quad Bayer 4x4 repeating tile has each color in 2x2 blocks, so a
/// 5x5 window (radius 2) guarantees at least 4 samples of every color for
/// any position in the tile.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let mut rgb = [0.0f32; 3];
            let mut count = [0u32; 3];

            let y_lo = y.saturating_sub(2);
            let y_hi = (y + 2).min(height - 1);
            let x_lo = x.saturating_sub(2);
            let x_hi = (x + 2).min(width - 1);

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
    fn solid_color_reconstruction() {
        for cfa in &[
            CfaPattern::quad_bayer_rggb(),
            CfaPattern::quad_bayer_bggr(),
            CfaPattern::quad_bayer_grbg(),
            CfaPattern::quad_bayer_gbrg(),
        ] {
            let w = 16;
            let h = 16;
            let input = vec![0.5f32; w * h];
            let mut output = vec![0.0f32; 3 * w * h];

            demosaic(&input, w, h, cfa, &mut output);

            for y in 2..h - 2 {
                for x in 2..w - 2 {
                    let idx = y * w + x;
                    for c in 0..3 {
                        let v = output[c * w * h + idx];
                        assert!(
                            (v - 0.5).abs() < 1e-6,
                            "ch {c} at ({y},{x}) = {v}, expected 0.5"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn known_pixels_preserved() {
        let cfa = CfaPattern::quad_bayer_rggb();
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
        demosaic(&input, w, h, &cfa, &mut output);

        // Interior pixels should reconstruct reasonably well
        for y in 2..h - 2 {
            for x in 2..w - 2 {
                let idx = y * w + x;
                let r = output[idx];
                let g = output[w * h + idx];
                let b = output[2 * w * h + idx];
                assert!((r - 0.8).abs() < 0.15, "R at ({y},{x}) = {r}");
                assert!((g - 0.5).abs() < 0.15, "G at ({y},{x}) = {g}");
                assert!((b - 0.3).abs() < 0.15, "B at ({y},{x}) = {b}");
            }
        }
    }
}
