use crate::bayer;
use crate::cfa::CfaPattern;
use crate::error::DemosaicError;
use crate::BayerAlgorithm;

/// 2x2-bin a single-channel image by averaging each 2x2 block of pixels.
///
/// Reduces dimensions by half in each axis.
///
/// # Arguments
/// - `input`: row-major single-channel data, length = `width * height`
/// - `width`, `height`: input dimensions, must both be even
/// - `output`: pre-allocated buffer, length = `(width / 2) * (height / 2)`
pub fn bin2x2(
    input: &[f32],
    width: usize,
    height: usize,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    if width % 2 != 0 || height % 2 != 0 {
        return Err(DemosaicError::DimensionsNotEven { width, height });
    }
    let npix = width * height;
    if input.len() != npix {
        return Err(DemosaicError::InputSizeMismatch { expected: npix, got: input.len() });
    }
    let out_w = width / 2;
    let out_h = height / 2;
    let out_npix = out_w * out_h;
    if output.len() != out_npix {
        return Err(DemosaicError::OutputSizeMismatch { expected: out_npix, got: output.len() });
    }

    for by in 0..out_h {
        for bx in 0..out_w {
            let r = 2 * by;
            let c = 2 * bx;
            let sum = input[r * width + c]
                + input[r * width + c + 1]
                + input[(r + 1) * width + c]
                + input[(r + 1) * width + c + 1];
            output[by * out_w + bx] = sum * 0.25;
        }
    }

    Ok(())
}

/// Demosaic a Quad Bayer CFA image via 2x2 binning + standard Bayer demosaic.
///
/// Produces a full-color image at half the input resolution in each dimension.
///
/// # Arguments
/// - `input`: Quad Bayer CFA data, row-major, length = `width * height`
/// - `width`, `height`: input dimensions (must be even)
/// - `cfa`: a 4x4 Quad Bayer CFA pattern (e.g. `CfaPattern::quad_bayer_rggb()`)
/// - `algorithm`: which Bayer demosaic to apply after binning
/// - `output`: pre-allocated buffer, planar CHW, length = `3 * (width/2) * (height/2)`
pub fn demosaic_quad_binned(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    algorithm: BayerAlgorithm,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    if !cfa.is_quad_bayer() {
        return Err(DemosaicError::UnsupportedAlgorithm {
            algorithm: "quad_binned",
            cfa: "non-Quad-Bayer",
        });
    }

    if width % 2 != 0 || height % 2 != 0 {
        return Err(DemosaicError::DimensionsNotEven { width, height });
    }

    let npix = width * height;
    if input.len() != npix {
        return Err(DemosaicError::InputSizeMismatch { expected: npix, got: input.len() });
    }

    let out_w = width / 2;
    let out_h = height / 2;
    let out_npix = out_w * out_h;
    if output.len() != 3 * out_npix {
        return Err(DemosaicError::OutputSizeMismatch {
            expected: 3 * out_npix,
            got: output.len(),
        });
    }

    // Step 1: Bin the quad Bayer input into a standard Bayer image.
    let mut binned = alloc::vec![0.0f32; out_npix];
    bin2x2(input, width, height, &mut binned)?;

    // Step 2: Derive the standard Bayer CFA pattern.
    let bayer_cfa = cfa.quad_to_bayer()
        .expect("is_quad_bayer() already validated");

    // Step 3: Dispatch to existing Bayer demosaic.
    match algorithm {
        BayerAlgorithm::Bilinear => bayer::bilinear(&binned, out_w, out_h, &bayer_cfa, output),
        BayerAlgorithm::Mhc => bayer::mhc(&binned, out_w, out_h, &bayer_cfa, output),
        BayerAlgorithm::Ppg => bayer::ppg(&binned, out_w, out_h, &bayer_cfa, output),
        BayerAlgorithm::Ahd => bayer::ahd(&binned, out_w, out_h, &bayer_cfa, output),
        BayerAlgorithm::Vng => bayer::vng(&binned, out_w, out_h, &bayer_cfa, output),
    }

    Ok(())
}

/// Demosaic a Quad Bayer CFA image via binning to interleaved RGB output.
///
/// Same as [`demosaic_quad_binned`], but outputs interleaved HWC layout.
pub fn demosaic_quad_binned_interleaved(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    algorithm: BayerAlgorithm,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    let out_w = width / 2;
    let out_h = height / 2;
    let out_npix = out_w * out_h;
    let mut planar = alloc::vec![0.0f32; 3 * out_npix];
    demosaic_quad_binned(input, width, height, cfa, algorithm, &mut planar)?;
    crate::planar_to_interleaved(&planar, output);
    Ok(())
}

#[cfg(test)]
mod tests {
    use alloc::vec;
    use super::*;
    use crate::cfa::Channel;

    #[test]
    fn bin2x2_basic() {
        let input = [1.0, 2.0, 3.0, 4.0,
                     5.0, 6.0, 7.0, 8.0,
                     9.0, 10.0, 11.0, 12.0,
                     13.0, 14.0, 15.0, 16.0f32];
        let mut output = [0.0f32; 4];
        bin2x2(&input, 4, 4, &mut output).unwrap();
        assert_eq!(output[0], (1.0 + 2.0 + 5.0 + 6.0) / 4.0);
        assert_eq!(output[1], (3.0 + 4.0 + 7.0 + 8.0) / 4.0);
        assert_eq!(output[2], (9.0 + 10.0 + 13.0 + 14.0) / 4.0);
        assert_eq!(output[3], (11.0 + 12.0 + 15.0 + 16.0) / 4.0);
    }

    #[test]
    fn bin2x2_odd_dimensions() {
        let input = [0.0f32; 15];
        let mut output = [0.0f32; 4];
        assert!(matches!(
            bin2x2(&input, 5, 3, &mut output),
            Err(DemosaicError::DimensionsNotEven { .. })
        ));
    }

    #[test]
    fn quad_binned_solid_color() {
        for cfa in &[
            CfaPattern::quad_bayer_rggb(),
            CfaPattern::quad_bayer_bggr(),
            CfaPattern::quad_bayer_grbg(),
            CfaPattern::quad_bayer_gbrg(),
        ] {
            let (w, h) = (64, 64);
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

            let ow = w / 2;
            let oh = h / 2;
            let out_npix = ow * oh;
            let mut output = vec![0.0f32; 3 * out_npix];

            for alg in &[
                BayerAlgorithm::Bilinear,
                BayerAlgorithm::Mhc,
                BayerAlgorithm::Ppg,
                BayerAlgorithm::Ahd,
            ] {
                demosaic_quad_binned(&input, w, h, cfa, *alg, &mut output).unwrap();

                let border = 4;
                for y in border..oh - border {
                    for x in border..ow - border {
                        let idx = y * ow + x;
                        let r = output[idx];
                        let g = output[out_npix + idx];
                        let b = output[2 * out_npix + idx];
                        assert!((r - 0.8).abs() < 0.05,
                            "{alg:?}: R at ({y},{x}) = {r}");
                        assert!((g - 0.5).abs() < 0.05,
                            "{alg:?}: G at ({y},{x}) = {g}");
                        assert!((b - 0.3).abs() < 0.05,
                            "{alg:?}: B at ({y},{x}) = {b}");
                    }
                }
            }
        }
    }

    #[test]
    fn quad_binned_output_dimensions() {
        let cfa = CfaPattern::quad_bayer_rggb();
        let (w, h) = (64, 64);
        let input = vec![0.5f32; w * h];
        let ow = w / 2;
        let oh = h / 2;

        // Correct size should succeed
        let mut output = vec![0.0f32; 3 * ow * oh];
        assert!(demosaic_quad_binned(&input, w, h, &cfa, BayerAlgorithm::Bilinear, &mut output).is_ok());

        // Wrong size should fail
        let mut wrong = vec![0.0f32; 3 * w * h];
        assert!(matches!(
            demosaic_quad_binned(&input, w, h, &cfa, BayerAlgorithm::Bilinear, &mut wrong),
            Err(DemosaicError::OutputSizeMismatch { .. })
        ));
    }
}
