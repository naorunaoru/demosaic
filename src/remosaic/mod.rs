mod bilinear_impl;

use crate::cfa::CfaPattern;
use crate::error::DemosaicError;

/// Convert Quad Bayer raw data to Standard Bayer raw data at the same resolution.
///
/// Uses bilinear interpolation (inverse-distance-squared weighted) for pixels
/// where the Quad Bayer and Standard Bayer patterns disagree on the color.
/// Pixels where both patterns agree pass through unchanged (37.5% of pixels).
///
/// After remosaicing, `output` contains standard Bayer data that can be passed
/// to [`crate::demosaic`] with the corresponding standard Bayer [`CfaPattern`].
///
/// # Arguments
/// - `input`: single-channel Quad Bayer CFA data, row-major, length = `width * height`
/// - `width`, `height`: image dimensions
/// - `quad_cfa`: the Quad Bayer CFA pattern (4x4)
/// - `output`: pre-allocated output buffer, length = `width * height`
pub fn remosaic(
    input: &[f32],
    width: usize,
    height: usize,
    quad_cfa: &CfaPattern,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    let npix = width * height;
    if input.len() != npix {
        return Err(DemosaicError::InputSizeMismatch { expected: npix, got: input.len() });
    }
    if output.len() != npix {
        return Err(DemosaicError::OutputSizeMismatch { expected: npix, got: output.len() });
    }
    if !quad_cfa.is_quad_bayer() {
        return Err(DemosaicError::UnsupportedAlgorithm {
            algorithm: "remosaic",
            cfa: "non-Quad-Bayer",
        });
    }

    let bayer_cfa = quad_cfa.quad_to_bayer()
        .expect("is_quad_bayer() already validated");

    bilinear_impl::remosaic(input, width, height, quad_cfa, &bayer_cfa, output);

    Ok(())
}
