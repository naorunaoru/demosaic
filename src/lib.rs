//! Demosaicing algorithms for Bayer, Quad Bayer, and X-Trans CFA sensors.
//!
//! Takes single-channel CFA (color filter array) sensor data and reconstructs
//! full-color RGB images. Supports 2x2 Bayer patterns (RGGB, BGGR, GRBG,
//! GBRG), 4x4 Quad Bayer patterns, and 6x6 Fujifilm X-Trans patterns.
//!
//! # Algorithms
//!
//! **Bayer:**
//! - [`Bilinear`](Algorithm::Bilinear) — simple neighbor averaging
//! - [`Mhc`](Algorithm::Mhc) — Malvar-He-Cutler gradient-corrected interpolation
//! - [`Ppg`](Algorithm::Ppg) — Patterned Pixel Grouping
//! - [`Ahd`](Algorithm::Ahd) — Adaptive Homogeneity-Directed
//!
//! **Quad Bayer:**
//! - [`Bilinear`](Algorithm::Bilinear) — 5x5 neighbor averaging
//! - [`QuadPpg`](Algorithm::QuadPpg) — gradient-based direct demosaicing
//! - [`demosaic_quad_binned`] — 2x2 binning + any Bayer algorithm (half resolution)
//! - [`remosaic`] — convert to standard Bayer, then use any Bayer algorithm
//!
//! **X-Trans:**
//! - [`Bilinear`](Algorithm::Bilinear) — simple neighbor averaging
//! - [`Markesteijn1`](Algorithm::Markesteijn1) — 1-pass, 4 directions
//! - [`Markesteijn3`](Algorithm::Markesteijn3) — 3-pass with median refinement
//! - [`Dht`](Algorithm::Dht) — Directional Homogeneity Test
//!
//! # Example
//!
//! ```
//! use demosaic::{demosaic, Algorithm, CfaPattern};
//!
//! let width = 4;
//! let height = 4;
//! let cfa = CfaPattern::bayer_rggb();
//! let input = vec![0.5f32; width * height];
//! let mut output = vec![0.0f32; 3 * width * height];
//!
//! demosaic(&input, width, height, &cfa, Algorithm::Bilinear, &mut output).unwrap();
//!
//! // output is planar CHW: [R plane, G plane, B plane]
//! let r_plane = &output[..width * height];
//! let g_plane = &output[width * height..2 * width * height];
//! let b_plane = &output[2 * width * height..];
//! ```

#![no_std]
#![warn(missing_docs)]

extern crate alloc;

mod cfa;
mod error;
mod bayer;
mod quad_bayer;
mod remosaic;
mod xtrans;
mod lab;

use core::fmt;

pub use cfa::{CfaPattern, Channel};
pub use error::DemosaicError;

/// Bayer-only demosaicing algorithm, used as inner algorithm for Quad Bayer binning.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BayerAlgorithm {
    /// Bilinear interpolation.
    Bilinear,
    /// Malvar-He-Cutler gradient-corrected interpolation.
    Mhc,
    /// Patterned Pixel Grouping.
    Ppg,
    /// Adaptive Homogeneity-Directed.
    Ahd,
}

impl fmt::Display for BayerAlgorithm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Bilinear => f.write_str("Bilinear"),
            Self::Mhc => f.write_str("MHC"),
            Self::Ppg => f.write_str("PPG"),
            Self::Ahd => f.write_str("AHD"),
        }
    }
}

/// Demosaicing algorithm selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Algorithm {
    /// Bilinear interpolation. Supported for both Bayer and X-Trans.
    Bilinear,
    /// Malvar-He-Cutler gradient-corrected linear interpolation. Bayer only.
    Mhc,
    /// Patterned Pixel Grouping. Bayer only.
    Ppg,
    /// Adaptive Homogeneity-Directed. Bayer only.
    Ahd,
    /// Markesteijn 1-pass (4 directions). X-Trans only.
    Markesteijn1,
    /// Markesteijn 3-pass (4 directions + 2 median refinements). X-Trans only.
    Markesteijn3,
    /// Directional Homogeneity Test (2 directions). X-Trans only.
    Dht,
    /// Quad-PPG: gradient-based direct demosaicing. Quad Bayer only.
    QuadPpg,
}

impl fmt::Display for Algorithm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Bilinear => f.write_str("Bilinear"),
            Self::Mhc => f.write_str("MHC"),
            Self::Ppg => f.write_str("PPG"),
            Self::Ahd => f.write_str("AHD"),
            Self::Markesteijn1 => f.write_str("Markesteijn (1-pass)"),
            Self::Markesteijn3 => f.write_str("Markesteijn (3-pass)"),
            Self::Dht => f.write_str("DHT"),
            Self::QuadPpg => f.write_str("Quad-PPG"),
        }
    }
}

/// Demosaic a single-channel CFA image to 3-channel RGB.
///
/// # Arguments
/// - `input`: single-channel CFA data, row-major, length = `width * height`
/// - `width`, `height`: image dimensions
/// - `cfa`: CFA pattern descriptor (Bayer or X-Trans, already shifted if needed)
/// - `algorithm`: which demosaicing method to use
/// - `output`: pre-allocated buffer for planar CHW output, length = `3 * width * height`
///   Layout: `[R plane (w*h), G plane (w*h), B plane (w*h)]`
pub fn demosaic(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    algorithm: Algorithm,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    let npix = width * height;
    if input.len() != npix {
        return Err(DemosaicError::InputSizeMismatch { expected: npix, got: input.len() });
    }
    if output.len() != 3 * npix {
        return Err(DemosaicError::OutputSizeMismatch { expected: 3 * npix, got: output.len() });
    }

    if cfa.is_bayer() {
        match algorithm {
            Algorithm::Bilinear => bayer::bilinear(input, width, height, cfa, output),
            Algorithm::Mhc => bayer::mhc(input, width, height, cfa, output),
            Algorithm::Ppg => bayer::ppg(input, width, height, cfa, output),
            Algorithm::Ahd => bayer::ahd(input, width, height, cfa, output),
            _ => {
                return Err(DemosaicError::UnsupportedAlgorithm {
                    algorithm: match algorithm {
                        Algorithm::Markesteijn1 => "Markesteijn1",
                        Algorithm::Markesteijn3 => "Markesteijn3",
                        Algorithm::Dht => "DHT",
                        Algorithm::QuadPpg => "Quad-PPG",
                        _ => unreachable!(),
                    },
                    cfa: "Bayer",
                });
            }
        }
    } else if cfa.is_quad_bayer() {
        match algorithm {
            Algorithm::Bilinear => quad_bayer::bilinear(input, width, height, cfa, output),
            Algorithm::QuadPpg => quad_bayer::qppg(input, width, height, cfa, output),
            _ => {
                return Err(DemosaicError::UnsupportedAlgorithm {
                    algorithm: match algorithm {
                        Algorithm::Mhc => "MHC",
                        Algorithm::Ppg => "PPG",
                        Algorithm::Ahd => "AHD",
                        Algorithm::Markesteijn1 => "Markesteijn1",
                        Algorithm::Markesteijn3 => "Markesteijn3",
                        Algorithm::Dht => "DHT",
                        _ => unreachable!(),
                    },
                    cfa: "Quad Bayer",
                });
            }
        }
    } else {
        match algorithm {
            Algorithm::Bilinear => xtrans::bilinear(input, width, height, cfa, output),
            Algorithm::Markesteijn1 | Algorithm::Markesteijn3 => {
                if width < 64 || height < 64 {
                    return Err(DemosaicError::ImageTooSmall {
                        min_width: 64,
                        min_height: 64,
                    });
                }
                if algorithm == Algorithm::Markesteijn3 {
                    xtrans::markesteijn3(input, width, height, cfa, output);
                } else {
                    xtrans::markesteijn1(input, width, height, cfa, output);
                }
            }
            Algorithm::Dht => {
                xtrans::dht(input, width, height, cfa, output);
            }
            _ => {
                return Err(DemosaicError::UnsupportedAlgorithm {
                    algorithm: match algorithm {
                        Algorithm::Mhc => "MHC",
                        Algorithm::Ppg => "PPG",
                        Algorithm::Ahd => "AHD",
                        Algorithm::QuadPpg => "Quad-PPG",
                        _ => unreachable!(),
                    },
                    cfa: "X-Trans",
                });
            }
        }
    }

    Ok(())
}

/// Demosaic a single-channel CFA image to interleaved 3-channel RGB.
///
/// Same as [`demosaic`], but outputs interleaved HWC layout: `[R,G,B, R,G,B, ...]`.
/// This is the format expected by most image crates (e.g. `image`, `png`).
///
/// Internally demosaics to a temporary planar buffer, then scatters into
/// `output` as interleaved.
pub fn demosaic_interleaved(
    input: &[f32],
    width: usize,
    height: usize,
    cfa: &CfaPattern,
    algorithm: Algorithm,
    output: &mut [f32],
) -> Result<(), DemosaicError> {
    let npix = width * height;
    let mut planar = alloc::vec![0.0f32; 3 * npix];
    demosaic(input, width, height, cfa, algorithm, &mut planar)?;
    planar_to_interleaved(&planar, output);
    Ok(())
}

/// Convert planar CHW layout to interleaved HWC between two buffers.
///
/// - `planar`: `[R0..Rn, G0..Gn, B0..Bn]`, length = `3 * width * height`
/// - `interleaved`: `[R0,G0,B0, R1,G1,B1, ...]`, length = `3 * width * height`
///
/// # Panics
///
/// Panics if `planar` and `interleaved` have different lengths or length is not
/// divisible by 3.
pub fn planar_to_interleaved(planar: &[f32], interleaved: &mut [f32]) {
    let len = planar.len();
    assert_eq!(len, interleaved.len());
    assert_eq!(len % 3, 0);
    let npix = len / 3;
    for i in 0..npix {
        interleaved[3 * i] = planar[i];
        interleaved[3 * i + 1] = planar[npix + i];
        interleaved[3 * i + 2] = planar[2 * npix + i];
    }
}

/// Convert interleaved HWC layout to planar CHW between two buffers.
///
/// - `interleaved`: `[R0,G0,B0, R1,G1,B1, ...]`, length = `3 * width * height`
/// - `planar`: `[R0..Rn, G0..Gn, B0..Bn]`, length = `3 * width * height`
///
/// # Panics
///
/// Panics if `interleaved` and `planar` have different lengths or length is not
/// divisible by 3.
pub fn interleaved_to_planar(interleaved: &[f32], planar: &mut [f32]) {
    let len = interleaved.len();
    assert_eq!(len, planar.len());
    assert_eq!(len % 3, 0);
    let npix = len / 3;
    for i in 0..npix {
        planar[i] = interleaved[3 * i];
        planar[npix + i] = interleaved[3 * i + 1];
        planar[2 * npix + i] = interleaved[3 * i + 2];
    }
}

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
    quad_bayer::binning::bin2x2(input, width, height, output)
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
    quad_bayer::binning::demosaic_quad_binned(input, width, height, cfa, algorithm, output)
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
    quad_bayer::binning::demosaic_quad_binned_interleaved(input, width, height, cfa, algorithm, output)
}

/// Convert Quad Bayer raw data to Standard Bayer raw data at the same resolution.
///
/// Uses bilinear interpolation for pixels where the Quad Bayer and Standard
/// Bayer patterns disagree on the color. Pixels where both patterns agree
/// pass through unchanged.
///
/// After remosaicing, `output` contains standard Bayer data that can be passed
/// to [`demosaic`] with the corresponding standard Bayer [`CfaPattern`]
/// (obtainable via [`CfaPattern::quad_to_bayer`]).
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
    remosaic::remosaic(input, width, height, quad_cfa, output)
}
