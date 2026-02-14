//! Demosaicing algorithms for Bayer and X-Trans CFA sensors.
//!
//! Takes single-channel CFA (color filter array) sensor data and reconstructs
//! full-color RGB images. Supports both 2x2 Bayer patterns (RGGB, BGGR, GRBG,
//! GBRG) and 6x6 Fujifilm X-Trans patterns.
//!
//! # Algorithms
//!
//! **Bayer:**
//! - [`Bilinear`](Algorithm::Bilinear) — simple neighbor averaging
//! - [`Mhc`](Algorithm::Mhc) — Malvar-He-Cutler gradient-corrected interpolation
//! - [`Ppg`](Algorithm::Ppg) — Patterned Pixel Grouping
//! - [`Ahd`](Algorithm::Ahd) — Adaptive Homogeneity-Directed
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

extern crate alloc;

mod cfa;
mod error;
mod bayer;
mod xtrans;
mod lab;

use core::fmt;

pub use cfa::{CfaPattern, Channel};
pub use error::DemosaicError;

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
                        _ => unreachable!(),
                    },
                    cfa: "Bayer",
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
                        _ => unreachable!(),
                    },
                    cfa: "X-Trans",
                });
            }
        }
    }

    Ok(())
}

