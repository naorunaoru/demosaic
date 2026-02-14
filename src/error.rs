use core::fmt;

/// Errors returned by demosaicing operations.
#[derive(Debug)]
pub enum DemosaicError {
    /// Input buffer length doesn't match width * height.
    InputSizeMismatch {
        /// Expected number of elements.
        expected: usize,
        /// Actual number of elements.
        got: usize,
    },
    /// Output buffer length doesn't match 3 * width * height.
    OutputSizeMismatch {
        /// Expected number of elements.
        expected: usize,
        /// Actual number of elements.
        got: usize,
    },
    /// Algorithm not supported for this CFA pattern.
    UnsupportedAlgorithm {
        /// Name of the algorithm.
        algorithm: &'static str,
        /// CFA type (e.g. "Bayer", "X-Trans").
        cfa: &'static str,
    },
    /// Image too small for the chosen algorithm.
    ImageTooSmall {
        /// Minimum required width.
        min_width: usize,
        /// Minimum required height.
        min_height: usize,
    },
}

impl fmt::Display for DemosaicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InputSizeMismatch { expected, got } => {
                write!(f, "input buffer: expected {expected} elements, got {got}")
            }
            Self::OutputSizeMismatch { expected, got } => {
                write!(f, "output buffer: expected {expected} elements, got {got}")
            }
            Self::UnsupportedAlgorithm { algorithm, cfa } => {
                write!(f, "algorithm '{algorithm}' not supported for {cfa} CFA")
            }
            Self::ImageTooSmall { min_width, min_height } => {
                write!(f, "image too small: minimum {min_width}x{min_height}")
            }
        }
    }
}

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "std")]
impl std::error::Error for DemosaicError {}
