use core::fmt;

/// Color channel in a CFA pattern.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Channel {
    /// Red channel.
    Red = 0,
    /// Green channel.
    Green = 1,
    /// Blue channel.
    Blue = 2,
}

impl fmt::Display for Channel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Red => f.write_str("R"),
            Self::Green => f.write_str("G"),
            Self::Blue => f.write_str("B"),
        }
    }
}

/// CFA pattern descriptor for both 2x2 Bayer and 6x6 X-Trans sensors.
///
/// Uses a fixed-size array internally — no heap allocation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CfaPattern {
    pattern: [Channel; 36],
    width: usize,
    height: usize,
}

use Channel::*;

/// Standard X-Trans 6x6 pattern used by Fujifilm ILC cameras (X-T*, X-Pro*, X-H*, X-E*, X-S*).
const XTRANS_DEFAULT: [Channel; 36] = [
    Red,   Blue,  Green, Blue,  Red,   Green,
    Green, Green, Red,   Green, Green, Blue,
    Green, Green, Blue,  Green, Green, Red,
    Blue,  Red,   Green, Red,   Blue,  Green,
    Green, Green, Blue,  Green, Green, Red,
    Green, Green, Red,   Green, Green, Blue,
];

impl CfaPattern {
    /// Create a Bayer RGGB pattern.
    pub fn bayer_rggb() -> Self {
        Self::bayer([Red, Green, Green, Blue])
    }

    /// Create a Bayer BGGR pattern.
    pub fn bayer_bggr() -> Self {
        Self::bayer([Blue, Green, Green, Red])
    }

    /// Create a Bayer GRBG pattern.
    pub fn bayer_grbg() -> Self {
        Self::bayer([Green, Red, Blue, Green])
    }

    /// Create a Bayer GBRG pattern.
    pub fn bayer_gbrg() -> Self {
        Self::bayer([Green, Blue, Red, Green])
    }

    fn bayer(pat: [Channel; 4]) -> Self {
        let mut pattern = [Green; 36];
        pattern[0] = pat[0];
        pattern[1] = pat[1];
        pattern[2] = pat[2];
        pattern[3] = pat[3];
        Self { pattern, width: 2, height: 2 }
    }

    /// Create a custom X-Trans 6x6 pattern.
    pub fn xtrans(pattern: [Channel; 36]) -> Self {
        Self { pattern, width: 6, height: 6 }
    }

    /// Standard Fujifilm X-Trans ILC pattern.
    pub fn xtrans_default() -> Self {
        Self { pattern: XTRANS_DEFAULT, width: 6, height: 6 }
    }

    /// Return a shifted view of this pattern.
    ///
    /// The shift (dy, dx) follows the additive convention:
    /// `shifted.color_at(row, col) == self.color_at(row + dy, col + dx)`
    pub fn shift(&self, dy: usize, dx: usize) -> Self {
        let mut shifted = [Green; 36];
        for y in 0..self.height {
            for x in 0..self.width {
                shifted[y * self.width + x] =
                    self.pattern[((y + dy) % self.height) * self.width + ((x + dx) % self.width)];
            }
        }
        Self { pattern: shifted, width: self.width, height: self.height }
    }

    /// Return the color channel at the given row and column (wraps modulo pattern size).
    #[inline]
    pub fn color_at(&self, row: usize, col: usize) -> Channel {
        self.pattern[(row % self.height) * self.width + (col % self.width)]
    }

    /// Pattern width: 2 for Bayer, 6 for X-Trans.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Pattern height: 2 for Bayer, 6 for X-Trans.
    pub fn height(&self) -> usize {
        self.height
    }

    /// Returns `true` if this is a 2x2 Bayer pattern.
    pub fn is_bayer(&self) -> bool {
        self.width == 2
    }

    /// Returns `true` if this is a 6x6 X-Trans pattern.
    pub fn is_xtrans(&self) -> bool {
        self.width == 6
    }
}

impl fmt::Display for CfaPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_bayer() {
            for i in 0..4 {
                write!(f, "{}", self.pattern[i])?;
            }
            Ok(())
        } else {
            write!(f, "X-Trans 6x6")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bayer_rggb_pattern() {
        let cfa = CfaPattern::bayer_rggb();
        assert_eq!(cfa.color_at(0, 0), Red);
        assert_eq!(cfa.color_at(0, 1), Green);
        assert_eq!(cfa.color_at(1, 0), Green);
        assert_eq!(cfa.color_at(1, 1), Blue);
        // Verify tiling
        assert_eq!(cfa.color_at(2, 2), Red);
        assert_eq!(cfa.color_at(3, 3), Blue);
    }

    #[test]
    fn xtrans_default_pattern() {
        let cfa = CfaPattern::xtrans_default();
        assert_eq!(cfa.width(), 6);
        assert_eq!(cfa.height(), 6);
        // Top-left corner: R B G B R G
        assert_eq!(cfa.color_at(0, 0), Red);
        assert_eq!(cfa.color_at(0, 1), Blue);
        assert_eq!(cfa.color_at(0, 2), Green);
        // Verify tiling
        assert_eq!(cfa.color_at(6, 0), Red);
        assert_eq!(cfa.color_at(0, 6), Red);
    }

    #[test]
    fn xtrans_shift() {
        let cfa = CfaPattern::xtrans_default();
        let shifted = cfa.shift(1, 1);
        // shifted.color_at(0,0) == cfa.color_at(1,1)
        assert_eq!(shifted.color_at(0, 0), cfa.color_at(1, 1));
        assert_eq!(shifted.color_at(2, 3), cfa.color_at(3, 4));
    }

    #[test]
    fn xtrans_every_3x3_has_all_colors() {
        let cfa = CfaPattern::xtrans_default();
        for by in 0..6 {
            for bx in 0..6 {
                let mut has = [false; 3];
                for y in 0..3 {
                    for x in 0..3 {
                        has[cfa.color_at(by + y, bx + x) as usize] = true;
                    }
                }
                assert!(has[0] && has[1] && has[2],
                    "3x3 block at ({by},{bx}) missing a color");
            }
        }
    }
}
