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

    /// Create a Quad Bayer RGGB pattern (4x4 repeating tile).
    ///
    /// ```text
    /// R R G G
    /// R R G G
    /// G G B B
    /// G G B B
    /// ```
    pub fn quad_bayer_rggb() -> Self {
        Self::quad_bayer([Red, Green, Green, Blue])
    }

    /// Create a Quad Bayer BGGR pattern (4x4 repeating tile).
    ///
    /// ```text
    /// B B G G
    /// B B G G
    /// G G R R
    /// G G R R
    /// ```
    pub fn quad_bayer_bggr() -> Self {
        Self::quad_bayer([Blue, Green, Green, Red])
    }

    /// Create a Quad Bayer GRBG pattern (4x4 repeating tile).
    ///
    /// ```text
    /// G G R R
    /// G G R R
    /// B B G G
    /// B B G G
    /// ```
    pub fn quad_bayer_grbg() -> Self {
        Self::quad_bayer([Green, Red, Blue, Green])
    }

    /// Create a Quad Bayer GBRG pattern (4x4 repeating tile).
    ///
    /// ```text
    /// G G B B
    /// G G B B
    /// R R G G
    /// R R G G
    /// ```
    pub fn quad_bayer_gbrg() -> Self {
        Self::quad_bayer([Green, Blue, Red, Green])
    }

    fn quad_bayer(base: [Channel; 4]) -> Self {
        let mut pattern = [Green; 36];
        for r in 0..2 {
            for c in 0..2 {
                let ch = base[r * 2 + c];
                let r0 = 2 * r;
                let c0 = 2 * c;
                pattern[r0 * 4 + c0] = ch;
                pattern[r0 * 4 + c0 + 1] = ch;
                pattern[(r0 + 1) * 4 + c0] = ch;
                pattern[(r0 + 1) * 4 + c0 + 1] = ch;
            }
        }
        Self { pattern, width: 4, height: 4 }
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

    /// Returns `true` if this is a 4x4 Quad Bayer pattern.
    pub fn is_quad_bayer(&self) -> bool {
        self.width == 4 && self.height == 4
    }

    /// Returns `true` if this is a 6x6 X-Trans pattern.
    pub fn is_xtrans(&self) -> bool {
        self.width == 6
    }

    /// For a Quad Bayer pattern, derive the corresponding 2x2 Bayer pattern
    /// that results from 2x2 binning.
    ///
    /// Returns `None` if this is not a Quad Bayer pattern.
    pub fn quad_to_bayer(&self) -> Option<CfaPattern> {
        if !self.is_quad_bayer() {
            return None;
        }
        Some(Self::bayer([
            self.color_at(0, 0),
            self.color_at(0, 2),
            self.color_at(2, 0),
            self.color_at(2, 2),
        ]))
    }
}

impl fmt::Display for CfaPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_bayer() {
            for i in 0..4 {
                write!(f, "{}", self.pattern[i])?;
            }
            Ok(())
        } else if self.is_quad_bayer() {
            write!(f, "Quad Bayer ")?;
            // Show the 2x2 base pattern (top-left of each quadrant).
            write!(f, "{}{}{}{}", self.pattern[0], self.pattern[2], self.pattern[8], self.pattern[10])
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
    fn quad_bayer_rggb_pattern() {
        let cfa = CfaPattern::quad_bayer_rggb();
        assert_eq!(cfa.width(), 4);
        assert_eq!(cfa.height(), 4);
        assert!(cfa.is_quad_bayer());
        assert!(!cfa.is_bayer());
        assert!(!cfa.is_xtrans());

        // Row 0: R R G G
        assert_eq!(cfa.color_at(0, 0), Red);
        assert_eq!(cfa.color_at(0, 1), Red);
        assert_eq!(cfa.color_at(0, 2), Green);
        assert_eq!(cfa.color_at(0, 3), Green);
        // Row 1: R R G G
        assert_eq!(cfa.color_at(1, 0), Red);
        assert_eq!(cfa.color_at(1, 1), Red);
        assert_eq!(cfa.color_at(1, 2), Green);
        assert_eq!(cfa.color_at(1, 3), Green);
        // Row 2: G G B B
        assert_eq!(cfa.color_at(2, 0), Green);
        assert_eq!(cfa.color_at(2, 1), Green);
        assert_eq!(cfa.color_at(2, 2), Blue);
        assert_eq!(cfa.color_at(2, 3), Blue);
        // Row 3: G G B B
        assert_eq!(cfa.color_at(3, 0), Green);
        assert_eq!(cfa.color_at(3, 1), Green);
        assert_eq!(cfa.color_at(3, 2), Blue);
        assert_eq!(cfa.color_at(3, 3), Blue);

        // Tiling
        assert_eq!(cfa.color_at(4, 0), Red);
        assert_eq!(cfa.color_at(0, 4), Red);
        assert_eq!(cfa.color_at(6, 6), Blue);
    }

    #[test]
    fn quad_bayer_all_variants() {
        let rggb = CfaPattern::quad_bayer_rggb();
        let bggr = CfaPattern::quad_bayer_bggr();
        let grbg = CfaPattern::quad_bayer_grbg();
        let gbrg = CfaPattern::quad_bayer_gbrg();

        // Top-left 2x2 block color
        assert_eq!(rggb.color_at(0, 0), Red);
        assert_eq!(bggr.color_at(0, 0), Blue);
        assert_eq!(grbg.color_at(0, 0), Green);
        assert_eq!(gbrg.color_at(0, 0), Green);

        // Bottom-right 2x2 block color
        assert_eq!(rggb.color_at(2, 2), Blue);
        assert_eq!(bggr.color_at(2, 2), Red);
        assert_eq!(grbg.color_at(2, 2), Green);
        assert_eq!(gbrg.color_at(2, 2), Green);

        // All should be quad Bayer
        for cfa in &[&rggb, &bggr, &grbg, &gbrg] {
            assert!(cfa.is_quad_bayer());
        }
    }

    #[test]
    fn quad_bayer_shift() {
        let cfa = CfaPattern::quad_bayer_rggb();
        let shifted = cfa.shift(1, 1);
        // shifted.color_at(0,0) == cfa.color_at(1,1)
        assert_eq!(shifted.color_at(0, 0), cfa.color_at(1, 1));
        assert_eq!(shifted.color_at(2, 3), cfa.color_at(3, 0));
    }

    #[test]
    fn quad_to_bayer_mapping() {
        assert_eq!(
            CfaPattern::quad_bayer_rggb().quad_to_bayer().unwrap(),
            CfaPattern::bayer_rggb()
        );
        assert_eq!(
            CfaPattern::quad_bayer_bggr().quad_to_bayer().unwrap(),
            CfaPattern::bayer_bggr()
        );
        assert_eq!(
            CfaPattern::quad_bayer_grbg().quad_to_bayer().unwrap(),
            CfaPattern::bayer_grbg()
        );
        assert_eq!(
            CfaPattern::quad_bayer_gbrg().quad_to_bayer().unwrap(),
            CfaPattern::bayer_gbrg()
        );
    }

    #[test]
    fn quad_to_bayer_returns_none_for_non_quad() {
        assert!(CfaPattern::bayer_rggb().quad_to_bayer().is_none());
        assert!(CfaPattern::xtrans_default().quad_to_bayer().is_none());
    }

    #[test]
    fn quad_bayer_color_count_per_tile() {
        // Each 4x4 tile should have 4R, 8G, 4B
        let cfa = CfaPattern::quad_bayer_rggb();
        let mut counts = [0u32; 3];
        for r in 0..4 {
            for c in 0..4 {
                counts[cfa.color_at(r, c) as usize] += 1;
            }
        }
        assert_eq!(counts[Red as usize], 4);
        assert_eq!(counts[Green as usize], 8);
        assert_eq!(counts[Blue as usize], 4);
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
