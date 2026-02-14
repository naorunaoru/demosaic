use crate::cfa::{CfaPattern, Channel};

/// Malvar-He-Cutler (2004) gradient-corrected linear demosaicing.
///
/// Uses 5×5 convolution kernels with integer coefficients to interpolate
/// missing color channels. Each CFA position has specific kernels that
/// incorporate gradient correction from the known channel, producing much
/// better results than simple bilinear averaging.
///
/// Reference: H.S. Malvar, L. He, R. Cutler, "High-Quality Linear
/// Interpolation for Demosaicing of Bayer-Patterned Color Images", ICASSP 2004.
pub fn demosaic(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    let npix = width * height;

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let ch = cfa.color_at(y, x);

            let r_off = 0;
            let g_off = npix;
            let b_off = 2 * npix;

            // Known channel: copy directly
            output[ch as usize * npix + idx] = input[idx];

            match ch {
                Channel::Red => {
                    output[g_off + idx] =
                        apply_kernel(input, width, height, y, x, &KERNEL_G_AT_RB, 8.0);
                    output[b_off + idx] =
                        apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_OPPOSITE, 16.0);
                }
                Channel::Blue => {
                    output[g_off + idx] =
                        apply_kernel(input, width, height, y, x, &KERNEL_G_AT_RB, 8.0);
                    output[r_off + idx] =
                        apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_OPPOSITE, 16.0);
                }
                Channel::Green => {
                    // Determine row type: is the adjacent pixel in this row R or B?
                    let adj_x = x ^ 1;
                    let adj_ch = cfa.color_at(y, adj_x);

                    if adj_ch == Channel::Red {
                        // Green in red row: R has same-row neighbors, B has same-column
                        output[r_off + idx] =
                            apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_G_SAME_ROW, 16.0);
                        output[b_off + idx] =
                            apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_G_SAME_COL, 16.0);
                    } else {
                        // Green in blue row: B has same-row neighbors, R has same-column
                        output[b_off + idx] =
                            apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_G_SAME_ROW, 16.0);
                        output[r_off + idx] =
                            apply_kernel(input, width, height, y, x, &KERNEL_RB_AT_G_SAME_COL, 16.0);
                    }
                }
            }
        }
    }
}

#[inline]
fn apply_kernel(
    input: &[f32],
    width: usize,
    height: usize,
    y: usize,
    x: usize,
    kernel: &[i32; 25],
    divisor: f32,
) -> f32 {
    let mut sum = 0.0f32;
    let mut k = 0;
    for ky in -2i32..=2 {
        let ny = (y as i32 + ky).clamp(0, height as i32 - 1) as usize;
        for kx in -2i32..=2 {
            let nx = (x as i32 + kx).clamp(0, width as i32 - 1) as usize;
            sum += kernel[k] as f32 * input[ny * width + nx];
            k += 1;
        }
    }
    sum / divisor
}

// Green at Red or Blue location.
// Divisor: 8
#[rustfmt::skip]
const KERNEL_G_AT_RB: [i32; 25] = [
     0,  0, -1,  0,  0,
     0,  0,  2,  0,  0,
    -1,  2,  4,  2, -1,
     0,  0,  2,  0,  0,
     0,  0, -1,  0,  0,
];

// R/B at Green location where the same color has neighbors in the same ROW.
// e.g. Red at Green-in-Red-row, or Blue at Green-in-Blue-row.
// Divisor: 16
#[rustfmt::skip]
const KERNEL_RB_AT_G_SAME_ROW: [i32; 25] = [
     0,  0,  1,  0,  0,
     0, -2,  0, -2,  0,
    -2,  8, 10,  8, -2,
     0, -2,  0, -2,  0,
     0,  0,  1,  0,  0,
];

// R/B at Green location where the same color has neighbors in the same COLUMN.
// Transpose of KERNEL_RB_AT_G_SAME_ROW.
// Divisor: 16
#[rustfmt::skip]
const KERNEL_RB_AT_G_SAME_COL: [i32; 25] = [
     0,  0, -2,  0,  0,
     0, -2,  8, -2,  0,
     1,  0, 10,  0,  1,
     0, -2,  8, -2,  0,
     0,  0, -2,  0,  0,
];

// R at B location or B at R location (diagonal interpolation).
// Divisor: 16
#[rustfmt::skip]
const KERNEL_RB_AT_OPPOSITE: [i32; 25] = [
     0,  0, -3,  0,  0,
     0,  4,  0,  4,  0,
    -3,  0, 12,  0, -3,
     0,  4,  0,  4,  0,
     0,  0, -3,  0,  0,
];

#[cfg(test)]
mod tests {
    use alloc::vec;
    use super::*;

    #[test]
    fn known_channel_preserved() {
        for cfa in &[
            CfaPattern::bayer_rggb(),
            CfaPattern::bayer_bggr(),
            CfaPattern::bayer_grbg(),
            CfaPattern::bayer_gbrg(),
        ] {
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
            demosaic(&input, w, h, cfa, &mut output);

            for y in 0..h {
                for x in 0..w {
                    let idx = y * w + x;
                    let ch = cfa.color_at(y, x) as usize;
                    assert_eq!(
                        output[ch * w * h + idx], input[idx],
                        "known channel mismatch at ({y},{x}) ch={ch}"
                    );
                }
            }
        }
    }

    #[test]
    fn uniform_color_reconstructed() {
        for cfa in &[
            CfaPattern::bayer_rggb(),
            CfaPattern::bayer_bggr(),
            CfaPattern::bayer_grbg(),
            CfaPattern::bayer_gbrg(),
        ] {
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
            demosaic(&input, w, h, cfa, &mut output);

            // Interior pixels (skip 3-pixel border)
            for y in 3..h - 3 {
                for x in 3..w - 3 {
                    let idx = y * w + x;
                    let r = output[idx];
                    let g = output[w * h + idx];
                    let b = output[2 * w * h + idx];
                    assert!((r - 0.8).abs() < 0.05, "R at ({y},{x}) = {r}");
                    assert!((g - 0.5).abs() < 0.05, "G at ({y},{x}) = {g}");
                    assert!((b - 0.3).abs() < 0.05, "B at ({y},{x}) = {b}");
                }
            }
        }
    }
}
