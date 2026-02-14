mod bilinear_impl;
mod qppg_impl;
pub(crate) mod binning;

use crate::CfaPattern;

pub fn bilinear(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    bilinear_impl::demosaic(input, width, height, cfa, output);
}

pub fn qppg(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    qppg_impl::demosaic(input, width, height, cfa, output);
}
