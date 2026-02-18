mod bilinear_impl;
mod mhc_impl;
mod ppg_impl;
mod ahd_impl;
mod vng_impl;

use crate::CfaPattern;

pub fn bilinear(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    bilinear_impl::demosaic(input, width, height, cfa, output);
}

pub fn mhc(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    mhc_impl::demosaic(input, width, height, cfa, output);
}

pub fn ppg(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    ppg_impl::demosaic(input, width, height, cfa, output);
}

pub fn ahd(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    ahd_impl::demosaic(input, width, height, cfa, output);
}

pub fn vng(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    vng_impl::demosaic(input, width, height, cfa, output);
}
