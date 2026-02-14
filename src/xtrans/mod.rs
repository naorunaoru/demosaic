mod bilinear_impl;
mod markesteijn_impl;
mod dht_impl;

use crate::CfaPattern;

pub fn bilinear(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    bilinear_impl::demosaic(input, width, height, cfa, output);
}

pub fn markesteijn1(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    markesteijn_impl::demosaic(input, width, height, cfa, output);
}

pub fn markesteijn3(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    markesteijn_impl::demosaic3(input, width, height, cfa, output);
}

pub fn dht(input: &[f32], width: usize, height: usize, cfa: &CfaPattern, output: &mut [f32]) {
    dht_impl::demosaic(input, width, height, cfa, output);
}
