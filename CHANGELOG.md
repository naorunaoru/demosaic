# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0]

### Added

- VNG (Variable Number of Gradients) demosaicing algorithm for Bayer patterns.
  8-directional gradient analysis with adaptive threshold selection.

## [0.2.0]

### Added

- Quad Bayer (4x4) CFA pattern support (`CfaPattern::quad_bayer_rggb()`, etc.).
- Quad-PPG direct demosaicing algorithm for Quad Bayer sensors.
- Bilinear demosaicing for Quad Bayer sensors (5x5 window).
- 2x2 binning + Bayer demosaic pipeline (`demosaic_quad_binned`).
- Remosaic: Quad Bayer to standard Bayer conversion (`remosaic`).
- `CfaPattern::quad_to_bayer()` for deriving the standard Bayer pattern from a Quad Bayer pattern.
- Interleaved (HWC) output via `demosaic_interleaved` and `demosaic_quad_binned_interleaved`.
- `planar_to_interleaved` and `interleaved_to_planar` conversion utilities.
- `bin2x2` standalone 2x2 binning utility.
- `BayerAlgorithm` enum for selecting Bayer-only algorithms in binning pipelines.

## [0.1.0]

### Added

- Bayer (2x2) CFA pattern support (RGGB, BGGR, GRBG, GBRG) with `CfaPattern::shift()` for crop offsets.
- X-Trans (6x6) CFA pattern support with standard Fujifilm ILC pattern and custom patterns.
- Bayer demosaicing: Bilinear, MHC (Malvar-He-Cutler), PPG (Patterned Pixel Grouping), AHD (Adaptive Homogeneity-Directed).
- X-Trans demosaicing: Bilinear, Markesteijn 1-pass, Markesteijn 3-pass, DHT (Directional Homogeneity Test).
- `no_std` support (requires `alloc`).
- `std` feature (default) for `std::error::Error` implementation.

[Unreleased]: https://github.com/naorunaoru/demosaic/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/naorunaoru/demosaic/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/naorunaoru/demosaic/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/naorunaoru/demosaic/releases/tag/v0.1.0
