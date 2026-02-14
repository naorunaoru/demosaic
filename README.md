# demosaic

[![crates.io](https://img.shields.io/crates/v/demosaic.svg)](https://crates.io/crates/demosaic)
[![docs.rs](https://docs.rs/demosaic/badge.svg)](https://docs.rs/demosaic)
[![CI](https://github.com/naorunaoru/demosaic/actions/workflows/ci.yml/badge.svg)](https://github.com/naorunaoru/demosaic/actions/workflows/ci.yml)

Demosaicing algorithms for Bayer, Quad Bayer, and X-Trans CFA sensors. Pure Rust, `no_std` compatible. *Not* blazingly fast, just moderately. Use at your own pace.

## Algorithms

| Algorithm | Bayer | Quad Bayer | X-Trans | Quality | Speed |
|-----------|:-----:|:----------:|:-------:|---------|-------|
| Bilinear | x | x | x | Low | Fastest |
| MHC | x | | | Medium | Fast |
| PPG | x | | | Medium-high | Medium |
| AHD | x | | | High | Slow |
| Quad-PPG | | x | | Medium-high | Medium |
| Markesteijn (1-pass) | | | x | Medium-high | Medium |
| Markesteijn (3-pass) | | | x | High | Slow |
| DHT | | | x | Medium-high | Medium |

Quad Bayer sensors also support two additional processing strategies:

| Strategy | Output resolution | Description |
|----------|-------------------|-------------|
| Binning + Bayer | Half (W/2 x H/2) | Average 2x2 same-color blocks, then apply any Bayer algorithm |
| Remosaic + Bayer | Full (W x H) | Convert to standard Bayer, then apply any Bayer algorithm |

## Usage

```rust
use demosaic::{demosaic, Algorithm, CfaPattern};

let width = 4032;
let height = 3024;
let cfa = CfaPattern::bayer_rggb();
let input: Vec<f32> = load_raw_sensor_data(); // single-channel, row-major

let mut output = vec![0.0f32; 3 * width * height];
demosaic(&input, width, height, &cfa, Algorithm::Ppg, &mut output).unwrap();

// output is planar CHW: [R plane (w*h), G plane (w*h), B plane (w*h)]
```

For interleaved RGB output (`[R,G,B, R,G,B, ...]`), use `demosaic_interleaved`:

```rust
use demosaic::demosaic_interleaved;

demosaic_interleaved(&input, width, height, &cfa, Algorithm::Ppg, &mut output).unwrap();

// output is interleaved HWC: [R,G,B, R,G,B, ...]
```

Standalone conversion utilities are also available: `planar_to_interleaved` and `interleaved_to_planar`.

Input and output are `f32` slices. The crate does not handle raw file parsing or color
space conversion.

### CFA patterns

```rust
// Bayer (2x2)
let cfa = CfaPattern::bayer_rggb();  // also bayer_bggr, bayer_grbg, bayer_gbrg

// Quad Bayer (4x4)
let cfa = CfaPattern::quad_bayer_rggb();  // also quad_bayer_bggr, quad_bayer_grbg, quad_bayer_gbrg

// X-Trans (6x6)
let cfa = CfaPattern::xtrans_default();  // standard Fujifilm ILC pattern
let cfa = CfaPattern::xtrans(pattern);   // custom 6x6 pattern

// Shift for top-left crop offset
let cfa = CfaPattern::bayer_rggb().shift(dy, dx);
```

### Quad Bayer

Quad Bayer sensors (Sony IMX586, Samsung ISOCELL GN1, etc.) use a 4x4 repeating pattern where each standard Bayer pixel is replaced by a 2x2 same-color block:

```
R R G G
R R G G
G G B B
G G B B
```

There are four ways to process Quad Bayer data, depending on your quality/resolution needs:

**Direct demosaic** (full resolution):
```rust
use demosaic::{demosaic, Algorithm, CfaPattern};

let cfa = CfaPattern::quad_bayer_rggb();
let mut output = vec![0.0f32; 3 * width * height];

// Bilinear (fast, low quality)
demosaic(&input, width, height, &cfa, Algorithm::Bilinear, &mut output).unwrap();

// Quad-PPG (gradient-based, higher quality)
demosaic(&input, width, height, &cfa, Algorithm::QuadPpg, &mut output).unwrap();
```

**Binning** (half resolution, best SNR):
```rust
use demosaic::{demosaic_quad_binned, BayerAlgorithm, CfaPattern};

let cfa = CfaPattern::quad_bayer_rggb();
let (ow, oh) = (width / 2, height / 2);
let mut output = vec![0.0f32; 3 * ow * oh];

demosaic_quad_binned(&input, width, height, &cfa, BayerAlgorithm::Ppg, &mut output).unwrap();
```

**Remosaic** (full resolution, leverages all existing Bayer algorithms):
```rust
use demosaic::{demosaic, remosaic, Algorithm, CfaPattern};

let cfa = CfaPattern::quad_bayer_rggb();
let bayer_cfa = cfa.quad_to_bayer().unwrap();

// Step 1: convert Quad Bayer to standard Bayer
let mut bayer_data = vec![0.0f32; width * height];
remosaic(&input, width, height, &cfa, &mut bayer_data).unwrap();

// Step 2: demosaic with any Bayer algorithm
let mut output = vec![0.0f32; 3 * width * height];
demosaic(&bayer_data, width, height, &bayer_cfa, Algorithm::Ahd, &mut output).unwrap();
```

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | Yes | Implements `std::error::Error` for `DemosaicError` |

Disable default features for `no_std` (requires `alloc`):

```toml
[dependencies]
demosaic = { version = "0.1", default-features = false }
```

## License

MIT OR Apache-2.0
