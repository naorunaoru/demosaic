# demosaic

Demosaicing algorithms for Bayer and X-Trans CFA sensors. Pure Rust, `no_std` compatible. *Not* blazingly fast, just moderately. Use at your own pace. 

## Algorithms

| Algorithm | Bayer | X-Trans | Quality | Speed |
|-----------|:-----:|:-------:|---------|-------|
| Bilinear | x | x | Low | Fastest |
| MHC | x | | Medium | Fast |
| PPG | x | | Medium-high | Medium |
| AHD | x | | High | Slow |
| Markesteijn (1-pass) | | x | Medium-high | Medium |
| Markesteijn (3-pass) | | x | High | Slow |
| DHT | | x | Medium-high | Medium |

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

// X-Trans (6x6)
let cfa = CfaPattern::xtrans_default();  // standard Fujifilm ILC pattern
let cfa = CfaPattern::xtrans(pattern);   // custom 6x6 pattern

// Shift for top-left crop offset
let cfa = CfaPattern::bayer_rggb().shift(dy, dx);
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
