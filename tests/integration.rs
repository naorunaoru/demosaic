use demosaic::{
    demosaic, demosaic_interleaved, demosaic_quad_binned, interleaved_to_planar,
    planar_to_interleaved, remosaic, Algorithm, BayerAlgorithm, CfaPattern, DemosaicError,
};

/// Helper: create CFA input where each pixel gets a value based on its filter color.
fn synthetic_input(width: usize, height: usize, cfa: &CfaPattern, rgb: [f32; 3]) -> Vec<f32> {
    let mut input = vec![0.0f32; width * height];
    for y in 0..height {
        for x in 0..width {
            input[y * width + x] = rgb[cfa.color_at(y, x) as usize];
        }
    }
    input
}

// ---------------------------------------------------------------------------
// Solid-color reconstruction: if every R pixel = every G pixel = every B pixel,
// all algorithms should reconstruct uniform planes (within tolerance).
// ---------------------------------------------------------------------------

fn assert_solid_reconstruction(
    algorithm: Algorithm,
    cfa: &CfaPattern,
    width: usize,
    height: usize,
    tolerance: f32,
) {
    let value = 0.5;
    let input = synthetic_input(width, height, cfa, [value; 3]);
    let mut output = vec![0.0f32; 3 * width * height];

    demosaic(&input, width, height, cfa, algorithm, &mut output).unwrap();

    let npix = width * height;
    // Skip a border of 8 pixels to avoid edge effects from algorithms that
    // don't fully interpolate the border region.
    let border = 8;
    for y in border..height.saturating_sub(border) {
        for x in border..width.saturating_sub(border) {
            let idx = y * width + x;
            for c in 0..3 {
                let v = output[c * npix + idx];
                assert!(
                    (v - value).abs() < tolerance,
                    "{:?} on {:?} CFA: pixel ({y},{x}) channel {c} = {v}, expected ~{value} (tol {tolerance})",
                    algorithm,
                    if cfa.is_bayer() { "Bayer" } else if cfa.is_quad_bayer() { "Quad Bayer" } else { "X-Trans" },
                );
            }
        }
    }
}

#[test]
fn solid_bayer_bilinear() {
    for cfa in &[
        CfaPattern::bayer_rggb(),
        CfaPattern::bayer_bggr(),
        CfaPattern::bayer_grbg(),
        CfaPattern::bayer_gbrg(),
    ] {
        assert_solid_reconstruction(Algorithm::Bilinear, cfa, 64, 64, 1e-6);
    }
}

#[test]
fn solid_bayer_mhc() {
    assert_solid_reconstruction(Algorithm::Mhc, &CfaPattern::bayer_rggb(), 64, 64, 1e-4);
}

#[test]
fn solid_bayer_ppg() {
    assert_solid_reconstruction(Algorithm::Ppg, &CfaPattern::bayer_rggb(), 64, 64, 1e-4);
}

#[test]
fn solid_bayer_ahd() {
    assert_solid_reconstruction(Algorithm::Ahd, &CfaPattern::bayer_rggb(), 64, 64, 1e-4);
}

#[test]
fn solid_bayer_vng() {
    assert_solid_reconstruction(Algorithm::Vng, &CfaPattern::bayer_rggb(), 64, 64, 1e-4);
}

#[test]
fn solid_xtrans_bilinear() {
    assert_solid_reconstruction(
        Algorithm::Bilinear,
        &CfaPattern::xtrans_default(),
        64,
        64,
        1e-6,
    );
}

#[test]
fn solid_xtrans_markesteijn1() {
    assert_solid_reconstruction(
        Algorithm::Markesteijn1,
        &CfaPattern::xtrans_default(),
        128,
        128,
        1e-2,
    );
}

#[test]
fn solid_xtrans_markesteijn3() {
    assert_solid_reconstruction(
        Algorithm::Markesteijn3,
        &CfaPattern::xtrans_default(),
        128,
        128,
        1e-2,
    );
}

#[test]
fn solid_xtrans_dht() {
    assert_solid_reconstruction(
        Algorithm::Dht,
        &CfaPattern::xtrans_default(),
        128,
        128,
        1e-2,
    );
}

// ---------------------------------------------------------------------------
// Color separation: given distinct per-channel values, the output planes
// should recover them (at least in the interior).
// ---------------------------------------------------------------------------

#[test]
fn color_separation_bayer() {
    let cfa = CfaPattern::bayer_rggb();
    let (w, h) = (64, 64);
    let rgb = [0.8, 0.5, 0.2];
    let input = synthetic_input(w, h, &cfa, rgb);
    let mut output = vec![0.0f32; 3 * w * h];

    demosaic(&input, w, h, &cfa, Algorithm::Bilinear, &mut output).unwrap();

    let npix = w * h;
    let border = 4;
    for y in border..h - border {
        for x in border..w - border {
            let idx = y * w + x;
            for c in 0..3 {
                let v = output[c * npix + idx];
                assert!(
                    (v - rgb[c]).abs() < 1e-4,
                    "pixel ({y},{x}) channel {c} = {v}, expected {}",
                    rgb[c]
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Error conditions
// ---------------------------------------------------------------------------

#[test]
fn error_input_size_mismatch() {
    let cfa = CfaPattern::bayer_rggb();
    let mut output = vec![0.0f32; 3 * 4 * 4];
    let result = demosaic(&[0.0; 10], 4, 4, &cfa, Algorithm::Bilinear, &mut output);
    assert!(matches!(
        result,
        Err(DemosaicError::InputSizeMismatch { expected: 16, got: 10 })
    ));
}

#[test]
fn error_output_size_mismatch() {
    let cfa = CfaPattern::bayer_rggb();
    let input = vec![0.0f32; 16];
    let mut output = vec![0.0f32; 10];
    let result = demosaic(&input, 4, 4, &cfa, Algorithm::Bilinear, &mut output);
    assert!(matches!(
        result,
        Err(DemosaicError::OutputSizeMismatch { expected: 48, got: 10 })
    ));
}

#[test]
fn error_unsupported_algorithm_bayer() {
    let cfa = CfaPattern::bayer_rggb();
    let input = vec![0.0f32; 16];
    let mut output = vec![0.0f32; 48];
    let result = demosaic(&input, 4, 4, &cfa, Algorithm::Markesteijn1, &mut output);
    assert!(matches!(
        result,
        Err(DemosaicError::UnsupportedAlgorithm { .. })
    ));
}

#[test]
fn error_unsupported_algorithm_xtrans() {
    let cfa = CfaPattern::xtrans_default();
    let (w, h) = (64, 64);
    let input = vec![0.0f32; w * h];
    let mut output = vec![0.0f32; 3 * w * h];
    let result = demosaic(&input, w, h, &cfa, Algorithm::Ppg, &mut output);
    assert!(matches!(
        result,
        Err(DemosaicError::UnsupportedAlgorithm { .. })
    ));
}

#[test]
fn error_xtrans_image_too_small() {
    let cfa = CfaPattern::xtrans_default();
    let (w, h) = (12, 12);
    let input = vec![0.0f32; w * h];
    let mut output = vec![0.0f32; 3 * w * h];
    let result = demosaic(&input, w, h, &cfa, Algorithm::Markesteijn1, &mut output);
    assert!(matches!(
        result,
        Err(DemosaicError::ImageTooSmall { .. })
    ));
}

// ---------------------------------------------------------------------------
// Interleaved output
// ---------------------------------------------------------------------------

#[test]
fn interleaved_matches_planar_bayer() {
    let cfa = CfaPattern::bayer_rggb();
    let (w, h) = (64, 64);
    let npix = w * h;
    let input = synthetic_input(w, h, &cfa, [0.8, 0.5, 0.2]);

    let mut planar = vec![0.0f32; 3 * npix];
    let mut interleaved = vec![0.0f32; 3 * npix];

    demosaic(&input, w, h, &cfa, Algorithm::Bilinear, &mut planar).unwrap();
    demosaic_interleaved(&input, w, h, &cfa, Algorithm::Bilinear, &mut interleaved).unwrap();

    for i in 0..npix {
        assert_eq!(interleaved[3 * i], planar[i], "R mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 1], planar[npix + i], "G mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 2], planar[2 * npix + i], "B mismatch at pixel {i}");
    }
}

#[test]
fn interleaved_matches_planar_xtrans() {
    let cfa = CfaPattern::xtrans_default();
    let (w, h) = (128, 128);
    let npix = w * h;
    let input = synthetic_input(w, h, &cfa, [0.7, 0.5, 0.3]);

    let mut planar = vec![0.0f32; 3 * npix];
    let mut interleaved = vec![0.0f32; 3 * npix];

    demosaic(&input, w, h, &cfa, Algorithm::Markesteijn1, &mut planar).unwrap();
    demosaic_interleaved(&input, w, h, &cfa, Algorithm::Markesteijn1, &mut interleaved).unwrap();

    for i in 0..npix {
        assert_eq!(interleaved[3 * i], planar[i], "R mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 1], planar[npix + i], "G mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 2], planar[2 * npix + i], "B mismatch at pixel {i}");
    }
}

// ---------------------------------------------------------------------------
// Layout conversion round-trips
// ---------------------------------------------------------------------------

#[test]
fn planar_interleaved_round_trip() {
    let npix = 100;
    let mut planar = vec![0.0f32; 3 * npix];
    for i in 0..3 * npix {
        planar[i] = i as f32;
    }

    let mut interleaved = vec![0.0f32; 3 * npix];
    planar_to_interleaved(&planar, &mut interleaved);

    // Verify interleaved layout.
    for i in 0..npix {
        assert_eq!(interleaved[3 * i], i as f32); // R
        assert_eq!(interleaved[3 * i + 1], (npix + i) as f32); // G
        assert_eq!(interleaved[3 * i + 2], (2 * npix + i) as f32); // B
    }

    // Round-trip back to planar.
    let mut back = vec![0.0f32; 3 * npix];
    interleaved_to_planar(&interleaved, &mut back);
    assert_eq!(planar, back);
}

// ---------------------------------------------------------------------------
// Quad Bayer: solid-color reconstruction
// ---------------------------------------------------------------------------

fn quad_bayer_cfas() -> [CfaPattern; 4] {
    [
        CfaPattern::quad_bayer_rggb(),
        CfaPattern::quad_bayer_bggr(),
        CfaPattern::quad_bayer_grbg(),
        CfaPattern::quad_bayer_gbrg(),
    ]
}

#[test]
fn solid_quad_bayer_bilinear() {
    for cfa in &quad_bayer_cfas() {
        assert_solid_reconstruction(Algorithm::Bilinear, cfa, 64, 64, 1e-6);
    }
}

#[test]
fn solid_quad_bayer_qppg() {
    for cfa in &quad_bayer_cfas() {
        assert_solid_reconstruction(Algorithm::QuadPpg, cfa, 64, 64, 1e-4);
    }
}

// ---------------------------------------------------------------------------
// Quad Bayer: color separation
// ---------------------------------------------------------------------------

#[test]
fn color_separation_quad_bayer() {
    let cfa = CfaPattern::quad_bayer_rggb();
    let (w, h) = (64, 64);
    let rgb = [0.8, 0.5, 0.2];
    let input = synthetic_input(w, h, &cfa, rgb);
    let mut output = vec![0.0f32; 3 * w * h];

    demosaic(&input, w, h, &cfa, Algorithm::Bilinear, &mut output).unwrap();

    let npix = w * h;
    let border = 4;
    for y in border..h - border {
        for x in border..w - border {
            let idx = y * w + x;
            for c in 0..3 {
                let v = output[c * npix + idx];
                assert!(
                    (v - rgb[c]).abs() < 0.05,
                    "pixel ({y},{x}) channel {c} = {v}, expected {}",
                    rgb[c]
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Quad Bayer: binning
// ---------------------------------------------------------------------------

#[test]
fn solid_quad_bayer_binned() {
    for cfa in &quad_bayer_cfas() {
        let (w, h) = (64, 64);
        let input = synthetic_input(w, h, cfa, [0.5; 3]);
        let ow = w / 2;
        let oh = h / 2;
        let out_npix = ow * oh;

        for alg in &[
            BayerAlgorithm::Bilinear,
            BayerAlgorithm::Mhc,
            BayerAlgorithm::Ppg,
            BayerAlgorithm::Ahd,
            BayerAlgorithm::Vng,
        ] {
            let mut output = vec![0.0f32; 3 * out_npix];
            demosaic_quad_binned(&input, w, h, cfa, *alg, &mut output).unwrap();

            let border = 4;
            for y in border..oh - border {
                for x in border..ow - border {
                    let idx = y * ow + x;
                    for c in 0..3 {
                        let v = output[c * out_npix + idx];
                        assert!(
                            (v - 0.5).abs() < 1e-4,
                            "{alg:?}: ch {c} at ({y},{x}) = {v}, expected ~0.5"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn quad_binned_dimensions() {
    let cfa = CfaPattern::quad_bayer_rggb();
    let (w, h) = (64, 64);
    let input = vec![0.5f32; w * h];
    let ow = w / 2;
    let oh = h / 2;

    // Correct output size
    let mut output = vec![0.0f32; 3 * ow * oh];
    assert!(demosaic_quad_binned(&input, w, h, &cfa, BayerAlgorithm::Bilinear, &mut output).is_ok());

    // Wrong output size
    let mut wrong = vec![0.0f32; 3 * w * h];
    assert!(matches!(
        demosaic_quad_binned(&input, w, h, &cfa, BayerAlgorithm::Bilinear, &mut wrong),
        Err(DemosaicError::OutputSizeMismatch { .. })
    ));
}

// ---------------------------------------------------------------------------
// Quad Bayer: remosaic
// ---------------------------------------------------------------------------

#[test]
fn remosaic_passthrough_preserved() {
    for cfa in &quad_bayer_cfas() {
        let (w, h) = (32, 32);
        let input = synthetic_input(w, h, cfa, [0.8, 0.5, 0.3]);
        let bayer_cfa = cfa.quad_to_bayer().unwrap();

        let mut output = vec![0.0f32; w * h];
        remosaic(&input, w, h, cfa, &mut output).unwrap();

        // Pass-through pixels should be exact
        for y in 0..h {
            for x in 0..w {
                let idx = y * w + x;
                if cfa.color_at(y, x) == bayer_cfa.color_at(y, x) {
                    assert_eq!(output[idx], input[idx],
                        "pass-through pixel ({y},{x}) changed");
                }
            }
        }
    }
}

#[test]
fn remosaic_then_demosaic() {
    let cfa = CfaPattern::quad_bayer_rggb();
    let (w, h) = (64, 64);
    let rgb = [0.8, 0.5, 0.3];
    let input = synthetic_input(w, h, &cfa, rgb);
    let bayer_cfa = cfa.quad_to_bayer().unwrap();

    // Remosaic
    let mut bayer_data = vec![0.0f32; w * h];
    remosaic(&input, w, h, &cfa, &mut bayer_data).unwrap();

    // Demosaic with PPG
    let mut output = vec![0.0f32; 3 * w * h];
    demosaic(&bayer_data, w, h, &bayer_cfa, Algorithm::Ppg, &mut output).unwrap();

    let npix = w * h;
    let border = 8;
    for y in border..h - border {
        for x in border..w - border {
            let idx = y * w + x;
            for c in 0..3 {
                let v = output[c * npix + idx];
                assert!(
                    (v - rgb[c]).abs() < 0.1,
                    "pixel ({y},{x}) channel {c} = {v}, expected {}",
                    rgb[c]
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Quad Bayer: error conditions
// ---------------------------------------------------------------------------

#[test]
fn error_unsupported_algorithm_quad_bayer() {
    let cfa = CfaPattern::quad_bayer_rggb();
    let (w, h) = (64, 64);
    let input = vec![0.0f32; w * h];
    let mut output = vec![0.0f32; 3 * w * h];

    // Bayer-only algorithms should be rejected
    for alg in &[Algorithm::Mhc, Algorithm::Ppg, Algorithm::Ahd] {
        assert!(matches!(
            demosaic(&input, w, h, &cfa, *alg, &mut output),
            Err(DemosaicError::UnsupportedAlgorithm { .. })
        ));
    }

    // X-Trans-only algorithms should be rejected
    for alg in &[Algorithm::Markesteijn1, Algorithm::Markesteijn3, Algorithm::Dht] {
        assert!(matches!(
            demosaic(&input, w, h, &cfa, *alg, &mut output),
            Err(DemosaicError::UnsupportedAlgorithm { .. })
        ));
    }
}

// ---------------------------------------------------------------------------
// Quad Bayer: interleaved output
// ---------------------------------------------------------------------------

#[test]
fn interleaved_matches_planar_quad_bayer() {
    let cfa = CfaPattern::quad_bayer_rggb();
    let (w, h) = (64, 64);
    let npix = w * h;
    let input = synthetic_input(w, h, &cfa, [0.7, 0.5, 0.3]);

    let mut planar = vec![0.0f32; 3 * npix];
    let mut interleaved = vec![0.0f32; 3 * npix];

    demosaic(&input, w, h, &cfa, Algorithm::Bilinear, &mut planar).unwrap();
    demosaic_interleaved(&input, w, h, &cfa, Algorithm::Bilinear, &mut interleaved).unwrap();

    for i in 0..npix {
        assert_eq!(interleaved[3 * i], planar[i], "R mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 1], planar[npix + i], "G mismatch at pixel {i}");
        assert_eq!(interleaved[3 * i + 2], planar[2 * npix + i], "B mismatch at pixel {i}");
    }
}
