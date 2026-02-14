use criterion::{criterion_group, criterion_main, Criterion};
use demosaic::{demosaic, Algorithm, CfaPattern, Channel};

fn synthetic_input(width: usize, height: usize, cfa: &CfaPattern) -> Vec<f32> {
    let mut input = vec![0.0f32; width * height];
    for y in 0..height {
        for x in 0..width {
            input[y * width + x] = match cfa.color_at(y, x) {
                Channel::Red => 0.7,
                Channel::Green => 0.5,
                Channel::Blue => 0.3,
            };
        }
    }
    input
}

fn bench_bayer(c: &mut Criterion) {
    let cfa = CfaPattern::bayer_rggb();
    let (w, h) = (4032, 3024);
    let input = synthetic_input(w, h, &cfa);
    let mut output = vec![0.0f32; 3 * w * h];

    let mut group = c.benchmark_group("bayer_4032x3024");
    for algo in [Algorithm::Bilinear, Algorithm::Mhc, Algorithm::Ppg, Algorithm::Ahd] {
        group.bench_function(format!("{algo}"), |b| {
            b.iter(|| demosaic(&input, w, h, &cfa, algo, &mut output).unwrap());
        });
    }
    group.finish();
}

fn bench_xtrans(c: &mut Criterion) {
    let cfa = CfaPattern::xtrans_default();
    let (w, h) = (4032, 3024);
    let input = synthetic_input(w, h, &cfa);
    let mut output = vec![0.0f32; 3 * w * h];

    let mut group = c.benchmark_group("xtrans_4032x3024");
    group.sample_size(10);
    for algo in [Algorithm::Bilinear, Algorithm::Markesteijn1, Algorithm::Markesteijn3, Algorithm::Dht] {
        group.bench_function(format!("{algo}"), |b| {
            b.iter(|| demosaic(&input, w, h, &cfa, algo, &mut output).unwrap());
        });
    }
    group.finish();
}

criterion_group!(benches, bench_bayer, bench_xtrans);
criterion_main!(benches);
