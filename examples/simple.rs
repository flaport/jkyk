use anyhow::Result;
use std::fs::File;
use std::io::Write;

const M: isize = 160;
const N: isize = 96;
const P: isize = 64;
const C: isize = 3;
const F: isize = 2;
const Q: isize = 96;

const SX: isize = M / 2;
const SY: isize = N / 2;
const SZ: isize = P / 2;

fn ramped_sin(omega: f32, width: f32, delay: f32, q: usize) -> Vec<f32> {
    let mut result = Vec::with_capacity(q);
    for i in 0..q {
        let t = i as f32 * omega;
        let value = ((1.0 + (t / width - delay).tanh()) / 2.0) * t.sin();
        result.push(value);
    }
    result
}

#[inline]
fn ix(m: isize, n: isize, p: isize, c: isize, f: isize) -> usize {
    ((((f * M + (m + M) % M) * N + (n + N) % N) * P + (p + P) % P) * C + c) as usize
}

fn run(eh: &mut [f32], st: &[f32]) {
    let sc = 0.99 / (3.0_f32).sqrt();
    let mut c1: isize;
    let mut c2: isize;
    let mut f: isize;
    let mut g: isize;
    let mut s: isize;
    let mut j0: isize;
    let mut j1: isize;
    let mut j2: isize;

    for i in 0..(2 * Q) {
        f = i % 2;
        g = 1 - f;
        s = 1 - (2 * f);

        for m in 0..M {
            for n in 0..N {
                for p in 0..P {
                    for c0 in 0..C {
                        c1 = (c0 + 1) % 3;
                        c2 = (c0 + 2) % 3;
                        j0 = s * ((c0 == 0) as isize);
                        j1 = s * ((c0 == 1) as isize);
                        j2 = s * ((c0 == 2) as isize);
                        eh[ix(m, n, p, c0, f)] += sc
                            * (eh[ix(m, n, p, c2, g)]
                                - eh[ix(m - j2, n - j0, p - j1, c2, g)]
                                - eh[ix(m, n, p, c1, g)]
                                + eh[ix(m - j1, n - j2, p - j0, c1, g)]);
                    }
                }
            }
        }
        if f == 0 {
            eh[ix(SX, SY, SZ, 2, 0)] += st[i as usize / 2];
        }
    }
}

fn main() -> Result<()> {
    let mut eh = vec![0.0_f32; (M * N * P * C * F) as usize];
    let st = ramped_sin(0.3, 5.0, 3.0, Q as usize);

    run(&mut eh, &st);

    // Write output to binary file
    let mut file = File::create("output.bin")?;

    // Convert f32 slice to bytes and write
    let byte_data: Vec<u8> = eh[..(2 * M * N * P * 3) as usize]
        .iter()
        .flat_map(|&f| f.to_le_bytes())
        .collect();

    file.write_all(&byte_data)?;

    Ok(())
}
