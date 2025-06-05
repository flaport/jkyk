use anyhow::Result;
use std::fs::File;
use std::io::Write;

const M: isize = 160;
const N: isize = 96;
const P: isize = 64;
const C: isize = 3;
const F: isize = 2;
const Q: isize = 96;
const W: isize = 8; // tile width
const H: isize = 8; // tile height
const W2: isize = W / 2;
const W32: isize = 3 * W / 2;

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

#[inline]
fn wix(m: isize, n: isize, p: isize, c: isize, f: isize) -> usize {
    ((((f * W32 + (m + W32) % W32) * W32 + (n + W32) % W32) * W32 + (p + W32) % W32) * C + c)
        as usize
}

fn main() -> Result<()> {
    let sc = 0.99 / (3.0_f32).sqrt();
    let mut eh = vec![0.0_f32; (M * N * P * C * F) as usize];
    let st = ramped_sin(0.3, 5.0, 3.0, Q as usize);

    let oms = M / W; // offset idxs in x/m direction
    let ons = N / W; // offset idxs in y/n direction
    let ops = P / W; // offset idxs in z/p direction

    let mut c1: isize;
    let mut c2: isize;
    let mut f: isize;
    let mut g: isize;
    let mut s: isize;
    let mut j0: isize;
    let mut j1: isize;
    let mut j2: isize;
    let mut fast = [0_f32; (3 * W / 2 * 3 * W / 2 * 3 * W / 2 * C * F) as usize];

    for mut i in 0..(2 * Q / H) {
        i *= H;
        for mut om in 0..oms {
            om *= W; // om: offset in x/m direction
            for mut on in 0..ons {
                on *= W; // on: offset in y/n direction
                for mut op in 0..ops {
                    op *= W; // op: offset in z/p direction
                    fill_fast(&mut fast, &eh, om, on, op);
                    for h in 0..H {
                        f = h % 2;
                        g = 1 - f;
                        s = (2 * f) - 1; // sign: -1 for E [f=0], +1 for H [f=1]
                        let m0 = W2 - (h + 1) / 2;
                        let m1 = W32 - (h + 1) / 2;
                        let n0 = W2 - (h + 1) / 2;
                        let n1 = W32 - (h + 1) / 2;
                        let p0 = W2 - (h + 1) / 2;
                        let p1 = W32 - (h + 1) / 2;
                        for m in m0..m1 {
                            for n in n0..n1 {
                                for p in p0..p1 {
                                    for c0 in 0..C {
                                        c1 = (c0 + 1) % 3;
                                        c2 = (c0 + 2) % 3;
                                        j0 = s * ((c0 == 0) as isize); // negative for E [f=0], positive for H [f=1]
                                        j1 = s * ((c0 == 1) as isize); // negative for E [f=0], positive for H [f=1]
                                        j2 = s * ((c0 == 2) as isize); // negative for E [f=0], positive for H [f=1]
                                        fast[wix(m, n, p, c0, f)] += sc
                                            * (fast[wix(m, n, p, c2, g)]
                                                - fast[wix(m + j2, n + j0, p + j1, c2, g)]
                                                - fast[wix(m, n, p, c1, g)]
                                                + fast[wix(m + j1, n + j2, p + j0, c1, g)]);
                                    }
                                    if f == 0
                                        && m == SX - om + W2
                                        && n == SY - on + W2
                                        && p == SZ - op + W2
                                    {
                                        fast[wix(
                                            SX - om + W2,
                                            SY - on + W2,
                                            SZ - op + W2,
                                            2,
                                            0,
                                        )] += st[(i + h) as usize / 2];
                                    }
                                }
                            }
                        }
                    }
                    extract_fast(&mut eh, &fast, om, on, op);
                }
            }
        }
    }

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

fn fill_fast(fast: &mut [f32], eh: &[f32], om: isize, on: isize, op: isize) {
    fast.fill(0.0);
    for f in 0..F {
        for (i, m) in (om - W2..om + W).enumerate() {
            if m < 0 || m >= M {
                continue;
            }
            for (j, n) in (on - W2..on + W).enumerate() {
                if n < 0 || n >= N {
                    continue;
                }
                for (k, p) in (op - W2..op + W).enumerate() {
                    if p < 0 || p >= P {
                        continue;
                    }
                    for c0 in 0..C {
                        fast[wix(i as isize, j as isize, k as isize, c0, f)] =
                            eh[ix(m, n, p, c0, f)];
                    }
                }
            }
        }
    }
}

fn extract_fast(eh: &mut [f32], fast: &[f32], om: isize, on: isize, op: isize) {
    for f in 0..F {
        for (i, m) in (om - W2..om + W).enumerate() {
            if m < 0 || m >= M {
                continue;
            }
            for (j, n) in (on - W2..on + W).enumerate() {
                if n < 0 || n >= N {
                    continue;
                }
                for (k, p) in (op - W2..op + W).enumerate() {
                    if p < 0 || p >= P {
                        continue;
                    }
                    for c0 in 0..C {
                        eh[ix(m, n, p, c0, f)] =
                            fast[wix(i as isize, j as isize, k as isize, c0, f)];
                    }
                }
            }
        }
    }
}
