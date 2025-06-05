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
const H: isize = 4; // tile height

const SX: isize = M / 2 + 1;
const SY: isize = N / 2 + 1;
const SZ: isize = P / 2 + 1;
const MV: isize = 2; // 1 mountain + 1 valley = 2

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
    ((((f * W + (m + W) % W) * W + (n + W) % W) * W + (p + W) % W) * C + c) as usize
}

fn main() -> Result<()> {
    let sc = 0.99 / (3.0_f32).sqrt();
    let mut eh = vec![0.0_f32; (M * N * P * C * F) as usize];
    eh[ix(SX, SY, SZ, 2, 0)] = 1.0;
    //let st = ramped_sin(0.3, 5.0, 3.0, Q as usize);

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
    let mut m0: isize = 0;
    let mut m1: isize;
    let mut n0: isize = 0;
    let mut n1: isize;
    let mut p0: isize = 0;
    let mut p1: isize;
    let mut fast = [0_f32; (W * W * W * 3 * 2) as usize];

    for _ in 0..(2 * Q / H) {
        //i *= H;
        for mvm in 0..MV {
            // mvm: 0: mountain, 1: valley (m-direction)
            for mut om in 0..oms {
                om *= W; // om: offset in x/m direction
                for mut on in 0..ons {
                    on *= W; // on: offset in y/n direction
                    for mut op in 0..ops {
                        op *= W; // op: offset in z/p direction
                        fill_fast(&mut fast, &eh, om, on, op);
                        if fast.iter().any(|&x| x > 0.0) {
                            println!("{:?}", &fast);
                        }
                        for h in 0..H {
                            f = h % 2;
                            g = 1 - f;
                            s = 1 - (2 * f);
                            m0 = if mvm == 0 {
                                0 + h / 2
                            } else {
                                W - 1 - (h + 1) / 2
                            };
                            m1 = if mvm == 0 {
                                W - 1 - (h + 1) / 2
                            } else {
                                W + h / 2
                            };
                            n0 = -(h + 1) / 2;
                            n1 = W - (h + 1) / 2;
                            p0 = -(h + 1) / 2;
                            p1 = W - (h + 1) / 2;
                            for m in m0..m1 {
                                for n in n0..n1 {
                                    for p in p0..p1 {
                                        for c0 in 0..C {
                                            c1 = (c0 + 1) % 3;
                                            c2 = (c0 + 2) % 3;
                                            j0 = s * ((c0 == 0) as isize);
                                            j1 = s * ((c0 == 1) as isize);
                                            j2 = s * ((c0 == 2) as isize);
                                            fast[wix(m, n, p, c0, f)] += sc
                                                * (fast[wix(m, n, p, c2, g)]
                                                    - fast[wix(m - j2, n - j0, p - j1, c2, g)]
                                                    - fast[wix(m, n, p, c1, g)]
                                                    + fast[wix(m - j1, n - j2, p - j0, c1, g)]);
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
    for f in 0..2 {
        for m in 0..W {
            for n in 0..W {
                for p in 0..W {
                    for c in 0..C {
                        fast[wix(m, n, p, c, f)] = eh[ix(om + m, on + n, op + p, c, f)];
                    }
                }
            }
        }
    }
}

fn extract_fast(eh: &mut [f32], fast: &[f32], om: isize, on: isize, op: isize) {
    for f in 0..2 {
        for m in 0..W {
            for n in 0..W {
                for p in 0..W {
                    for c in 0..C {
                        eh[ix(om + m, on + n, op + p, c, f)] = fast[wix(m, n, p, c, f)];
                    }
                }
            }
        }
    }
}
