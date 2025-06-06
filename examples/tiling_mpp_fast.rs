use anyhow::Result;
use std::fs::File;
use std::io::Write;

const M: isize = 160;
const N: isize = 96;
const P: isize = 64;
const C: isize = 3;
const F: isize = 2;
const W: isize = 8;
const H: isize = W - 2;
const Q: isize = 16 * H;
const W2: isize = W / 2;
const W32: isize = 3 * W2;

const SX: isize = M / 2;
const SY: isize = N / 2;
const SZ: isize = P / 2;
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
    ((((f * W + m) * W32 + n) * W32 + p) * C + c) as usize // no bounds wrapping needed
}

fn main() -> Result<()> {
    let sc = 0.99 / (3.0_f32).sqrt();
    let mut eh = vec![0.0_f32; (M * N * P * C * F) as usize];
    let st = ramped_sin(0.3, 5.0, 3.0, Q as usize);

    let oms = M / W; // offset idxs in x/m direction
    let ons = N / W + 1; // offset idxs in y/n direction (one extra bc no periodic boundaries)
    let ops = P / W + 1; // offset idxs in z/p direction (one extra bc no periodic boundaries)

    let mut c1: isize;
    let mut c2: isize;
    let mut f: isize;
    let mut g: isize;
    let mut s: isize;
    let mut j0: isize;
    let mut j1: isize;
    let mut j2: isize;

    // NOTE: when W == H, the fast array dimension size for mountain valley dimensions
    // should be of size W + 1. in that case also the wix() function should be adjusted
    // to use W+1 instead as well.
    let mut fast = [0_f32; (W * W32 * W32 * C * F) as usize];

    for mut i in 0..(2 * Q / H) {
        i *= H;
        for mvm in 0..MV {
            // mvm: 0: mountain, 1: valley (m-direction)
            for mut om in 0..oms {
                om *= W; // om: offset in x/m direction
                for mut on in 0..ons {
                    on *= W; // on: offset in y/n direction
                    for mut op in 0..ops {
                        op *= W; // op: offset in z/p direction
                        fill_fast(&mut fast, &eh, om, mvm, on, op);
                        for h in 0..H {
                            f = h % 2;
                            g = 1 - f;
                            s = (2 * f) - 1; // sign: -1 for E [f=0], +1 for H [f=1]

                            // global coords
                            let gm0 = om + (1 - mvm) * (h / 2 + 1) + mvm * (W - (h + 1) / 2);
                            let gm1 = om + (1 - mvm) * (W - (h + 1) / 2) + mvm * (W + 1 + h / 2);
                            let gn0 = on + W - (h + 1) / 2;
                            let gn1 = on + 2 * W - (h + 1) / 2;
                            let gp0 = op + W - (h + 1) / 2;
                            let gp1 = op + 2 * W - (h + 1) / 2;

                            // local coords
                            let m0 = (1 - mvm) * (h / 2 + 1) + mvm * (W2 - (h + 1) / 2);
                            let m1 = (1 - mvm) * (W - (h + 1) / 2) + mvm * (W2 + 1 + h / 2);
                            let n0 = W2 - (h + 1) / 2;
                            let n1 = W32 - (h + 1) / 2;
                            let p0 = W2 - (h + 1) / 2;
                            let p1 = W32 - (h + 1) / 2;

                            // update tile
                            for (m, gm) in (m0..m1).zip(gm0..gm1) {
                                for (n, gn) in (n0..n1).zip(gn0..gn1) {
                                    for (p, gp) in (p0..p1).zip(gp0..gp1) {
                                        for c0 in 0..C {
                                            c1 = (c0 + 1) % 3;
                                            c2 = (c0 + 2) % 3;
                                            j0 = s * ((c0 == 0) as isize); // positive for E [f=0], negative for H [f=1]
                                            j1 = s * ((c0 == 1) as isize); // positive for E [f=0], negative for H [f=1]
                                            j2 = s * ((c0 == 2) as isize); // positive for E [f=0], negative for H [f=1]
                                            eh[ix(m, n, p, c0, f)] += sc
                                                * (eh[ix(m, n, p, c2, g)]
                                                    - eh[ix(m + j2, n + j0, p + j1, c2, g)]
                                                    - eh[ix(m, n, p, c1, g)]
                                                    + eh[ix(m + j1, n + j2, p + j0, c1, g)]);
                                        }
                                        if f == 0 && gm == SX && gn == SY && gp == SZ {
                                            fast[wix(m, n, p, 2, 0)] += st[(i + h) as usize / 2];
                                        }
                                    }
                                }
                            }
                        }
                        extract_fast(&mut eh, &fast, om, mvm, on, op);
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

fn fill_fast(fast: &mut [f32], eh: &[f32], om: isize, mvm: isize, on: isize, op: isize) {
    fast.fill(0.0);
    let dwm = mvm * W2;
    for f in 0..F {
        for (i, m) in (om + dwm..om + W + dwm).enumerate() {
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

fn extract_fast(eh: &mut [f32], fast: &[f32], om: isize, mvm: isize, on: isize, op: isize) {
    let dwm = mvm * W2;
    for f in 0..F {
        for (i, m) in (om + dwm..om + W + dwm).enumerate() {
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
