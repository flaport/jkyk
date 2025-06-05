use anyhow::Result;
use std::fs::File;
use std::io::Write;

const M: usize = 160;
const N: usize = 100;
const P: usize = 60;
const C: usize = 3;
const F: usize = 2;
const Q: usize = 100;

const Sx: usize = M / 2;
const Sy: usize = N / 2;
const Sz: usize = P / 2;

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
fn ix(m: usize, dm: isize, n: usize, dn: isize, p: usize, dp: isize, c: usize, f: usize) -> usize {
    let mm = ((m as isize + dm + M as isize) % (M as isize)) as usize;
    let nn = ((n as isize + dn + N as isize) % (N as isize)) as usize;
    let pp = ((p as isize + dp + P as isize) % (P as isize)) as usize;
    return (((f * M + mm) * N + nn) * P + pp) * C + c;
}

fn main() -> Result<()> {
    let Sc = 0.99 / (3.0_f32).sqrt();
    let mut EH = vec![0.0_f32; M * N * P * C * F];
    let St = ramped_sin(0.3, 5.0, 3.0, Q);

    let mut c1: usize;
    let mut c2: usize;
    let mut f: isize;
    let mut g: isize;
    let mut s: isize;
    let mut j0: isize;
    let mut j1: isize;
    let mut j2: isize;

    for i in 0..(2 * Q) {
        f = (i % 2) as isize;
        g = 1 - f;
        s = 1 - (2 * f as isize);

        for m in 0..M {
            for n in 0..N {
                for p in 0..P {
                    for c0 in 0..C {
                        c1 = (c0 + 1) % 3;
                        c2 = (c0 + 2) % 3;
                        j0 = s * ((c0 == 0) as isize);
                        j1 = s * ((c0 == 1) as isize);
                        j2 = s * ((c0 == 2) as isize);
                        EH[ix(m, -f, n, -f, p, -f, c0, f as usize)] += Sc
                            * (EH[ix(m, -f, n, -f, p, -f, c2, g as usize)]
                                - EH[ix(m, -j2 - f, n, -j0 - f, p, -j1 - f, c2, g as usize)]
                                - EH[ix(m, -f, n, -f, p, -f, c1, g as usize)]
                                + EH[ix(m, -j1 - f, n, -j2 - f, p, -j0 - f, c1, g as usize)]);
                    }
                }
            }
        }
        if f == 0 {
            EH[ix(Sx, 0, Sy, 0, Sz, 0, 2, 0)] += St[i / 2];
        }
    }

    // Write output to binary file
    let mut file = File::create("output.bin")?;

    // Convert f32 slice to bytes and write
    let byte_data: Vec<u8> = EH[..2 * M * N * P * 3]
        .iter()
        .flat_map(|&f| f.to_le_bytes())
        .collect();

    file.write_all(&byte_data)?;

    Ok(())
}
