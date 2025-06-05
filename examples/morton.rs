use anyhow::Result;
use std::fs::File;
use std::io::Write;

const M: isize = 128;
const N: isize = 128;
const P: isize = 128;
const C: isize = 3;
const F: isize = 2;
const Q: isize = 100;

const SX: isize = M / 2;
const SY: isize = N / 2;
const SZ: isize = P / 2;

fn main() -> Result<()> {
    let sc = 0.99 / (3.0_f32).sqrt();
    let mut eh = vec![0.0_f32; (M * N * P * C * F) as usize];
    let st = ramped_sin(0.3, 5.0, 3.0, Q as usize);

    #[inline]
    fn ix(m: isize, n: isize, p: isize, c: isize, f: isize) -> usize {
        let mm = (m + M) % M;
        let nn = (n + N) % N;
        let pp = (p + P) % P;
        //((((f * M + mm) * N + nn) * P + pp) * C + c) as usize
        (f * M * N * P * C) as usize + morton_3d(mm, nn, pp) * (C as usize) + (c as usize)
    }

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
                        eh[ix(m - f, n - f, p - f, c0, f)] += sc
                            * (eh[ix(m - f, n - f, p - f, c2, g)]
                                - eh[ix(m - j2 - f, n - j0 - f, p - j1 - f, c2, g)]
                                - eh[ix(m - f, n - f, p - f, c1, g)]
                                + eh[ix(m - j1 - f, n - j2 - f, p - j0 - f, c1, g)]);
                    }
                }
            }
        }
        if f == 0 {
            eh[ix(SX, SY, SZ, 2, 0)] += st[i as usize / 2];
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

fn ramped_sin(omega: f32, width: f32, delay: f32, q: usize) -> Vec<f32> {
    let mut result = Vec::with_capacity(q);
    for i in 0..q {
        let t = i as f32 * omega;
        let value = ((1.0 + (t / width - delay).tanh()) / 2.0) * t.sin();
        result.push(value);
    }
    result
}

fn part_bits(n: usize) -> usize {
    let mut n = n & 0x3ff; // Keep only lower 10 bits (sufficient for our dimensions)
    n = (n | (n << 16)) & 0x30000ff;
    n = (n | (n << 8)) & 0x300f00f;
    n = (n | (n << 4)) & 0x30c30c3;
    n = (n | (n << 2)) & 0x9249249;
    n
}

fn morton_3d(x: isize, y: isize, z: isize) -> usize {
    let x = part_bits(x as usize);
    let y = part_bits(y as usize);
    let z = part_bits(z as usize);
    z | (y << 1) | (x << 2)
}
