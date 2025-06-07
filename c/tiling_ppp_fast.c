#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const long M = 160;
const long N = 96;
const long P = 64;
const long C = 3;
const long F = 2;
const long Q = 96;
const long W = 8; // tile width
const long H = 8; // tile height
const long W2 = W / 2;
const long W32 = 3 * W / 2;

const long SX = M / 2;
const long SY = N / 2;
const long SZ = P / 2;

void ramped_sin(float omega, float width, float delay, int q, float* result) {
    for (int i = 0; i < q; i++) {
        float t = i * omega;
        float value = ((1.0f + tanhf(t / width - delay)) / 2.0f) * sinf(t);
        result[i] = value;
    }
}

static inline size_t ix(long m, long n, long p, long c, long f) {
    return (size_t)((((f * M + (m + M) % M) * N + (n + N) % N) * P + (p + P) % P) * C + c);
}

static inline size_t wix(long m, long n, long p, long c, long f) {
    return (size_t)((((f * W32 + (m + W32) % W32) * W32 + (n + W32) % W32) * W32 + (p + W32) % W32) * C + c);
}

void fill_fast(float* fast, float* eh, long om, long on, long op) {
    memset(fast, 0, (size_t)(3 * W / 2 * 3 * W / 2 * 3 * W / 2 * C * F) * sizeof(float));
    for (long f = 0; f < F; f++) {
        long i = 0;
        for (long m = om - W2; m < om + W; m++) {
            if (m < 0 || m >= M) {
                i++;
                continue;
            }
            long j = 0;
            for (long n = on - W2; n < on + W; n++) {
                if (n < 0 || n >= N) {
                    j++;
                    continue;
                }
                long k = 0;
                for (long p = op - W2; p < op + W; p++) {
                    if (p < 0 || p >= P) {
                        k++;
                        continue;
                    }
                    for (long c0 = 0; c0 < C; c0++) {
                        fast[wix(i, j, k, c0, f)] = eh[ix(m, n, p, c0, f)];
                    }
                    k++;
                }
                j++;
            }
            i++;
        }
    }
}

void extract_fast(float* eh, float* fast, long om, long on, long op) {
    for (long f = 0; f < F; f++) {
        long i = 0;
        for (long m = om - W2; m < om + W; m++) {
            if (m < 0 || m >= M) {
                i++;
                continue;
            }
            long j = 0;
            for (long n = on - W2; n < on + W; n++) {
                if (n < 0 || n >= N) {
                    j++;
                    continue;
                }
                long k = 0;
                for (long p = op - W2; p < op + W; p++) {
                    if (p < 0 || p >= P) {
                        k++;
                        continue;
                    }
                    for (long c0 = 0; c0 < C; c0++) {
                        eh[ix(m, n, p, c0, f)] = fast[wix(i, j, k, c0, f)];
                    }
                    k++;
                }
                j++;
            }
            i++;
        }
    }
}

void run(float* eh, float* st) {
    float sc = 0.99f / sqrtf(3.0f);
    long oms = M / W; // offset idxs in x/m direction
    long ons = N / W; // offset idxs in y/n direction
    long ops = P / W; // offset idxs in z/p direction

    long c1;
    long c2;
    long f;
    long g;
    long s;
    long j0;
    long j1;
    long j2;
    float fast[3 * W / 2 * 3 * W / 2 * 3 * W / 2 * C * F];

    for (long i = 0; i < (2 * Q / H); i++) {
        long i_times_H = i * H;
        for (long om = 0; om < oms; om++) {
            long om_times_W = om * W; // om: offset in x/m direction
            for (long on = 0; on < ons; on++) {
                long on_times_W = on * W; // on: offset in y/n direction
                for (long op = 0; op < ops; op++) {
                    long op_times_W = op * W; // op: offset in z/p direction
                    fill_fast(fast, eh, om_times_W, on_times_W, op_times_W);
                    for (long h = 0; h < H; h++) {
                        f = h % 2;
                        g = 1 - f;
                        s = (2 * f) - 1; // sign: -1 for E [f=0], +1 for H [f=1]
                        long m0 = W2 - (h + 1) / 2;
                        long m1 = W32 - (h + 1) / 2;
                        long n0 = W2 - (h + 1) / 2;
                        long n1 = W32 - (h + 1) / 2;
                        long p0 = W2 - (h + 1) / 2;
                        long p1 = W32 - (h + 1) / 2;
                        for (long m = m0; m < m1; m++) {
                            for (long n = n0; n < n1; n++) {
                                for (long p = p0; p < p1; p++) {
                                    for (long c0 = 0; c0 < C; c0++) {
                                        c1 = (c0 + 1) % 3;
                                        c2 = (c0 + 2) % 3;
                                        j0 = s * (long)(c0 == 0); // negative for E [f=0], positive for H [f=1]
                                        j1 = s * (long)(c0 == 1); // negative for E [f=0], positive for H [f=1]
                                        j2 = s * (long)(c0 == 2); // negative for E [f=0], positive for H [f=1]
                                        fast[wix(m, n, p, c0, f)] += sc
                                            * (fast[wix(m, n, p, c2, g)]
                                                - fast[wix(m + j2, n + j0, p + j1, c2, g)]
                                                - fast[wix(m, n, p, c1, g)]
                                                + fast[wix(m + j1, n + j2, p + j0, c1, g)]);
                                    }
                                    if (f == 0
                                        && m == SX - om_times_W + W2
                                        && n == SY - on_times_W + W2
                                        && p == SZ - op_times_W + W2) {
                                        fast[wix(
                                            SX - om_times_W + W2,
                                            SY - on_times_W + W2,
                                            SZ - op_times_W + W2,
                                            2,
                                            0)] += st[(i_times_H + h) / 2];
                                    }
                                }
                            }
                        }
                    }
                    extract_fast(eh, fast, om_times_W, on_times_W, op_times_W);
                }
            }
        }
    }
}

int main() {
    float* eh = calloc((size_t)(M * N * P * C * F), sizeof(float));
    float* st = calloc((size_t)Q , sizeof(float));
    ramped_sin(0.3f, 5.0f, 3.0f, (int)Q, st);

    run(eh, st);

    // Write output to binary file
    FILE* file = fopen("output.bin", "wb");
    if (!file) {
        fprintf(stderr, "Error opening output file\n");
        free(eh);
        free(st);
        return 1;
    }

    // Write f32 data as bytes
    size_t elements_to_write = (size_t)(2 * M * N * P * 3);
    fwrite(eh, sizeof(float), elements_to_write, file);

    fclose(file);
    free(eh);
    free(st);

    return 0;
}

