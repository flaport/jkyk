#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// Compiler optimization hints
#pragma GCC optimize("O3,unroll-loops,inline-functions")
#pragma GCC target("native")

const long M = 160;
const long N = 96;
const long P = 64;
const long C = 3;
const long F = 2;
const long Q = 96;
const long W = 8;
const long H = 6;
const long W2 = W / 2;
const long W32 = 3 * W2;

const long SX = M / 2;
const long SY = N / 2;
const long SZ = P / 2;
const long MV = 2; // 1 mountain + 1 valley = 2

void ramped_sin(float omega, float width, float delay, int q, float *result)
{
    for (int i = 0; i < q; i++)
    {
        float t = i * omega;
        float value = ((1.0f + tanhf(t / width - delay)) / 2.0f) * sinf(t);
        result[i] = value;
    }
}

static inline size_t ix(long m, long n, long p, long c, long f)
{
    // Optimize modulo operations - these are expensive
    long m_wrapped = (m < 0) ? (m + M) : ((m >= M) ? (m - M) : m);
    long n_wrapped = (n < 0) ? (n + N) : ((n >= N) ? (n - N) : n);
    long p_wrapped = (p < 0) ? (p + P) : ((p >= P) ? (p - P) : p);

    return (size_t)((((f * M + m_wrapped) * N + n_wrapped) * P + p_wrapped) * C + c);
}

static inline size_t wix(long m, long n, long p, long c, long f)
{
    // Optimize modulo operations
    long m_wrapped = (m < 0) ? (m + W) : ((m >= W) ? (m - W) : m);
    long n_wrapped = (n < 0) ? (n + W32) : ((n >= W32) ? (n - W32) : n);
    long p_wrapped = (p < 0) ? (p + W32) : ((p >= W32) ? (p - W32) : p);

    return (size_t)((((f * W + m_wrapped) * W32 + n_wrapped) * W32 + p_wrapped) * C + c);
}

void fill_fast(float *restrict fast, float *restrict eh, long om, long on, long op, long mvm);
void extract_fast(float *restrict eh, float *restrict fast, long om, long on, long op, long mvm);

int main()
{
    float sc = 0.99f / sqrtf(3.0f);
    float *eh = calloc((size_t)(M * N * P * C * F), sizeof(float));
    float *st = malloc((size_t)Q * sizeof(float));
    ramped_sin(0.3f, 5.0f, 3.0f, (int)Q, st);

    long oms = M / W;     // offset idxs in x/m direction
    long ons = N / W + 1; // offset idxs in y/n direction (one extra bc no periodic boundaries)
    long ops = P / W + 1; // offset idxs in z/p direction (one extra bc no periodic boundaries)

    for (long i = 0; i < (2 * Q / H); i++)
    {
        long i_times_H = i * H;
        for (long mvm = 0; mvm < MV; mvm++)
        {
            // Pre-allocate arrays for better memory management
            const size_t fast_size = W * W32 * W32 * C * F;

#pragma omp parallel
            {
                // Thread-private allocation
                float *fast = malloc(fast_size * sizeof(float));

#pragma omp for
                for (long om = 0; om < oms; om++)
                {
                    long c1;
                    long c2;
                    long f;
                    long g;
                    long s;
                    long j0;
                    long j1;
                    long j2;
                    long om_times_W = om * W; // om: offset in x/m direction
                    for (long on = 0; on < ons; on++)
                    {
                        long on_times_W = on * W; // on: offset in y/n direction
                        for (long op = 0; op < ops; op++)
                        {
                            long op_times_W = op * W; // op: offset in z/p direction
                            fill_fast(fast, eh, om_times_W, on_times_W, op_times_W, mvm);
                            for (long h = 0; h < H; h++)
                            {
                                f = h % 2;
                                g = 1 - f;
                                s = (2 * f) - 1; // sign: -1 for E [f=0], +1 for H [f=1]

                                // global coords
                                long gm0 = om_times_W + (1 - mvm) * (h / 2 + 1) + mvm * (W - (h + 1) / 2);
                                long gm1 = om_times_W + (1 - mvm) * (W - (h + 1) / 2) + mvm * (W + 1 + h / 2);
                                long gn0 = on_times_W - (h + 1) / 2;
                                long gn1 = on_times_W + W - (h + 1) / 2;
                                long gp0 = op_times_W - (h + 1) / 2;
                                long gp1 = op_times_W + W - (h + 1) / 2;

                                // local coords
                                long m0 = (1 - mvm) * (h / 2 + 1) + mvm * (W2 - (h + 1) / 2);
                                long m1 = (1 - mvm) * (W - (h + 1) / 2) + mvm * (W2 + 1 + h / 2);
                                long n0 = W2 - (h + 1) / 2;
                                long n1 = W32 - (h + 1) / 2;
                                long p0 = W2 - (h + 1) / 2;
                                long p1 = W32 - (h + 1) / 2;

                                long gm = gm0;
                                for (long m = m0; m < m1; m++, gm++)
                                {
                                    long gn = gn0;
                                    for (long n = n0; n < n1; n++, gn++)
                                    {
                                        long gp = gp0;
                                        for (long p = p0; p < p1; p++, gp++)
                                        {
                                            // Pre-calculate common indices for better cache performance
                                            size_t base_idx = wix(m, n, p, 0, f);
                                            size_t base_idx_g = wix(m, n, p, 0, g);

                                            for (long c0 = 0; c0 < C; c0++)
                                            {
                                                c1 = (c0 + 1) % 3;
                                                c2 = (c0 + 2) % 3;
                                                j0 = s * (long)(c0 == 0); // negative for E [f=0], positive for H [f=1]
                                                j1 = s * (long)(c0 == 1); // negative for E [f=0], positive for H [f=1]
                                                j2 = s * (long)(c0 == 2); // negative for E [f=0], positive for H [f=1]

                                                // Use pre-calculated base indices and add component offset
                                                size_t idx_c0_f = base_idx + c0;
                                                size_t idx_c2_g = base_idx_g + c2;
                                                size_t idx_c1_g = base_idx_g + c1;

                                                fast[idx_c0_f] += sc * (fast[idx_c2_g] -
                                                                        fast[wix(m + j2, n + j0, p + j1, c2, g)] -
                                                                        fast[idx_c1_g] +
                                                                        fast[wix(m + j1, n + j2, p + j0, c1, g)]);
                                            }
                                            if (f == 0 && gm == SX && gn == SY && gp == SZ)
                                            {
                                                fast[wix(m, n, p, 2, 0)] += st[(i_times_H + h) / 2];
                                            }
                                        }
                                    }
                                }
                            }
                            extract_fast(eh, fast, om_times_W, on_times_W, op_times_W, mvm);
                        }
                    }
                }

                // Free thread-private memory
                free(fast);
            }
        }
    }

    // Write output to binary file
    FILE *file = fopen("output.bin", "wb");
    if (!file)
    {
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

void fill_fast(float *restrict fast, float *restrict eh, long om, long on, long op, long mvm)
{
    const size_t fast_size = (size_t)(W * W32 * W32 * C * F);
    memset(fast, 0, fast_size * sizeof(float));

    const long dwm = mvm * W2;
    const long om_dwm = om + dwm;
    const long on_w2 = on - W2;
    const long op_w2 = op - W2;

    for (long f = 0; f < F; f++)
    {
        for (long i = 0, m = om_dwm; i < W; i++, m++)
        {
            for (long j = 0, n = on_w2; j < W32; j++, n++)
            {
                if (n < 0 || n >= N)
                    continue;

                for (long k = 0, p = op_w2; k < W32; k++, p++)
                {
                    if (p < 0 || p >= P)
                        continue;

                    // Batch copy all components at once for better cache utilization
                    const size_t fast_base = wix(i, j, k, 0, f);
                    const size_t eh_base = ix(m, n, p, 0, f);

                    fast[fast_base] = eh[eh_base];         // c0 = 0
                    fast[fast_base + 1] = eh[eh_base + 1]; // c0 = 1
                    fast[fast_base + 2] = eh[eh_base + 2]; // c0 = 2
                }
            }
        }
    }
}

void extract_fast(float *restrict eh, float *restrict fast, long om, long on, long op, long mvm)
{
    const long dwm = mvm * W2;
    const long om_dwm = om + dwm;
    const long on_w2 = on - W2;
    const long op_w2 = op - W2;

    for (long f = 0; f < F; f++)
    {
        for (long i = 0, m = om_dwm; i < W; i++, m++)
        {
            for (long j = 0, n = on_w2; j < W32; j++, n++)
            {
                if (n < 0 || n >= N)
                    continue;

                for (long k = 0, p = op_w2; k < W32; k++, p++)
                {
                    if (p < 0 || p >= P)
                        continue;

                    // Batch copy all components at once for better cache utilization
                    const size_t fast_base = wix(i, j, k, 0, f);
                    const size_t eh_base = ix(m, n, p, 0, f);

                    eh[eh_base] = fast[fast_base];         // c0 = 0
                    eh[eh_base + 1] = fast[fast_base + 1]; // c0 = 1
                    eh[eh_base + 2] = fast[fast_base + 2]; // c0 = 2
                }
            }
        }
    }
}
