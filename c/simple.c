#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 160
#define N 100
#define P 60
#define Q 100

const int Sx = M / 2;
const int Sy = N / 2;
const int Sz = P / 2;

void ramped_sin(float omega, float width, float delay, size_t q, float *result) {
    float t;
    for (size_t i = 0; i < q; ++i) {
        t = i * omega;
        result[i] = ((1.0 + tanh(t / width - delay)) / 2.0) * sin(t);
    }
}


int main() {
    const float Sc = 0.99 / sqrt(3.0);
    float St[Q];
    ramped_sin(0.3, 5.0, 3.0, Q, St);
    float *EH = (float *)calloc(2 * M * N * P * 3, sizeof(float));

    size_t q, r, t, c1, c2, m_q, n_q, p_q, m_j2_q, n_j0_q, p_j1_q, m_j1_q, n_j2_q, p_j0_q;
    int s, j0, j1, j2;

    for (size_t i = 0; i < 2 * Q; ++i) {
        q = i % 2;
        s = 1 - (2 * q);
        r = (1 - q) * M * N * P * 3;
        t = q * M * N * P * 3;

        for (size_t m = 0; m < M; ++m) {
            for (size_t n = 0; n < N; ++n) {
                for (size_t p = 0; p < P; ++p) {
                    for (size_t c0 = 0; c0 < 3; ++c0) {
                        c1 = (c0 + 1) % 3;
                        c2 = (c0 + 2) % 3;
                        j0 = s * (c0 == 0);
                        j1 = s * (c0 == 1);
                        j2 = s * (c0 == 2);
                        m_q = ((M + m - q) % M) * N * P * 3;
                        n_q = ((N + n - q) % N) * P * 3;
                        p_q = ((P + p - q) % P) * 3;
                        m_j2_q = ((M + m - j2 - q) % M) * N * P * 3;
                        n_j0_q = ((N + n - j0 - q) % N) * P * 3;
                        p_j1_q = ((P + p - j1 - q) % P) * 3;
                        m_j1_q = ((M + m - j1 - q) % M) * N * P * 3;
                        n_j2_q = ((N + n - j2 - q) % N) * P * 3;
                        p_j0_q = ((P + p - j0 - q) % P) * 3;
                        EH[t + m_q + n_q + p_q + c0] += Sc
                            * (EH[r + m_q + n_q + p_q + c2]
                                - EH[r + m_j2_q + n_j0_q + p_j1_q + c2]
                                - EH[r + m_q + n_q + p_q + c1]
                                + EH[r + m_j1_q + n_j2_q + p_j0_q + c1]);
                    }
                }
            }
        }
        if (q == 0) {
            EH[t + (size_t)(Sx * N * P * 3) + (size_t)(Sy * P * 3) + (size_t)(Sz * 3) + 2] += St[i / 2];
        }
    }

    FILE *file = fopen("output.bin", "wb");
    fwrite(EH, sizeof(float), 2 * M * N * P * 3, file);
    fclose(file);

    free(EH);

    return 0;
}

