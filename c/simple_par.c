#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <omp.h>

#define M 160
#define N 96
#define P 64
#define C 3
#define F 2
#define Q 96

const int Sx = M / 2;
const int Sy = N / 2;
const int Sz = P / 2;

void ramped_sin(float omega, float width, float delay, size_t q,
                float *result) {
  float t;
  for (size_t i = 0; i < q; ++i) {
    t = i * omega;
    result[i] = ((1.0 + tanh(t / width - delay)) / 2.0) * sin(t);
  }
}

static inline size_t ix(int m, int n, int p, int c, int f) {
  return (((f * M + (m + M) % M) * N + (n + N) % N) * P + (p + P) % P) * C + c;
}

int main() {
  const float Sc = 0.99 / sqrt(3.0);

  float *EH = (float *)calloc(M * N * P * C * F, sizeof(float));

  float St[Q];
  ramped_sin(0.3, 5.0, 3.0, Q, St);

  size_t f, g;
  int s, m, n, p, c0, c1, c2, j0, j1, j2;

  for (size_t i = 0; i < 2 * Q; ++i) {
    f = i % 2;
    g = 1 - f;
    s = 1 - (2 * f);

    #pragma omp parallel for
    for (m = 0; m < M; ++m) {
      for (n = 0; n < N; ++n) {
        for (p = 0; p < P; ++p) {
          for (c0 = 0; c0 < 3; ++c0) {
            c1 = (c0 + 1) % 3;
            c2 = (c0 + 2) % 3;
            j0 = s * (c0 == 0);
            j1 = s * (c0 == 1);
            j2 = s * (c0 == 2);
            EH[ix(m, n, p, c0, f)] +=
                Sc * (EH[ix(m, n, p, c2, g)] -
                      EH[ix(m - j2, n - j0, p - j1, c2, g)] -
                      EH[ix(m, n, p, c1, g)] +
                      EH[ix(m - j1, n - j2, p - j0, c1, g)]);
          }
        }
      }
    }
    if (f == 0) {
      EH[ix(Sx, Sy, Sz, 2, 0)] += St[i / 2];
    }
  }

  FILE *file = fopen("output.bin", "wb");
  fwrite(EH, sizeof(float), 2 * M * N * P * 3, file);
  fclose(file);

  free(EH);

  return 0;
}
