#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define M 160
#define N 100
#define P 60
#define C 3
#define F 2
#define Q 100

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

// Helper function to separate bits for Morton encoding
static inline uint32_t part_bits(uint32_t n) {
  n &= 0x3ff; // Keep only lower 10 bits (sufficient for our dimensions)
  n = (n | (n << 16)) & 0x30000ff;
  n = (n | (n << 8)) & 0x300f00f;
  n = (n | (n << 4)) & 0x30c30c3;
  n = (n | (n << 2)) & 0x9249249;
  return n;
}

// Morton 3D encoding - interleaves bits of x, y, z coordinates
static inline uint32_t morton3D(int x, int y, int z) {
  return part_bits(x) | (part_bits(y) << 1) | (part_bits(z) << 2);
}

static inline size_t ix(int m, int n, int p, int c, int f) {
  // Apply modulo with proper handling of negative indices
  int mod_m = (m + M) % M;
  int mod_n = (n + N) % N;
  int mod_p = (p + P) % P;

  return f * M * N * P * C + C * morton3D(mod_m, mod_n, mod_p) + c;
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

    for (m = 0; m < M; ++m) {
      for (n = 0; n < N; ++n) {
        for (p = 0; p < P; ++p) {
          for (c0 = 0; c0 < 3; ++c0) {
            c1 = (c0 + 1) % 3;
            c2 = (c0 + 2) % 3;
            j0 = s * (c0 == 0);
            j1 = s * (c0 == 1);
            j2 = s * (c0 == 2);
            EH[ix(m - f, n - f, p - f, c0, f)] +=
                Sc * (EH[ix(m - f, n - f, p - f, c2, g)] -
                      EH[ix(m - j2 - f, n - j0 - f, p - j1 - f, c2, g)] -
                      EH[ix(m - f, n - f, p - f, c1, g)] +
                      EH[ix(m - j1 - f, n - j2 - f, p - j0 - f, c1, g)]);
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
