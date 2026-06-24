#include "Wigner.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace cascade_enrich {

namespace {

// Exact factorial cache (n! as double; exact for n up to ~22, ample here).
double fact(int n) {
  static std::vector<double> cache{1.0};
  if (n < 0)
    return 0.0;
  while (static_cast<int>(cache.size()) <= n)
    cache.push_back(cache.back() * static_cast<double>(cache.size()));
  return cache[static_cast<std::size_t>(n)];
}

// Triangle condition on 2*j arguments: |j1-j2| <= j3 <= j1+j2 and j1+j2+j3
// integer (i.e. two_j1+two_j2+two_j3 even).
bool triangle(int a, int b, int c) {
  if ((a + b + c) % 2 != 0)
    return false;
  if (c < std::abs(a - b) || c > a + b)
    return false;
  return true;
}

// Triangle coefficient Delta(abc) = sqrt[(a+b-c)!(a-b+c)!(-a+b+c)!/(a+b+c+1)!],
// arguments are 2*j.
double delta_coeff(int a, int b, int c) {
  const double num = fact((a + b - c) / 2) * fact((a - b + c) / 2) *
                     fact((-a + b + c) / 2);
  const double den = fact((a + b + c) / 2 + 1);
  return std::sqrt(num / den);
}

int phase(int n) { return (n % 2 == 0) ? 1 : -1; }

} // namespace

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0)
    return 0.0;
  if (!triangle(j1, j2, j3))
    return 0.0;
  if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3)
    return 0.0;
  if (((j1 + m1) & 1) || ((j2 + m2) & 1) || ((j3 + m3) & 1))
    return 0.0;  // ji+mi must be integer

  const double pref =
      phase((j1 - j2 - m3) / 2) * delta_coeff(j1, j2, j3) *
      std::sqrt(fact((j1 + m1) / 2) * fact((j1 - m1) / 2) *
                fact((j2 + m2) / 2) * fact((j2 - m2) / 2) *
                fact((j3 + m3) / 2) * fact((j3 - m3) / 2));

  // Sum over integer k where every factorial argument is non-negative.
  const int kmin = std::max(0, std::max((j2 - j3 - m1) / 2, (j1 - j3 + m2) / 2));
  const int kmax = std::min((j1 + j2 - j3) / 2,
                            std::min((j1 - m1) / 2, (j2 + m2) / 2));

  double sum = 0.0;
  for (int k = kmin; k <= kmax; ++k) {
    const double den = fact(k) * fact((j1 + j2 - j3) / 2 - k) *
                       fact((j1 - m1) / 2 - k) * fact((j2 + m2) / 2 - k) *
                       fact((j3 - j2 + m1) / 2 + k) *
                       fact((j3 - j1 - m2) / 2 + k);
    sum += phase(k) / den;
  }
  return pref * sum;
}

double wigner_6j(int a, int b, int c, int d, int e, int f) {
  // Four triads must each satisfy the triangle condition.
  if (!triangle(a, b, c) || !triangle(a, e, f) ||
      !triangle(d, b, f) || !triangle(d, e, c))
    return 0.0;

  const double pre = delta_coeff(a, b, c) * delta_coeff(a, e, f) *
                     delta_coeff(d, b, f) * delta_coeff(d, e, c);

  // Sums of the four triads and the three "square" sums (all in 2*j units / 2).
  const int s1 = (a + b + c) / 2;
  const int s2 = (a + e + f) / 2;
  const int s3 = (d + b + f) / 2;
  const int s4 = (d + e + c) / 2;
  const int t1 = (a + b + d + e) / 2;
  const int t2 = (b + c + e + f) / 2;
  const int t3 = (a + c + d + f) / 2;

  const int kmin = std::max(std::max(s1, s2), std::max(s3, s4));
  const int kmax = std::min(t1, std::min(t2, t3));

  double sum = 0.0;
  for (int k = kmin; k <= kmax; ++k) {
    const double den = fact(k - s1) * fact(k - s2) * fact(k - s3) *
                       fact(k - s4) * fact(t1 - k) * fact(t2 - k) *
                       fact(t3 - k);
    sum += phase(k) * fact(k + 1) / den;
  }
  return pre * sum;
}

} // namespace cascade_enrich
