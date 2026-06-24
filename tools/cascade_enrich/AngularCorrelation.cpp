#include "AngularCorrelation.h"

#include "Wigner.h"

#include <cmath>

namespace cascade_enrich {

double F_coefficient(int k, int L, int Lp, int two_j_far, int two_j_inter) {
  if (k == 0)
    return 1.0;  // F_0 normalization
  // F_k(L L' j_far j) = (-1)^(j_far + j - 1)
  //    * sqrt[(2L+1)(2L'+1)(2k+1)(2j+1)]
  //    * (L L' k / 1 -1 0)             [3j]
  //    * { L L' k / j j j_far }        [6j]
  const double threej = wigner_3j(2 * L, 2 * Lp, 2 * k, 2, -2, 0);
  if (threej == 0.0)
    return 0.0;
  const double sixj =
      wigner_6j(2 * L, 2 * Lp, 2 * k, two_j_inter, two_j_inter, two_j_far);
  if (sixj == 0.0)
    return 0.0;

  // Phase (-1)^(j_far + j_inter - 1); (two_j_far+two_j_inter)/2 is integer for a
  // valid gamma cascade (their difference is bounded by an integer multipole).
  const int phase_exp = (two_j_far + two_j_inter) / 2 - 1;
  const double phase = (phase_exp % 2 == 0) ? 1.0 : -1.0;

  const double norm = std::sqrt(static_cast<double>((2 * L + 1) * (2 * Lp + 1) *
                                                    (2 * k + 1) *
                                                    (two_j_inter + 1)));
  return phase * norm * threej * sixj;
}

double A_coefficient(int k, int L, int Lp, double delta,
                     int two_j_far, int two_j_inter) {
  if (k == 0)
    return 1.0;
  const double f_LL = F_coefficient(k, L, L, two_j_far, two_j_inter);
  const double f_LLp = F_coefficient(k, L, Lp, two_j_far, two_j_inter);
  const double f_LpLp = F_coefficient(k, Lp, Lp, two_j_far, two_j_inter);
  const double num = f_LL + 2.0 * delta * f_LLp + delta * delta * f_LpLp;
  return num / (1.0 + delta * delta);
}

AngularCorrCoeffs angular_correlation(int two_j0,
                                      int L1, int L1p, double delta1,
                                      int two_j1_inter,
                                      int L2, int L2p, double delta2,
                                      int two_j2) {
  AngularCorrCoeffs c;
  if (two_j0 < 0 || two_j1_inter < 0 || two_j2 < 0)
    return c;  // unknown spin -> isotropic
  if (L1 <= 0 || L1p <= 0 || L2 <= 0 || L2p <= 0)
    return c;
  // A spin-0 intermediate cannot carry an alignment: correlation is isotropic.
  if (two_j1_inter == 0) {
    c.valid = true;
    return c;
  }

  c.a2 = A_coefficient(2, L1, L1p, delta1, two_j0, two_j1_inter) *
         A_coefficient(2, L2, L2p, delta2, two_j2, two_j1_inter);
  c.a4 = A_coefficient(4, L1, L1p, delta1, two_j0, two_j1_inter) *
         A_coefficient(4, L2, L2p, delta2, two_j2, two_j1_inter);
  c.valid = true;
  return c;
}

} // namespace cascade_enrich
