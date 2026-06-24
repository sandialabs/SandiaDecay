// Tests for the Wigner symbols and gamma-gamma angular-correlation coefficients.
// Returns non-zero on failure.

#include "AngularCorrelation.h"
#include "Wigner.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace cascade_enrich;

namespace {
int g_failures = 0;

void check_close(double got, double want, double tol, const std::string& msg) {
  if (std::fabs(got - want) > tol) {
    std::cerr << "FAIL: " << msg << " (got " << got << ", want " << want
              << " +/- " << tol << ")\n";
    ++g_failures;
  } else {
    std::cout << "ok:   " << msg << " = " << got << "\n";
  }
}
} // namespace

int main() {
  // --- Wigner 3j known values (arguments are 2*j, 2*m) ---
  // (1 1 0 / 0 0 0) = -1/sqrt(3)
  check_close(wigner_3j(2, 2, 0, 0, 0, 0), -1.0 / std::sqrt(3.0), 1e-12,
              "3j(1 1 0;0 0 0)");
  // (1 1 2 / 0 0 0) = sqrt(2/15)
  check_close(wigner_3j(2, 2, 4, 0, 0, 0), std::sqrt(2.0 / 15.0), 1e-12,
              "3j(1 1 2;0 0 0)");
  // (2 2 0 / 0 0 0) = 1/sqrt(5)
  check_close(wigner_3j(4, 4, 0, 0, 0, 0), 1.0 / std::sqrt(5.0), 1e-12,
              "3j(2 2 0;0 0 0)");

  // --- Wigner 6j known values ---
  // {1 1 1 / 1 1 1} = 1/6
  check_close(wigner_6j(2, 2, 2, 2, 2, 2), 1.0 / 6.0, 1e-12, "6j(1 1 1;1 1 1)");
  // {1 1 2 / 1 1 2} = -1/30  (standard table value)
  check_close(wigner_6j(2, 2, 4, 2, 2, 4), 1.0 / 30.0, 1e-12, "6j(1 1 2;1 1 2)");
  // {2 2 2 / 2 2 2} = -3/70
  check_close(wigner_6j(4, 4, 4, 4, 4, 4), -3.0 / 70.0, 1e-12, "6j(2 2 2;2 2 2)");

  // --- Co-60 angular correlation: 4(E2)2(E2)0 cascade ---
  // Literature: a2 = 0.1020, a4 = 0.0091.
  // J0=4 (twoJ=8), intermediate J1=2 (twoJ=4), J2=0; both pure E2 (L=2, d=0).
  {
    AngularCorrCoeffs c = angular_correlation(/*2J0*/ 8,
                                              /*L1*/ 2, /*L1p*/ 2, /*d1*/ 0.0,
                                              /*2J1*/ 4,
                                              /*L2*/ 2, /*L2p*/ 2, /*d2*/ 0.0,
                                              /*2J2*/ 0);
    check_close(c.a2, 0.1020, 0.001, "Co-60 a2 (4(E2)2(E2)0)");
    check_close(c.a4, 0.0091, 0.0005, "Co-60 a4 (4(E2)2(E2)0)");
  }

  // --- A 0(E1)1(E1)0 cascade: classic dipole-dipole, a2 = 0.5, a4 = 0. ---
  // Pure E1 (L=1), J0=0, intermediate J1=1, J2=0.
  {
    AngularCorrCoeffs c = angular_correlation(0, 1, 1, 0.0, 2, 1, 1, 0.0, 0);
    check_close(c.a2, 0.5, 1e-6, "0(E1)1(E1)0 a2 = 1/2");
    check_close(c.a4, 0.0, 1e-9, "0(E1)1(E1)0 a4 = 0");
  }

  // --- Spin-0 intermediate -> isotropic ---
  {
    AngularCorrCoeffs c = angular_correlation(2, 1, 1, 0.0, 0, 1, 1, 0.0, 2);
    check_close(c.a2, 0.0, 1e-12, "J_inter=0 -> a2=0");
    if (!c.valid) { std::cerr << "FAIL: J_inter=0 should be valid\n"; ++g_failures; }
  }

  if (g_failures == 0) {
    std::cout << "\nALL ANGULAR-CORRELATION TESTS PASSED\n";
    return 0;
  }
  std::cerr << "\n" << g_failures << " test(s) FAILED\n";
  return 1;
}
