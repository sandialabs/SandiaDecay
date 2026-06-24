#ifndef CASCADE_ENRICH_ANGULARCORRELATION_H
#define CASCADE_ENRICH_ANGULARCORRELATION_H

// Gamma-gamma directional correlation coefficients a2/a4 for a two-step cascade
//
//     J0  --(gamma1: L1/L1', mixing d1)-->  J1  --(gamma2: L2/L2', d2)-->  J2
//
// through a shared intermediate level J1, computed from the standard
// F-coefficient (Racah) formalism. The correlation function is
//
//     W(theta) = 1 + a2*P2(cos theta) + a4*P4(cos theta),
//
// with a_k = A_k(gamma1) * A_k(gamma2), and per-transition
//
//     A_k = [F_k(L,L) + 2 d F_k(L,L') + d^2 F_k(L',L')] / (1 + d^2),
//
// where each F_k references the transition's "far" level and the shared
// intermediate level. Spins are passed as 2*J (integers).

namespace cascade_enrich {

/// F_k(L, Lp; j_far, j_inter) coefficient. Angular momenta L,Lp,k are integers;
/// the two spins are passed as 2*J. j_inter is the shared intermediate level.
double F_coefficient(int k, int L, int Lp, int two_j_far, int two_j_inter);

/// Per-transition orientation coefficient A_k (mixing-weighted), referencing the
/// transition's far level and the shared intermediate level.
double A_coefficient(int k, int L, int Lp, double mixing_delta,
                     int two_j_far, int two_j_inter);

struct AngularCorrCoeffs {
  double a2 = 0.0;
  double a4 = 0.0;
  bool valid = false;  ///< false if any spin is unknown / inputs degenerate
};

/// a2/a4 for the J0->J1->J2 cascade. Pass 2*J for each spin; (L,Lp,delta) per
/// transition. Returns valid=false when a spin is unknown (negative) or the
/// geometry is degenerate (e.g. J0 or J2 == 0 gives a known nonzero result; a
/// spin-0 intermediate gives isotropic a2=a4=0, valid=true).
AngularCorrCoeffs angular_correlation(int two_j0,
                                      int L1, int L1p, double delta1,
                                      int two_j1_inter,
                                      int L2, int L2p, double delta2,
                                      int two_j2);

} // namespace cascade_enrich

#endif
