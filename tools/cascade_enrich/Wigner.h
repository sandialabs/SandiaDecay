#ifndef CASCADE_ENRICH_WIGNER_H
#define CASCADE_ENRICH_WIGNER_H

// Exact Wigner 3j / 6j symbols via the Racah single-sum formulas.
//
// All angular momenta and projections are passed as TWICE their value (integers)
// so half-integer spins are exact: e.g. j = 3/2 is passed as 3, m = -1/2 as -1.

namespace cascade_enrich {

/// Wigner 3j symbol (j1 j2 j3 / m1 m2 m3), arguments are 2*j and 2*m.
double wigner_3j(int two_j1, int two_j2, int two_j3,
                 int two_m1, int two_m2, int two_m3);

/// Wigner 6j symbol { j1 j2 j3 / j4 j5 j6 }, arguments are 2*j.
double wigner_6j(int two_j1, int two_j2, int two_j3,
                 int two_j4, int two_j5, int two_j6);

} // namespace cascade_enrich

#endif
