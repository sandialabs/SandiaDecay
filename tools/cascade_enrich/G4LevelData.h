#ifndef CASCADE_ENRICH_G4LEVELDATA_H
#define CASCADE_ENRICH_G4LEVELDATA_H

// Parser for GEANT4 PhotonEvaporation level data (z{Z}.a{A} ASCII files).
//
// These files describe the nuclear de-excitation scheme of a daughter nucleus:
// a series of levels, each followed by its gamma transitions to lower levels.
// We extract, per transition, the multipole orders + mixing ratio (for
// gamma-gamma angular correlations) and the total + subshell internal-conversion
// coefficients (for rigorous x-ray cascade coincidences), and per level the spin
// (2J) needed by the F-coefficient formalism.
//
// File format (PhotonEvaporation6.1.2 README):
//   Level line   : index floatFlag Eexc(keV) halflife(s) JPi nGammas
//   Transition   : daughterIdx Egamma(keV) relIntensity multipoleCode mixing
//                  alphaTotal [10 normalized subshell ICC fractions, iff alpha>0]
//   JPi == 99    => unknown; sign of JPi gives parity.
//   multipoleCode: single digit d in 1..9 => E0,E1,M1,E2,M2,E3,M3,E4,M4 (order
//                  L = d/2); 100*Nx+Ny => mixed, lower multipole Nx + higher Ny.

#include <string>
#include <vector>

namespace cascade_enrich {

struct G4Level {
  int index = 0;
  double excitation_keV = 0.0;
  double halflife_s = -1.0;   ///< -1 => stable ground state
  int twoJ = -1;              ///< 2*J; -1 if unknown (JPi == 99)
  int parity = 0;             ///< +1 / -1; 0 if unknown
  bool j_known = false;
};

struct G4Transition {
  int from_level = -1;        ///< emitting (upper) level index
  int to_level = -1;          ///< daughter (lower) level index
  double gamma_keV = 0.0;
  double rel_intensity = 0.0; ///< relative gamma emission intensity within the level
  int multipole_code = 0;     ///< raw PhotonEvaporation code (0 => unknown)
  int L0 = 0;                 ///< lower multipole angular-momentum order
  int Lp = 0;                 ///< higher multipole order (== L0 when pure)
  double mixing_delta = 0.0;  ///< multipole mixing ratio delta
  double alpha_total = 0.0;   ///< total internal conversion coefficient I_ce/I_gamma
  double icc_shell[10] = {0,0,0,0,0,0,0,0,0,0}; ///< normalized fractions: K,L1-3,M1-5,outer
  bool has_icc = false;

  /// Absolute K-shell ICC = alpha_total * f_K.
  double alpha_K() const { return alpha_total * icc_shell[0]; }
  /// Absolute total L-shell ICC = alpha_total * (f_L1+f_L2+f_L3).
  double alpha_L() const {
    return alpha_total * (icc_shell[1] + icc_shell[2] + icc_shell[3]);
  }
  /// Absolute per-subshell L ICC. sub = 1,2,3 for L1,L2,L3.
  double alpha_Lsub(int sub) const { return alpha_total * icc_shell[sub]; }
};

struct G4LevelScheme {
  int Z = 0, A = 0;
  std::vector<G4Level> levels;
  std::vector<G4Transition> transitions;
  bool loaded = false;

  const G4Level* level(int index) const;
};

/// Decode a PhotonEvaporation multipolarity code into lower/higher angular
/// momentum orders. Returns false (L0=Lp=0) for an unknown code (0).
bool decode_multipolarity(int code, int& L0, int& Lp);

/// Load PhotonEvaporation z{Z}.a{A}. `data_dir` is the PhotonEvaporation6.1.2
/// directory (e.g. $G4LEVELGAMMADATA). On a missing/unreadable file returns a
/// scheme with loaded == false.
G4LevelScheme load_photon_evaporation(int Z, int A, const std::string& data_dir);

} // namespace cascade_enrich

#endif
