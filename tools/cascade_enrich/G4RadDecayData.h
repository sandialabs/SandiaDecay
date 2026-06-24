#ifndef CASCADE_ENRICH_G4RADDECAYDATA_H
#define CASCADE_ENRICH_G4RADDECAYDATA_H

// Parser for GEANT4 RadioactiveDecay data (z{Z}.a{A} ASCII files), specifically
// the electron-capture shell fractions needed to attribute K/L vacancies (and
// hence x-ray emission) to each daughter level fed by EC.
//
// File format (RadioactiveDecay6.1.2 README_RDM):
//   P <excitation> <flag> <halflife>         parent state header
//   <Mode> <daughterEx> <intensity>          summary record (3 tokens; daughterEx=0)
//   <Mode> <daughterEx> <flag> <intensity> <Q>   per-level record (5 tokens)
//   EC modes: KshellEC / LshellEC / MshellEC / NshellEC. Per-level intensities
//   are percentages of the summary; we normalize the K/L/M/N split per level.

#include <string>
#include <vector>

namespace cascade_enrich {

struct ECShellFractions {
  double daughter_ex_keV = 0.0;
  double fK = 0.0, fL = 0.0, fM = 0.0, fN = 0.0;  ///< normalized shell fractions (sum 1)
  double total_intensity = 0.0;  ///< sum of raw K+L+M+N intensities for this level
};

struct G4ECData {
  int Z = 0, A = 0;
  std::vector<ECShellFractions> per_level;  ///< ground-state parent, one per fed daughter level
  ECShellFractions summary;                 ///< overall K/L/M/N split (fallback)
  bool has_summary = false;
  bool loaded = false;

  /// EC shell fractions for the daughter level nearest `ex_keV` (within `tol`);
  /// falls back to the summary split when no per-level record matches.
  const ECShellFractions* for_level(double ex_keV, double tol = 0.6) const;
};

/// Load RadioactiveDecay z{Z}.a{A} (the PARENT nuclide's file; EC fractions live
/// there keyed by daughter excitation). Only the ground-state parent decay is
/// captured. Missing file => loaded == false.
G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir);

} // namespace cascade_enrich

#endif
