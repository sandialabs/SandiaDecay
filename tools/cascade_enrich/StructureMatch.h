#ifndef CASCADE_ENRICH_STRUCTUREMATCH_H
#define CASCADE_ENRICH_STRUCTUREMATCH_H

// Structure-aware matching of SandiaDecay gammas to GEANT4 level-scheme
// transitions.
//
// Energy-only matching is ambiguous in dense schemes (many transitions within a
// fraction of a keV). This matcher additionally uses the SandiaDecay coincidence
// graph: two coincident gammas must lie in the same de-excitation cascade, so
// their G4 transitions should be *adjacent* (share a level). We collect all
// candidates within a wide energy window, then iteratively pick each gamma's
// candidate that best agrees with its already-chosen coincidence partners
// (iterated conditional modes), and finally accept only matches that are clean,
// structurally supported, or uniquely close -- so the result is never worse than
// the energy-only baseline.

#include "G4LevelData.h"

#include <vector>

namespace cascade_enrich {

struct SdGamma {
  double energy = 0.0;
  std::vector<int> partners;  ///< indices (into the gamma list) of coincidence partners
};

struct GammaMatch {
  int g4 = -1;                 ///< index into G4LevelScheme::transitions; -1 = unmatched
  double gap = 0.0;            ///< signed E_sd - E_g4 of the chosen match
  bool accepted = false;
  bool ambiguous = false;      ///< 2+ candidates within the tight tolerance
  bool structure_used = false; ///< accepted on structural support (>=1 chaining partner)
  bool structure_changed = false; ///< final pick differs from the nearest-energy candidate
};

struct MatchParams {
  double wide_tol = 2.5;    ///< candidate-collection window (keV)
  double tight_tol = 0.6;   ///< "clean unique" window (keV)
  double widen_tol = 1.0;   ///< unique-within-this accepted even without structure
  double accept_tol = 1.5;  ///< hard cap on |gap| for any accepted match
  double energy_sigma = 0.7;///< gaussian energy-score width (keV)
  double struct_bonus = 0.5;///< score added per chaining coincidence partner
  int iterations = 6;
};

/// Two transitions are adjacent (share a level in a descending cascade) if one's
/// daughter level is the other's parent level.
inline bool adjacent(const G4Transition& a, const G4Transition& b) {
  return a.to_level == b.from_level || b.to_level == a.from_level;
}

/// Match each SandiaDecay gamma to a G4 transition using energy + coincidence
/// structure. Returns one GammaMatch per input gamma (same order).
std::vector<GammaMatch> match_gammas(const G4LevelScheme& scheme,
                                     const std::vector<SdGamma>& gammas,
                                     const MatchParams& p = {});

} // namespace cascade_enrich

#endif
