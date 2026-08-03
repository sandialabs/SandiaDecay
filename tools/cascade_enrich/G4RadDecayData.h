#ifndef CASCADE_ENRICH_G4RADDECAYDATA_H
#define CASCADE_ENRICH_G4RADDECAYDATA_H

// Parser for GEANT4 RadioactiveDecay data (z{Z}.a{A} ASCII files), specifically
// the electron-capture shell fractions needed to attribute K/L vacancies (and
// hence x-ray emission) to each daughter level fed by EC, plus the real
// per-level beta-plus records needed to identify the level fed by a positron.
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
  double q_keV = 0.0;            ///< EC Q-value to this level (SandiaDecay keys EC on this)
  double fK = 0.0, fL = 0.0, fM = 0.0, fN = 0.0;  ///< normalized shell fractions (sum 1)
  double total_intensity = 0.0;  ///< sum of raw K+L+M+N intensities for this level
};

struct BetaPlusRecord {
  double daughter_ex_keV = 0.0;
  std::string floating_flag;
  double q_keV = 0.0;           ///< EC Q value; beta-plus endpoint is Q - 2*m_e
  double intensity = 0.0;       ///< raw GEANT4 per-level intensity (percent)
};

/** One five-column RadioactiveDecay alpha or beta-minus daughter-feed record.
    `intensity_percent` is absolute percent of parent decays; divide by
    `100*DecayModeLevelFeeds::branch_fraction` for probability conditional on
    the selected decay mode. The floating flag is part of the daughter-level
    identity and must not be discarded: distinct fixed and floating levels may
    have exactly the same excitation energy. */
struct DaughterLevelFeedRecord {
  double daughter_ex_keV = 0.0;
  std::string floating_flag;
  double intensity_percent = 0.0;
  double q_keV = 0.0;
  std::string forbiddenness;  ///< raw optional RDM beta-shape token
};

struct DecayModeLevelFeeds {
  double branch_fraction = 0.0;  ///< selected P-block mode summary, fraction of unity
  int summary_count = 0;         ///< must be exactly one before feeds are deliverable
  std::vector<DaughterLevelFeedRecord> records;
};

enum class ParentBlockSelectionMethod {
  None,
  GroundExcitation,
  MassExcitation,
  HalfLife,
  SingleBlockHalfLifeDrift,
  OrdinalFallback
};

struct ParentStateHint {
  int isomer = 0;
  double excitation_keV = -1.0;  ///< 0 for ground, mass-derived for isomer
  double half_life_s = -1.0;     ///< negative/nonfinite when unavailable
  double ground_half_life_s = -1.0;  ///< rejects a ground block with similar lifetime
  bool half_life_identity_unique = true;  ///< false when another same-Z/A state shares it
};

const char* parent_block_selection_name(ParentBlockSelectionMethod method);

struct G4ECData {
  int Z = 0, A = 0;
  int parent_isomer = 0;                    ///< requested Sandia isomer number
  int parent_block = -1;                    ///< actual selected P-block ordinal
  double parent_ex_keV = 0.0;
  double parent_half_life_s = 0.0;
  ParentBlockSelectionMethod selection_method = ParentBlockSelectionMethod::None;
  std::vector<ECShellFractions> per_level;  ///< selected parent, one per fed daughter level
  std::vector<BetaPlusRecord> beta_plus;    ///< real 5-token BetaPlus records only
  DecayModeLevelFeeds alpha;
  DecayModeLevelFeeds beta_minus;
  ECShellFractions summary;                 ///< overall K/L/M/N split (fallback)
  bool has_summary = false;
  bool loaded = false;                      ///< requested P block was found

  /// Exact per-level EC shell fractions for the daughter level nearest
  /// `ex_keV` (within `tol`). Does not fall back to the summary: callers must
  /// not manufacture a fed level from a summary record.
  const ECShellFractions* for_level(double ex_keV, double tol = 0.6) const;

  /// Exact per-level EC shell fractions for the level whose EC Q-value is
  /// nearest `q_keV` (SandiaDecay's <electronCapture energy=...> is the Q to
  /// the level). Does not fall back to the summary.
  const ECShellFractions* for_q(double q_keV, double tol = 1.0,
                                bool* ambiguous = nullptr) const;

  /// Shell fractions usable for vacancy yields. This may fall back to the
  /// summary split, but the returned summary has no daughter-level identity.
  const ECShellFractions* shells_for_q(double q_keV, double tol = 1.0,
                                       bool* ambiguous = nullptr) const;

  /// Real per-level BetaPlus record nearest an EC-Q value. Summary BetaPlus
  /// records are deliberately excluded because they carry no fed level.
  const BetaPlusRecord* beta_plus_for_q(double q_keV, double tol = 1.0,
                                        bool* ambiguous = nullptr) const;
};

/// Load RadioactiveDecay z{Z}.a{A} (the PARENT nuclide's file; EC fractions live
/// there keyed by daughter excitation). Parent P blocks are selected from
/// physical state hints, because GEANT4 omits stable ground-state blocks and a
/// Sandia isomer number is therefore not generally a P-block ordinal. All
/// states use half-life identity first; excitation corroborates or breaks ties.
/// A conservative single-block lifetime-drift fallback is explicitly reported.
/// Ordinal is a last fallback only when all physical hints are absent.
G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir,
                                 const ParentStateHint& parent);

/// Compatibility overload for parser clients without physical state metadata.
/// Non-ground selection will be reported as OrdinalFallback.
G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir,
                                 int parent_isomer = 0);

} // namespace cascade_enrich

#endif
