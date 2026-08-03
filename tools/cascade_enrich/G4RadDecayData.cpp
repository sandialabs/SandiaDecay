#include "G4RadDecayData.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>

namespace cascade_enrich {

const char* parent_block_selection_name(ParentBlockSelectionMethod method) {
  switch (method) {
    case ParentBlockSelectionMethod::GroundExcitation: return "ground-excitation";
    case ParentBlockSelectionMethod::MassExcitation: return "mass-excitation";
    case ParentBlockSelectionMethod::HalfLife: return "half-life";
    case ParentBlockSelectionMethod::SingleBlockHalfLifeDrift: return "single-block-half-life-drift";
    case ParentBlockSelectionMethod::OrdinalFallback: return "ordinal-fallback";
    case ParentBlockSelectionMethod::None: return "no-match";
  }
  return "no-match";
}

const ECShellFractions* G4ECData::for_level(double ex_keV, double tol) const {
  const ECShellFractions* best = nullptr;
  double best_d = tol;
  for (const ECShellFractions& f : per_level) {
    const double d = std::fabs(f.daughter_ex_keV - ex_keV);
    if (d <= best_d) {
      best_d = d;
      best = &f;
    }
  }
  return best;
}

const ECShellFractions* G4ECData::for_q(double q_keV, double tol,
                                        bool* ambiguous) const {
  if (ambiguous) *ambiguous = false;
  const ECShellFractions* best = nullptr;
  double best_d = std::numeric_limits<double>::infinity();
  double second_d = std::numeric_limits<double>::infinity();
  for (const ECShellFractions& f : per_level) {
    const double d = std::fabs(f.q_keV - q_keV);
    if (d < best_d) {
      second_d = best_d;
      best_d = d;
      best = &f;
    } else if (d < second_d) {
      second_d = d;
    }
  }
  if (!best || best_d > tol)
    return nullptr;
  constexpr double kNearestMarginKeV = 0.5;
  if (second_d <= tol && second_d - best_d <= kNearestMarginKeV) {
    if (ambiguous) *ambiguous = true;
    return nullptr;
  }
  return best;
}

const ECShellFractions* G4ECData::shells_for_q(double q_keV, double tol,
                                               bool* ambiguous) const {
  if (const ECShellFractions* exact = for_q(q_keV, tol, ambiguous))
    return exact;
  return has_summary ? &summary : nullptr;
}

const BetaPlusRecord* G4ECData::beta_plus_for_q(double q_keV, double tol,
                                                bool* ambiguous) const {
  if (ambiguous) *ambiguous = false;
  const BetaPlusRecord* best = nullptr;
  double best_d = std::numeric_limits<double>::infinity();
  double second_d = std::numeric_limits<double>::infinity();
  for (const BetaPlusRecord& b : beta_plus) {
    const double d = std::fabs(b.q_keV - q_keV);
    if (d < best_d) {
      second_d = best_d;
      best_d = d;
      best = &b;
    } else if (d < second_d) {
      second_d = d;
    }
  }
  if (!best || best_d > tol)
    return nullptr;
  constexpr double kNearestMarginKeV = 0.5;
  if (second_d <= tol && second_d - best_d <= kNearestMarginKeV) {
    if (ambiguous) *ambiguous = true;
    return nullptr;
  }
  return best;
}

namespace {

std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> t;
  std::istringstream is(line);
  std::string tok;
  while (is >> tok)
    t.push_back(tok);
  return t;
}

// Returns 'K'/'L'/'M'/'N' for an EC mode token like "KshellEC", else 0.
char ec_shell(const std::string& mode) {
  if (mode.size() < 7)
    return 0;
  if (mode.compare(1, 6, "shellE") == 0 && mode.back() == 'C') {
    const char s = mode[0];
    if (s == 'K' || s == 'L' || s == 'M' || s == 'N')
      return s;
  }
  return 0;
}

struct Accum { double K = 0, L = 0, M = 0, N = 0; };

struct ParentBlock {
  double excitation_keV = 0.0;
  double half_life_s = 0.0;
  std::vector<std::vector<std::string>> records;
};

bool usable_half_life(double value) {
  return std::isfinite(value) && value > 0.0;
}

bool half_life_agrees(double a, double b) {
  // The two files are independently rounded exports of essentially the same
  // evaluated data. Five percent admits normal rounding/evaluation drift while
  // preventing a merely similar-lived ground state from masquerading as an
  // isomer.
  return usable_half_life(a) && usable_half_life(b) &&
         std::fabs(std::log(a / b)) <= std::log(1.05);
}

bool half_life_compatible_with_mass_match(double a, double b) {
  // A sub-5-keV mass match already supplies strong state identity. Use a wider
  // lifetime consistency gate here for independently rounded evaluations (for
  // example Er-149m is 8.9 s vs 9.6 s), while still excluding prompt states.
  return usable_half_life(a) && usable_half_life(b) &&
         std::fabs(std::log(a / b)) <= std::log(1.25);
}

double half_life_error(double a, double b) {
  return (usable_half_life(a) && usable_half_life(b))
      ? std::fabs(std::log(a / b))
      : std::numeric_limits<double>::infinity();
}

void add(Accum& a, char shell, double v) {
  switch (shell) {
    case 'K': a.K += v; break;
    case 'L': a.L += v; break;
    case 'M': a.M += v; break;
    case 'N': a.N += v; break;
    default: break;
  }
}

ECShellFractions normalize(double ex, double q, const Accum& a) {
  ECShellFractions f;
  f.daughter_ex_keV = ex;
  f.q_keV = q;
  const double tot = a.K + a.L + a.M + a.N;
  f.total_intensity = tot;
  if (tot > 0.0) {
    f.fK = a.K / tot;
    f.fL = a.L / tot;
    f.fM = a.M / tot;
    f.fN = a.N / tot;
  }
  return f;
}

} // namespace

G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir,
                                 const ParentStateHint& parent) {
  G4ECData data;
  data.Z = Zparent;
  data.A = Aparent;
  data.parent_isomer = parent.isomer;

  if (parent.isomer < 0)
    return data;

  std::string path = data_dir;
  if (!path.empty() && path.back() != '/')
    path += '/';
  path += "z" + std::to_string(Zparent) + ".a" + std::to_string(Aparent);

  std::ifstream in(path);
  if (!in)
    return data;  // loaded == false

  // Read every P block first. Stable ground states are often omitted entirely
  // (Ba-137 begins with its 661.659-keV isomer), so isomerNumber cannot select a
  // block until excitation/half-life metadata has been considered.
  std::vector<ParentBlock> blocks;
  ParentBlock* current = nullptr;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    const std::vector<std::string> t = tokenize(line);
    if (t.empty())
      continue;
    if (t[0] == "P") {
      ParentBlock b;
      b.excitation_keV = (t.size() > 1) ? std::atof(t[1].c_str()) : 0.0;
      b.half_life_s = (t.size() > 3) ? std::atof(t[3].c_str()) : 0.0;
      blocks.push_back(b);
      current = &blocks.back();
    } else if (current) {
      current->records.push_back(t);
    }
  }

  constexpr double kExcitationToleranceKeV = 5.0;
  const bool have_exc = std::isfinite(parent.excitation_keV) &&
                        parent.excitation_keV >= 0.0;
  const bool have_half = usable_half_life(parent.half_life_s);
  const auto is_known_ground_block = [&](std::size_t i) {
    return parent.isomer > 0 && std::fabs(blocks[i].excitation_keV) <= 0.1 &&
           half_life_error(blocks[i].half_life_s, parent.ground_half_life_s) <
               half_life_error(blocks[i].half_life_s, parent.half_life_s);
  };
  const auto nearest_excitation = [&](const std::vector<int>& candidates,
                                      bool require_close) {
    int best_index = -1;
    double best = require_close ? kExcitationToleranceKeV
                                : std::numeric_limits<double>::infinity();
    for (const int i : candidates) {
      const double d = std::fabs(blocks[i].excitation_keV - parent.excitation_keV);
      if (d < best) { best = d; best_index = i; }
    }
    return best_index;
  };

  int selected = -1;

  // Half-life is the state identity for both ground and isomer parents. GEANT4
  // may shift the long-lived state to ex=0 and put the physical ground state at
  // a nonzero P excitation (Eu-152, Y-85, Np-236), so ex=0 is not a ground-state
  // discriminator. Excitation is used only to resolve half-life ambiguity.
  if (have_half) {
    std::vector<int> half_candidates;
    double best_half = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < blocks.size(); ++i) {
      if (is_known_ground_block(i) ||
          !half_life_agrees(blocks[i].half_life_s, parent.half_life_s))
        continue;
      const double e = half_life_error(blocks[i].half_life_s, parent.half_life_s);
      if (e + 1e-12 < best_half) {
        best_half = e;
        half_candidates.assign(1, static_cast<int>(i));
      } else if (std::fabs(e - best_half) <= 1e-12) {
        half_candidates.push_back(static_cast<int>(i));
      }
    }

    if (!half_candidates.empty()) {
      if (!parent.half_life_identity_unique) {
        // Another Sandia state has the same lifetime. Require an independent
        // excitation match; otherwise both states would silently claim one P
        // block.
        if (have_exc)
          selected = nearest_excitation(half_candidates, true);
      } else if (half_candidates.size() == 1) {
        selected = half_candidates.front();
      } else if (have_exc) {
        selected = nearest_excitation(half_candidates, true);
      }
      if (selected >= 0)
        data.selection_method = ParentBlockSelectionMethod::HalfLife;
    }
  }

  // If strict lifetime identity misses because evaluations drift, a precise
  // excitation match may recover the state, but only with a broad lifetime
  // consistency check so a prompt state cannot be selected.
  if (selected < 0 && have_exc) {
    std::vector<int> excitation_candidates;
    for (std::size_t i = 0; i < blocks.size(); ++i) {
      const double d = std::fabs(blocks[i].excitation_keV - parent.excitation_keV);
      if (d <= kExcitationToleranceKeV && !is_known_ground_block(i) &&
          (!have_half || half_life_compatible_with_mass_match(
                             blocks[i].half_life_s, parent.half_life_s)))
        excitation_candidates.push_back(static_cast<int>(i));
    }
    selected = nearest_excitation(excitation_candidates, true);
    if (selected >= 0)
      data.selection_method = (parent.isomer == 0)
          ? ParentBlockSelectionMethod::GroundExcitation
          : ParentBlockSelectionMethod::MassExcitation;
  }

  // Final physical fallback: a file with one P block and modest lifetime drift.
  // Do not use it for a state whose lifetime is shared with a sibling or when
  // that sole block is better identified as the known ground state.
  if (selected < 0 && blocks.size() == 1 && have_half &&
      parent.half_life_identity_unique && !is_known_ground_block(0) &&
      half_life_compatible_with_mass_match(blocks[0].half_life_s,
                                           parent.half_life_s)) {
    selected = 0;
    data.selection_method = ParentBlockSelectionMethod::SingleBlockHalfLifeDrift;
  }

  // Auditable compatibility fallback only when no physical hint exists.
  if (selected < 0 && !have_exc && !have_half &&
      parent.isomer < static_cast<int>(blocks.size())) {
    selected = parent.isomer;
    data.selection_method = ParentBlockSelectionMethod::OrdinalFallback;
  }

  if (selected < 0)
    return data;

  data.parent_block = selected;
  data.parent_ex_keV = blocks[selected].excitation_keV;
  data.parent_half_life_s = blocks[selected].half_life_s;

  // Group selected-block per-level (5-token) EC records by rounded daughter
  // excitation, with 3-token summary records separately. Alpha/BetaMinus are
  // also retained verbatim as transition-level daughter-feed evidence; unlike
  // endpoint matching, those records already carry the exact daughter
  // excitation and floating-level family.
  std::map<long, Accum> per_level_acc;       // key = round(daughter_ex*1000)
  std::map<long, double> per_level_ex;       // key -> representative ex
  std::map<long, double> per_level_q;        // key -> representative Q
  Accum summary_acc;
  for (const std::vector<std::string>& t : blocks[selected].records) {
    DecayModeLevelFeeds* mode_feeds = nullptr;
    if (t[0] == "Alpha")
      mode_feeds = &data.alpha;
    else if (t[0] == "BetaMinus")
      mode_feeds = &data.beta_minus;
    if (mode_feeds) {
      if (t.size() == 3) {
        mode_feeds->branch_fraction = std::atof(t[2].c_str());
        ++mode_feeds->summary_count;
      } else if (t.size() >= 5) {
        DaughterLevelFeedRecord r;
        r.daughter_ex_keV = std::atof(t[1].c_str());
        r.floating_flag = t[2];
        r.intensity_percent = std::atof(t[3].c_str());
        r.q_keV = std::atof(t[4].c_str());
        if (t.size() >= 6)
          r.forbiddenness = t[5];
        mode_feeds->records.push_back(r);
      }
      continue;
    }

    if (t[0] == "BetaPlus" && t.size() >= 5) {
      BetaPlusRecord b;
      b.daughter_ex_keV = std::atof(t[1].c_str());
      b.floating_flag = t[2];
      b.intensity = std::atof(t[3].c_str());
      b.q_keV = std::atof(t[4].c_str());
      data.beta_plus.push_back(b);
      continue;
    }

    const char shell = ec_shell(t[0]);
    if (!shell)
      continue;

    if (t.size() == 3) {
      // Summary record: Mode daughterEx(=0) intensity(fraction)
      add(summary_acc, shell, std::atof(t[2].c_str()));
    } else if (t.size() >= 5) {
      // Per-level record: Mode daughterEx flag intensity(%) Q
      const double ex = std::atof(t[1].c_str());
      const double inten = std::atof(t[3].c_str());
      const double q = std::atof(t[4].c_str());
      const long key = std::lround(ex * 1000.0);
      add(per_level_acc[key], shell, inten);
      per_level_ex[key] = ex;
      per_level_q[key] = q;
    }
  }

  for (const auto& kv : per_level_acc)
    data.per_level.push_back(
        normalize(per_level_ex[kv.first], per_level_q[kv.first], kv.second));

  if (summary_acc.K + summary_acc.L + summary_acc.M + summary_acc.N > 0.0) {
    data.summary = normalize(0.0, 0.0, summary_acc);
    data.has_summary = true;
  }

  data.loaded = true;
  return data;
}

G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir,
                                 int parent_isomer) {
  ParentStateHint parent;
  parent.isomer = parent_isomer;
  return load_ec_shell_fractions(Zparent, Aparent, data_dir, parent);
}

} // namespace cascade_enrich
