#include "StructureMatch.h"

#include <algorithm>
#include <cmath>

namespace cascade_enrich {

namespace {
double energy_score(double gap, double sigma) {
  const double x = gap / sigma;
  return std::exp(-x * x);
}
} // namespace

std::vector<GammaMatch> match_gammas(const G4LevelScheme& scheme,
                                     const std::vector<SdGamma>& gammas,
                                     const MatchParams& p) {
  const std::size_t n = gammas.size();
  std::vector<GammaMatch> out(n);
  if (!scheme.loaded)
    return out;

  // Candidate G4 transition indices within the wide window, per gamma, sorted by
  // ascending |gap| (so candidates[g][0] is the nearest-energy choice).
  std::vector<std::vector<int>> cands(n);
  std::vector<int> nearest(n, -1);
  std::vector<int> assign(n, -1);
  for (std::size_t g = 0; g < n; ++g) {
    std::vector<std::pair<double,int>> c;
    for (std::size_t t = 0; t < scheme.transitions.size(); ++t) {
      const double d = std::fabs(scheme.transitions[t].gamma_keV - gammas[g].energy);
      if (d <= p.wide_tol) c.emplace_back(d, static_cast<int>(t));
    }
    std::sort(c.begin(), c.end());
    for (auto& pr : c) cands[g].push_back(pr.second);
    if (!cands[g].empty()) {
      nearest[g] = cands[g][0];
      assign[g] = cands[g][0];  // initialize at nearest energy
    }
  }

  // Iterated conditional modes: re-pick each gamma's candidate to maximize
  // energy score + structural agreement with partners' current picks.
  for (int it = 0; it < p.iterations; ++it) {
    bool changed = false;
    for (std::size_t g = 0; g < n; ++g) {
      if (cands[g].empty()) continue;
      double best_s = -1e18; int best = assign[g];
      for (int c : cands[g]) {
        const G4Transition& tc = scheme.transitions[c];
        const double gap = gammas[g].energy - tc.gamma_keV;
        double s = energy_score(gap, p.energy_sigma);
        for (int pidx : gammas[g].partners) {
          if (pidx < 0 || pidx >= static_cast<int>(n)) continue;
          const int pa = assign[pidx];
          if (pa < 0) continue;
          if (adjacent(tc, scheme.transitions[pa])) s += p.struct_bonus;
        }
        if (s > best_s) { best_s = s; best = c; }
      }
      if (best != assign[g]) { assign[g] = best; changed = true; }
    }
    if (!changed) break;
  }

  // Acceptance + diagnostics.
  for (std::size_t g = 0; g < n; ++g) {
    GammaMatch& m = out[g];
    if (assign[g] < 0) continue;
    const G4Transition& tc = scheme.transitions[assign[g]];
    const double gap = gammas[g].energy - tc.gamma_keV;
    m.g4 = assign[g];
    m.gap = gap;

    // Ambiguity: count candidates within the tight tolerance.
    int n_tight = 0, n_widen = 0;
    for (int c : cands[g]) {
      const double d = std::fabs(scheme.transitions[c].gamma_keV - gammas[g].energy);
      if (d <= p.tight_tol) ++n_tight;
      if (d <= p.widen_tol) ++n_widen;
    }
    m.ambiguous = (n_tight >= 2);

    // Structural support: a coincidence partner whose own (assigned) transition
    // is adjacent to this one.
    bool supported = false;
    for (int pidx : gammas[g].partners) {
      if (pidx < 0 || pidx >= static_cast<int>(n)) continue;
      const int pa = assign[pidx];
      if (pa >= 0 && adjacent(tc, scheme.transitions[pa])) { supported = true; break; }
    }

    const double agap = std::fabs(gap);
    const bool clean  = (agap <= p.tight_tol);
    const bool widen_unique = (agap <= p.widen_tol && n_widen == 1);
    const bool struct_ok = (agap <= p.accept_tol && supported);

    m.accepted = clean || widen_unique || struct_ok;
    if (m.accepted) {
      m.structure_used = struct_ok && !clean && !widen_unique;
      m.structure_changed = (assign[g] != nearest[g]);
    } else {
      m.g4 = -1;  // reject -> behaves as un-enriched (no worse than baseline)
    }
  }
  return out;
}

} // namespace cascade_enrich
