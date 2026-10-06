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

  // Per-gamma effective tolerances: max(fixed floor, rel_tol * energy).
  auto eff = [&](double floor, double E) { return std::max(floor, p.rel_tol * E); };
  std::vector<double> tol_wide(n), tol_tight(n), tol_widen(n), tol_accept(n), sig(n);
  for (std::size_t g = 0; g < n; ++g) {
    const double E = gammas[g].energy;
    tol_wide[g]   = eff(p.wide_tol, E);
    tol_tight[g]  = eff(p.tight_tol, E);
    tol_widen[g]  = eff(p.widen_tol, E);
    tol_accept[g] = eff(p.accept_tol, E);
    sig[g]        = eff(p.energy_sigma, E);
  }

  // Candidate G4 transition indices within the wide window, per gamma, sorted by
  // ascending |gap| (so candidates[g][0] is the nearest-energy choice).
  std::vector<std::vector<int>> cands(n);
  std::vector<int> nearest(n, -1);
  std::vector<int> assign(n, -1);
  for (std::size_t g = 0; g < n; ++g) {
    if (gammas[g].reviewed_unplaced)
      continue;  // no candidates: it can neither match nor collide with a placed line
    std::vector<std::pair<double,int>> c;
    for (std::size_t t = 0; t < scheme.transitions.size(); ++t) {
      const double d = std::fabs(scheme.transitions[t].gamma_keV - gammas[g].energy);
      if (d <= tol_wide[g]) c.emplace_back(d, static_cast<int>(t));
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
        double s = energy_score(gap, sig[g]);
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
      if (d <= tol_tight[g]) ++n_tight;
      if (d <= tol_widen[g]) ++n_widen;
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
    const bool clean  = (agap <= tol_tight[g]);
    const bool widen_unique = (agap <= tol_widen[g] && n_widen == 1);
    const bool struct_ok = (agap <= tol_accept[g] && supported);

    m.accepted = clean || widen_unique || struct_ok;
    if (m.accepted) {
      m.structure_used = struct_ok && !clean && !widen_unique;
      m.structure_changed = (assign[g] != nearest[g]);
    } else {
      m.g4 = -1;  // reject -> behaves as un-enriched (no worse than baseline)
    }
  }

  // Enforce the graph invariant required by every downstream consumer: one
  // Sandia gamma may own a G4 transition edge, never two. A full global solve
  // would still need nuclear-data evidence to distinguish close candidate
  // edges with materially different ICCs, so conservatively reject every
  // member of an accepted-edge collision. In particular, do not keep the first
  // record and call the remainder duplicates merely because G4 has one edge.
  std::vector<std::vector<std::size_t>> owners(scheme.transitions.size());
  for (std::size_t g = 0; g < n; ++g)
    if (out[g].accepted && out[g].g4 >= 0)
      owners[static_cast<std::size_t>(out[g].g4)].push_back(g);
  for (const std::vector<std::size_t>& collision : owners) {
    if (collision.size() < 2) continue;
    for (const std::size_t g : collision) {
      out[g].accepted = false;
      out[g].collision_rejected = true;
      out[g].g4 = -1;
    }
  }
  return out;
}

} // namespace cascade_enrich
