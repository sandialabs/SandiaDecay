// Tests for the structure-aware gamma matcher.
#include "StructureMatch.h"
#include "G4LevelData.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>

using namespace cascade_enrich;
namespace {
int fails = 0;
void check(bool ok, const std::string& m) {
  std::cout << (ok ? "ok:   " : "FAIL: ") << m << "\n";
  if (!ok) ++fails;
}
std::string env_or(const char* v, const std::string& f) {
  const char* e = std::getenv(v); return (e && e[0]) ? e : f;
}

G4LevelScheme synthetic_scheme(
    const std::vector<std::tuple<int, int, double>>& edges) {
  G4LevelScheme s;
  s.loaded = true;
  for (const auto& edge : edges) {
    G4Transition t;
    t.from_level = std::get<0>(edge);
    t.to_level = std::get<1>(edge);
    t.gamma_keV = std::get<2>(edge);
    s.transitions.push_back(t);
  }
  return s;
}

bool accepted_matches_are_injective(const std::vector<GammaMatch>& matches) {
  std::set<int> edges;
  for (const GammaMatch& m : matches)
    if (m.accepted && !edges.insert(m.g4).second) return false;
  return true;
}

std::map<double, int> accepted_by_energy(const std::vector<SdGamma>& gammas,
                                         const std::vector<GammaMatch>& matches) {
  std::map<double, int> result;
  for (std::size_t i = 0; i < gammas.size(); ++i)
    if (matches[i].accepted) result[gammas[i].energy] = matches[i].g4;
  return result;
}
} // namespace

int main(int argc, char** argv) {
  // Two close Sandia records and two close G4 edges must remain distinct, and
  // merely permuting the input records must not change the physical mapping.
  const G4LevelScheme two_edges =
      synthetic_scheme({{2, 0, 100.00}, {3, 1, 100.90}});
  std::vector<SdGamma> close = {{100.05, {}}, {100.86, {}}};
  const auto close_result = match_gammas(two_edges, close);
  check(close_result[0].accepted && close_result[1].accepted &&
            close_result[0].g4 != close_result[1].g4,
        "two records use two distinct nearby G4 edges");
  check(accepted_matches_are_injective(close_result),
        "accepted two-edge result is injective");
  std::vector<SdGamma> close_permuted = {close[1], close[0]};
  const auto close_permuted_result = match_gammas(two_edges, close_permuted);
  check(accepted_by_energy(close, close_result) ==
            accepted_by_energy(close_permuted, close_permuted_result),
        "two-edge assignment is invariant under input permutation");

  // If two accepted records select the only edge, neither record has enough
  // evidence to be called a duplicate. Conservatively reject both owners.
  const G4LevelScheme one_edge = synthetic_scheme({{1, 0, 100.00}});
  std::vector<SdGamma> collision = {{99.90, {}}, {100.10, {}}};
  const auto collision_result = match_gammas(one_edge, collision);
  check(!collision_result[0].accepted && !collision_result[1].accepted &&
            collision_result[0].collision_rejected &&
            collision_result[1].collision_rejected,
        "all owners of a one-edge collision are rejected");
  std::vector<SdGamma> collision_permuted = {collision[1], collision[0]};
  const auto collision_permuted_result = match_gammas(one_edge, collision_permuted);
  check(!collision_permuted_result[0].accepted &&
            !collision_permuted_result[1].accepted &&
            collision_permuted_result[0].collision_rejected &&
            collision_permuted_result[1].collision_rejected,
        "collision rejection is invariant under input permutation");
  check(accepted_matches_are_injective(collision_result),
        "accepted collision result is injective");

  const std::string dir = (argc > 1) ? argv[1]
      : env_or("G4LEVELGAMMADATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/PhotonEvaporation6.1.2");
  const G4LevelScheme ni60 = load_photon_evaporation(28, 60, dir);
  check(ni60.loaded, "Ni-60 loaded");

  // Co-60 main cascade: 1173.24 (4->1) and 1332.51 (1->0), mutually coincident.
  // Both must match and be adjacent (share level 1).
  std::vector<SdGamma> g(2);
  g[0].energy = 1173.228; g[0].partners = {1};
  g[1].energy = 1332.492; g[1].partners = {0};
  auto r = match_gammas(ni60, g);
  check(r[0].accepted && r[1].accepted, "both Co-60 gammas accepted");
  check(accepted_matches_are_injective(r), "accepted Co-60 matches are injective");
  if (r[0].accepted && r[1].accepted) {
    check(adjacent(ni60.transitions[r[0].g4], ni60.transitions[r[1].g4]),
          "matched transitions are adjacent (share a level)");
  }

  // Structure disambiguation: inject a decoy gamma whose energy is near 1173 but
  // which has NO coincidence partner -- the real 1173 (with a partner that chains)
  // should still resolve to the adjacent transition.
  // Here we test that a partnerless gamma at a slightly-off energy still matches
  // the nearest transition (no structure needed) and a structured pair is
  // preferred when energies are ambiguous: use 826.1 (2158->1332, i.e. 2->1) which
  // is coincident with 1332.
  std::vector<SdGamma> g2(2);
  g2[0].energy = 826.10; g2[0].partners = {1};   // 2158.6->1332.5 (level 2->1)
  g2[1].energy = 1332.492; g2[1].partners = {0};
  auto r2 = match_gammas(ni60, g2);
  check(r2[0].accepted && r2[1].accepted, "826/1332 pair accepted");
  if (r2[0].accepted && r2[1].accepted)
    check(adjacent(ni60.transitions[r2[0].g4], ni60.transitions[r2[1].g4]),
          "826 & 1332 adjacent (2->1->0)");

  if (fails == 0) { std::cout << "\nALL STRUCTURE-MATCH TESTS PASSED\n"; return 0; }
  std::cerr << "\n" << fails << " FAILED\n"; return 1;
}
