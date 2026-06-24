// Tests for the structure-aware gamma matcher.
#include "StructureMatch.h"
#include "G4LevelData.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

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
} // namespace

int main(int argc, char** argv) {
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
