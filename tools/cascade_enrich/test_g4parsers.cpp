// Standalone sanity tests for the GEANT4 level / decay parsers used by the
// cascade-detail enrichment. Returns non-zero on failure.
//
// Data directories are taken from $G4LEVELGAMMADATA / $G4RADIOACTIVEDATA, or the
// 1st/2nd command-line arguments, or a built-in default install path.

#include "G4LevelData.h"
#include "G4RadDecayData.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace cascade_enrich;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_failures;
  } else {
    std::cout << "ok:   " << msg << "\n";
  }
}

void check_close(double got, double want, double tol, const std::string& msg) {
  const bool ok = std::fabs(got - want) <= tol;
  if (!ok)
    std::cerr << "FAIL: " << msg << " (got " << got << ", want " << want
              << " +/- " << tol << ")\n";
  else
    std::cout << "ok:   " << msg << " = " << got << "\n";
  if (!ok)
    ++g_failures;
}

std::string env_or(const char* var, const std::string& fallback) {
  const char* v = std::getenv(var);
  return (v && v[0]) ? std::string(v) : fallback;
}

const G4Transition* find_gamma(const G4LevelScheme& s, double e, double tol) {
  const G4Transition* best = nullptr;
  double bd = tol;
  for (const G4Transition& t : s.transitions) {
    const double d = std::fabs(t.gamma_keV - e);
    if (d <= bd) { bd = d; best = &t; }
  }
  return best;
}

} // namespace

int main(int argc, char** argv) {
  const std::string photoneva = (argc > 1)
      ? argv[1]
      : env_or("G4LEVELGAMMADATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/PhotonEvaporation6.1.2");
  const std::string raddecay = (argc > 2)
      ? argv[2]
      : env_or("G4RADIOACTIVEDATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/RadioactiveDecay6.1.2");

  // --- Multipolarity decode ---
  {
    int L0, Lp;
    check(decode_multipolarity(4, L0, Lp) && L0 == 2 && Lp == 2, "decode E2 (4) -> L=2");
    check(decode_multipolarity(304, L0, Lp) && L0 == 1 && Lp == 2, "decode M1+E2 (304) -> 1,2");
    check(decode_multipolarity(407, L0, Lp) && L0 == 2 && Lp == 3, "decode E2+M3 (407) -> 2,3");
    check(decode_multipolarity(708, L0, Lp) && L0 == 3 && Lp == 4, "decode M3+E4 (708) -> 3,4");
    check(!decode_multipolarity(0, L0, Lp), "decode unknown (0) -> false");
  }

  // --- Ni-60 (z28.a60): Co-60 daughter level scheme ---
  {
    const G4LevelScheme ni60 = load_photon_evaporation(28, 60, photoneva);
    check(ni60.loaded, "load Ni-60 PhotonEvaporation");
    if (ni60.loaded) {
      // Level spins: ground 0+, first excited 2+ (1332.5), the 2505.75 level 4+.
      const G4Level* l0 = ni60.level(0);
      const G4Level* l1 = ni60.level(1);
      const G4Level* l4 = ni60.level(4);
      check(l0 && l0->twoJ == 0, "Ni-60 ground twoJ=0");
      check(l1 && l1->twoJ == 4 && std::fabs(l1->excitation_keV - 1332.514) < 0.1,
            "Ni-60 level1 = 1332.5 keV, twoJ=4 (2+)");
      check(l4 && l4->twoJ == 8 && std::fabs(l4->excitation_keV - 2505.753) < 0.1,
            "Ni-60 level4 = 2505.75 keV, twoJ=8 (4+)");

      const G4Transition* g1332 = find_gamma(ni60, 1332.514, 0.2);
      const G4Transition* g1173 = find_gamma(ni60, 1173.239, 0.2);
      check(g1332 && g1332->from_level == 1 && g1332->to_level == 0,
            "1332 keV transition 1->0");
      check(g1332 && g1332->L0 == 2, "1332 keV is E2 (L=2)");
      check(g1173 && g1173->from_level == 4 && g1173->to_level == 1,
            "1173 keV transition 4->1");
      check(g1173 && g1173->L0 == 2, "1173 keV lower multipole L=2");
      if (g1173)
        check_close(g1173->alpha_total, 0.0001722, 1e-6, "1173 keV alpha_total");
      if (g1173)
        check_close(g1173->icc_shell[0], 0.8993, 0.01, "1173 keV f_K");
    }
  }

  // --- Ba-133 EC shell fractions (z56.a133) ---
  {
    const G4ECData ba133 = load_ec_shell_fractions(56, 133, raddecay);
    check(ba133.loaded, "load Ba-133 RadioactiveDecay");
    if (ba133.loaded) {
      check(ba133.has_summary, "Ba-133 has summary EC split");
      // Summary: K 0.68638, L 0.23978, M 0.07384 (NshellEC absent) -> normalized.
      check_close(ba133.summary.fK, 0.6864, 0.005, "Ba-133 summary f_K");
      // Per-level: daughter 437.0113 keV (the 81 keV gamma's parent level):
      // K 57.206, L 21.399, M 6.6243 -> f_K = 57.206 / 85.23 = 0.671.
      const ECShellFractions* f = ba133.for_level(437.011, 0.6);
      check(f != nullptr, "Ba-133 EC level 437.01 keV found");
      if (f)
        check_close(f->fK, 0.671, 0.01, "Ba-133 437 keV level f_K");
    }
  }

  // --- I-125 EC shell fractions (z53.a125) ---
  {
    const G4ECData i125 = load_ec_shell_fractions(53, 125, raddecay);
    check(i125.loaded, "load I-125 RadioactiveDecay");
    if (i125.loaded) {
      const ECShellFractions* f = i125.for_level(35.4925, 0.6);
      check(f != nullptr, "I-125 EC level 35.49 keV found");
      if (f)
        check_close(f->fK, 0.8011, 0.005, "I-125 35.49 keV level f_K");
    }
  }

  if (g_failures == 0) {
    std::cout << "\nALL PARSER TESTS PASSED\n";
    return 0;
  }
  std::cerr << "\n" << g_failures << " parser test(s) FAILED\n";
  return 1;
}
