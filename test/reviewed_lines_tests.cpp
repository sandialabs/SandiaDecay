// Regression tests for individually reviewed decay lines (2026-10) that are not part of the U-235 / U-238 reviews;
// each value is checked against ENSDF and DDEP, see the commit messages.  Run once per delivered XML variant.
//
// Usage: reviewed_lines_tests <sandia.decay[.min|.nocoinc.min].xml>

#include "SandiaDecay.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int g_failures = 0;

void check(bool value, const std::string& message) {
  std::cout << (value ? "ok: " : "FAIL: ") << message << "\n";
  if (!value) ++g_failures;
}

bool near(double value, double expected, double rel_tol) {
  return std::fabs(value - expected) <= rel_tol * std::fabs(expected);
}

const SandiaDecay::Transition* find_transition(const SandiaDecay::SandiaDecayDataBase& db,
                                               const std::string& parent, const std::string& child,
                                               SandiaDecay::DecayMode mode) {
  const SandiaDecay::Nuclide* p = db.nuclide(parent);
  if (!p) return nullptr;
  for (const SandiaDecay::Transition* tr : p->decaysToChildren)
    if (tr->child && tr->child->symbol == child && tr->mode == mode)
      return tr;
  return nullptr;
}

// Sum of intensities of particles of `type` with energy in [lo_keV, hi_keV].
double sum_in(const SandiaDecay::Transition* tr, SandiaDecay::ProductType type, double lo_keV, double hi_keV,
              int* count = nullptr) {
  double s = 0.0;
  int n = 0;
  for (const SandiaDecay::RadParticle& p : tr->products)
    if (p.type == type && p.energy >= lo_keV && p.energy <= hi_keV) {
      s += p.intensity;
      ++n;
    }
  if (count) *count = n;
  return s;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <sandia.decay.xml>\n";
    return 2;
  }
  const std::string file = argv[1];

  try {
    SandiaDecay::SandiaDecayDataBase db(file);

    // --- Sr-89 / Y-89m (2026-10-19): ENSDF NDS 114, 1 (2013) and DDEP (2004): 909 keV is emitted only by the
    //     Y-89m IT, 0.00956(6) per 100 Sr-89 decays; beta- to Y-89m 0.00964 %. ---
    const SandiaDecay::Transition* sr_gs = find_transition(db, "Sr89", "Y89", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* sr_m = find_transition(db, "Sr89", "Y89m", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* y_it = find_transition(db, "Y89m", "Y89", SandiaDecay::IsometricTransitionDecay);
    check(sr_gs && sr_m && y_it, "Sr89 -> Y89, Sr89 -> Y89m and Y89m -> Y89 transitions exist");
    if (sr_gs && sr_m && y_it) {
      int n = 0;
      check(sum_in(sr_gs, SandiaDecay::GammaParticle, 0.0, 1.0E6, &n) == 0.0 && n == 0,
            "Sr89 -> Y89 ground-state feed emits no gamma");
      check(sum_in(sr_m, SandiaDecay::GammaParticle, 0.0, 1.0E6, &n) == 0.0 && n == 0,
            "Sr89 -> Y89m feed emits no gamma");
      check(near(sr_m->branchRatio, 9.64E-05, 1e-6), "Sr89 -> Y89m branch ratio 9.64E-05");
      check(near(sum_in(y_it, SandiaDecay::GammaParticle, 908.9, 909.0, &n), 0.9916, 1e-6) && n == 1,
            "Y89m IT: 908.96 keV gamma 0.9916 per decay");

      // Full chain: per Sr-89 decay, with Y-89m (16 s) in secular equilibrium.
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity(db.nuclide("Sr89"), 1.0);
      double i909 = 0.0;
      for (const SandiaDecay::EnergyRatePair& g : mix.gammas(3600.0, SandiaDecay::NuclideMixture::OrderByEnergy, true))
        if (std::fabs(g.energy - 908.96) < 0.5) i909 += g.numPerSecond;
      i909 /= mix.activity(3600.0, db.nuclide("Sr89"));
      check(near(i909, 9.56E-05, 0.01), "Sr89 chain: 909 keV yield 9.56E-05 per Sr-89 decay (+-1 %), got "
                                            + std::to_string(i909));
    }
  } catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << "\n";
    return 1;
  }

  std::cout << (g_failures ? "FAILED: " : "PASSED: ") << g_failures << " failure(s)\n";
  return g_failures ? 1 : 0;
}
