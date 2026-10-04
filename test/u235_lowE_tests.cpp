// Regression tests for the reviewed low-energy U-235 / U-238 data (2026-10): line energies and intensities checked
// against ENSDF (Singh, Tuli and Browne, NDS 185, 560 (2022)), DDEP and LANL micro-calorimeter and IDB planar HPGe
// spectra of U standards; see the commit message.  Run once per delivered XML variant.
//
// Usage: u235_lowE_tests <sandia.decay[.min|.nocoinc.min].xml>

#include "SandiaDecay.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
    const SandiaDecay::SandiaDecayDataBase db(file);
    const SandiaDecay::Transition* u5 = find_transition(db, "U235", "Th231", SandiaDecay::AlphaDecay);
    const SandiaDecay::Transition* u8 = find_transition(db, "U238", "Th234", SandiaDecay::AlphaDecay);
    check(u5 && u8, file + ": U-235 and U-238 alpha transitions present");
    if (!(u5 && u8)) return 1;
    using SandiaDecay::GammaParticle;
    using SandiaDecay::XrayParticle;

    // --- U-235 gamma energies (ENSDF; micro-cal 109.1781(45) and 163.341(11) keV) ---
    int n = 0;
    check(sum_in(u5, GammaParticle, 108.94, 109.04, &n) == 0.0 && n == 0, "U235: no gamma at 108.99 keV");
    check(near(sum_in(u5, GammaParticle, 109.170, 109.182, &n), 0.0154, 1e-6) && n == 1,
          "U235: 109.176 keV gamma, 1.54 per 100");
    check(sum_in(u5, GammaParticle, 163.37, 163.39, &n) == 0.0 && n == 0, "U235: no gamma at 163.38 keV");
    check(near(sum_in(u5, GammaParticle, 163.352, 163.362, &n), 0.0508, 1e-6) && n == 1,
          "U235: 163.357 keV gamma, 5.08 per 100");
    // --- U-235 gamma intensities (ENSDF) ---
    check(near(sum_in(u5, GammaParticle, 172.95, 173.05), 1.144e-4, 1e-3), "U235: 173.0 keV 0.01144 per 100");
    check(near(sum_in(u5, GammaParticle, 198.92, 198.935, &n), 4.233e-4, 1e-3) && n == 1,
          "U235: 198.927 keV 0.04233 per 100");

    // --- Th K x-rays from U-235 alpha: ENSDF group totals, K'b1 (K-M) and K'b2 (K-N, K-O) ---
    check(near(sum_in(u5, XrayParticle, 93.30, 93.40), 0.04839, 1e-3), "U235: Th Ka1 4.839 per 100");
    check(near(sum_in(u5, XrayParticle, 89.90, 90.00), 0.029914, 1e-3), "U235: Th Ka2 2.991 per 100");
    check(near(sum_in(u5, XrayParticle, 104.5, 106.6), 0.017219, 1e-3), "U235: Th K'b1 (K-M) 1.722 per 100");
    check(near(sum_in(u5, XrayParticle, 108.3, 110.0), 0.0058027, 1e-3), "U235: Th K'b2 (K-N, K-O) 0.580 per 100");

    // --- Th K x-rays from U-238 alpha: same split (ENSDF K'b1/Ka1 0.356, K'b2/Ka1 0.120) ---
    check(near(sum_in(u8, XrayParticle, 104.5, 106.6), 3.9806e-6, 1e-3), "U238: Th K'b1 3.98e-6 per decay");
    check(near(sum_in(u8, XrayParticle, 108.3, 110.0), 1.3415e-6, 1e-3), "U238: Th K'b2 1.34e-6 per decay");
  } catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << "\n";
    return 1;
  }
  std::cout << (g_failures ? "FAILED" : "PASSED") << " (" << g_failures << " failures)\n";
  return g_failures ? 1 : 0;
}
