// Regression tests for the Th, Pa and U K x-ray structure (2026-10 review): decay x-rays of every Th/Pa/U daughter
// and the Pb/Th/Pa/U fluorescence blocks.  Energies are Deslattes et al. 2003 (NIST SRD 128) theory, xraylib for
// K-O and K-P; sub-line rates xraylib (K-L, K-M2/M3, K-N2/N3), EADL2023 (K-M4, K-M5, K-N4,5, K-O, K-P) and xraydb
// (K-L1); each decay keeps its evaluated K-alpha, K'b1 and K'b2 group totals.  LANL micro-calorimeter spectra
// resolve K-N2/K-N3, K-M4/K-M5 and K-O at these energies.  Run once per delivered XML variant.
//
// Usage: xray_k_structure_tests <sandia.decay[.min|.nocoinc.min].xml>

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

// Intensity of the single x-ray within 1 eV of `energy` (-1 if there is not exactly one).
double xray_at(const SandiaDecay::Transition* tr, double energy) {
  double intensity = -1.0;
  int n = 0;
  for (const SandiaDecay::RadParticle& p : tr->products)
    if (p.type == SandiaDecay::XrayParticle && std::fabs(p.energy - energy) < 0.001) {
      intensity = p.intensity;
      ++n;
    }
  return (n == 1) ? intensity : -1.0;
}

double xray_sum(const SandiaDecay::Transition* tr, double lo_keV, double hi_keV) {
  double s = 0.0;
  for (const SandiaDecay::RadParticle& p : tr->products)
    if (p.type == SandiaDecay::XrayParticle && p.energy >= lo_keV && p.energy <= hi_keV)
      s += p.intensity;
  return s;
}

double element_xray_at(const SandiaDecay::Element* el, double energy) {
  for (const SandiaDecay::EnergyIntensityPair& x : el->xrays)
    if (std::fabs(x.energy - energy) < 0.001)
      return x.intensity;
  return -1.0;
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
    const SandiaDecay::Transition* pa4 = find_transition(db, "Pa234", "U234", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* th1 = find_transition(db, "Th231", "Pa231", SandiaDecay::BetaDecay);
    check(u5 && pa4 && th1, file + ": U-235 alpha, Pa-234 and Th-231 beta- transitions present");
    if (!(u5 && pa4 && th1)) return 1;

    // --- Th K x-rays from U-235 alpha: split K'b1 (K-M) and K'b2 (K-N, K-O, K-P) sub-lines ---
    struct Line { double energy, intensity; const char* label; };
    const std::vector<Line> th_k = {
      {104.8172, 0.0057017, "K-M2"}, {105.6019, 0.011118, "K-M3"}, {106.1564, 0.00018525, "K-M4"},
      {106.3149, 0.00021399, "K-M5"}, {108.4812, 0.0014846, "K-N2"}, {108.6809, 0.0030021, "K-N3"},
      {108.9546, 0.00011803, "K-N4,5"}, {109.4535, 0.0010291, "K-O2,3"}, {109.6305, 0.0001689, "K-P2,3"}};
    for (const Line& l : th_k)
      check(near(xray_at(u5, l.energy), l.intensity, 1e-3), std::string("U235: Th ") + l.label + " at "
            + std::to_string(l.energy) + " keV");
    check(xray_at(u5, 108.583) < 0.0, "U235: no merged Th K-N2,3 line at 108.583 keV");

    // --- U K x-rays from Pa-234 beta-: K'b2 split (old merged K-N2,3 at 114.445 keV) ---
    check(near(xray_at(pa4, 114.3292), 0.0060245, 1e-3), "Pa234: U K-N2 at 114.3292 keV");
    check(near(xray_at(pa4, 114.5562), 0.012277, 1e-3), "Pa234: U K-N3 at 114.5562 keV");
    check(near(xray_at(pa4, 115.3897), 0.0042768, 1e-3), "Pa234: U K-O2,3 at 115.3897 keV");
    check(xray_at(pa4, 114.445) < 0.0, "Pa234: no merged U K-N2,3 line at 114.445 keV");
    check(near(xray_sum(pa4, 114.0, 116.0), 0.017540 + 0.006263, 1e-3), "Pa234: U K'b2 total unchanged");
    check(near(xray_sum(pa4, 98.3, 98.5), 0.1946, 1e-3), "Pa234: U K-L3 unchanged");

    // --- Pa K x-rays from Th-231 beta- ---
    check(near(xray_at(th1, 112.3915), 0.00014918, 1e-3), "Th231: Pa K-O2,3 at 112.3915 keV");
    check(near(xray_sum(th1, 111.0, 113.0), 0.00062 + 0.000216, 1e-3), "Th231: Pa K'b2 total unchanged");

    // --- No transition with a Th, Pa or U daughter keeps the old merged K-N2,3 energies ---
    int bad = 0, n_checked = 0;
    for (const SandiaDecay::Nuclide* nuc : db.nuclides()) {
      for (const SandiaDecay::Transition* tr : nuc->decaysToChildren) {
        if (!tr || !tr->child) continue;
        const short z = tr->child->atomicNumber;
        const double old_kn23 = (z == 90) ? 108.583 : (z == 91) ? 111.486 : (z == 92) ? 114.445 : -1.0;
        if (old_kn23 < 0.0) continue;
        ++n_checked;
        for (const SandiaDecay::RadParticle& p : tr->products)
          if (p.type == SandiaDecay::XrayParticle && std::fabs(p.energy - old_kn23) < 0.001) ++bad;
      }
    }
    check(n_checked > 40 && bad == 0, "no Th/Pa/U daughter has an x-ray at the old merged K-N2,3 energy ("
          + std::to_string(bad) + " found)");

    // --- Fluorescence (element) blocks: 12 K lines, relative to K-L3 = 1 ---
    if (db.xmlContainedElementalXRayInfo()) {
      const SandiaDecay::Element* u = db.element("U");
      const SandiaDecay::Element* th = db.element("Th");
      const SandiaDecay::Element* pb = db.element("Pb");
      check(u && th && pb, "U, Th and Pb elements present");
      if (u && th && pb) {
        check(u->xrays.size() == 12, "U: 12 fluorescence K x-rays");
        check(near(element_xray_at(u, 98.4336), 1.0, 1e-6), "U: K-L3 = 1");
        check(near(element_xray_at(u, 94.6531), 0.6248, 1e-3), "U: K-L2 0.6248");
        check(near(element_xray_at(u, 114.5562), 0.061923, 1e-3), "U: K-N3 0.061923");
        check(near(element_xray_at(u, 115.3897), 0.021571, 1e-3), "U: K-O2,3 0.021571");
        check(near(element_xray_at(th, 109.4535), 0.020803, 1e-3), "Th: K-O2,3 0.020803");
        check(near(element_xray_at(pb, 74.9702), 1.0, 1e-6), "Pb: K-L3 at 74.9702 keV");
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << "\n";
    return 1;
  }
  std::cout << (g_failures ? "FAILED" : "PASSED") << " (" << g_failures << " failures)\n";
  return g_failures ? 1 : 0;
}
