// Regression tests for the reviewed A=234 decay data of the U-238 chain
// (U-238 alpha, Th-234, Pa-234m beta-minus and IT, Pa-234), audited against
// S. Ota, NDS 207, 351 (2026) and JRC Detective uranium spectra; see the
// commit message.  Run once per delivered XML variant;
// coincidence checks are skipped for a file without coincidence data.
//
// Usage: u238_chain_tests <sandia.decay[.min|.nocoinc.min].xml>

#include "SandiaDecay.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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

// Gammas of a transition within 0.05 keV of an energy.
std::vector<const SandiaDecay::RadParticle*> gammas_at(const SandiaDecay::Transition* tr,
                                                        double energy_keV) {
  std::vector<const SandiaDecay::RadParticle*> out;
  if (!tr) return out;
  for (const SandiaDecay::RadParticle& p : tr->products)
    if (p.type == SandiaDecay::GammaParticle && std::fabs(p.energy - energy_keV) <= 0.05)
      out.push_back(&p);
  return out;
}

double sum_intensity(const std::vector<const SandiaDecay::RadParticle*>& gs) {
  double s = 0.0;
  for (const SandiaDecay::RadParticle* g : gs) s += g->intensity;
  return s;
}

// The single gamma at `energy_keV` whose intensity is closest to `intensity`.
const SandiaDecay::RadParticle* gamma_like(const SandiaDecay::Transition* tr, double energy_keV,
                                           double intensity) {
  const SandiaDecay::RadParticle* best = nullptr;
  for (const SandiaDecay::RadParticle* g : gammas_at(tr, energy_keV))
    if (!best || std::fabs(g->intensity - intensity) < std::fabs(best->intensity - intensity))
      best = g;
  return best;
}

// P(partner | gamma) for partner = the product at `partner` in the same transition.
double coincidence(const SandiaDecay::Transition* tr, const SandiaDecay::RadParticle* g,
                   const SandiaDecay::RadParticle* partner) {
  if (!tr || !g || !partner) return -1.0;
  const size_t index = static_cast<size_t>(partner - &tr->products[0]);
  double p = 0.0;
  for (const auto& c : g->coincidences)
    if (c.first == index) p += c.second;
  return p;
}

// Photons per U-238 decay near an energy, for U-238 aged 20 years.
double chain_yield(const SandiaDecay::SandiaDecayDataBase& db, double energy_keV, double half_width) {
  SandiaDecay::NuclideMixture mix;
  const SandiaDecay::Nuclide* u238 = db.nuclide("U238");
  mix.addAgedNuclideByActivity(u238, 1.0E6, 20.0 * 365.25 * 86400.0);
  const double act = mix.activity(0.0, u238);
  double s = 0.0;
  for (const SandiaDecay::EnergyRatePair& e : mix.gammas(0.0, SandiaDecay::NuclideMixture::OrderByEnergy, true))
    if (std::fabs(e.energy - energy_keV) <= half_width) s += e.numPerSecond;
  return s / act;
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
    using SandiaDecay::DecayMode;
    const SandiaDecay::Transition* pa = find_transition(db, "Pa234", "U234", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* pam = find_transition(db, "Pa234m", "U234", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* it = find_transition(db, "Pa234m", "Pa234", SandiaDecay::IsometricTransitionDecay);
    const SandiaDecay::Transition* u8 = find_transition(db, "U238", "Th234", SandiaDecay::AlphaDecay);
    check(pa && pam && it && u8, file + ": A=234 transitions present");
    if (!(pa && pam && it && u8)) return 1;

    bool has_coinc = false;
    for (const SandiaDecay::RadParticle& p : pa->products)
      has_coinc = has_coinc || !p.coincidences.empty();

    // --- Pa-234 photon intensities (ENSDF RI x the library normalization 1.03) ---
    const auto g880 = gammas_at(pa, 880.5);
    check(g880.size() == 2, "Pa234 880.5 keV has two placements");
    check(near(sum_intensity(g880), 0.10403, 2e-3), "Pa234 880.5 keV total 0.10403");
    const SandiaDecay::RadParticle* e2_880 = gamma_like(pa, 880.5, 0.0618);
    const SandiaDecay::RadParticle* e1_880 = gamma_like(pa, 880.5, 0.04223);
    check(e2_880 && e1_880 && e2_880 != e1_880 && near(e2_880->intensity, 0.0618, 2e-3) &&
              near(e1_880->intensity, 0.04223, 2e-3),
          "Pa234 880.5 keV E2 (1023.92) 0.0618 and E1 (1023.787) 0.04223");
    check(near(sum_intensity(gammas_at(pa, 883.24)), 0.09579, 2e-3), "Pa234 883.24 keV 0.09579");
    check(near(sum_intensity(gammas_at(pa, 980.3)), 0.03, 2e-3), "Pa234 980.3 keV doublet total 0.03");
    check(gammas_at(pa, 513.4).size() == 2 && near(sum_intensity(gammas_at(pa, 513.4)), 0.007519, 2e-3),
          "Pa234 513.4 keV: two placements, total 0.007519");
    check(!gammas_at(pa, 1393.9).empty() && gammas_at(pa, 1393.57).empty(), "Pa234 1393.9 keV (not 1393.57)");
    check(near(sum_intensity(gammas_at(pa, 926.72)), 0.09597, 2e-3), "Pa234 926.72 keV 0.09597 (restored)");
    check(near(sum_intensity(gammas_at(pa, 946.0)), 0.1411, 2e-3), "Pa234 946.0 keV 0.1411 (restored)");
    check(gammas_at(pa, 926.0).empty(), "Pa234 926.0 keV not emitted (absent in data; folded into 926.72)");
    check(gammas_at(pa, 331.4).size() == 1 && near(sum_intensity(gammas_at(pa, 331.4)), 0.000721, 2e-3),
          "Pa234 331.4 keV counted once, as in ENSDF");
    // Transition flux must not be stored as photon intensity: I(1+alpha) <= 1.
    for (const SandiaDecay::Transition* tr : {pa, pam})
      for (const SandiaDecay::RadParticle& p : tr->products)
        if (p.type == SandiaDecay::GammaParticle && !p.e0_verified && p.icc_total > 0 && p.icc_total < 1.0e4 &&
            p.energy < 120.0f)
          check(p.intensity * (1.0 + p.icc_total) <= 1.0,
                tr->parent->symbol + " " + std::to_string(p.energy) + " keV: I(1+alpha) <= 1");
    check(near(sum_intensity(gammas_at(pa, 34.3)), 3.399e-5, 1e-2), "Pa234 34.3 keV 3.4e-5 (was transition flux)");
    check(near(sum_intensity(gammas_at(pa, 75.0)), 7.035e-6, 1e-2), "Pa234 75.0 keV 7.0e-6 (was transition flux)");
    check(near(sum_intensity(gammas_at(pam, 43.5)), 2.0e-5, 1e-2), "Pa234m 43.5 keV 2.0e-5 (was transition flux)");

    // --- Line ownership: each decay carries its own lines ---
    check(!gammas_at(pam, 880.9).empty() && !gammas_at(pam, 883.22).empty(), "Pa234m direct 880.9/883.22 keV present");
    // Pa-234m's own 926/946 keV lines, no longer holding 0.0016 x Pa-234's; both owners are scaled alike
    // (x1.071, x1.054) to keep the measured chain totals.  926.61 uses the Pa-234 branching of the shared
    // 926.7 keV level (0.935 x 883.22 keV).
    check(near(sum_intensity(gammas_at(pam, 945.9)), 1.066e-4, 1e-2), "Pa234m 945.9 keV holds only its own line");
    check(near(sum_intensity(gammas_at(pam, 926.61)), 3.546e-5, 1e-2), "Pa234m 926.61 keV holds only its own line");
    check(!gammas_at(pam, 781.75).empty() && gammas_at(pam, 781.37).empty(), "Pa234m 781.75 keV (not 781.37)");
    check(near(sum_intensity(gammas_at(pam, 739.95)), 1.096e-4, 1e-2) && gammas_at(pam, 740.10).empty(),
          "Pa234m 739.95 keV (energy confirmed in spectra) at the ENSDF intensity");
    for (const auto& e0 : {std::make_pair(pa, 810.0), std::make_pair(find_transition(db, "U234", "Th230", SandiaDecay::AlphaDecay), 634.9)}) {
      const auto g = gammas_at(e0.first, e0.second);
      check(g.size() == 1 && g[0]->e0_verified, "pure E0 " + std::to_string(e0.second) + " keV is e0Verified");
    }
    bool u8_high = false;
    for (const SandiaDecay::RadParticle& p : u8->products)
      u8_high = u8_high || (p.type == SandiaDecay::GammaParticle && p.energy > 200.0f);
    check(!u8_high, "U238 alpha carries no U-234 transitions (> 200 keV)");
    for (double e : {994.93, 1124.93, 1548.12, 1591.65, 2001.0})
      check(!gammas_at(pam, e).empty(), "Pa234m carries " + std::to_string(e) + " keV (moved from U238 alpha)");
    check(near(sum_intensity(gammas_at(it, 73.92)), 1.0 / 11.6, 2e-3), "Pa234m IT 73.92 keV = 1/(1+10.6)");

    // --- Chain yields per U-238 decay at 20 y: what RelActCalcAuto sees ---
    check(near(chain_yield(db, 1001.0, 0.2), 8.37e-3, 5e-3), "chain 1001 keV");
    check(near(chain_yield(db, 766.4, 0.1), 2.941e-3, 5e-3), "chain 766.4 keV");
    check(near(chain_yield(db, 258.26, 0.1), 7.288e-4, 5e-3), "chain 258.3 keV");
    check(near(chain_yield(db, 880.7, 0.5), 2.06e-4, 1e-2), "chain 880.5+880.9 keV");
    check(near(chain_yield(db, 883.23, 0.1), 1.886e-4, 1e-2), "chain 883.2 keV");
    check(near(chain_yield(db, 926.65, 0.2), 1.890e-4, 1e-2), "chain 926.6-926.7 keV (measured 1.89e-4)");
    check(near(chain_yield(db, 945.95, 0.2), 3.322e-4, 1e-2), "chain 945.9-946.0 keV (measured 3.32e-4)");
    check(chain_yield(db, 43.5, 0.05) < 1e-4, "chain 43.5 keV is not a 1% phantom line");

    // --- Decay x-rays: K lines within the K-vacancy budget, physical Ka2/Ka1 ---
    // X-rays are selected by energy (the API carries no line assignment).
    auto xray = [](const SandiaDecay::Transition* tr, double lo_keV, double hi_keV) {
      double s = 0.0;
      for (const SandiaDecay::RadParticle& p : tr->products)
        if (p.type == SandiaDecay::XrayParticle && p.energy >= lo_keV && p.energy <= hi_keV) s += p.intensity;
      return s;
    };
    const SandiaDecay::Transition* th = find_transition(db, "Th234", "Pa234m", SandiaDecay::BetaDecay);
    check(near(xray(u8, 93.3, 93.4), 1.119e-5, 1e-2), "U238 alpha Th Ka1 1.119e-5 (was 25x the K-vacancy budget)");
    check(th && near(xray(th, 95.8, 95.9), 1.4e-4, 1e-2) && xray(th, 92.25, 92.32) < xray(th, 95.8, 95.9),
          "Th234 Pa Ka1 1.4e-4 and Ka2 < Ka1");
    check(near(xray(pam, 98.4, 98.5), 1.087e-4, 1e-2), "Pa234m U Ka1 1.087e-4 (was 21x)");
    check(near(xray(it, 5.0, 25.0), 0.3377, 1e-2), "Pa234m IT Pa L x-rays 0.338 per IT");

    // --- Coincidences (P(partner | gamma), following-direction) ---
    if (has_coinc) {
      const SandiaDecay::RadParticle* g99 = gamma_like(pa, 99.86, 1.0);
      const SandiaDecay::RadParticle* g43 = gamma_like(pa, 43.49, 1.0);
      for (double e : {103.77, 738.0}) {   // feed 1023.787 (3-): E1 880.5 only
        const SandiaDecay::RadParticle* f = gamma_like(pa, e, 1.0);
        check(near(coincidence(pa, f, e1_880), 0.2866, 0.05) && coincidence(pa, f, e2_880) < 1e-3,
              "Pa234 " + std::to_string(e) + " keV feeds 1023.787: E1 880.5 at ~0.287, no E2");
      }
      for (double e : {397.7, 472.3}) {    // feed 1023.92 (4+): E2 880.5 only
        const SandiaDecay::RadParticle* f = gamma_like(pa, e, 1.0);
        check(near(coincidence(pa, f, e2_880), 0.815, 0.05) && coincidence(pa, f, e1_880) < 1e-3,
              "Pa234 " + std::to_string(e) + " keV feeds 1023.92: E2 880.5 at ~0.815, no E1");
      }
      for (const SandiaDecay::RadParticle* g : {e1_880, e2_880}) {
        check(near(coincidence(pa, g, g99), 0.0693, 0.03), "Pa234 880.5 keV -> 99.86 keV ~0.069");
        check(near(coincidence(pa, g, g43), 0.0014, 0.05), "Pa234 880.5 keV -> 43.49 keV ~0.0014");
      }
      // The 809.9 keV 0+ level decays mostly by the 810 keV E0, so its feeders see 766.4 keV only ~27%
      // of the time; treating the E0's sentinel ICC as a photon multiplier would drive this to ~0.
      const SandiaDecay::RadParticle* g766 = gamma_like(pa, 766.4, 1.0);
      for (double e : {316.7, 427.4})
        check(near(coincidence(pa, gamma_like(pa, e, 1.0), g766), 0.2737, 0.05),
              "Pa234 " + std::to_string(e) + " keV feeds 809.9: 766.4 keV ~0.27");
      // Final level found though the evaluated energies miss it by 0.63 keV.
      check(near(coincidence(pa, gamma_like(pa, 1354.6, 1.0), gamma_like(pa, 152.71, 1.0)), 0.318, 0.05),
            "Pa234 1354.6 keV -> 152.71 keV ~0.32");
      const SandiaDecay::RadParticle* m258 = gamma_like(pam, 258.26, 1.0);
      const SandiaDecay::RadParticle* m743 = gamma_like(pam, 742.77, 1.0);
      check(near(coincidence(pam, m258, m743), 0.619, 0.05), "Pa234m 258.26 keV -> 742.77 keV ~0.62");
      // Pa-234m's 810.4 keV level: the 811 keV E0 to ground takes ~70% of its decays.
      check(near(coincidence(pam, gamma_like(pam, 691.08, 1.0), gamma_like(pam, 766.37, 1.0)), 0.290, 0.05),
            "Pa234m 691.08 keV -> 766.37 keV ~0.29");
    }
  } catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << "\n";
    return 1;
  }
  std::cout << (g_failures ? "FAILED " : "PASSED ") << file << " (" << g_failures << " failures)\n";
  return g_failures ? 1 : 0;
}
