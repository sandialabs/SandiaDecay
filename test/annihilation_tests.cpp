// Tests that NuclideMixture adds the 511 keV annihilation entry only when the decay chain has a positron
// emitter, and that, when it does, the entry is present at every age (so returned energies are age-independent).
// Run once per delivered XML variant.
//
// Usage: annihilation_tests <sandia.decay[.min|.nocoinc.min].xml>

#include "SandiaDecay.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int g_failures = 0;

const double sm_annih_energy = 510.998910;

void check(bool value, const std::string& message) {
  std::cout << (value ? "ok: " : "FAIL: ") << message << "\n";
  if (!value) ++g_failures;
}

// Number of entries at exactly the annihilation energy, and their summed rate (or count).
int num_annih(const std::vector<SandiaDecay::EnergyRatePair>& lines, double* rate = nullptr) {
  int n = 0;
  if (rate) *rate = 0.0;
  for (const SandiaDecay::EnergyRatePair& l : lines)
    if (std::fabs(l.energy - sm_annih_energy) < 1.0E-6) {
      ++n;
      if (rate) *rate += l.numPerSecond;
    }
  return n;
}

int num_annih(const std::vector<SandiaDecay::EnergyCountPair>& lines, double* count = nullptr) {
  int n = 0;
  if (count) *count = 0.0;
  for (const SandiaDecay::EnergyCountPair& l : lines)
    if (std::fabs(l.energy - sm_annih_energy) < 1.0E-6) {
      ++n;
      if (count) *count += l.count;
    }
  return n;
}

// Checks every annihilation-aware API for `symbol` (1 Bq at t=0) has `expected_entries` annihilation entries,
// at age 0 and at one half-life; returns the 511 rate per parent decay at one half-life.
double check_nuclide(const SandiaDecay::SandiaDecayDataBase& db, const std::string& symbol, int expected_entries) {
  const SandiaDecay::Nuclide* nuc = db.nuclide(symbol);
  check(nuc != nullptr, symbol + " is in the database");
  if (!nuc) return 0.0;

  using Order = SandiaDecay::NuclideMixture;
  SandiaDecay::NuclideMixture mix;
  mix.addNuclideByActivity(nuc, 1.0);

  const double hl = nuc->halfLife;
  double rate = 0.0;
  for (const double age : {0.0, hl}) {
    const std::string when = symbol + (age == 0.0 ? " at age 0" : " at one half-life");
    const std::vector<SandiaDecay::EnergyRatePair> gams = mix.gammas(age, Order::OrderByEnergy, true);
    const std::vector<SandiaDecay::EnergyRatePair> phots = mix.photons(age, Order::OrderByEnergy);
    check(num_annih(gams, &rate) == expected_entries,
          when + ": gammas(...,true) has " + std::to_string(expected_entries) + " annihilation entries");
    check(num_annih(phots) == expected_entries,
          when + ": photons() has " + std::to_string(expected_entries) + " annihilation entries");
    check(num_annih(mix.gammas(age, Order::OrderByEnergy, false)) == 0,
          when + ": gammas(...,false) has no annihilation entry");

    const double dt = 0.01 * hl;
    check(num_annih(mix.decayGammasInInterval(age, dt, true, Order::OrderByEnergy, 5)) == expected_entries,
          when + ": decayGammasInInterval(...,true) has " + std::to_string(expected_entries)
              + " annihilation entries");
    check(num_annih(mix.decayPhotonsInInterval(age, dt, Order::OrderByEnergy, 5)) == expected_entries,
          when + ": decayPhotonsInInterval() has " + std::to_string(expected_entries) + " annihilation entries");

    // Same energies at every age, and the interval API agrees with the instantaneous one.
    check(gams.size() == mix.gammas(0.0, Order::OrderByEnergy, true).size(),
          when + ": gammas(...,true) size is age-independent");
    check(phots.size() == mix.decayPhotonsInInterval(age, dt, Order::OrderByEnergy, 5).size(),
          when + ": photons() and decayPhotonsInInterval() return the same number of lines");
  }

  return rate / mix.activity(hl, nuc);
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

    // No positron emitter anywhere in the chain: no annihilation entry at all (previously a zero-rate 511
    //  entry was always added, which callers had to filter out).
    for (const char* sym : {"H3", "C14", "P32", "Sr90", "Cs137", "Co60", "Am241", "Po210"})
      check_nuclide(db, sym, 0);

    // Ru106 -> Rh106 emits a real 511.86 keV gamma, but no positrons: that line is kept, and no annihilation
    //  entry is added.
    check_nuclide(db, "Ru106", 0);
    {
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity(db.nuclide("Ru106"), 1.0);
      int n_real = 0;
      for (const SandiaDecay::EnergyRatePair& g :
           mix.gammas(86400.0, SandiaDecay::NuclideMixture::OrderByEnergy, true))
        if (std::fabs(g.energy - 511.86) < 0.1 && g.numPerSecond > 0.0) ++n_real;
      check(n_real == 1, "Ru106 keeps its real 511.86 keV gamma");
    }

    // Positron emitters: one annihilation entry, two photons per positron.
    const double na22 = check_nuclide(db, "Na22", 1);
    check(std::fabs(na22 - 2.0 * 0.903) < 0.02, "Na22: ~1.81 annihilation photons per decay, got "
                                                    + std::to_string(na22));
    check_nuclide(db, "Cu64", 1);

    // Ge68 is pure EC; its Ga68 daughter is the positron emitter.  At age 0 Ga68 has no activity, but the
    //  entry is still present (at zero rate) so the energies do not depend on age.
    {
      check_nuclide(db, "Ge68", 1);
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity(db.nuclide("Ge68"), 1.0);
      double rate0 = -1.0;
      num_annih(mix.gammas(0.0, SandiaDecay::NuclideMixture::OrderByEnergy, true), &rate0);
      check(rate0 == 0.0, "Ge68 at age 0: annihilation entry has zero rate");

      const double age = 10.0 * 3600.0;  // Ga68 (68 min) in equilibrium
      double rate_eq = 0.0;
      num_annih(mix.gammas(age, SandiaDecay::NuclideMixture::OrderByEnergy, true), &rate_eq);
      rate_eq /= mix.activity(age, db.nuclide("Ge68"));
      check(rate_eq > 1.7 && rate_eq < 1.85, "Ge68 in equilibrium: ~1.78 annihilation photons per decay, got "
                                                 + std::to_string(rate_eq));
    }
  } catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << "\n";
    return 1;
  }

  std::cout << (g_failures ? "FAILED: " : "PASSED: ") << g_failures << " failure(s)\n";
  return g_failures ? 1 : 0;
}
