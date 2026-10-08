#include "SandiaDecay.h"
#include "rapidxml/rapidxml.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <set>
#include <vector>

namespace {

void check(bool value, const std::string& message) {
  if (!value)
    throw std::runtime_error(message);
  std::cout << "ok: " << message << "\n";
}

struct TestParticle : SandiaDecay::RadParticle {
  explicit TestParticle(const rapidxml::xml_node<char>* node)
      : SandiaDecay::RadParticle(node) {}
};

SandiaDecay::RadParticle parse(char* text) {
  rapidxml::xml_document<char> doc;
  doc.parse<0>(text);
  return TestParticle(doc.first_node());
}

SandiaDecay::Transition parse_transition(char* text) {
  rapidxml::xml_document<char> doc;
  doc.parse<0>(text);
  const std::vector<const SandiaDecay::Nuclide*> nuclides;
  return SandiaDecay::Transition(doc.first_node(), nuclides);
}

bool transition_parse_throws(char* text) {
  try {
    (void)parse_transition(text);
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

const SandiaDecay::Transition* find_transition(
    const SandiaDecay::SandiaDecayDataBase& db, const std::string& parent,
    const std::string& child, SandiaDecay::DecayMode mode) {
  const SandiaDecay::Nuclide* p = db.nuclide(parent);
  if (!p) return nullptr;
  for (const SandiaDecay::Transition* tr : p->decaysToChildren)
    if (tr->child && tr->child->symbol == child && tr->mode == mode)
      return tr;
  return nullptr;
}

const SandiaDecay::RadParticle* find_gamma(
    const SandiaDecay::SandiaDecayDataBase& db, const std::string& parent,
    const std::string& child, double energy_keV) {
  const SandiaDecay::Nuclide* p = db.nuclide(parent);
  if (!p) return nullptr;
  for (const SandiaDecay::Transition* tr : p->decaysToChildren) {
    if (!tr->child || tr->child->symbol != child) continue;
    for (const SandiaDecay::RadParticle& particle : tr->products)
      if (particle.type == SandiaDecay::GammaParticle &&
          std::fabs(particle.energy - energy_keV) <= 0.05)
        return &particle;
  }
  return nullptr;
}

int verified_e0_count(const SandiaDecay::SandiaDecayDataBase& db) {
  int count = 0;
  for (const SandiaDecay::Nuclide* parent : db.nuclides())
    for (const SandiaDecay::Transition* tr : parent->decaysToChildren)
      for (const SandiaDecay::RadParticle& particle : tr->products)
        if (particle.e0_verified) ++count;
  return count;
}

double bounded_single_gamma_yield(const SandiaDecay::RadParticle& gamma) {
  const double alpha = std::max(0.0, static_cast<double>(gamma.icc_total));
  const bool sentinel_e0 = alpha >= 1.0e4 && gamma.energy > 100.0f;
  if (gamma.e0_verified || sentinel_e0) return 0.0;
  return std::min(std::max(0.0, static_cast<double>(gamma.intensity)),
                  1.0 / (1.0 + alpha));
}

double rate_near(const std::vector<SandiaDecay::EnergyRatePair>& values,
                 double energy_keV) {
  double result = 0.0;
  for (const SandiaDecay::EnergyRatePair& value : values)
    if (std::fabs(value.energy - energy_keV) < 0.05)
      result += value.numPerSecond;
  return result;
}

double count_near(const std::vector<SandiaDecay::EnergyCountPair>& values,
                  double energy_keV) {
  double result = 0.0;
  for (const SandiaDecay::EnergyCountPair& value : values)
    if (std::fabs(value.energy - energy_keV) < 0.05)
      result += value.count;
  return result;
}

void check_database_provenance(const std::string& filename) {
  const SandiaDecay::SandiaDecayDataBase db(filename);
  const SandiaDecay::RadParticle* n16 = find_gamma(db, "N16", "O16", 6048.2);
  const SandiaDecay::RadParticle* ra232 = find_gamma(db, "Ra232", "Ac232", 373.3);
  const SandiaDecay::RadParticle* sr78 = find_gamma(db, "Sr78", "Rb78", 189.8);
  const SandiaDecay::RadParticle* nb98 = find_gamma(db, "Nb98", "Mo98", 734.6);
  const SandiaDecay::RadParticle* u235 = find_gamma(db, "U235", "Th231", 19.59);
  const SandiaDecay::RadParticle* eu152_e0 =
      find_gamma(db, "Eu152", "Gd152", 615.4);
  const SandiaDecay::RadParticle* eu152m_615 =
      find_gamma(db, "Eu152m", "Gd152", 615.44);
  check(n16 != nullptr, filename + " contains N16->O16 6048.2-keV record");
  check(ra232 != nullptr, filename + " contains Ra232->Ac232 373.3-keV record");
  check(sr78 != nullptr, filename + " contains Sr78->Rb78 189.8-keV record");
  check(nb98 != nullptr, filename + " contains Nb98->Mo98 734.6-keV record");
  check(u235 != nullptr, filename + " contains U235->Th231 19.59-keV record");
  check(eu152_e0 != nullptr, filename + " contains Eu152->Gd152 615.4-keV record");
  check(n16 && n16->e0_verified,
        filename + " retains verified N16 E0");
  check(ra232 && !ra232->e0_verified,
        filename + " does not promote advisory Ac232 mult=1");
  check(sr78 && !sr78->e0_verified,
        filename + " does not promote advisory Rb78 mult=1");
  check(ra232 && ra232->multipole == 1 && ra232->icc_total == 0.0f &&
            std::fabs(bounded_single_gamma_yield(*ra232) - 0.62) < 1e-5,
        filename + " keeps Ac232 advisory mult=1 as a photon");
  check(sr78 && sr78->multipole == 1 && sr78->icc_total == 0.0f &&
            std::fabs(bounded_single_gamma_yield(*sr78) - 0.16) < 1e-5,
        filename + " keeps Rb78 advisory mult=1 as a photon");
  check(n16 && bounded_single_gamma_yield(*n16) == 0.0,
        filename + " verified N16 E0 has zero single-gamma yield");
  const bool eu152_topology_ok = eu152_e0 &&
      (filename.find("nocoinc") != std::string::npos
           ? (eu152_e0->from_level == -1 && eu152_e0->to_level == -1)
           : (eu152_e0->from_level == 2 && eu152_e0->to_level == 0));
  check(eu152_e0 && eu152_e0->e0_verified && eu152_e0->multipole == 1 &&
            eu152_topology_ok &&
            eu152_e0->icc_total > 1.0e20f &&
            bounded_single_gamma_yield(*eu152_e0) == 0.0,
        filename + " retains verified Eu152 615.4-keV E0 physics/topology");
  check(eu152m_615 && !eu152m_615->e0_verified,
        filename + " does not apply Eu152 E0 provenance to Eu152m");
  // Reviewed U-238-chain pure E0 transitions (ENSDF NDS 207, 2026 and NDS 197, 2024).
  const bool nocoinc = filename.find("nocoinc") != std::string::npos;
  for (const auto& e0 : {std::make_tuple("Pa234", "U234", 810.0f), std::make_tuple("U234", "Th230", 634.9f)}) {
    const SandiaDecay::RadParticle* g =
        find_gamma(db, std::get<0>(e0), std::get<1>(e0), std::get<2>(e0));
    const bool topology_ok = g && (nocoinc ? (g->from_level == -1 && g->to_level == -1)
                                           : (g->from_level == 7 && g->to_level == 0));
    check(g && g->e0_verified && topology_ok && bounded_single_gamma_yield(*g) == 0.0,
          filename + " marks reviewed " + std::get<0>(e0) + "->" + std::get<1>(e0) + " " +
              std::to_string(std::get<2>(e0)) + "-keV pure E0");
  }
  check(nb98 && nb98->icc_total > 1.0e20f && nb98->icc_K > 1.0e20f &&
            nb98->multipole == 1 && bounded_single_gamma_yield(*nb98) == 0.0,
        filename + " Nb98 sentinel E0 has zero single-gamma yield");
  check(u235 && u235->icc_total > 114.0f && u235->icc_K == 0.0f &&
            u235->multipole == 3 &&
            std::fabs(u235->intensity - 0.00533f) < 1e-8f,
        filename + " applies reviewed U235 19.59-keV intensity");   // 185.7-keV level balance (2026-10-07)
  if (u235) {
    const double occurrence = u235->intensity * (1.0 + u235->icc_total);
    check(std::fabs(occurrence - 0.616681) < 2e-6 && occurrence < 1.0,
          filename + " U235 19.59-keV transition occurrence is physical");
    check(std::fabs(bounded_single_gamma_yield(*u235) - u235->intensity) < 1e-9,
          filename + " U235 direct and fallback photon yields close");
  }

  const SandiaDecay::Transition* eu152_sm = find_transition(
      db, "Eu152", "Sm152", SandiaDecay::ElectronCaptureDecay);
  check(eu152_sm != nullptr, filename + " contains Eu152->Sm152 EC branch");
  if (eu152_sm) {
    double weak_parent_mass = 0.0;
    double ec_k_parent_mass = 0.0;
    double positron_parent_mass = 0.0;
    double sm_xray_parent_mass = 0.0;
    for (const SandiaDecay::RadParticle& particle : eu152_sm->products) {
      if (particle.type == SandiaDecay::CaptureElectronParticle) {
        weak_parent_mass += eu152_sm->branchRatio * particle.intensity;
        ec_k_parent_mass += eu152_sm->branchRatio * particle.intensity *
                            particle.ec_K;
      } else if (particle.type == SandiaDecay::PositronParticle) {
        weak_parent_mass += eu152_sm->branchRatio * particle.intensity;
        positron_parent_mass += eu152_sm->branchRatio * particle.intensity;
      } else if (particle.type == SandiaDecay::XrayParticle) {
        sm_xray_parent_mass += eu152_sm->branchRatio * particle.intensity;
      }
    }
    check(std::fabs(weak_parent_mass - 0.720282) < 2e-7,
          filename + " Eu152 weak products are multiplied by branch ratio once");
    const double ec_k_tolerance =
        filename.find(".min.xml") == std::string::npos ? 2e-7 : 3e-5;
    check(std::fabs(ec_k_parent_mass - 0.600803823) < ec_k_tolerance,
          filename + " Eu152 EC K-vacancy source preserves parent-absolute mass");
    check(std::fabs(positron_parent_mass - 0.000138) < 2e-8,
          filename + " Eu152 positron source preserves parent-absolute mass");
    check(weak_parent_mass < eu152_sm->branchRatio,
          filename + " Eu152 keeps the genuine unrepresented weak remainder");
    check(std::fabs(sm_xray_parent_mass - 0.7379236) < 2e-7,
          filename + " Eu152 explicit Sm K x-rays are multiplied by branch ratio once");

    // Exercise the public API that originally exposed the duplicate branch
    // factor, not merely the transition's stored values.
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(db.nuclide("Eu152"), 1.0);
    const std::vector<SandiaDecay::EnergyRatePair> xrays =
        mixture.xrays(0.0, SandiaDecay::NuclideMixture::OrderByEnergy);
    const double sm_energies[] =
        {39.097, 39.522, 40.118, 45.294, 45.414, 45.741, 46.578, 46.705};
    double api_sm_xray_rate = 0.0;
    for (const SandiaDecay::EnergyRatePair& xray : xrays)
      for (double energy : sm_energies)
        if (std::fabs(xray.energy - energy) < 0.001)
          api_sm_xray_rate += xray.numPerSecond;
    check(std::fabs(api_sm_xray_rate - 0.7379236) < 2e-7,
          filename + " NuclideMixture returns absolute Eu152 Sm K-x-ray yield");

    // The raw E0 record must remain available to cascade consumers, but every
    // public emitted-photon wrapper shares decayParticle() and must suppress it.
    check(rate_near(mixture.decayParticle(
                        0.0, SandiaDecay::GammaParticle,
                        SandiaDecay::NuclideMixture::OrderByEnergy),
                    615.4) == 0.0,
          filename + " generic particle API suppresses Eu152 verified E0 photon");
    check(rate_near(mixture.gammas(
                        0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false),
                    615.4) == 0.0,
          filename + " gamma API suppresses Eu152 verified E0 photon");
    check(rate_near(mixture.photons(
                        0.0, SandiaDecay::NuclideMixture::OrderByEnergy),
                    615.4) == 0.0,
          filename + " photon API suppresses Eu152 verified E0 photon");
    check(count_near(mixture.decayParticlesInInterval(
                         0.0, 1.0, SandiaDecay::GammaParticle,
                         SandiaDecay::NuclideMixture::OrderByEnergy),
                     615.4) == 0.0,
          filename + " generic interval API suppresses Eu152 verified E0 photon");
    check(count_near(mixture.decayGammasInInterval(
                         0.0, 1.0, false,
                         SandiaDecay::NuclideMixture::OrderByEnergy),
                     615.4) == 0.0,
          filename + " gamma interval API suppresses Eu152 verified E0 photon");
    check(count_near(mixture.decayPhotonsInInterval(
                         0.0, 1.0,
                         SandiaDecay::NuclideMixture::OrderByEnergy),
                     615.4) == 0.0,
          filename + " photon interval API suppresses Eu152 verified E0 photon");
  }

  SandiaDecay::NuclideMixture n16_mix;
  n16_mix.addNuclideByActivity(db.nuclide("N16"), 1.0);
  check(rate_near(n16_mix.gammas(
                      0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false),
                  6048.2) == 0.0,
        filename + " gamma API suppresses N16 verified E0 photon");

  // Multipolarity code 1 alone is advisory.  These measured lines have no
  // verified-E0 provenance and must remain in the public gamma stream.
  SandiaDecay::NuclideMixture ra232_mix;
  ra232_mix.addNuclideByActivity(db.nuclide("Ra232"), 1.0);
  check(rate_near(ra232_mix.gammas(
                      0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false),
                  373.3) > 0.0,
        filename + " advisory Ra232 mult=1 line remains a photon");
  SandiaDecay::NuclideMixture sr78_mix;
  sr78_mix.addNuclideByActivity(db.nuclide("Sr78"), 1.0);
  check(rate_near(sr78_mix.gammas(
                      0.0, SandiaDecay::NuclideMixture::OrderByEnergy, false),
                  189.8) > 0.0,
        filename + " advisory Sr78 mult=1 line remains a photon");

  const SandiaDecay::Transition* u235_alpha = find_transition(
      db, "U235", "Th231", SandiaDecay::AlphaDecay);
  check(u235_alpha != nullptr, filename + " contains U235->Th231 alpha branch");
  if (u235_alpha) {
    const double alpha_energy[] = {
        3976.0, 4080.6, 4152.5, 4216.0, 4219.0, 4248.0, 4266.2,
        4273.0, 4282.6, 4323.0, 4323.2, 4357.0, 4364.6, 4377.0,
        4395.7, 4415.1, 4438.3, 4502.9, 4556.5, 4597.7};
    const double alpha_probability[] = {
        0.00006915006016, 0.0001718901495, 0.002933902552,
        0.06025905243, 0.001215101057, 0.0006915006016,
        0.002272101977, 0.0003259902836, 0.001086600945,
        0.001649701435, 0.03536503077, 0.002074501805,
        0.1871001628, 0.003457503008, 0.5708804967,
        0.03052502656, 0.002420202106, 0.012645011,
        0.03773603283, 0.047121041};
    const double alpha_hindrance[] = {
        12.0, 41.8, 11.5, 1.94, 117.0, 320.0, 150.0, 1200.0,
        410.0, 580.0, 26.0, 870.0, 10.4, 740.0, 6.02, 160.0,
        3100.0, 1850.0, 1570.0, 2550.0};
    const int alpha_level[] = {
        49, 37, 30, 26, 22, 21, 18, 17, 16, 13,
        14, 10, 8, 6, 5, 4, 3, 2, 1, 0};
    std::size_t alpha_index = 0;
    double alpha_sum = 0.0;
    bool omitted_rows_absent = true;
    for (const SandiaDecay::RadParticle& particle : u235_alpha->products) {
      if (particle.type != SandiaDecay::AlphaParticle) continue;
      omitted_rows_absent &= std::fabs(particle.energy - 4042.0) > 0.1 &&
                             std::fabs(particle.energy - 4212.0) > 0.1 &&
                             std::fabs(particle.energy - 4296.0) > 0.1;
      check(alpha_index < 20 &&
                std::fabs(particle.energy - alpha_energy[alpha_index]) < 0.001 &&
                std::fabs(particle.intensity - alpha_probability[alpha_index]) < 2e-7 &&
                std::fabs(particle.hindrance - alpha_hindrance[alpha_index]) <
                    std::max(1e-5, alpha_hindrance[alpha_index] * 1e-6),
            filename + " U235 projected alpha energy/probability/hindrance entry " +
                std::to_string(alpha_index));
      alpha_sum += particle.intensity;
      ++alpha_index;
    }
    check(alpha_index == 20 && std::fabs(alpha_sum - 1.0) < 2e-7,
          filename + " U235 projected alpha products contain 20 entries and sum to one");
    check(omitted_rows_absent,
          filename + " U235 omits documented unmapped 4042/4212/4296-keV rows");

    const int feed_level[] = {
        0, 1, 2, 3, 4, 5, 6, 8, 10, 13,
        14, 16, 17, 18, 21, 22, 26, 30, 37, 49};
    const double feed_probability[] = {
        0.047121041, 0.03773603283, 0.012645011,
        0.002420202106, 0.03052502656, 0.5708804967,
        0.003457503008, 0.1871001628, 0.002074501805,
        0.001649701435, 0.03536503077, 0.001086600945,
        0.0003259902836, 0.002272101977, 0.0006915006016,
        0.001215101057, 0.06025905243, 0.002933902552,
        0.0001718901495, 0.00006915006016};
    const double feed_q[] = {
        4678.1, 4636.148, 4581.971, 4515.99, 4492.384,
        4472.79, 4456.703, 4441.201, 4437.221, 4402.675,
        4400.5, 4361.019, 4353.187, 4347.5, 4326.589,
        4300.527, 4290.272, 4225.908, 4147.87, 4044.054};
    check(u235_alpha->level_feeds.size() == 20,
          filename + " U235 projected feed law contains 20 entries");
    double feed_sum = 0.0;
    for (std::size_t i = 0; i < u235_alpha->level_feeds.size() && i < 20; ++i) {
      const SandiaDecay::LevelFeed& feed = u235_alpha->level_feeds[i];
      check(feed.level == feed_level[i] &&
                std::fabs(feed.probability - feed_probability[i]) < 2e-7 &&
                std::fabs(feed.q_keV - feed_q[i]) < 0.001 &&
                feed.forbiddenness == SandiaDecay::NoForbiddenness,
            filename + " U235 exact projected level/Q/feed entry " +
                std::to_string(i));
      feed_sum += feed.probability;
    }
    check(std::fabs(feed_sum - 1.0) < 2e-7,
          filename + " U235 projected level feeds sum to one");
    for (std::size_t i = 0; i < 20; ++i) {
      const auto feed_it = std::find_if(
          u235_alpha->level_feeds.begin(), u235_alpha->level_feeds.end(),
          [&](const SandiaDecay::LevelFeed& feed) {
            return feed.level == alpha_level[i];
          });
      check(feed_it != u235_alpha->level_feeds.end() &&
                std::fabs(feed_it->probability - alpha_probability[i]) < 2e-7,
            filename + " U235 alpha product probability equals mapped feed " +
                std::to_string(i));
    }
  }

  const SandiaDecay::Transition* th234 = find_transition(
      db, "Th234", "Pa234m", SandiaDecay::BetaDecay);
  check(th234 != nullptr, filename + " contains Th234->Pa234m beta branch");
  if (th234) {
    // 92.38/92.80 keV moved to 92.362/92.787 keV by the 2026-10-05 data review (micro-calorimeter against the
    //  Th Ka1 / U Ka2 x-rays)
    const double energies[] = {62.86, 73.92, 83.3, 87.02, 92.362, 92.787,
                               103.59, 103.71, 108.0, 112.8, 132.9, 184.8, 63.29};   // 103.35 -> 103.59 (2026-10-06)
    const double old_i[] = {0.0002112, 0.0, 0.000792, 0.0001914, 0.028116,
                            0.027715, 0.00004158, 0.0, 0.0001056, 0.002772,
                            0.0, 0.000132, 0.0485};
    // 92.38, 92.80 and 112.8 keV were then raised x1.16 by the 2026-10-05 data review (measured spectra of certified
    //  U3O8 standards); the rest keep the 25/33 (ENSDF NR) correction.
    const auto expected_i = [&](std::size_t i) {
      if (energies[i] == 92.362) return 0.02471;
      if (energies[i] == 92.787) return 0.02436;
      if (energies[i] == 112.8) return 0.002436;
      if (energies[i] == 108.0) return 1.0e-5;   // 2026-10-06 data review (micro-calorimeter)
      return old_i[i] * 25.0 / 33.0;
    };
    std::size_t gamma_index = 0;
    std::vector<double> beta_intensities;
    for (const SandiaDecay::RadParticle& particle : th234->products) {
      if (particle.type == SandiaDecay::GammaParticle) {
        check(gamma_index < 13 &&
                  std::fabs(particle.energy - energies[gamma_index]) < 0.001 &&
                  std::fabs(particle.intensity - expected_i(gamma_index)) < 2e-8,
              filename + " retains exact corrected Th234 gamma vector entry " +
                  std::to_string(gamma_index));
        ++gamma_index;
      } else if (particle.type == SandiaDecay::BetaParticle) {
        beta_intensities.push_back(particle.intensity);
      }
    }
    check(gamma_index == 13, filename + " has all 13 corrected Th234 gammas");
    const double expected_beta[] = {0.029, 0.00021, 0.192, 0.703};
    check(beta_intensities.size() == 4,
          filename + " retains four Th234 beta products");
    for (std::size_t i = 0; i < beta_intensities.size() && i < 4; ++i)
      check(std::fabs(beta_intensities[i] - expected_beta[i]) < 1e-8,
            filename + " leaves Th234 beta intensity unchanged " +
                std::to_string(i));
  }

  const SandiaDecay::Transition* u238 = find_transition(
      db, "U238", "Th234", SandiaDecay::AlphaDecay);
  check(u238 != nullptr, filename + " contains U238->Th234 alpha branch");
  if (u238) {
    const double alpha_energy[] = {4038.0, 4151.0, 4198.0};
    const double alpha_probability[] = {0.00068, 0.2292, 0.77012};
    const double alpha_hindrance[] = {45.0, 1.4, 1.0};
    std::size_t alpha_index = 0;
    double alpha_sum = 0.0;
    for (const SandiaDecay::RadParticle& particle : u238->products) {
      if (particle.type != SandiaDecay::AlphaParticle) continue;
      check(alpha_index < 3 &&
                std::fabs(particle.energy - alpha_energy[alpha_index]) < 1e-6 &&
                std::fabs(particle.intensity - alpha_probability[alpha_index]) < 1e-7 &&
                std::fabs(particle.hindrance - alpha_hindrance[alpha_index]) < 1e-6,
            filename + " U238 reviewed alpha energy/intensity/hindrance entry " +
                std::to_string(alpha_index));
      alpha_sum += particle.intensity;
      ++alpha_index;
    }
    check(alpha_index == 3 && std::fabs(alpha_sum - 1.0) < 1e-7,
          filename + " U238 corrected alpha products sum to one");

    const double feed_probability[] = {0.77012, 0.2292, 0.00068};
    const double feed_q[] = {4269.9, 4220.35, 4106.85};
    check(u238->level_feeds.size() == 3,
          filename + " U238 retains exact three-level alpha feed law");
    double feed_sum = 0.0;
    for (std::size_t i = 0; i < u238->level_feeds.size() && i < 3; ++i) {
      const SandiaDecay::LevelFeed& feed = u238->level_feeds[i];
      check(feed.level == static_cast<int>(i) &&
                std::fabs(feed.probability - feed_probability[i]) < 1e-7 &&
                std::fabs(feed.q_keV - feed_q[i]) < 0.001 &&
                feed.forbiddenness == SandiaDecay::NoForbiddenness,
            filename + " U238 reviewed feed probability/Q/topology entry " +
                std::to_string(i));
      feed_sum += feed.probability;
      check(std::fabs(feed.probability - alpha_probability[2 - i]) < 1e-7,
            filename + " U238 alpha product and daughter feed agree " +
                std::to_string(i));
    }
    check(std::fabs(feed_sum - 1.0) < 1e-7,
          filename + " U238 corrected level-feed probabilities sum to one");
  }

  bool topology_well_formed = true;
  bool topology_injective = true;
  for (const SandiaDecay::Nuclide* parent : db.nuclides()) {
    for (const SandiaDecay::Transition* tr : parent->decaysToChildren) {
      std::set<std::pair<int, int>> claimed_edges;
      for (const SandiaDecay::RadParticle& particle : tr->products) {
        if (particle.type != SandiaDecay::GammaParticle) continue;
        const bool has_from = particle.from_level >= 0;
        const bool has_to = particle.to_level >= 0;
        if (has_from != has_to) topology_well_formed = false;
        if (has_from && !claimed_edges.insert(
                            {particle.from_level, particle.to_level}).second)
          topology_injective = false;
      }
    }
  }
  check(topology_well_formed,
        filename + " every matched gamma has both level endpoints");
  check(topology_injective,
        filename + " accepted gamma-to-level-edge mapping is injective per transition");

  const SandiaDecay::Transition* pa234 = find_transition(
      db, "Pa234", "U234", SandiaDecay::BetaDecay);
  check(pa234 != nullptr, filename + " contains Pa234->U234 beta branch");
  if (pa234) {
    int n880 = 0, n980 = 0;
    bool pair880_unmatched = true, pair980_unmatched = true;
    for (const SandiaDecay::RadParticle& particle : pa234->products) {
      if (particle.type != SandiaDecay::GammaParticle) continue;
      if (std::fabs(particle.energy - 880.5) < 0.01) {
        ++n880;
        pair880_unmatched &= particle.from_level < 0 && particle.to_level < 0;
      }
      if (std::fabs(particle.energy - 980.3) < 0.01) {
        ++n980;
        pair980_unmatched &= particle.from_level < 0 && particle.to_level < 0;
      }
    }
    check(n880 == 2 && pair880_unmatched,
          filename + " keeps both unequal-intensity Pa234 880.5-keV records unmatched");
    check(n980 == 2 && pair980_unmatched,
          filename + " keeps both unequal-intensity Pa234 980.3-keV records unmatched");
  }

  long feed_transitions = 0;
  bool feeds_unique_and_normalized = true;
  for (const SandiaDecay::Nuclide* parent : db.nuclides()) {
    for (const SandiaDecay::Transition* tr : parent->decaysToChildren) {
      if (tr->level_feeds.empty()) continue;
      ++feed_transitions;
      double sum = 0.0;
      std::set<int> levels;
      for (const SandiaDecay::LevelFeed& feed : tr->level_feeds) {
        sum += feed.probability;
        levels.insert(feed.level);
      }
      feeds_unique_and_normalized &=
          levels.size() == tr->level_feeds.size() && std::fabs(sum - 1.0) < 2e-4;
    }
  }
  check(feeds_unique_and_normalized,
        filename + " every level-feed law is unique and normalized");
  check(feed_transitions > 100,
        filename + " retains database-wide direct level-feed laws");
  {
    const SandiaDecay::Transition* ufeed = find_transition(
        db, "U235", "Th231", SandiaDecay::AlphaDecay);
    const SandiaDecay::Transition* tfeed = find_transition(
        db, "Th234", "Pa234m", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* pmfeed = find_transition(
        db, "Pa234m", "U234", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* pfeed = find_transition(
        db, "Pa234", "U234", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* th231feed = find_transition(
        db, "Th231", "Pa231", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* pbfeed = find_transition(
        db, "Pb214", "Bi214", SandiaDecay::BetaDecay);
    const SandiaDecay::Transition* bifeed = find_transition(
        db, "Bi214", "Po214", SandiaDecay::BetaDecay);
    check(ufeed && ufeed->level_feeds.size() == 20,
          filename + " U235 has complete 20-level alpha feed law");
    check(tfeed && tfeed->level_feeds.size() == 4,
          filename + " Th234 has complete four-level beta feed law");
    check(pmfeed && pmfeed->level_feeds.size() == 23,
          filename + " Pa234m has complete 23-level beta feed law");
    check(pfeed && pfeed->level_feeds.size() == 55,
          filename + " Pa234 has complete 55-level beta feed law");
    check(th231feed && th231feed->level_feeds.size() == 14,
          filename + " Th231 has complete 14-level beta feed law");
    check(pbfeed && pbfeed->level_feeds.size() == 6,
          filename + " Pb214 has complete six-level beta feed law");
    check(bifeed && bifeed->level_feeds.size() == 70,
          filename + " Bi214 has complete 70-level beta feed law");
    if (ufeed && ufeed->level_feeds.size() == 20) {
      check(ufeed->level_feeds[0].level == 0 &&
                ufeed->level_feeds[0].q_keV > 4678.0f &&
                std::fabs(ufeed->level_feeds[0].probability - 0.04712104f) < 2e-7f,
            filename + " U235 ground feed probability and total alpha Q");
      check(ufeed->level_feeds[5].level == 5 &&
                std::fabs(ufeed->level_feeds[5].probability - 0.5708805f) < 2e-7f,
            filename + " U235 dominant 205.31-keV feed probability");
    }
    if (tfeed && tfeed->level_feeds.size() == 4)
      check(tfeed->level_feeds[0].level == 2 &&
                std::fabs(tfeed->level_feeds[0].q_keV - 200.08f) < 0.001f &&
                std::fabs(tfeed->level_feeds[0].probability - 0.7807808f) < 2e-7f,
            filename + " Th234 +X feed maps to level 2 with beta endpoint");
  }
  check(verified_e0_count(db) == 4,
        filename + " has the expected verified-E0 count");
}

}  // namespace

int main() {
  try {
    char positron_long[] =
        "<positron energy=\"901.55\" intensity=\"0.22768\" lvl=\"1\"/>";
    const SandiaDecay::RadParticle pl = parse(positron_long);
    check(pl.type == SandiaDecay::PositronParticle,
          "long positron particle type");
    check(pl.fed_level == 1, "long positron lvl parsed");

    char positron_short[] = "<p e=\"901.55\" i=\"0.22768\" lv=\"7\"/>";
    const SandiaDecay::RadParticle ps = parse(positron_short);
    check(ps.type == SandiaDecay::PositronParticle,
          "short positron particle type");
    check(ps.fed_level == 7, "short positron lv parsed");

    char capture_long[] =
        "<electronCapture energy=\"1923.55\" intensity=\"0.76304\" "
        "ecK=\"0.8732\" lvl=\"1\"/>";
    const SandiaDecay::RadParticle ec = parse(capture_long);
    check(ec.type == SandiaDecay::CaptureElectronParticle,
          "long electron-capture particle type");
    check(ec.fed_level == 1, "long electron-capture lvl parsed");
    check(ec.ec_K > 0.8731f && ec.ec_K < 0.8733f,
          "long electron-capture shell fraction parsed");

    char no_level[] = "<positron energy=\"100\" intensity=\"0.1\"/>";
    check(parse(no_level).fed_level == -1,
          "missing positron level remains explicit");

    char verified_long[] =
        "<gamma energy=\"6048.2\" intensity=\"0.0014\" mult=\"1\" "
        "e0Verified=\"1\"/>";
    check(parse(verified_long).e0_verified,
          "long e0Verified provenance flag parsed");
    char verified_short[] = "<g e=\"6048.2\" i=\"0.0014\" e0=\"1\"/>";
    check(parse(verified_short).e0_verified,
          "short e0 provenance flag parsed");
    char advisory_only[] =
        "<gamma energy=\"373.3\" intensity=\"0.62\" mult=\"1\"/>";
    check(!parse(advisory_only).e0_verified,
          "advisory mult=1 alone is not verified E0 provenance");
    char e0_zero[] = "<gamma energy=\"1\" e0Verified=\"0\"/>";
    char e0_negative[] = "<gamma energy=\"1\" e0Verified=\"-1\"/>";
    char e0_two[] = "<gamma energy=\"1\" e0Verified=\"2\"/>";
    char e0_text[] = "<gamma energy=\"1\" e0Verified=\"yes\"/>";
    char e0_junk[] = "<gamma energy=\"1\" e0Verified=\"1junk\"/>";
    check(!parse(e0_zero).e0_verified && !parse(e0_negative).e0_verified &&
              !parse(e0_two).e0_verified && !parse(e0_text).e0_verified &&
              !parse(e0_junk).e0_verified,
          "noncanonical long E0 provenance tokens fail closed");
    char short_e0_zero[] = "<g e=\"1\" e0=\"0\"/>";
    char short_e0_negative[] = "<g e=\"1\" e0=\"-1\"/>";
    char short_e0_two[] = "<g e=\"1\" e0=\"2\"/>";
    char short_e0_text[] = "<g e=\"1\" e0=\"yes\"/>";
    char short_e0_junk[] = "<g e=\"1\" e0=\"1junk\"/>";
    check(!parse(short_e0_zero).e0_verified &&
              !parse(short_e0_negative).e0_verified &&
              !parse(short_e0_two).e0_verified &&
              !parse(short_e0_text).e0_verified &&
              !parse(short_e0_junk).e0_verified,
          "noncanonical short E0 provenance tokens fail closed");

    char feed_long[] =
        "<transition mode=\"a\" branchRatio=\"1\">"
        "<levelFeed lvl=\"0\" probability=\"0.25\" q=\"10\"/>"
        "<levelFeed lvl=\"3\" probability=\"0.75\" q=\"7\"/>"
        "<alpha energy=\"1\" intensity=\"1\"/>"
        "</transition>";
    const SandiaDecay::Transition long_feed = parse_transition(feed_long);
    check(long_feed.level_feeds.size() == 2 && long_feed.products.size() == 1 &&
              long_feed.level_feeds[1].level == 3,
          "long levelFeed children parse separately from particles");
    char feed_short[] =
        "<t m=\"b-\" br=\"1\"><lf l=\"0\" p=\"0.4\" q=\"10\"/>"
        "<lf l=\"7\" p=\"0.6\" q=\"3\" f=\"1u\"/>"
        "<b e=\"1\" i=\"1\"/></t>";
    const SandiaDecay::Transition short_feed = parse_transition(feed_short);
    check(short_feed.level_feeds.size() == 2 && short_feed.products.size() == 1 &&
              short_feed.level_feeds[1].level == 7 &&
              short_feed.level_feeds[1].forbiddenness ==
                  SandiaDecay::FirstUniqueForbidden,
          "short lf children parse endpoint/shape separately from particles");
    char duplicate_feed[] =
        "<transition mode=\"a\"><levelFeed lvl=\"0\" probability=\"0.5\" q=\"10\"/>"
        "<levelFeed lvl=\"0\" probability=\"0.5\" q=\"10\"/></transition>";
    check(transition_parse_throws(duplicate_feed),
          "duplicate daughter-level feed records fail closed");
    char malformed_feed[] =
        "<transition mode=\"b-\"><levelFeed lvl=\"0junk\" probability=\"1\" "
        "q=\"10\"/></transition>";
    check(transition_parse_throws(malformed_feed),
          "noncanonical level-feed index fails closed");
    char malformed_feed_float[] =
        "<transition mode=\"b-\"><levelFeed lvl=\"0\" probability=\"1junk\" "
        "q=\"10\"/></transition>";
    check(transition_parse_throws(malformed_feed_float),
          "noncanonical level-feed floating value fails closed");
    char unknown_shape[] =
        "<transition mode=\"b-\"><levelFeed lvl=\"0\" probability=\"1\" q=\"10\" "
        "forbiddenness=\"futureShape\"/></transition>";
    check(transition_parse_throws(unknown_shape),
          "unknown level-feed beta shape fails closed");
    char unsupported_mode[] =
        "<transition mode=\"ec\"><levelFeed lvl=\"0\" probability=\"1\" "
        "q=\"10\"/></transition>";
    check(transition_parse_throws(unsupported_mode),
          "level-feed metadata on an unsupported decay mode fails closed");

    const std::string source = SANDIA_DECAY_SOURCE_DIR;
    check_database_provenance(source + "/sandia.decay.xml");
    check_database_provenance(source + "/sandia.decay.min.xml");
    check_database_provenance(source + "/sandia.decay.nocoinc.min.xml");
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
