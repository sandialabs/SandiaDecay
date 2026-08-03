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
      check(l0 && l0->floating_flag == "-",
            "Ni-60 fixed-level floating flag retained");
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

  // --- GEANT4 multipolarity is advisory, not evaluated E0 provenance ---
  // PhotonEvaporation labels each of these records with code 1. The current
  // adopted ENSDF tables, however, retain the Ac-232 and Rb-78 lines as measured
  // gammas with parenthesized/ambiguous level spins and no adopted M/ICC value.
  // GEANT4 also stores code 1 independently of alpha (zero for those two;
  // finite for Gd-154). Keep all raw fields intact so a downstream consumer
  // cannot mistake the code alone for a photon veto or electron source.
  {
    const G4LevelScheme ac232 = load_photon_evaporation(89, 232, photoneva);
    const G4Transition* g373 = find_gamma(ac232, 373.3, 0.05);
    check(ac232.loaded && g373, "load Ac-232 373.3-keV advisory mult=1 record");
    check(g373 && g373->multipole_code == 1,
          "Ac-232 373.3-keV raw GEANT4 multipole code retained");
    check(g373 && g373->alpha_total == 0.0,
          "Ac-232 373.3-keV mult=1 candidate retains zero ICC");
    check(g373 && std::fabs(g373->rel_intensity - 62.0) < 1e-12,
          "Ac-232 373.3-keV mult=1 candidate retains gamma intensity");

    const G4LevelScheme rb78 = load_photon_evaporation(37, 78, photoneva);
    const G4Transition* g190 = find_gamma(rb78, 189.54, 0.05);
    check(rb78.loaded && g190, "load Rb-78 189.54-keV advisory mult=1 record");
    check(g190 && g190->multipole_code == 1 && g190->alpha_total == 0.0,
          "Rb-78 mult=1 candidate independently retains zero ICC");
    check(g190 && std::fabs(g190->rel_intensity - 50.0) < 1e-12,
          "Rb-78 mult=1 candidate retains gamma intensity");

    const G4LevelScheme gd154 = load_photon_evaporation(64, 154, photoneva);
    const G4Transition* g603 = find_gamma(gd154, 602.6683, 0.01);
    check(gd154.loaded && g603, "load Gd-154 finite-ICC mult=1 record");
    check(g603 && g603->multipole_code == 1,
          "Gd-154 602.668-keV raw GEANT4 multipole code retained");
    if (g603)
      check_close(g603->alpha_total, 0.038, 1e-12,
                  "Gd-154 mult=1 candidate independently retains finite ICC");
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

  // --- Parent P-block selection + real per-level beta-plus records (Rb-82) ---
  {
    ParentStateHint ground;
    ground.isomer = 0;
    ground.excitation_keV = 0.0;
    ground.half_life_s = 75.45;
    ParentStateHint isomer;
    isomer.isomer = 1;
    isomer.excitation_keV =
        (81.918282673 - 81.918208598) * 931494.10242;
    isomer.half_life_s = 23299.2;
    ParentStateHint absent;
    absent.isomer = 99;
    absent.excitation_keV = 9999.0;
    absent.half_life_s = 123.0;
    const G4ECData rb82 = load_ec_shell_fractions(37, 82, raddecay, ground);
    const G4ECData rb82m = load_ec_shell_fractions(37, 82, raddecay, isomer);
    const G4ECData missing = load_ec_shell_fractions(37, 82, raddecay, absent);
    check(rb82.loaded, "load Rb-82 ground-state P block");
    check(rb82m.loaded, "load Rb-82m isomer P block");
    check(!missing.loaded, "missing Rb-82 P block is explicit no-match");
    if (rb82.loaded) {
      check_close(rb82.parent_ex_keV, 0.0, 1e-6, "Rb-82 parent excitation");
      check(rb82.selection_method == ParentBlockSelectionMethod::HalfLife,
            "Rb-82 selected by half-life identity");
      const BetaPlusRecord* b = rb82.beta_plus_for_q(4403.0, 0.1);
      check(b != nullptr, "Rb-82 real BetaPlus Q=4403 record found");
      if (b) {
        check_close(b->daughter_ex_keV, 0.0, 1e-6,
                    "Rb-82 Q=4403 feeds Kr-82 ground level");
        check_close(b->intensity, 81.777, 0.001,
                    "Rb-82 ground-level BetaPlus intensity retained");
      }
      check(rb82.for_q(123.0, 0.01) == nullptr,
            "per-level EC Q miss does not return summary");
      check(rb82.shells_for_q(123.0, 0.01) == &rb82.summary,
            "EC shell lookup may explicitly fall back to summary");
    }
    if (rb82m.loaded) {
      check_close(rb82m.parent_ex_keV, 69.0, 0.01,
                  "Rb-82m selects second P block");
      check(rb82m.parent_block == 1 &&
                rb82m.selection_method == ParentBlockSelectionMethod::HalfLife,
            "Rb-82m half-life match is ordinal-aligned block 1");
      const BetaPlusRecord* b = rb82m.beta_plus_for_q(1823.16, 0.1);
      check(b != nullptr, "Rb-82m BetaPlus record comes from isomer block");
      if (b)
        check_close(b->daughter_ex_keV, 2648.36, 0.01,
                    "Rb-82m Q=1823.16 daughter level");
    }
  }

  // Stable ground-state blocks may be absent: isomer number is not a P-block
  // ordinal. Ba-137m is block 0, and Ba-133m is block 2 after a prompt state.
  {
    ParentStateHint ba137m_hint;
    ba137m_hint.isomer = 1;
    ba137m_hint.excitation_keV =
        (136.906537704 - 136.905827384) * 931494.10242;
    ba137m_hint.half_life_s = 153.12;
    const G4ECData ba137m = load_ec_shell_fractions(56, 137, raddecay, ba137m_hint);
    check(ba137m.loaded, "load Ba-137m physical parent block");
    if (ba137m.loaded) {
      check(ba137m.parent_block == 0, "Ba-137m selects nonordinal block 0");
      check_close(ba137m.parent_ex_keV, 661.659, 0.01,
                  "Ba-137m parent excitation 661.659 keV");
      check(ba137m.selection_method == ParentBlockSelectionMethod::HalfLife,
            "Ba-137m selected by half-life identity");
    }

    ParentStateHint ba133m_hint;
    ba133m_hint.isomer = 1;
    ba133m_hint.excitation_keV =
        (132.906316936 - 132.90600749) * 931494.10242;
    ba133m_hint.half_life_s = 140040.0;
    const G4ECData ba133m = load_ec_shell_fractions(56, 133, raddecay, ba133m_hint);
    check(ba133m.loaded && ba133m.parent_block == 2,
          "Ba-133m skips prompt ordinal block and selects 288-keV block");

    ParentStateHint ba137_ground_hint;
    ba137_ground_hint.isomer = 0;
    ba137_ground_hint.excitation_keV = 0.0;
    const G4ECData ba137_ground =
        load_ec_shell_fractions(56, 137, raddecay, ba137_ground_hint);
    check(!ba137_ground.loaded,
          "Ba-137 ground is explicit no-match when G4 omits ground P block");

    ParentStateHint stale_mass = ba137m_hint;
    stale_mass.excitation_keV = 2349.1;
    const G4ECData by_half = load_ec_shell_fractions(56, 137, raddecay, stale_mass);
    check(by_half.loaded && by_half.parent_block == 0 &&
              by_half.selection_method == ParentBlockSelectionMethod::HalfLife,
          "half-life fallback rejects nearby prompt state with stale mass hint");

    ParentStateHint hg191m_hint;
    hg191m_hint.isomer = 1;
    hg191m_hint.half_life_s = 3048.0;
    hg191m_hint.ground_half_life_s = 2940.0;
    const G4ECData hg191m = load_ec_shell_fractions(80, 191, raddecay, hg191m_hint);
    check(!hg191m.loaded,
          "half-life fallback does not relabel near-lived Hg-191 ground block");

    ParentStateHint ambiguous_half;
    ambiguous_half.isomer = 1;
    ambiguous_half.half_life_s = 23299.2;
    ambiguous_half.half_life_identity_unique = false;
    const G4ECData rb82m_ambiguous =
        load_ec_shell_fractions(37, 82, raddecay, ambiguous_half);
    check(!rb82m_ambiguous.loaded,
          "ambiguous sibling-isomer half-life gives explicit no-match");
  }

  // Several GEANT4 files shift the long-lived ground state to a nonzero P
  // excitation and the isomer to P ex=0. Half-life, not P excitation or block
  // ordinal, must lead state identity. Equal-lifetime siblings are suppressed
  // unless excitation independently identifies one.
  {
    const auto hint = [](int isomer, double ex, double half, double ground_half,
                         bool unique = true) {
      ParentStateHint h;
      h.isomer = isomer;
      h.excitation_keV = ex;
      h.half_life_s = half;
      h.ground_half_life_s = ground_half;
      h.half_life_identity_unique = unique;
      return h;
    };
    const auto block = [&](int Z, int A, const ParentStateHint& h) {
      return load_ec_shell_fractions(Z, A, raddecay, h);
    };

    const G4ECData eu_g = block(63, 152, hint(0, 0.0, 4.27195e8, 4.27195e8));
    const G4ECData eu_m = block(63, 152, hint(1, 45.5994, 33521.8, 4.27195e8));
    check(eu_g.loaded && eu_g.parent_block == 1,
          "Eu-152 ground selects nonzero-excitation P block by half-life");
    check(eu_m.loaded && eu_m.parent_block == 0,
          "Eu-152m selects ex=0 P block by half-life");

    const G4ECData y_g = block(39, 85, hint(0, 0.0, 9648.0, 9648.0));
    const G4ECData y_m = block(39, 85, hint(1, 19.7998, 17496.0, 9648.0));
    check(y_g.loaded && y_g.parent_block == 1,
          "Y-85 ground selects 19.68-keV P block by half-life");
    check(y_m.loaded && y_m.parent_block == 0,
          "Y-85m selects ex=0 P block by half-life");

    const G4ECData np_g = block(93, 236, hint(0, 0.0, 4.82831e12, 4.82831e12));
    const G4ECData np_m = block(93, 236, hint(1, -1.0, 81000.0, 4.82831e12));
    check(np_g.loaded && np_g.parent_block == 1,
          "Np-236 ground tolerates evaluated half-life drift to P block 1");
    check(np_m.loaded && np_m.parent_block == 0,
          "Np-236m selects ex=0 P block by half-life");

    const G4ECData sc_g = block(21, 44, hint(0, 0.0, 14292.0, 14292.0));
    const G4ECData sc_m = block(21, 44, hint(1, 270.958, 210996.0, 14292.0));
    check(sc_g.loaded && sc_g.parent_block == 4,
          "Sc-44 ground selects 271.241-keV P block by half-life");
    check(sc_m.loaded && sc_m.parent_block == 0,
          "Sc-44m selects ex=0 P block despite mass-derived excitation");

    const G4ECData rh_g = block(45, 100, hint(0, 0.0, 74880.0, 74880.0));
    const G4ECData rh_m = block(45, 100, hint(1, 107.60, 276.0, 74880.0, false));
    const G4ECData rh_m2 = block(45, 100, hint(2, 105.00, 276.0, 74880.0, false));
    check(rh_g.loaded && rh_g.parent_block == 3,
          "Rh-100 ground selects 107.59-keV P block by half-life");
    check(!rh_m.loaded && !rh_m2.loaded,
          "Rh-100 equal-lifetime isomers do not both claim P block 0");

    const G4ECData ta_g = block(73, 178, hint(0, 0.0, 8496.0, 8496.0));
    const G4ECData ta_m = block(73, 178, hint(1, -1.0, 558.6, 8496.0, false));
    const G4ECData ta_m2 = block(73, 178, hint(2, -1.0, 558.6, 8496.0, false));
    check(!ta_g.loaded, "Ta-178 ground does not claim the 558.6-s P block");
    check(!ta_m.loaded && !ta_m2.loaded,
          "Ta-178 equal-lifetime isomers leave ambiguous P block unmatched");
  }

  // Near-tied Q candidates cannot safely supply daughter-level identity.
  {
    G4ECData qdata;
    qdata.per_level.resize(2);
    qdata.per_level[0].q_keV = 100.0;
    qdata.per_level[1].q_keV = 100.4;
    bool ambiguous = false;
    check(qdata.for_q(100.0, 5.0, &ambiguous) == nullptr && ambiguous,
          "EC Q near tie is suppressed and audited");
    qdata.per_level[1].q_keV = 100.6;
    ambiguous = false;
    check(qdata.for_q(100.0, 5.0, &ambiguous) == &qdata.per_level[0] && !ambiguous,
          "EC Q candidate beyond ambiguity margin keeps nearest match");

    qdata.beta_plus.resize(2);
    qdata.beta_plus[0].q_keV = 200.0;
    qdata.beta_plus[1].q_keV = 200.4;
    ambiguous = false;
    check(qdata.beta_plus_for_q(200.0, 5.0, &ambiguous) == nullptr && ambiguous,
          "beta-plus Q near tie is suppressed and audited");
    qdata.beta_plus[1].q_keV = 200.6;
    ambiguous = false;
    check(qdata.beta_plus_for_q(200.0, 5.0, &ambiguous) == &qdata.beta_plus[0] &&
              !ambiguous,
          "beta-plus Q candidate beyond ambiguity margin keeps nearest match");
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

  // --- Exact alpha/beta-minus daughter-feed records + floating identity ---
  {
    const G4LevelScheme pa234 = load_photon_evaporation(91, 234, photoneva);
    const G4Level* fixed_7392 = pa234.level(1);
    const G4Level* floating_7392 = pa234.level(2);
    check(pa234.loaded && fixed_7392 && floating_7392,
          "load Pa-234 duplicate-energy fixed/floating levels");
    check(fixed_7392 && floating_7392 &&
              std::fabs(fixed_7392->excitation_keV - 73.92) < 1e-9 &&
              std::fabs(floating_7392->excitation_keV - 73.92) < 1e-9 &&
              fixed_7392->floating_flag == "-" &&
              floating_7392->floating_flag == "+X",
          "Pa-234 73.92-keV floating flag is exact level identity");

    ParentStateHint th_hint;
    th_hint.isomer = 0;
    th_hint.excitation_keV = 0.0;
    th_hint.half_life_s = 2082240.0;
    const G4ECData th234 =
        load_ec_shell_fractions(90, 234, raddecay, th_hint);
    check(th234.loaded &&
              th234.selection_method == ParentBlockSelectionMethod::HalfLife,
          "Th-234 feed parent block selected physically");
    check(th234.beta_minus.summary_count == 1,
          "Th-234 has exactly one beta-minus mode summary");
    check_close(th234.beta_minus.branch_fraction, 1.0, 1e-12,
                "Th-234 beta-minus mode branch");
    check(th234.beta_minus.records.size() == 4,
          "Th-234 retains all four beta-minus daughter feeds");
    if (th234.beta_minus.records.size() == 4) {
      const DaughterLevelFeedRecord& main = th234.beta_minus.records[0];
      check_close(main.daughter_ex_keV, 73.92, 1e-9,
                  "Th-234 dominant feed excitation");
      check(main.floating_flag == "+X",
            "Th-234 dominant feed retains +X identity");
      check_close(main.intensity_percent, 78.078, 1e-9,
                  "Th-234 dominant feed percent");
    }

    ParentStateHint u235_hint;
    u235_hint.isomer = 0;
    u235_hint.excitation_keV = 0.0;
    u235_hint.half_life_s = 2.221655e16;
    const G4ECData u235 =
        load_ec_shell_fractions(92, 235, raddecay, u235_hint);
    check(u235.loaded && u235.alpha.summary_count == 1 &&
              u235.alpha.records.size() == 20,
          "U-235 retains complete 20-record alpha feed law");
    if (u235.alpha.records.size() == 20) {
      check_close(u235.alpha.records[0].intensity_percent, 4.7121, 1e-9,
                  "U-235 alpha ground-state feed percent");
      check_close(u235.alpha.records[5].daughter_ex_keV, 205.3101, 1e-9,
                  "U-235 dominant alpha daughter excitation");
      check_close(u235.alpha.records[5].intensity_percent, 57.088, 1e-9,
                  "U-235 dominant alpha feed percent");
    }

    ParentStateHint na28_hint;
    na28_hint.isomer = 0;
    na28_hint.excitation_keV = 0.0;
    const G4ECData na28 =
        load_ec_shell_fractions(11, 28, raddecay, na28_hint);
    check(na28.loaded && na28.beta_minus.summary_count == 1,
          "Na-28 beta-minus parent block and summary load");
    bool found_unique_first = false;
    for (const DaughterLevelFeedRecord& record : na28.beta_minus.records)
      if (std::fabs(record.daughter_ex_keV - 5171.3) < 0.01 &&
          record.forbiddenness == "uniqueFirstForbidden")
        found_unique_first = true;
    check(found_unique_first,
          "Na-28 retains raw unique-first-forbidden beta shape token");
  }

  if (g_failures == 0) {
    std::cout << "\nALL PARSER TESTS PASSED\n";
    return 0;
  }
  std::cerr << "\n" << g_failures << " parser test(s) FAILED\n";
  return 1;
}
