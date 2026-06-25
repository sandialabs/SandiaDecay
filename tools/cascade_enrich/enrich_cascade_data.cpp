/* SandiaDecay cascade-detail enrichment tool.
 *
 * Reads sandia.decay.xml and, using the GEANT4 11.4.0 PhotonEvaporation and
 * RadioactiveDecay level data, injects the extra per-transition / per-particle
 * data needed for high-fidelity true-coincidence-summing:
 *
 *   <gamma ... iccTotal= iccK= iccL1= iccL2= iccL3= mult= delta=>  IC + multipole
 *   <coincidentGamma ... a2= a4=>                    gamma-gamma angular correlation
 *   <electronCapture ... ecK= ecL=>                  EC shell fractions (K, L)
 *
 * a2/a4 are computed from the nuclear spin sequence + multipole mixing via the
 * F-coefficient (Racah) formalism (see AngularCorrelation). All injected
 * attributes are optional; a database without them parses identically (the
 * mcDetEff adapter / InterSpec fall back to the transition-group x-ray model and
 * isotropic emission).
 *
 * Usage: enrich_cascade_data <in.xml> <out.xml> [photonEvapDir] [radDecayDir]
 *   dirs default to $G4LEVELGAMMADATA / $G4RADIOACTIVEDATA.
 */

#include "AngularCorrelation.h"
#include "G4LevelData.h"
#include "G4RadDecayData.h"
#include "StructureMatch.h"

#include <set>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

using namespace cascade_enrich;
using rapidxml::xml_node;
using rapidxml::xml_attribute;

namespace {

// Space-only gate (NOT a correctness guard): the 10 keV x-ray cut is enforced
// downstream by the per-subshell L fluorescence table (sub-10-keV lines removed
// -> empty subshells for light elements) and by the engine skipping <10 keV
// lines, so no sub-10-keV L x-ray can ever be emitted regardless of this gate.
// Its only purpose is to avoid writing dead iccL1/L2/L3 + ecL attributes onto
// the many low-Z IC gammas whose L x-rays are all cut. It is set to the element
// table's actual onset (Z=69, Tm: the lowest Z keeping any L line >= 10 keV after
// the cut -- the L1->N Lgamma line) so it can never over-cut data the energy cut
// would keep. Regenerate the onset from tools/extract_g4_xs --fluor-l if the cut
// changes.
constexpr int kLxrayMinZ = 69;

std::string env_or(const char* var, const std::string& fallback) {
  const char* v = std::getenv(var);
  return (v && v[0]) ? std::string(v) : fallback;
}

double attr_double(xml_node<char>* node, const char* name, bool& present) {
  xml_attribute<char>* a = node->first_attribute(name);
  present = (a != nullptr);
  return present ? std::atof(a->value()) : 0.0;
}

// Append (or overwrite) a double attribute, allocating strings in the document.
void set_attr(rapidxml::xml_document<char>& doc, xml_node<char>* node,
              const char* name, double value) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6G", value);
  const char* val = doc.allocate_string(buf);
  xml_attribute<char>* existing = node->first_attribute(name);
  if (existing) {
    existing->value(val);
  } else {
    node->append_attribute(doc.allocate_attribute(doc.allocate_string(name), val));
  }
}

void set_attr_int(rapidxml::xml_document<char>& doc, xml_node<char>* node,
                  const char* name, int value) {
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%d", value);
  const char* val = doc.allocate_string(buf);
  if (node->first_attribute(name))
    node->first_attribute(name)->value(val);
  else
    node->append_attribute(doc.allocate_attribute(doc.allocate_string(name), val));
}

long zakey(int Z, int A) { return (static_cast<long>(Z) << 16) | A; }

// Index of the level in `s` nearest `ex_keV` (within tol), else -1.
int level_index_for_energy(const G4LevelScheme& s, double ex_keV, double tol = 1.5) {
  int best = -1;
  double bd = tol;
  for (const G4Level& l : s.levels) {
    const double d = std::fabs(l.excitation_keV - ex_keV);
    if (d <= bd) { bd = d; best = l.index; }
  }
  return best;
}

// Nearest G4 transition to energy e (any distance); for documenting misses.
const G4Transition* nearest_g4(const G4LevelScheme& s, double e, double& gap) {
  const G4Transition* best = nullptr;
  double bd = 1e18;
  for (const G4Transition& t : s.transitions) {
    const double d = std::fabs(t.gamma_keV - e);
    if (d < bd) { bd = d; best = &t; }
  }
  gap = best ? (e - best->gamma_keV) : 0.0;
  return best;
}

int level_twoJ(const G4LevelScheme& s, int idx, bool& known) {
  const G4Level* l = s.level(idx);
  if (!l || !l->j_known) { known = false; return -1; }
  known = true;
  return l->twoJ;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <in.xml> <out.xml> [photonEvapDir] [radDecayDir]\n";
    return 1;
  }
  const std::string in_path = argv[1];
  const std::string out_path = argv[2];
  const std::string photoneva = (argc > 3)
      ? argv[3]
      : env_or("G4LEVELGAMMADATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/PhotonEvaporation6.1.2");
  const std::string raddecay = (argc > 4)
      ? argv[4]
      : env_or("G4RADIOACTIVEDATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/RadioactiveDecay6.1.2");

  {
    std::ifstream test(out_path);
    if (test.is_open()) {
      std::cerr << "Output file '" << out_path << "' exists - not overwriting.\n";
      return 1;
    }
  }

  try {
    rapidxml::file<char> xmlfile(in_path.c_str());
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_full>(xmlfile.data());

    xml_node<char>* root = doc.first_node();
    if (root && root->first_attribute("version"))
      root = root->next_sibling();
    if (!root) { std::cerr << "No document root.\n"; return 1; }

    // Symbol -> (Z,A) from <nuclide> elements.
    std::map<std::string, std::pair<int,int>> sym_to_za;
    for (xml_node<char>* n = root->first_node("nuclide"); n;
         n = n->next_sibling("nuclide")) {
      xml_attribute<char>* sa = n->first_attribute("symbol");
      xml_attribute<char>* za = n->first_attribute("atomicNumber");
      xml_attribute<char>* ma = n->first_attribute("massNumber");
      if (sa && za && ma)
        sym_to_za[sa->value()] = {std::atoi(za->value()), std::atoi(ma->value())};
    }

    std::map<long, G4LevelScheme> level_cache;
    auto get_scheme = [&](int Z, int A) -> const G4LevelScheme& {
      const long k = zakey(Z, A);
      auto it = level_cache.find(k);
      if (it == level_cache.end())
        it = level_cache.emplace(k, load_photon_evaporation(Z, A, photoneva)).first;
      return it->second;
    };
    std::map<long, G4ECData> ec_cache;
    auto get_ec = [&](int Z, int A) -> const G4ECData& {
      const long k = zakey(Z, A);
      auto it = ec_cache.find(k);
      if (it == ec_cache.end())
        it = ec_cache.emplace(k, load_ec_shell_fractions(Z, A, raddecay)).first;
      return it->second;
    };

    long n_trans = 0, n_gamma = 0, n_accepted = 0, n_struct_changed = 0,
         n_struct_used = 0, n_ambiguous = 0, n_coinc = 0, n_coinc_a2a4 = 0,
         n_icc = 0, n_ec = 0, n_ec_matched = 0, n_doc = 0;

    for (xml_node<char>* tr = root->first_node("transition"); tr;
         tr = tr->next_sibling("transition")) {
      ++n_trans;
      xml_attribute<char>* pa = tr->first_attribute("parent");
      xml_attribute<char>* ca = tr->first_attribute("child");
      if (!pa || !ca) continue;
      auto cit = sym_to_za.find(ca->value());
      if (cit == sym_to_za.end()) continue;
      const int daughter_Z = cit->second.first;
      const G4LevelScheme& scheme = get_scheme(cit->second.first, cit->second.second);

      // Collect gammas (node, energy, intensity) and gamma-id -> index.
      struct GNode { xml_node<char>* node; double energy; double intensity; };
      std::vector<GNode> gs;
      std::map<std::string, int> id_to_idx;
      for (xml_node<char>* g = tr->first_node("gamma"); g;
           g = g->next_sibling("gamma")) {
        bool he; const double E = attr_double(g, "energy", he);
        if (!he) continue;
        bool hi; const double I = attr_double(g, "intensity", hi);
        xml_attribute<char>* ida = g->first_attribute("id");
        if (ida && ida->value_size())
          id_to_idx[ida->value()] = static_cast<int>(gs.size());
        gs.push_back({g, E, hi ? I : 0.0});
      }
      const std::size_t ng = gs.size();
      n_gamma += static_cast<long>(ng);
      if (ng == 0) continue;

      // Resolve coincidence partners (by id) and run the structure-aware match.
      std::vector<SdGamma> sg(ng);
      for (std::size_t i = 0; i < ng; ++i) {
        sg[i].energy = gs[i].energy;
        for (xml_node<char>* cg = gs[i].node->first_node("coincidentGamma"); cg;
             cg = cg->next_sibling("coincidentGamma")) {
          xml_attribute<char>* ida = cg->first_attribute("id");
          if (!ida) continue;
          auto it = id_to_idx.find(ida->value());
          if (it != id_to_idx.end()) sg[i].partners.push_back(it->second);
        }
      }
      const std::vector<GammaMatch> mm = match_gammas(scheme, sg);

      std::vector<int> from_lvl(ng, -1), to_lvl(ng, -1);
      std::vector<const G4Transition*> gtr(ng, nullptr);
      std::set<int> populated;
      for (std::size_t i = 0; i < ng; ++i) {
        if (mm[i].ambiguous) ++n_ambiguous;
        if (!mm[i].accepted) continue;
        ++n_accepted;
        if (mm[i].structure_changed) ++n_struct_changed;
        if (mm[i].structure_used) ++n_struct_used;
        gtr[i] = &scheme.transitions[mm[i].g4];
        from_lvl[i] = gtr[i]->from_level;
        to_lvl[i] = gtr[i]->to_level;
        populated.insert(gtr[i]->from_level);
        populated.insert(gtr[i]->to_level);
        // Internal-conversion coefficients (raw inputs for the vacancy-level
        // x-ray model). Written only when nonzero, on IC-active gammas.
        // K is written for all Z (K x-rays are >10 keV for every EC/IC-relevant
        // daughter). L is written per-subshell (iccL1/iccL2/iccL3) ONLY for
        // daughters whose L x-rays clear the 10 keV file cut (Z >= kLxrayMinZ);
        // below that the L x-rays are sub-10 keV and are not modelled, so the
        // (space-costly) per-subshell L coefficients are omitted entirely.
        if (gtr[i]->alpha_total > 0.0) {
          set_attr(doc, gs[i].node, "iccTotal", gtr[i]->alpha_total);
          set_attr(doc, gs[i].node, "iccK", gtr[i]->alpha_K());
          if (daughter_Z >= kLxrayMinZ) {
            set_attr(doc, gs[i].node, "iccL1", gtr[i]->alpha_Lsub(1));
            set_attr(doc, gs[i].node, "iccL2", gtr[i]->alpha_Lsub(2));
            set_attr(doc, gs[i].node, "iccL3", gtr[i]->alpha_Lsub(3));
          }
          ++n_icc;
        }
        // Level-scheme topology (G4 level indices) for the level-path cascade
        // realization: which daughter levels this gamma connects.
        set_attr_int(doc, gs[i].node, "fl", gtr[i]->from_level);
        set_attr_int(doc, gs[i].node, "tl", gtr[i]->to_level);
      }

      // a2/a4 on <coincidentGamma> for accepted, adjacent gamma pairs (lean: only
      // the derived angular-correlation coefficients are written).
      for (std::size_t i = 0; i < ng; ++i) {
        if (!gtr[i]) continue;
        for (xml_node<char>* cg = gs[i].node->first_node("coincidentGamma"); cg;
             cg = cg->next_sibling("coincidentGamma")) {
          ++n_coinc;
          xml_attribute<char>* ida = cg->first_attribute("id");
          if (!ida) continue;
          auto it = id_to_idx.find(ida->value());
          if (it == id_to_idx.end()) continue;
          const int j = it->second;
          if (j < 0 || !gtr[j]) continue;

          int twoJ_up = -1, twoJ_inter = -1, twoJ_low = -1;
          int Lu0 = 0, Lup = 0, Ll0 = 0, Llp = 0;
          double du = 0.0, dl = 0.0;
          bool adj = false, k1 = false, k2 = false, k3 = false;
          if (to_lvl[i] == from_lvl[j]) {
            adj = true;
            twoJ_up = level_twoJ(scheme, from_lvl[i], k1);
            twoJ_inter = level_twoJ(scheme, to_lvl[i], k2);
            twoJ_low = level_twoJ(scheme, to_lvl[j], k3);
            Lu0 = gtr[i]->L0; Lup = gtr[i]->Lp; du = gtr[i]->mixing_delta;
            Ll0 = gtr[j]->L0; Llp = gtr[j]->Lp; dl = gtr[j]->mixing_delta;
          } else if (to_lvl[j] == from_lvl[i]) {
            adj = true;
            twoJ_up = level_twoJ(scheme, from_lvl[j], k1);
            twoJ_inter = level_twoJ(scheme, to_lvl[j], k2);
            twoJ_low = level_twoJ(scheme, to_lvl[i], k3);
            Lu0 = gtr[j]->L0; Lup = gtr[j]->Lp; du = gtr[j]->mixing_delta;
            Ll0 = gtr[i]->L0; Llp = gtr[i]->Lp; dl = gtr[i]->mixing_delta;
          }
          if (!adj || !(k1 && k2 && k3)) continue;
          const AngularCorrCoeffs c = angular_correlation(
              twoJ_up, Lu0, Lup, du, twoJ_inter, Ll0, Llp, dl, twoJ_low);
          if (!c.valid) continue;
          set_attr(doc, cg, "a2", c.a2);
          set_attr(doc, cg, "a4", c.a4);
          ++n_coinc_a2a4;
        }
      }

      // Documentation comment: SandiaDecay gammas with no accepted G4 match, and
      // G4 transitions on populated cascade levels that no SD gamma claimed.
      std::string doc_sd, doc_g4;
      int nd = 0;
      for (std::size_t i = 0; i < ng; ++i) {
        if (mm[i].accepted) continue;
        if (nd++ >= 25) { doc_sd += " ..."; break; }
        const char* reason;
        if (std::fabs(gs[i].energy - 510.999) < 1.0) {
          reason = "annih511";
        } else {
          double gap = 0.0;
          const G4Transition* nt =
              scheme.loaded ? nearest_g4(scheme, gs[i].energy, gap) : nullptr;
          reason = (nt && std::fabs(gap) <= 3.0) ? "neartol" : "absent";
        }
        char b[96];
        std::snprintf(b, sizeof(b), " %.2f(I=%.3g,%s)",
                      gs[i].energy, gs[i].intensity, reason);
        doc_sd += b;
      }
      if (scheme.loaded) {
        int ng4 = 0;
        for (const G4Transition& t : scheme.transitions) {
          if (t.rel_intensity < 5.0 || !populated.count(t.from_level)) continue;
          bool claimed = false;
          for (std::size_t i = 0; i < ng; ++i)
            if (gtr[i] == &t) { claimed = true; break; }
          if (claimed) continue;
          if (ng4++ >= 25) { doc_g4 += " ..."; break; }
          char b[64];
          std::snprintf(b, sizeof(b), " %.2f(rI=%.3g)", t.gamma_keV, t.rel_intensity);
          doc_g4 += b;
        }
      }
      if (!doc_sd.empty() || !doc_g4.empty()) {
        std::string txt = " cascade-match:";
        if (!scheme.loaded) txt += " (no G4 level file for daughter)";
        if (!doc_sd.empty()) txt += " SD-unmatched[keV]:" + doc_sd + ";";
        if (!doc_g4.empty()) txt += " G4-unclaimed-on-populated[keV]:" + doc_g4 + ";";
        txt += " ";
        tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                          doc.allocate_string(txt.c_str())));
        ++n_doc;
      }

      // electronCapture shell fractions (raw input for the vacancy-level x-ray
      // model), matched by EC Q-value to the parent's RadioactiveDecay data.
      auto pit = sym_to_za.find(pa->value());
      if (pit != sym_to_za.end()) {
        const G4ECData& ec = get_ec(pit->second.first, pit->second.second);
        if (ec.loaded) {
          for (xml_node<char>* e = tr->first_node("electronCapture"); e;
               e = e->next_sibling("electronCapture")) {
            ++n_ec;
            bool he; const double q = attr_double(e, "energy", he);
            const ECShellFractions* f = ec.for_q(he ? q : 0.0, 1.5);
            if (!f) continue;
            // K vacancy fraction for all Z (K x-rays > 10 keV everywhere this
            // matters). L only for daughters whose L x-rays clear the 10 keV cut;
            // M/N x-rays are sub-10 keV for every element and are never modelled,
            // so ecM/ecN are dropped entirely.
            set_attr(doc, e, "ecK", f->fK);
            if (daughter_Z >= kLxrayMinZ)
              set_attr(doc, e, "ecL", f->fL);
            // Daughter level fed by this EC branch (index in the daughter scheme),
            // for the level-path cascade realization.
            const int lvl = scheme.loaded
                ? level_index_for_energy(scheme, f->daughter_ex_keV) : -1;
            if (lvl >= 0) set_attr_int(doc, e, "lvl", lvl);
            ++n_ec_matched;
          }
        }
      }
    }

    std::ofstream out(out_path, std::ios::out | std::ios::binary);
    if (!out) { std::cerr << "Cannot open output.\n"; return 1; }
    rapidxml::print<char>(out, doc, 0);

    const double pct = n_gamma ? 100.0 * n_accepted / n_gamma : 0.0;
    std::printf("Enriched %ld transitions (%ld documented with miss comments):\n",
                n_trans, n_doc);
    std::printf("  gammas: %ld; accepted G4 matches: %ld (%.1f%%)\n",
                n_gamma, n_accepted, pct);
    std::printf("  ambiguous (2+ within tol): %ld; structure changed the pick: %ld; "
                "accepted via structure only: %ld\n",
                n_ambiguous, n_struct_changed, n_struct_used);
    std::printf("  coincidences scanned: %ld; got a2/a4: %ld\n", n_coinc, n_coinc_a2a4);
    std::printf("  gammas with ICC written: %ld; electronCapture with shell fractions: %ld/%ld\n",
                n_icc, n_ec_matched, n_ec);
    std::printf("Wrote %s\n", out_path.c_str());
  } catch (std::exception& e) {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return 0;
}
