/* SandiaDecay cascade-detail enrichment tool.
 *
 * Reads sandia.decay.xml and, using the GEANT4 11.4.0 PhotonEvaporation and
 * RadioactiveDecay level data, injects the extra per-transition / per-particle
 * data needed for high-fidelity true-coincidence-summing:
 *
 *   <gamma ... iccTotal= iccK= iccL= mult= delta=>   internal conversion + multipole
 *   <coincidentGamma ... a2= a4=>                    gamma-gamma angular correlation
 *   <electronCapture ... ecK= ecL= ecM= ecN=>        EC shell fractions
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

const double kGammaMatchTolKeV = 0.6;
const double kECQMatchTolKeV = 1.5;

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

// Per-gamma resolution against the daughter G4 level scheme.
struct GammaInfo {
  xml_node<char>* node = nullptr;
  double energy = 0.0;
  int from_level = -1, to_level = -1;
  bool matched = false;
};

const G4Transition* match_gamma(const G4LevelScheme& s, double e) {
  const G4Transition* best = nullptr;
  double bd = kGammaMatchTolKeV;
  for (const G4Transition& t : s.transitions) {
    const double d = std::fabs(t.gamma_keV - e);
    if (d <= bd) { bd = d; best = &t; }
  }
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
    std::map<long, G4ECData> ec_cache;
    auto get_scheme = [&](int Z, int A) -> const G4LevelScheme& {
      const long k = zakey(Z, A);
      auto it = level_cache.find(k);
      if (it == level_cache.end())
        it = level_cache.emplace(k, load_photon_evaporation(Z, A, photoneva)).first;
      return it->second;
    };
    auto get_ec = [&](int Z, int A) -> const G4ECData& {
      const long k = zakey(Z, A);
      auto it = ec_cache.find(k);
      if (it == ec_cache.end())
        it = ec_cache.emplace(k, load_ec_shell_fractions(Z, A, raddecay)).first;
      return it->second;
    };

    long n_trans = 0, n_gamma = 0, n_gamma_matched = 0, n_icc = 0;
    long n_coinc = 0, n_coinc_a2a4 = 0, n_ec = 0, n_ec_matched = 0;

    for (xml_node<char>* tr = root->first_node("transition"); tr;
         tr = tr->next_sibling("transition")) {
      ++n_trans;
      xml_attribute<char>* pa = tr->first_attribute("parent");
      xml_attribute<char>* ca = tr->first_attribute("child");
      if (!pa || !ca) continue;
      auto pit = sym_to_za.find(pa->value());
      auto cit = sym_to_za.find(ca->value());
      if (pit == sym_to_za.end() || cit == sym_to_za.end()) continue;
      const int Zp = pit->second.first, Ap = pit->second.second;
      const int Zd = cit->second.first, Ad = cit->second.second;

      const G4LevelScheme& scheme = get_scheme(Zd, Ad);

      // Pass 1: gammas -> match to G4, write ICC, record level structure.
      std::map<std::string, GammaInfo*> id_map;  // gamma id -> info
      std::vector<GammaInfo> ginfos;
      ginfos.reserve(16);
      for (xml_node<char>* g = tr->first_node("gamma"); g;
           g = g->next_sibling("gamma")) {
        ++n_gamma;
        GammaInfo info;
        info.node = g;
        bool has_e = false;
        info.energy = attr_double(g, "energy", has_e);
        if (!has_e) continue;
        ginfos.push_back(info);
      }
      for (GammaInfo& info : ginfos) {
        if (scheme.loaded) {
          const G4Transition* t = match_gamma(scheme, info.energy);
          if (t) {
            info.matched = true;
            info.from_level = t->from_level;
            info.to_level = t->to_level;
            ++n_gamma_matched;
            set_attr(doc, info.node, "iccTotal", t->alpha_total);
            set_attr(doc, info.node, "iccK", t->alpha_K());
            set_attr(doc, info.node, "iccL", t->alpha_L());
            set_attr_int(doc, info.node, "mult", t->multipole_code);
            set_attr(doc, info.node, "delta", t->mixing_delta);
            ++n_icc;
          }
        }
        xml_attribute<char>* ida = info.node->first_attribute("id");
        if (ida && ida->value_size())
          id_map[ida->value()] = &info;
      }

      // Pass 2: gamma-gamma angular correlations on <coincidentGamma>.
      if (scheme.loaded) {
        for (GammaInfo& gi : ginfos) {
          if (!gi.matched) continue;
          for (xml_node<char>* cg = gi.node->first_node("coincidentGamma"); cg;
               cg = cg->next_sibling("coincidentGamma")) {
            ++n_coinc;
            xml_attribute<char>* ida = cg->first_attribute("id");
            if (!ida) continue;
            auto pit2 = id_map.find(ida->value());
            if (pit2 == id_map.end()) continue;
            const GammaInfo& pj = *pit2->second;
            if (!pj.matched) continue;

            // Determine the shared intermediate level (adjacent cascade only).
            int twoJ_far_up = -1, twoJ_inter = -1, twoJ_far_low = -1;
            int Lu0 = 0, Lup = 0, Ll0 = 0, Llp = 0;
            double du = 0.0, dl = 0.0;
            bool adjacent = false;
            const G4Transition* tg = match_gamma(scheme, gi.energy);
            const G4Transition* tp = match_gamma(scheme, pj.energy);
            if (!tg || !tp) continue;

            if (gi.to_level == pj.from_level) {
              // gi is upper (from_up -> inter), pj is lower (inter -> far_low)
              adjacent = true;
              bool k1, k2, k3;
              twoJ_far_up = level_twoJ(scheme, gi.from_level, k1);
              twoJ_inter = level_twoJ(scheme, gi.to_level, k2);
              twoJ_far_low = level_twoJ(scheme, pj.to_level, k3);
              Lu0 = tg->L0; Lup = tg->Lp; du = tg->mixing_delta;
              Ll0 = tp->L0; Llp = tp->Lp; dl = tp->mixing_delta;
              if (!(k1 && k2 && k3)) adjacent = false;
            } else if (pj.to_level == gi.from_level) {
              // pj is upper, gi is lower
              adjacent = true;
              bool k1, k2, k3;
              twoJ_far_up = level_twoJ(scheme, pj.from_level, k1);
              twoJ_inter = level_twoJ(scheme, pj.to_level, k2);
              twoJ_far_low = level_twoJ(scheme, gi.to_level, k3);
              Lu0 = tp->L0; Lup = tp->Lp; du = tp->mixing_delta;
              Ll0 = tg->L0; Llp = tg->Lp; dl = tg->mixing_delta;
              if (!(k1 && k2 && k3)) adjacent = false;
            }

            if (!adjacent) continue;
            AngularCorrCoeffs c = angular_correlation(
                twoJ_far_up, Lu0, Lup, du, twoJ_inter, Ll0, Llp, dl, twoJ_far_low);
            if (!c.valid) continue;
            set_attr(doc, cg, "a2", c.a2);
            set_attr(doc, cg, "a4", c.a4);
            ++n_coinc_a2a4;
          }
        }
      }

      // electronCapture shell fractions (matched by EC Q-value).
      const G4ECData& ec = get_ec(Zp, Ap);
      if (ec.loaded) {
        for (xml_node<char>* e = tr->first_node("electronCapture"); e;
             e = e->next_sibling("electronCapture")) {
          ++n_ec;
          bool has_e = false;
          const double q = attr_double(e, "energy", has_e);
          const ECShellFractions* f =
              has_e ? ec.for_q(q, kECQMatchTolKeV) : ec.for_q(0.0, kECQMatchTolKeV);
          if (!f) continue;
          set_attr(doc, e, "ecK", f->fK);
          set_attr(doc, e, "ecL", f->fL);
          set_attr(doc, e, "ecM", f->fM);
          set_attr(doc, e, "ecN", f->fN);
          ++n_ec_matched;
        }
      }
    }

    std::ofstream out(out_path, std::ios::out | std::ios::binary);
    if (!out) { std::cerr << "Cannot open output.\n"; return 1; }
    rapidxml::print<char>(out, doc, 0);

    std::printf("Enriched %ld transitions:\n", n_trans);
    std::printf("  gammas: %ld (%ld matched to G4, %ld ICC written)\n",
                n_gamma, n_gamma_matched, n_icc);
    std::printf("  coincidences: %ld (%ld got a2/a4)\n", n_coinc, n_coinc_a2a4);
    std::printf("  electronCapture: %ld (%ld got shell fractions)\n",
                n_ec, n_ec_matched);
    std::printf("Wrote %s\n", out_path.c_str());
  } catch (std::exception& e) {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return 0;
}
