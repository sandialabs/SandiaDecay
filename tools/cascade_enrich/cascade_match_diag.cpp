/* Diagnostic: characterize how SandiaDecay gammas match GEANT4 PhotonEvaporation
 * transitions, so we can see WHY matches miss and whether we can do better.
 *
 * For every <transition> with gammas, loads the daughter G4 level scheme and, for
 * each SandiaDecay gamma, finds the nearest G4 transition (regardless of
 * tolerance) and records the signed energy gap, the G4 transition's relative
 * intensity, and whether a second G4 candidate lies within the tolerance
 * (ambiguity). Also lists "significant" G4 transitions (rel_intensity >= thresh)
 * that no SandiaDecay gamma claimed.
 *
 * Usage: cascade_match_diag <in.xml> [photonEvapDir]
 *   -> writes sd_misses.csv, g4_unclaimed.csv, and a stdout summary.
 */

#include "G4LevelData.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

using namespace cascade_enrich;
using rapidxml::xml_node;
using rapidxml::xml_attribute;

namespace {
const double kTol = 0.6;          // current production tolerance (keV)
const double kG4SignificantI = 5; // rel_intensity threshold for "G4 has a gamma SD lacks"

std::string env_or(const char* var, const std::string& fb) {
  const char* v = std::getenv(var);
  return (v && v[0]) ? std::string(v) : fb;
}
double dattr(xml_node<char>* n, const char* name, bool& ok) {
  xml_attribute<char>* a = n->first_attribute(name);
  ok = (a != nullptr);
  return ok ? std::atof(a->value()) : 0.0;
}

// Nearest G4 transition to energy e; returns nullptr only if scheme empty.
const G4Transition* nearest(const G4LevelScheme& s, double e, double& gap,
                            int& n_within_tol) {
  const G4Transition* best = nullptr;
  double bd = 1e9;
  n_within_tol = 0;
  for (const G4Transition& t : s.transitions) {
    const double d = std::fabs(t.gamma_keV - e);
    if (d <= kTol) ++n_within_tol;
    if (d < bd) { bd = d; best = &t; }
  }
  gap = best ? (e - best->gamma_keV) : 0.0;  // signed: SD - G4
  return best;
}
long zakey(int Z, int A) { return (static_cast<long>(Z) << 16) | A; }
} // namespace

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "Usage: %s <in.xml> [photonEvapDir]\n", argv[0]); return 1; }
  const std::string photoneva = (argc > 2) ? argv[2]
      : env_or("G4LEVELGAMMADATA",
               "/Users/wcjohns/install/GEANT4-11.4.0/share/Geant4/data/PhotonEvaporation6.1.2");

  rapidxml::file<char> xmlfile(argv[1]);
  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_default>(xmlfile.data());
  xml_node<char>* root = doc.first_node();
  if (root && root->first_attribute("version")) root = root->next_sibling();

  std::map<std::string, std::pair<int,int>> za;
  for (xml_node<char>* n = root->first_node("nuclide"); n; n = n->next_sibling("nuclide")) {
    auto* s = n->first_attribute("symbol"); auto* z = n->first_attribute("atomicNumber");
    auto* m = n->first_attribute("massNumber");
    if (s && z && m) za[s->value()] = {std::atoi(z->value()), std::atoi(m->value())};
  }
  std::map<long, G4LevelScheme> cache;
  auto scheme = [&](int Z, int A) -> const G4LevelScheme& {
    long k = zakey(Z, A);
    auto it = cache.find(k);
    if (it == cache.end()) it = cache.emplace(k, load_photon_evaporation(Z, A, photoneva)).first;
    return it->second;
  };

  FILE* fsd = std::fopen("sd_misses.csv", "w");
  std::fprintf(fsd, "parent,child,E_sd,I_sd,nearestG4,gap,g4_relI,reason\n");
  FILE* fg4 = std::fopen("g4_unclaimed.csv", "w");
  std::fprintf(fg4, "parent,child,E_g4,relI,from,to\n");

  long n_g = 0, n_match = 0, n_nofile = 0, n_ambig = 0;
  double I_tot = 0, I_match = 0, I_miss = 0;
  // gap histogram of "near misses" (nearest exists but > tol)
  std::map<int, long> gap_hist;  // bin = round(gap*10) for |gap| in (tol, 3]

  for (xml_node<char>* tr = root->first_node("transition"); tr; tr = tr->next_sibling("transition")) {
    auto* pa = tr->first_attribute("parent"); auto* ca = tr->first_attribute("child");
    if (!pa || !ca) continue;
    auto pit = za.find(pa->value()), cit = za.find(ca->value());
    if (cit == za.end()) continue;
    const G4LevelScheme& s = scheme(cit->second.first, cit->second.second);

    std::vector<const G4Transition*> claimed;
    for (xml_node<char>* g = tr->first_node("gamma"); g; g = g->next_sibling("gamma")) {
      bool ok; double E = dattr(g, "energy", ok); if (!ok) continue;
      bool oki; double I = dattr(g, "intensity", oki); if (!oki) I = 0.0;
      ++n_g; I_tot += I;
      if (!s.loaded) {
        ++n_nofile; I_miss += I;
        std::fprintf(fsd, "%s,%s,%.3f,%.4g,,,,nofile\n", pa->value(), ca->value(), E, I);
        continue;
      }
      double gap; int nwt; const G4Transition* nt = nearest(s, E, gap, nwt);
      if (nt && std::fabs(gap) <= kTol) {
        ++n_match; I_match += I; claimed.push_back(nt);
        if (nwt > 1) ++n_ambig;
      } else {
        I_miss += I;
        const char* reason = nt ? "neartol" : "noscheme";
        std::fprintf(fsd, "%s,%s,%.3f,%.4g,%.3f,%.3f,%.4g,%s\n",
                     pa->value(), ca->value(), E, I,
                     nt ? nt->gamma_keV : 0.0, gap, nt ? nt->rel_intensity : 0.0, reason);
        if (nt && std::fabs(gap) <= 3.0) {
          int b = (int)std::lround(gap * 10.0);
          ++gap_hist[b];
        }
      }
    }
    // G4 transitions no SD gamma claimed, with significant rel_intensity.
    if (s.loaded) {
      for (const G4Transition& t : s.transitions) {
        if (t.rel_intensity < kG4SignificantI) continue;
        bool was_claimed = false;
        for (const G4Transition* c : claimed) if (c == &t) { was_claimed = true; break; }
        if (was_claimed) continue;
        // only report if NO sd gamma is anywhere near it (else it's just the
        // weaker partner of a matched line)
        std::fprintf(fg4, "%s,%s,%.3f,%.4g,%d,%d\n", pa->value(), ca->value(),
                     t.gamma_keV, t.rel_intensity, t.from_level, t.to_level);
      }
    }
  }
  std::fclose(fsd); std::fclose(fg4);

  std::printf("SandiaDecay gammas: %ld\n", n_g);
  std::printf("  matched (<=%.2f keV): %ld (%.1f%%), by intensity %.1f%%\n",
              kTol, n_match, 100.0 * n_match / n_g, 100.0 * I_match / I_tot);
  std::printf("  unmatched: %ld (%.1f%%), by intensity %.1f%%\n",
              n_g - n_match, 100.0 * (n_g - n_match) / n_g, 100.0 * I_miss / I_tot);
  std::printf("  of which no-G4-file: %ld\n", n_nofile);
  std::printf("  ambiguous matches (2+ G4 within tol): %ld\n", n_ambig);
  std::printf("near-miss signed-gap histogram (SD - G4, keV), |gap| in (%.1f, 3]:\n", kTol);
  for (auto& kv : gap_hist)
    std::printf("  %+5.1f : %ld\n", kv.first / 10.0, kv.second);
  std::printf("Wrote sd_misses.csv, g4_unclaimed.csv\n");
  return 0;
}
