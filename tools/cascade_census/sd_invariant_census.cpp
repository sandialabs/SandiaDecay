/* Database-wide invariant census for sandia.decay.xml and its minimized variants.
 *
 * WHY THIS EXISTS
 * ---------------
 * Several counts describing the coincidence data have been quoted across notes,
 * commit messages and TODO items with incompatible definitions -- "2192", "2134",
 * "3089", "3106", "10355" and "3172 of 131851" have all appeared, and some of them
 * measure genuinely different quantities. This tool assigns each one a single
 * normative definition, emits it with that definition in the CSV header, and lets
 * CTest assert the subset that must hold. Quote numbers from here, by ID.
 *
 * WHAT IT MEASURES
 * ----------------
 * Two passes over one database:
 *
 *   parsed pass  -- walks the SandiaDecayDataBase object model, i.e. exactly what a
 *                   consumer receives. This is the authoritative view for anything
 *                   a consumer can act on.
 *   raw pass     -- a streaming byte scan of the same file. Needed because the
 *                   parser silently drops coincidence links whose `id` matches no
 *                   product (SandiaDecay.cpp, Transition::set) and transitions whose
 *                   parent has no <nuclide> record. Those defects are invisible from
 *                   the object model by construction, so they can only be counted
 *                   against the file itself.
 *
 * ONE DATABASE PER INVOCATION. Loading full+min+nocoinc together costs roughly
 * 900 MB of rapidxml DOM; run it once per variant (as the existing xml_file_tests
 * targets do) and diff the CSVs to compare variants.
 *
 * USAGE
 *   sd_invariant_census <sandia.decay.xml> [--assert] [--csv <dir>] [--tag <name>]
 *
 *   --assert   exit non-zero if a MUST-HOLD invariant fails. Statistical counts
 *              (C1..C5) are reported but never fail the build: they describe known
 *              data defects whose remediation is tracked separately.
 *   --csv DIR  write machine-readable tables to DIR.
 *   --tag NAME label used in output and filenames (default: derived from filename).
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "rapidxml/rapidxml.hpp"

#include "SandiaDecay.h"

namespace {

// ---------------------------------------------------------------------------
// Normative tolerances. Changing any of these changes what the published counts
// mean, so they are named, defined once, and echoed into every CSV header.
// ---------------------------------------------------------------------------

/** C1: a directional coincidence value above this is not representable as a
    probability. The 1e-9 slack excludes float round-trip noise (values that are
    exactly 1.0 in the evaluation but land at 1.0000000000000002 after a
    decimal->float->decimal round trip). C1b reports the strict >1.0 count so the
    size of that noise band stays visible. */
constexpr double kOverUnityEps = 1e-9;

/** C2: relative slack on the Frechet upper bound, so that a pair which saturates
    the bound exactly is not counted as violating it. */
constexpr double kFrechetRelEps = 1e-6;

/** C5: relative slack on the "a transition happens at most once per decay" bound. */
constexpr double kOccurrenceEps = 1e-9;

/** Level-feed probabilities must sum to one. SandiaDecay's own parser already
    throws outside 2e-4; this reports the achieved margin. */
constexpr double kFeedSumTol = 2e-4;

/** Practical reporting thresholds, in photons per parent decay. Record counts
    alone overstate the importance of vanishingly rare records, so every defect is
    also reported weighted by branchRatio*intensity above these cuts. */
constexpr double kWeightThresholds[] = {1e-5, 1e-4, 1e-3};

struct Counter {
  long long n = 0;              // records matching
  long long n_ge[3] = {0,0,0};  // matching AND above each weight threshold
  double weight = 0.0;          // summed photons per parent decay
  void add(double w) {
    ++n;
    weight += w;
    for (int i = 0; i < 3; ++i)
      if (w >= kWeightThresholds[i]) ++n_ge[i];
  }
};

/** Bounded photon marginal: a gamma cannot be emitted more often than its
    transition occurs, and the transition emits a photon only 1/(1+alpha) of the
    time. This is the same ceiling CeeLo applies in bounded_gamma_intensity(); it
    holds whichever of the tabulated intensity/ICC pair is the wrong one. */
double bounded_intensity(const SandiaDecay::RadParticle &p) {
  const double raw = p.intensity;
  if (p.icc_total > 0.0) return std::min<double>(raw, 1.0 / (1.0 + p.icc_total));
  return raw;
}

bool is_gamma(const SandiaDecay::RadParticle &p) {
  return p.type == SandiaDecay::GammaParticle;
}

// ---------------------------------------------------------------------------
// Raw pass. Counts what is in the FILE, including records the parser discards,
// and attributes every discarded coincidence link to a cause.
//
// This runs BEFORE the database is loaded and frees its DOM before returning,
// so peak resident memory is max(DOM, database) rather than their sum.
//
// The parser's transition-retention rule is SandiaDecay.cpp: a transition is
// kept only when parent != NULL, branchRatio > 0, and the parent's half-life is
// not infinite. Critically, a <nuclide> with NO halfLife attribute defaults to
// infinity, so such a parent is treated as stable and loses every decay it has.
// ---------------------------------------------------------------------------
struct RawCounts {
  long long coincidence_elements = 0;   // <coincidentGamma .../> or <c .../>
  long long gamma_elements = 0;
  long long transition_elements = 0;
  long long nuclide_elements = 0;

  // transitions the parser will discard, by cause
  long long tr_dropped = 0;
  long long tr_drop_no_parent_record = 0;
  long long tr_drop_no_branchratio = 0;
  long long tr_drop_branchratio_le0 = 0;
  long long tr_drop_parent_no_halflife = 0;  // <nuclide> present but halfLife absent
  long long tr_drop_parent_stable = 0;       // halfLife="INF"

  // coincidence links lost, by cause
  long long links_in_dropped_transitions = 0;
  long long links_unresolved_id = 0;         // id matches no product in the transition

  bool ok = false;
};

RawCounts scan_raw(const std::string &path) {
  RawCounts rc;

  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in) return rc;
  const std::streamsize size = in.tellg();
  in.seekg(0);
  std::vector<char> buf(static_cast<size_t>(size) + 1);
  if (!in.read(buf.data(), size)) return rc;
  buf[static_cast<size_t>(size)] = '\0';

  rapidxml::xml_document<char> doc;
  try {
    doc.parse<rapidxml::parse_non_destructive>(buf.data());
  } catch (std::exception &) {
    return rc;
  }

  const rapidxml::xml_node<char> *root = doc.first_node();
  if (!root) return rc;

  // Minimized variants rename elements and attributes (mapping is documented in
  // the min file's own header comment). Note the deliberate collisions: element
  // <c> is coincidentGamma while attribute c is child, and element <p> is
  // positron while attribute p is parent. Element and attribute name spaces are
  // separate, so this is unambiguous -- but it is easy to get wrong.
  const bool shortnames = (root->first_node("t") != nullptr);
  const char *k_tr   = shortnames ? "t"  : "transition";
  const char *k_nuc  = shortnames ? "n"  : "nuclide";
  const char *k_coin = shortnames ? "c"  : "coincidentGamma";
  const char *k_gam  = shortnames ? "g"  : "gamma";
  const char *k_sym  = shortnames ? "s"  : "symbol";
  const char *k_hl   = shortnames ? "hl" : "halfLife";
  const char *k_br   = shortnames ? "br" : "branchRatio";
  const char *k_par  = shortnames ? "p"  : "parent";

  auto attr = [](const rapidxml::xml_node<char> *n, const char *k)
      -> const rapidxml::xml_attribute<char> * {
    return n->first_attribute(k, std::strlen(k));
  };
  auto sval = [](const rapidxml::xml_attribute<char> *a) {
    return a ? std::string(a->value(), a->value_size()) : std::string();
  };

  // Parent half-life state, mirroring the parser: absent attribute => infinite.
  enum class HL { Finite, Infinite, MissingAttribute };
  std::map<std::string, HL> parent_hl;
  for (const rapidxml::xml_node<char> *n = root->first_node(k_nuc); n;
       n = n->next_sibling(k_nuc)) {
    ++rc.nuclide_elements;
    const std::string sym = sval(attr(n, k_sym));
    if (sym.empty()) continue;
    const rapidxml::xml_attribute<char> *hl = attr(n, k_hl);
    if (!hl) { parent_hl[sym] = HL::MissingAttribute; continue; }
    // sscanf("%lf") accepts "INF" on most libc implementations and yields
    // infinity, so a successful parse is not by itself evidence of a finite
    // half-life -- the value has to be tested too. SandiaDecay reaches the same
    // state either way (parse failure also stores infinity).
    double v = 0.0;
    const std::string s = sval(hl);
    const bool parsed = (std::sscanf(s.c_str(), "%lf", &v) == 1);
    parent_hl[sym] = (parsed && !std::isinf(v)) ? HL::Finite : HL::Infinite;
  }

  for (const rapidxml::xml_node<char> *t = root->first_node(k_tr); t;
       t = t->next_sibling(k_tr)) {
    ++rc.transition_elements;

    const std::string par = sval(attr(t, k_par));
    const rapidxml::xml_attribute<char> *br_att = attr(t, k_br);
    const auto it = parent_hl.find(par);

    bool dropped = false;
    if (it == parent_hl.end()) { ++rc.tr_drop_no_parent_record; dropped = true; }
    else if (!br_att)          { ++rc.tr_drop_no_branchratio;   dropped = true; }
    else if (std::atof(sval(br_att).c_str()) <= 0.0)
                               { ++rc.tr_drop_branchratio_le0;  dropped = true; }
    else if (it->second == HL::MissingAttribute)
                               { ++rc.tr_drop_parent_no_halflife; dropped = true; }
    else if (it->second == HL::Infinite)
                               { ++rc.tr_drop_parent_stable;    dropped = true; }
    if (dropped) ++rc.tr_dropped;

    // Product `id`s available as coincidence-link targets in this transition.
    std::set<std::string> ids;
    for (const rapidxml::xml_node<char> *p = t->first_node(); p;
         p = p->next_sibling()) {
      const rapidxml::xml_attribute<char> *id = attr(p, "id");
      if (id) ids.insert(sval(id));
    }

    for (const rapidxml::xml_node<char> *p = t->first_node(); p;
         p = p->next_sibling()) {
      const std::string pname(p->name(), p->name_size());
      if (pname == k_gam) ++rc.gamma_elements;
      for (const rapidxml::xml_node<char> *c = p->first_node(k_coin); c;
           c = c->next_sibling(k_coin)) {
        ++rc.coincidence_elements;
        if (dropped) { ++rc.links_in_dropped_transitions; continue; }
        const rapidxml::xml_attribute<char> *id = attr(c, "id");
        if (!id || !ids.count(sval(id))) ++rc.links_unresolved_id;
      }
    }
  }

  rc.ok = true;
  return rc;  // doc and buf released here, before the database is loaded
}

std::string basename_of(const std::string &p) {
  const size_t s = p.find_last_of("/\\");
  return (s == std::string::npos) ? p : p.substr(s + 1);
}

}  // namespace


int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <sandia.decay.xml> [--assert] [--csv <dir>] [--tag <name>]\n";
    return 2;
  }

  const std::string xml_path = argv[1];
  bool do_assert = false;
  std::string csv_dir, tag, dump_tr;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--assert") do_assert = true;
    else if (a == "--csv" && i + 1 < argc) csv_dir = argv[++i];
    else if (a == "--tag" && i + 1 < argc) tag = argv[++i];
    else if (a == "--dump-transitions" && i + 1 < argc) dump_tr = argv[++i];
    else { std::cerr << "unknown argument: " << a << "\n"; return 2; }
  }
  if (tag.empty()) tag = basename_of(xml_path);

  // ---- raw pass first: cheap, and it tells us what the parser will drop ----
  const RawCounts raw = scan_raw(xml_path);
  if (!raw.ok) {
    std::cerr << "could not read " << xml_path << "\n";
    return 2;
  }

  SandiaDecay::SandiaDecayDataBase db;
  try {
    db.initialize(xml_path);
  } catch (std::exception &e) {
    std::cerr << "failed to parse " << xml_path << ": " << e.what() << "\n";
    return 2;
  }

  // ------------------------- counters -------------------------
  Counter c1_over_unity, c1b_strict, c2_frechet, c2b_frechet_bounded;
  Counter c2_frechet_lower, c5_occurrence, bayes_reverse_gt1;
  long long n_links = 0, n_gammas = 0, n_transitions = 0, n_products = 0;
  long long n_links_with_corr = 0, n_self_links = 0;
  long long c3_bidirectional = 0;

  // MUST-HOLD violations
  long long v_nonfinite = 0, v_negative_intensity = 0, v_index_out_of_range = 0;
  long long v_corr_out_of_range = 0, v_nondescending = 0, v_feed_sum = 0;
  long long v_feed_dup_level = 0, v_e0_with_photon = 0, v_exclusivity = 0;
  long long n_e0_verified = 0, v_e0_not_gamma = 0;
  double worst_feed_sum_err = 0.0;

  // stratification: per decay mode, and per parent nuclide for the worst offenders
  std::map<int, Counter> c2_by_mode;
  std::map<std::string, Counter> c1_by_parent;

  struct Row { std::string parent, child; int mode; double e_i, e_j, I_i, I_j, s_ij, implied, bound, w; };
  std::vector<Row> c1_rows, c2_rows;

  std::ofstream tr_dump;
  if (!dump_tr.empty()) { tr_dump.open(dump_tr); tr_dump << "parent,child,mode,branchRatio\n"; }

  const std::vector<const SandiaDecay::Nuclide *> &nucs = db.nuclides();
  for (const SandiaDecay::Nuclide *nuc : nucs) {
    if (!nuc) continue;
    for (const SandiaDecay::Transition *tr : nuc->decaysToChildren) {
      if (!tr) continue;
      ++n_transitions;
      const double br = tr->branchRatio;
      const std::string parent = nuc->symbol;
      const std::string child = tr->child ? tr->child->symbol : std::string("-");
      if (tr_dump.is_open())
        tr_dump << parent << "," << child << "," << (int)tr->mode << "," << br << "\n";

      // ---- level-feed law invariants ----
      if (!tr->level_feeds.empty()) {
        double sum = 0.0;
        std::set<int> seen;
        for (const SandiaDecay::LevelFeed &lf : tr->level_feeds) {
          sum += lf.probability;
          if (!seen.insert(lf.level).second) ++v_feed_dup_level;
        }
        const double err = std::fabs(sum - 1.0);
        worst_feed_sum_err = std::max(worst_feed_sum_err, err);
        if (err > kFeedSumTol) ++v_feed_sum;
      }

      n_products += static_cast<long long>(tr->products.size());

      for (size_t ip = 0; ip < tr->products.size(); ++ip) {
        const SandiaDecay::RadParticle &p = tr->products[ip];
        if (is_gamma(p)) ++n_gammas;

        if (!std::isfinite(p.energy) || !std::isfinite(p.intensity)
            || !std::isfinite(p.icc_total) || !std::isfinite(p.icc_K))
          ++v_nonfinite;
        if (p.intensity < 0.0) ++v_negative_intensity;

        // C5: a transition cannot occur more than once per parent decay, so
        // I_gamma * (1 + alpha) <= 1. A violation means the tabulated intensity
        // and ICC cannot both be right on the same basis.
        if (is_gamma(p) && p.icc_total > 0.0) {
          const double occ = static_cast<double>(p.intensity) * (1.0 + p.icc_total);
          if (occ > 1.0 + kOccurrenceEps) c5_occurrence.add(br * p.intensity);
        }

        // Verified pure-E0 records keep their historical stored intensity; the
        // photon is suppressed at emission time (SandiaDecay.cpp, the
        // e0_verified branch of the decay-particle assembly), not zeroed in the
        // product record. So a nonzero intensity here is expected, and this is
        // informational rather than a must-hold invariant. What IS asserted is
        // that the flag only ever lands on a gamma.
        if (p.e0_verified) {
          ++n_e0_verified;
          if (p.intensity > 0.0) ++v_e0_with_photon;
          if (!is_gamma(p)) ++v_e0_not_gamma;
        }

        // Descending topology: a gamma must go from a higher level to a lower one.
        if (p.from_level >= 0 && p.to_level >= 0 && p.to_level >= p.from_level)
          ++v_nondescending;

        const double Ii_raw = p.intensity;
        const double Ii_bnd = bounded_intensity(p);

        for (const SandiaDecay::CoincidenceData &cd : p.coincidences) {
          ++n_links;
          const double s = cd.second;

          if (!std::isfinite(s) || !std::isfinite(cd.a2) || !std::isfinite(cd.a4))
            ++v_nonfinite;
          if (cd.first >= tr->products.size()) { ++v_index_out_of_range; continue; }
          if (cd.first == ip) ++n_self_links;
          if (cd.a2 != 0.0f || cd.a4 != 0.0f) ++n_links_with_corr;
          // W(theta) = 1 + a2*P2 + a4*P4 must stay a sane angular distribution.
          if (std::fabs(cd.a2) > 2.0 || std::fabs(cd.a4) > 2.0) ++v_corr_out_of_range;

          const SandiaDecay::RadParticle &q = tr->products[cd.first];
          const double Ij_raw = q.intensity;
          const double Ij_bnd = bounded_intensity(q);
          const double w = br * Ii_raw;

          // C1 / C1b: a directional conditional above 1 is not a probability.
          if (s > 1.0) {
            c1b_strict.add(w);
            if (s > 1.0 + kOverUnityEps) {
              c1_over_unity.add(w);
              c1_by_parent[parent].add(w);
              if (c1_rows.size() < 20000)
                c1_rows.push_back({parent, child, (int)tr->mode, p.energy, q.energy,
                                   Ii_raw, Ij_raw, s, Ii_raw * s, std::min(Ii_raw, Ij_raw), w});
            }
          }

          // C2: Frechet upper bound. The joint P(i and j) implied by the stored
          // conditional cannot exceed min of the two marginals.
          const double implied = Ii_raw * s;
          const double ub_raw = std::min(Ii_raw, Ij_raw);
          if (implied > ub_raw * (1.0 + kFrechetRelEps)) {
            c2_frechet.add(w);
            c2_by_mode[(int)tr->mode].add(w);
            if (c2_rows.size() < 20000)
              c2_rows.push_back({parent, child, (int)tr->mode, p.energy, q.energy,
                                 Ii_raw, Ij_raw, s, implied, ub_raw, w});
          }
          // Frechet lower bound: P(i and j) >= max(0, Ii + Ij - 1).
          const double lb_raw = std::max(0.0, Ii_raw + Ij_raw - 1.0);
          if (implied < lb_raw * (1.0 - kFrechetRelEps)) c2_frechet_lower.add(w);

          // C2b: same test against the bounded marginals CeeLo actually uses.
          const double implied_b = Ii_bnd * s;
          if (implied_b > std::min(Ii_bnd, Ij_bnd) * (1.0 + kFrechetRelEps))
            c2b_frechet_bounded.add(br * Ii_bnd);

          // Bayes: the reverse conditional implied by this link. If it exceeds 1
          // the two directions cannot share a probability basis -- this is the
          // cleanest single proof that the stored values are not conditionals.
          if (Ij_raw > 0.0 && (implied / Ij_raw) > 1.0 + kOverUnityEps)
            bayes_reverse_gt1.add(w);

          // C3: is the reverse link also stored?
          for (const SandiaDecay::CoincidenceData &back : q.coincidences)
            if (back.first == ip) { ++c3_bidirectional; break; }
        }
      }
    }
  }

  // C4 is only observable against the file: the parser drops links whose id
  // matches no product, so they never reach the object model. The raw pass
  // attributes each lost link to a cause; this is the cross-check that the two
  // passes agree.
  const long long c4_dangling = raw.links_unresolved_id;
  const long long links_lost_total = raw.coincidence_elements - n_links;
  const long long links_lost_unaccounted =
      links_lost_total - raw.links_unresolved_id - raw.links_in_dropped_transitions;

  // ------------------------- report -------------------------
  auto pct = [](long long a, long long b) {
    return b ? (100.0 * static_cast<double>(a) / static_cast<double>(b)) : 0.0;
  };

  std::cout << "=== SandiaDecay invariant census: " << tag << " ===\n"
            << "file                 " << xml_path << "\n"
            << "transitions          " << n_transitions << "  (raw file: " << raw.transition_elements << ")\n"
            << "products             " << n_products << "\n"
            << "gammas               " << n_gammas << "  (raw file: " << raw.gamma_elements << ")\n"
            << "coincidence links    " << n_links << "  (raw file: " << raw.coincidence_elements << ")\n"
            << "  with a2/a4         " << n_links_with_corr
            << "  (" << pct(n_links_with_corr, n_links) << "%)\n"
            << "\n-- statistical counts (reported, never fatal) --\n"
            << "C1_over_unity        " << c1_over_unity.n
            << "   links with P > 1+" << kOverUnityEps << "\n"
            << "C1b_strict           " << c1b_strict.n
            << "   links with P > 1.0 exactly (float-noise band = "
            << (c1b_strict.n - c1_over_unity.n) << ")\n"
            << "C2_frechet           " << c2_frechet.n
            << "   I_i*P(j|i) > min(I_i,I_j), RAW intensity, per transition\n"
            << "C2_frechet_lower     " << c2_frechet_lower.n
            << "   I_i*P(j|i) < max(0, I_i+I_j-1)\n"
            << "C2b_frechet_bounded  " << c2b_frechet_bounded.n
            << "   same test on bounded marginals min(I, 1/(1+alpha))\n"
            << "C3_bidirectional     " << c3_bidirectional
            << "   ordered links whose reverse is also stored\n"
            << "C4_dangling          " << c4_dangling
            << "   links whose id matches no product in a RETAINED transition\n"
            << "C5_occurrence        " << c5_occurrence.n
            << "   gammas with I*(1+alpha) > 1\n"
            << "Bayes_reverse_gt1    " << bayes_reverse_gt1.n
            << "   links whose implied reverse conditional exceeds 1\n";

  std::cout << "\n-- weighted coverage (photons per parent decay) --\n";
  std::cout << "id                   n        >=1e-5   >=1e-4   >=1e-3   sum_w\n";
  auto line = [&](const char *id, const Counter &c) {
    std::printf("%-20s %-8lld %-8lld %-8lld %-8lld %.6g\n",
                id, c.n, c.n_ge[0], c.n_ge[1], c.n_ge[2], c.weight);
  };
  line("C1_over_unity", c1_over_unity);
  line("C2_frechet", c2_frechet);
  line("C2b_frechet_bounded", c2b_frechet_bounded);
  line("C5_occurrence", c5_occurrence);
  line("Bayes_reverse_gt1", bayes_reverse_gt1);

  std::cout << "\n-- records the parser discards (invisible to consumers) --\n"
            << "transitions dropped       " << raw.tr_dropped
            << "  of " << raw.transition_elements << " in file\n"
            << "  no <nuclide> for parent  " << raw.tr_drop_no_parent_record << "\n"
            << "  no branchRatio attribute " << raw.tr_drop_no_branchratio << "\n"
            << "  branchRatio <= 0         " << raw.tr_drop_branchratio_le0 << "\n"
            << "  parent halfLife ABSENT   " << raw.tr_drop_parent_no_halflife
            << "   (defaults to INF => parent treated as stable)\n"
            << "  parent halfLife = INF    " << raw.tr_drop_parent_stable << "\n"
            << "links lost total          " << links_lost_total << "\n"
            << "  in dropped transitions   " << raw.links_in_dropped_transitions << "\n"
            << "  unresolved id (C4)       " << raw.links_unresolved_id << "\n"
            << "  unaccounted              " << links_lost_unaccounted << "\n";

  std::cout << "\n-- MUST-HOLD invariants --\n"
            << "nonfinite values          " << v_nonfinite << "\n"
            << "negative intensities      " << v_negative_intensity << "\n"
            << "partner index out of range" << v_index_out_of_range << "\n"
            << "self links                " << n_self_links << "\n"
            << "|a2| or |a4| > 2          " << v_corr_out_of_range << "\n"
            << "non-descending fl->tl     " << v_nondescending << "\n"
            << "level-feed sum off > " << kFeedSumTol << " " << v_feed_sum
            << "   (worst |sum-1| = " << worst_feed_sum_err << ")\n"
            << "duplicate feed levels     " << v_feed_dup_level << "\n"
            << "verified E0 records       " << n_e0_verified
            << "   (with nonzero stored intensity: " << v_e0_with_photon
            << " -- expected; suppressed at emission)\n"
            << "e0_verified on a non-gamma" << v_e0_not_gamma << "\n";

  if (!csv_dir.empty()) {
    const std::string pre = csv_dir + "/census_" + tag;
    {
      std::ofstream o(pre + "_summary.csv");
      o << "# SandiaDecay invariant census, normative definitions\n"
        << "# file," << xml_path << "\n"
        << "# over_unity_eps," << kOverUnityEps
        << ",frechet_rel_eps," << kFrechetRelEps
        << ",occurrence_eps," << kOccurrenceEps << "\n"
        << "# C1_over_unity: directed links with P(j|i) > 1+over_unity_eps\n"
        << "# C1b_strict: directed links with P(j|i) > 1.0 exactly\n"
        << "# C2_frechet: I_i*P(j|i) > min(I_i,I_j)*(1+eps), RAW intensity, per transition\n"
        << "# C2b_frechet_bounded: as C2 with marginals bounded by min(I,1/(1+alpha))\n"
        << "# C3_bidirectional: ordered links whose reverse link is also present\n"
        << "# C4_dangling: links present in the file but dropped at parse (id unresolved)\n"
        << "# C5_occurrence: gammas with I*(1+alpha) > 1+eps\n"
        << "id,count,n_ge_1e-5,n_ge_1e-4,n_ge_1e-3,sum_weight\n";
      auto row = [&](const char *id, const Counter &c) {
        o << id << "," << c.n << "," << c.n_ge[0] << "," << c.n_ge[1] << ","
          << c.n_ge[2] << "," << c.weight << "\n";
      };
      row("C1_over_unity", c1_over_unity);
      row("C1b_strict", c1b_strict);
      row("C2_frechet", c2_frechet);
      row("C2_frechet_lower", c2_frechet_lower);
      row("C2b_frechet_bounded", c2b_frechet_bounded);
      row("C5_occurrence", c5_occurrence);
      row("Bayes_reverse_gt1", bayes_reverse_gt1);
      o << "C3_bidirectional," << c3_bidirectional << ",,,,\n";
      o << "C4_dangling," << c4_dangling << ",,,,\n";
      o << "links_lost_total," << links_lost_total << ",,,,\n";
      o << "links_in_dropped_transitions," << raw.links_in_dropped_transitions << ",,,,\n";
      o << "links_lost_unaccounted," << links_lost_unaccounted << ",,,,\n";
      o << "transitions_dropped," << raw.tr_dropped << ",,,,\n";
      o << "tr_drop_no_parent_record," << raw.tr_drop_no_parent_record << ",,,,\n";
      o << "tr_drop_no_branchratio," << raw.tr_drop_no_branchratio << ",,,,\n";
      o << "tr_drop_branchratio_le0," << raw.tr_drop_branchratio_le0 << ",,,,\n";
      o << "tr_drop_parent_halflife_absent," << raw.tr_drop_parent_no_halflife << ",,,,\n";
      o << "tr_drop_parent_stable," << raw.tr_drop_parent_stable << ",,,,\n";
      o << "n_links," << n_links << ",,,,\n";
      o << "n_links_raw_file," << raw.coincidence_elements << ",,,,\n";
      o << "n_links_with_a2a4," << n_links_with_corr << ",,,,\n";
      o << "n_gammas," << n_gammas << ",,,,\n";
      o << "n_transitions," << n_transitions << ",,,,\n";
    }
    auto dump_rows = [&](const std::string &path, const std::vector<Row> &rows) {
      std::ofstream o(path);
      o << "parent,child,mode,E_i_keV,E_j_keV,I_i,I_j,P_j_given_i,implied_joint,"
           "frechet_upper_bound,weight_photons_per_decay\n";
      for (const Row &r : rows)
        o << r.parent << "," << r.child << "," << r.mode << "," << r.e_i << ","
          << r.e_j << "," << r.I_i << "," << r.I_j << "," << r.s_ij << ","
          << r.implied << "," << r.bound << "," << r.w << "\n";
    };
    dump_rows(pre + "_C1_over_unity.csv", c1_rows);
    dump_rows(pre + "_C2_frechet.csv", c2_rows);
    {
      std::ofstream o(pre + "_C1_by_parent.csv");
      o << "parent,count,n_ge_1e-5,n_ge_1e-4,n_ge_1e-3,sum_weight\n";
      for (const auto &kv : c1_by_parent)
        o << kv.first << "," << kv.second.n << "," << kv.second.n_ge[0] << ","
          << kv.second.n_ge[1] << "," << kv.second.n_ge[2] << ","
          << kv.second.weight << "\n";
    }
    std::cout << "\nwrote CSVs to " << csv_dir << "\n";
  }

  if (do_assert) {
    const long long fatal = v_nonfinite + v_negative_intensity + v_index_out_of_range
                          + v_corr_out_of_range + v_nondescending + v_feed_sum
                          + v_feed_dup_level + v_e0_not_gamma + v_exclusivity
                          + n_self_links;
    if (fatal) {
      std::cerr << "\nFAILED: " << fatal << " must-hold invariant violation(s)\n";
      return 1;
    }
    std::cout << "\nall must-hold invariants pass\n";
  }
  return 0;
}
