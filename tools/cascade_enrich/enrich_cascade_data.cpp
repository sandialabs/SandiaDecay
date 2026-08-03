/* SandiaDecay cascade-detail enrichment tool.
 *
 * Reads sandia.decay.xml and, using the GEANT4 11.4.0 PhotonEvaporation and
 * RadioactiveDecay level data, injects the extra per-transition / per-particle
 * data needed for high-fidelity true-coincidence-summing:
 *
 *   <gamma ... iccTotal= iccK= iccL1= iccL2= iccL3= fl= tl= e0Verified=>
 *                                                           IC + level topology
 *                                                           + reviewed E0 provenance
 *   <coincidentGamma ... a2= a4=>                    gamma-gamma angular correlation
 *   <electronCapture ... ecK= ecL= lvl=>             EC shell fractions (K, L) + fed level
 *   <positron ... lvl=>                              beta-plus fed level
 *
 * a2/a4 are computed from the nuclear spin sequence + multipole mixing via the
 * F-coefficient (Racah) formalism (see AngularCorrelation). All injected
 * attributes are optional; a database without them parses identically (the
 * CeeLo adapter / InterSpec fall back to the transition-group x-ray model and
 * isotropic emission).
 *
 * The tool is IDEMPOTENT: strip_cascade_detail() first removes any prior cascade
 * attributes + glossary/cascade-match comments, so the committed sandia.decay.xml
 * is the single source of truth and re-running on the tool's own output is a
 * fixed point (verify with a double-run diff).
 *
 * A per-transition "cascade-match:" audit comment records SandiaDecay gammas the
 * structure-aware matcher (StructureMatch) could not map to a GEANT4 transition,
 * bucketed by a per-transition detection threshold (real gap / below-threshold /
 * I=0 / SD-duplicate), plus the GEANT4 branches SandiaDecay lacks.
 *
 * Usage: enrich_cascade_data <in.xml> <out.xml> [photonEvapDir] [radDecayDir]
 *   dirs default to $G4LEVELGAMMADATA / $G4RADIOACTIVEDATA.
 */

#include "AngularCorrelation.h"
#include "G4LevelData.h"
#include "G4RadDecayData.h"
#include "StructureMatch.h"

#include <set>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
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

// Reviewed primary-evaluation provenance for pure electric-monopole
// transitions. This table is deliberately independent of GEANT4's advisory
// multipole code. Each entry must identify the Sandia transition and its
// already-accepted G4 topology tightly enough to fail closed if either data set
// drifts. Additions require a primary evaluated source and a regression test.
struct VerifiedE0Override {
  const char* parent;
  const char* child;
  const char* mode;
  double branch_ratio;
  double sandia_energy_keV;
  double g4_energy_keV;
  int from_level;
  int to_level;
  double min_alpha_total;
};

constexpr VerifiedE0Override kVerifiedE0Overrides[] = {
  // NNDC ENSDF, O-16 Adopted Levels, Gammas, 1993Ti07 (full evaluation by
  // J. H. Kelley et al., evaluated Dec 1992, literature cutoff 31-Dec-1992):
  // the 6049.4-keV 0+ level to 0.0-keV 0+ ground-state transition is listed as
  // 6048.2(10) [E0], I(gamma+ce)=100.
  // https://www.nndc.bnl.gov/ensnds/16/O/adopted.pdf
  {"N16", "O16", "b-", 1.0, 6048.2, 6049.4, 1, 0, 0.0},
  // NNDC ENSDF, Eu-152 beta-minus decay to Gd-152 (NDS 111, 2013):
  // the 615.415-keV 0+ level has a pure 615.41-keV E0 transition to the
  // 0+ ground state, I(gamma+ce)=0.0376(12), with no photon observed.  The
  // ordinary energy matcher otherwise selects an unrelated exactly-615.4-keV
  // high-level E2 edge (149->107) instead of the physical 615.386-keV 2->0
  // sentinel.  Parent-state identity excludes Eu152m's separate 615.44 line.
  // https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=152eu+bM+decay+%2813.517+y%29&nucleus=152GD
  {"Eu152", "Gd152", "b-", 0.279, 615.4, 615.386, 2, 0, 1.0e4},
};

const VerifiedE0Override* verified_e0_identity(
    const char* parent, const char* child, const char* mode, double branch_ratio,
    double sandia_energy_keV) {
  constexpr double kSandiaEnergyToleranceKeV = 0.05;
  constexpr double kBranchTolerance = 1e-12;
  for (const VerifiedE0Override& item : kVerifiedE0Overrides) {
    if (std::string(parent) == item.parent && std::string(child) == item.child &&
        std::string(mode) == item.mode &&
        std::fabs(branch_ratio - item.branch_ratio) <= kBranchTolerance &&
        std::fabs(sandia_energy_keV - item.sandia_energy_keV) <=
            kSandiaEnergyToleranceKeV)
      return &item;
  }
  return nullptr;
}

bool matches_verified_e0_edge(const VerifiedE0Override& item,
                              const G4Transition& g4) {
  constexpr double kG4EnergyToleranceKeV = 0.05;
  if (std::fabs(g4.gamma_keV - item.g4_energy_keV) >
          kG4EnergyToleranceKeV ||
      g4.from_level != item.from_level || g4.to_level != item.to_level ||
      g4.multipole_code != 1)
    return false;
  if (item.min_alpha_total <= 0.0)
    return true;
  if (!g4.has_icc || !std::isfinite(g4.alpha_total) ||
      g4.alpha_total < item.min_alpha_total)
    return false;
  double shell_sum = 0.0;
  for (double shell : g4.icc_shell) {
    if (!std::isfinite(shell) || shell < 0.0)
      return false;
    shell_sum += shell;
  }
  return g4.icc_shell[0] > 0.0 && shell_sum > 0.999 && shell_sum < 1.001;
}

// Reviewed correction for a demonstrable percent-vs-fraction transcription in
// the base SandiaDecay gamma intensity. Like the E0 table above, every field is
// part of the identity gate: a changed parent branch, photon energy, G4 edge or
// ICC fails closed instead of modifying an unrelated line.
struct GammaIntensityOverride {
  const char* parent;
  const char* child;
  const char* mode;
  double branch_ratio;
  double sandia_energy_keV;
  double original_intensity;
  double corrected_intensity;
  double g4_energy_keV;
  double g4_icc_total;
  int from_level;
  int to_level;
};

constexpr GammaIntensityOverride kGammaIntensityOverrides[] = {
  // NNDC ENSDF 235U alpha decay evaluation, NDS 185, 560 (2022), evaluated
  // Aug-2022 with literature cutoff 31-Aug-2022. The 205.3098-keV Th-231 level
  // lists G 19.55(5), %IG=0.583(11), [M1], alpha=114.7(19). `%IG` is percent per
  // parent decay, hence the transition-conditional fraction is 0.00583, not
  // 0.61. The photon energy remains the existing SandiaDecay 19.59 keV.
  // https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=235u+a+decay+%287.04eP8+y%29&nucleus=231TH
  {"U235", "Th231", "a", 1.0, 19.59, 0.61, 0.00583,
   19.5941, 114.7, 5, 4},
  // Same evaluation (NNDC ENSDF 235U alpha decay, NDS 185, 560 (2022), cutoff
  // 31-Aug-2022). The 452.189-keV Th-231 level lists G 64.45(5), %IG=0.0017(14),
  // [M1+E2], alpha=52, to the 387.824-keV level. Per parent decay that is
  // 1.7e-5(1.4e-5); SandiaDecay carried 0.0065, a factor 382 high, unchanged
  // since the 2019 InterSpec import and with no recorded provenance. The
  // conversion coefficient already agreed exactly (52), and the G4 edge
  // (level 30 -> 26, 64.364 keV) matches the ENSDF level pair to 0.004 keV, so
  // only the photon intensity was wrong.
  //
  // The stated ENSDF uncertainty is 82% of the value, so this is a weakly
  // determined line -- but the correction is a factor 382, far outside it.
  // https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=235u+a+decay+%287.04eP8+y%29&nucleus=231TH
  {"U235", "Th231", "a", 1.0, 64.35, 0.0065, 1.7e-5,
   64.364, 52.0, 30, 26},
};

// REVIEWED REMOVAL, applied directly to the delivered XML (2026-08-03).
//
// This one is recorded here rather than mechanised. A removal cannot use the
// override table's "fires exactly once" gate, because the delivered XML is both
// this tool's input and its output: once the record is gone, a regeneration
// finds nothing to remove. Deleting it from the delivered file is sufficient --
// enrichment is idempotent on that file, and the intensity is data this tool
// never generates. Re-enriching an OLDER base (e.g. the upstream sandialabs
// file) would reintroduce it; that is a Stage-6 concern, where the ENSDF join
// replaces energy-proximity matching altogether and can drop unsupported
// records systematically instead of one at a time.
//
// Four U235 -> Th231 records total were found unsupported. 64.35 keV is an
// intensity correction and lives in the table above; the other three are
// deletions, recorded here:
//
// 26.55 keV, intensity 0.0054, alpha 192, occurrence 1.042/decay.
//   - NNDC ENSDF 235U alpha decay (NDS 185, 560 (2022), cutoff 31-Aug-2022)
//     lists NO gamma at 26.55 keV; below 45 keV it has only 19.55, 31.6, 34.7,
//     41.1 and 42.02 (plus an unplaced 15.784).
//   - The record carried intensity 0.0054 with alpha=192, i.e. an occurrence of
//     1.042 per decay: one record consuming an entire decay's worth of flux,
//     which by itself rejected the whole U-235 level scheme.
//   - It was matched to the G4 edge level 11 -> 6 (247.5867 -> 221.3966 keV),
//     whose energy is 26.19 keV -- 0.36 keV from the Sandia energy, a loose
//     match. That nuclear transition is real and heavily converted (alpha=192,
//     so only 1 photon per 193 transitions); what is unsupported is a
//     detectable PHOTON line at this energy and intensity.
//   - Microcalorimeter spectra examined by the maintainer show no evidence of
//     this line. The value is unchanged since the 2019 InterSpec import with no
//     recorded provenance.
//   If the topology edge is ever needed, it should return as a 26.19 keV record
//   with an evaluated intensity -- not as this line.
//
// 38.90 keV, intensity 0.0011, alpha 59.5, occurrence 0.0665/decay.
//   - ENSDF lists no gamma at this energy in the 235U alpha decay at all.
//   - GEANT4 has the corresponding edge (level 44 -> 40, 593.6173 -> 554.6503,
//     38.967 keV, alpha 59.5), so the transition exists; the photon intensity
//     does not come from any evaluation. Its iccL1/L2/L3 were all zero while
//     iccTotal was 59.5, which cannot be right at Z=90 and is a further sign
//     the record was not evaluated.
//   - Microcalorimeter spectra show no such line.
//
// 41.40 keV, intensity 0.0003, alpha 51, occurrence 0.0156/decay.
//   - Corresponds to ENSDF 41.100 keV (277.60 -> 236.893, cc=51.0), whose
//     intensity field is BLANK -- unmeasured, which is NOT the same as an
//     evaluated zero. GEANT4 independently carries the edge (level 14 -> 8,
//     40.701 keV) with relative intensity 0.
//   - Both are consistent with no observable photon, and microcalorimeter
//     spectra show none where one would be expected. SandiaDecay's 3e-4 has no
//     traceable source.
//
// REVIEWED, NOT CHANGED (2026-08-03). Recorded so the next reader does not
// re-derive the same analysis:
//
// U235 -> Th231, 182.504 keV (387.824 -> 205.3098). SandiaDecay carries 0.340
//   per 100 decays; ENSDF gives 0.411(15) -- a 17% gap at roughly 5 sigma, the
//   largest remaining amplitude disagreement in this branch. The maintainer's
//   microcalorimeter fit, normalized against the 185.71 keV line, lands near
//   92% of the ENSDF value (~0.378), i.e. between the two. No change was made:
//   0.378 is not an evaluated number, and adopting it would mean sourcing an
//   intensity to a local measurement rather than to an evaluation. Revisit if
//   that measurement is written up, and record it as measurement provenance
//   rather than as ENSDF.
//
// Together these four moved the branch's level-scheme feed sum from 2.37344 to
// 1.00130 -- i.e. to within 0.13% of the exact physical bound of 1, which is
// itself strong evidence that they were the defective records and that the
// remaining ~50 are mutually consistent. The branch is now ACCEPTED and uses
// the exact level walk instead of the bounded pairwise fallback.

const GammaIntensityOverride* gamma_intensity_override(
    const char* parent, const char* child, const char* mode, double branch_ratio,
    double sandia_energy_keV, const G4Transition& g4) {
  constexpr double kSandiaEnergyToleranceKeV = 0.005;
  constexpr double kG4EnergyToleranceKeV = 0.005;
  constexpr double kIccTolerance = 1.0e-6;
  constexpr double kBranchTolerance = 1.0e-12;
  for (const GammaIntensityOverride& item : kGammaIntensityOverrides) {
    if (std::string(parent) == item.parent && std::string(child) == item.child &&
        std::string(mode) == item.mode &&
        std::fabs(branch_ratio - item.branch_ratio) <= kBranchTolerance &&
        std::fabs(sandia_energy_keV - item.sandia_energy_keV) <=
            kSandiaEnergyToleranceKeV &&
        std::fabs(g4.gamma_keV - item.g4_energy_keV) <= kG4EnergyToleranceKeV &&
        std::fabs(g4.alpha_total - item.g4_icc_total) <= kIccTolerance &&
        g4.from_level == item.from_level && g4.to_level == item.to_level)
      return &item;
  }
  return nullptr;
}

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

void set_attr_precise(rapidxml::xml_document<char>& doc, xml_node<char>* node,
                      const char* name, double value) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.10G", value);
  const char* val = doc.allocate_string(buf);
  xml_attribute<char>* existing = node->first_attribute(name);
  if (existing)
    existing->value(val);
  else
    node->append_attribute(doc.allocate_attribute(doc.allocate_string(name), val));
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

void set_attr_string(rapidxml::xml_document<char>& doc, xml_node<char>* node,
                     const char* name, const char* value) {
  const char* val = doc.allocate_string(value);
  if (node->first_attribute(name))
    node->first_attribute(name)->value(val);
  else
    node->append_attribute(
        doc.allocate_attribute(doc.allocate_string(name), val));
}

// Remove an attribute by name if present.
void del_attr(xml_node<char>* node, const char* name) {
  if (xml_attribute<char>* a = node->first_attribute(name))
    node->remove_attribute(a);
}

// Remove every cascade-detail artifact this tool emits, so the enricher is
// IDEMPOTENT: running it on an already-enriched file yields the same result as
// running on the clean base. Attributes must be REMOVED (not merely overwritten):
// if a gamma matched on a previous run (got fl/tl) but is not accepted this run,
// an overwrite-only pass would leave stale fl/tl that silently mis-route cascade
// flux. Both comment kinds are also removed (root glossary + per-transition
// "cascade-match:" notes) so they are not duplicated. The attribute list must
// stay in sync with the set_attr/set_attr_int calls below and with
// minimize_xml_file.cpp's remove_cascade_detail_attribs().
void strip_cascade_detail(rapidxml::xml_document<char>& /*doc*/, xml_node<char>* root) {
  // Root-level glossary comment(s).
  std::vector<xml_node<char>*> root_comments;
  for (xml_node<char>* c = root->first_node(); c; c = c->next_sibling())
    if (c->type() == rapidxml::node_comment &&
        std::string(c->value(), c->value_size())
            .find("field glossary and cascade-match audit legend") != std::string::npos)
      root_comments.push_back(c);
  for (xml_node<char>* c : root_comments) root->remove_node(c);

  // Per-transition cascade attributes + per-transition cascade-match comments.
  for (xml_node<char>* tr = root->first_node("transition"); tr;
       tr = tr->next_sibling("transition")) {
    std::vector<xml_node<char>*> old_feeds;
    for (xml_node<char>* f = tr->first_node("levelFeed"); f;
         f = f->next_sibling("levelFeed"))
      old_feeds.push_back(f);
    for (xml_node<char>* f : old_feeds) tr->remove_node(f);
    for (xml_node<char>* g = tr->first_node("gamma"); g;
         g = g->next_sibling("gamma")) {
      for (const char* a : {"fl", "tl", "iccTotal", "iccK", "iccL1", "iccL2",
                            "iccL3", "mult", "delta", "e0Verified"})
        del_attr(g, a);
      for (xml_node<char>* cg = g->first_node("coincidentGamma"); cg;
           cg = cg->next_sibling("coincidentGamma")) {
        del_attr(cg, "a2");
        del_attr(cg, "a4");
      }
    }
    for (xml_node<char>* e = tr->first_node("electronCapture"); e;
         e = e->next_sibling("electronCapture"))
      for (const char* a : {"ecK", "ecL", "ecM", "ecN", "lvl"})
        del_attr(e, a);
    for (xml_node<char>* p = tr->first_node("positron"); p;
         p = p->next_sibling("positron"))
      del_attr(p, "lvl");
    std::vector<xml_node<char>*> tr_comments;
    for (xml_node<char>* c = tr->first_node(); c; c = c->next_sibling())
      if (c->type() == rapidxml::node_comment) {
        const std::string text(c->value(), c->value_size());
        if (text.find("cascade-match:") != std::string::npos ||
            text.find("cascade-feed:") != std::string::npos ||
            text.find("intensity-fix:") != std::string::npos)
          tr_comments.push_back(c);
      }
    for (xml_node<char>* c : tr_comments) tr->remove_node(c);
  }
}

// Human-readable header comment: what the abbreviated attributes mean and how to
// read the per-transition "cascade-match:" audit comments. Aimed at readers who
// are not nuclear-data experts. NOTE: XML comments may not contain a double
// hyphen, so this text uses "=" rule lines and only single hyphens.
const char* const kHeaderDoc =
"\n"
"==============================================================================\n"
" SandiaDecay database - field glossary and cascade-match audit legend\n"
"==============================================================================\n"
" Most attributes below are standard nuclear-decay quantities. The fields marked\n"
" [C] are \"cascade detail\", added from GEANT4 11.4.0 evaluated nuclear data\n"
" (PhotonEvaporation 6.1.2 and RadioactiveDecay 6.1.2, both derived from ENSDF)\n"
" to support true-coincidence (cascade) summing and the atomic x-ray emitted\n"
" after internal conversion or electron capture. An attribute is OMITTED (not\n"
" set to zero) when it is unknown or does not apply.\n"
"\n"
" Nuclide / transition:\n"
"   atomicNumber    Z (proton number).       massNumber   A (nucleon number).\n"
"   isomerNumber    0 = ground state, 1 = first metastable isomer, ...\n"
"   halfLife        half-life in seconds (\"INF\" = stable).\n"
"   branchRatio     fraction of parent decays that take this transition.\n"
"   productIntensityBasis  optional closed-enum provenance for reviewed\n"
"                   EC/beta-plus/x-ray product-intensity vectors.\n"
"                   transitionConditional means each\n"
"                   intensity is conditional on this transition, as required\n"
"                   by the SandiaDecay API; multiply by branchRatio once.\n"
"   mode            decay mode (e.g. b-, b+, ec, alpha, it).\n"
"   levelFeed       [C] complete daughter-level population law for an alpha or\n"
"                   beta-minus transition. lvl is the exact PhotonEvaporation\n"
"                   level index and probability is conditional on this decay\n"
"                   transition; q is the RDM final Q (the beta endpoint, but\n"
"                   total decay Q rather than emitted-alpha kinetic energy), and\n"
"                   forbiddenness carries the evaluated beta-spectrum shape\n"
"                   (absent = allowed). Entries include ground level 0 and sum\n"
"                   to one.\n"
"                   The entire list is omitted on any identity/closure failure.\n"
"\n"
" Emitted-particle energy / intensity (on beta, positron, electronCapture,\n"
" alpha, gamma, xray):\n"
"   energy          keV.\n"
"   intensity       emission probability conditional on this decay transition;\n"
"                   multiply by branchRatio for probability per parent decay.\n"
"   relintensity    (xray) intensity relative to that atom's strongest x-ray.\n"
"   logFT           log10 of the beta/EC comparative half-life (ft); a smaller\n"
"                   value means a more allowed (faster) beta/EC branch.\n"
"   forbiddenness   beta/EC forbiddenness class (\"\" = allowed, \"1U\" = first\n"
"                   forbidden unique, ...).\n"
"\n"
" [C] on <gamma> (internal conversion, IC, competes with gamma emission):\n"
"   iccTotal        total IC coefficient alpha = (IC electrons)/(emitted gammas).\n"
"   iccK            partial IC coefficient for a K-shell conversion.\n"
"   iccL1,iccL2,iccL3   partial IC coefficients for the L1/L2/L3 subshells,\n"
"                   written only for progeny whose L x-rays exceed 10 keV\n"
"                   (Z >= 69).\n"
"   mult            Advisory GEANT4 multipolarity code copied verbatim: 1..9 =\n"
"                   E0,E1,M1,E2,M2,E3,M3,E4,M4; 100*Nx+Ny = mixed. It is not\n"
"                   primary-evaluation provenance: PhotonEvaporation contains\n"
"                   mult=\"1\" candidates with zero or finite iccTotal attached\n"
"                   to tabulated gamma lines (for example Ac-232 373.3 keV and\n"
"                   Rb-78 189.54 keV). Do not use mult=\"1\" alone to suppress a\n"
"                   photon or synthesize an electron. GEANT4 itself chooses\n"
"                   gamma versus conversion from iccTotal. A pure-E0 decision\n"
"                   needs the energy-aware ICC-sentinel convention or separately\n"
"                   sourced evaluated provenance; this enricher supplies neither\n"
"                   claim through mult.\n"
"   delta           multipole mixing ratio, written when nonzero.\n"
"   e0Verified      1 only for a pure-E0 transition verified against a reviewed\n"
"                   primary nuclear-data evaluation. It is independent of mult\n"
"                   and absent otherwise. Curated records are O-16 6048.2 keV\n"
"                   [E0] (NNDC ENSDF 1993Ti07) and Gd-152 615.4 keV [E0]\n"
"                   following Eu-152 beta-minus decay (NDS 111, 2013).\n"
"   fl, tl          \"from level\" and \"to level\": indices of the nuclear levels\n"
"                   (in the GEANT4 progeny level scheme) that this gamma\n"
"                   connects. Together they give the cascade topology, i.e. which\n"
"                   gammas follow one another within a decay.\n"
"\n"
" [C] on <electronCapture> and <positron>:\n"
"   ecK, ecL        fraction of electron captures from the K / L atomic shell.\n"
"   lvl             index of the progeny level fed by this individual EC or\n"
"                   beta-plus branch. A positron is matched by endpoint + 2*m_e\n"
"                   to a real per-level GEANT4 BetaPlus Q record. Summary records\n"
"                   never supply lvl because they do not identify a level.\n"
"\n"
" [C] on <coincidentGamma> (\"id\" = the partner gamma emitted in the same decay):\n"
"   a2, a4          gamma-gamma angular-correlation coefficients in\n"
"                   W(theta) = 1 + a2*P2(cos theta) + a4*P4(cos theta). Non-zero\n"
"                   means the two gammas prefer a correlated relative angle, used\n"
"                   when summing coincident gammas; a2 = a4 = 0 means isotropic.\n"
"\n"
"==============================================================================\n"
" Per-transition \"cascade-match:\" audit comments.\n"
" The cascade detail above was produced by adding information from GEANT4 to the\n"
" existing SandiaDecay data: each SandiaDecay transition's gammas were matched to\n"
" the GEANT4 level scheme of the progeny nuclide. That matching is not always\n"
" complete, and the leftover cases fall into two kinds: GEANT4 had a transition\n"
" that no SandiaDecay gamma accounted for, or SandiaDecay had a gamma that GEANT4\n"
" had no matching transition for. Where either happened, a note is recorded on\n"
" the <transition> element. These notes do not change the values SandiaDecay\n"
" returns; they are kept for tracking and to guide possible future improvements.\n"
" A note may contain any of these sections:\n"
"\n"
"   SD-unmatched-above-threshold[keV]: E(I=i,reason)\n"
"       A SandiaDecay gamma (energy E keV, transition-conditional intensity i)\n"
"       ABOVE this transition's detection threshold T (see below) that was NOT\n"
"       matched to a GEANT4 nuclear transition, so it carries no cascade detail\n"
"       (no iccX / fl / tl / a2 / a4). It is still emitted normally; it is just\n"
"       treated as an isolated, isotropic gamma by the cascade model. These are\n"
"       the genuine, potentially-detectable gaps. reason is one of:\n"
"         absent    no GEANT4 transition within 3 keV of E (often a SandiaDecay\n"
"                   line, for example an x-ray, that GEANT4 does not model as a\n"
"                   nuclear transition; it is NOT force-matched).\n"
"         neartol   a GEANT4 transition lies within 3 keV but the structure-aware\n"
"                   matcher rejected it (ambiguous or inconsistent levels).\n"
"         edge-collision\n"
"                   two or more otherwise-accepted SandiaDecay records selected\n"
"                   the same GEANT4 edge. All owners were rejected: one G4 edge\n"
"                   is not evidence that either Sandia record is a duplicate.\n"
"         annih511  the 511 keV annihilation line (not a nuclear transition).\n"
"\n"
"   N SD gammas below detectable threshold T=t/decay: E E ...\n"
"       Unmatched gammas whose intensity is below the per-transition detection\n"
"       threshold T = max(1e-5/decay, 0.1% of the strongest gamma in the\n"
"       transition). These are too weak to detect in practice, so they are noted\n"
"       as a collapsed list rather than flagged as real gaps.\n"
"\n"
"   N SD gammas with I=0 (not emitted, suppressed)\n"
"       Unmatched gammas carrying zero intensity (not emitted at all); counted\n"
"       only, not listed.\n"
"\n"
"   G4-listed-not-in-SD[keV]: E(rI=r)\n"
"       A GEANT4 transition (energy E keV, RELATIVE intensity r on GEANT4's\n"
"       per-level 0 to 100 scale) on a level that SandiaDecay does populate, but\n"
"       which no SandiaDecay gamma claimed: a branch GEANT4 lists that SandiaDecay\n"
"       does not have as a gamma. Only branches with rI >= 5 are shown. This flags\n"
"       a SandiaDecay data gap (documented only; no gamma is added).\n"
"\n"
"   (no G4 level file for progeny)\n"
"       GEANT4 has no level file for this progeny, so no cascade detail could be\n"
"       derived for the transition.\n"
"\n"
" Each explicit list is truncated at 25 entries (\"...\"). These audit comments\n"
" appear only in this human-readable file; they are removed from the .min files.\n"
" A \"cascade-feed: omitted(reason)\" comment records why an alpha/beta-minus\n"
" transition has no levelFeed law. An \"intensity-fix:\" comment records a\n"
" narrowly gated primary-evaluation correction applied by this enricher.\n"
"==============================================================================\n";

// Insert kHeaderDoc as the first node of the document root, so the file opens
// with the glossary. Idempotency is provided by strip_cascade_detail(), which is
// called first to remove any prior glossary/cascade-match comments and cascade
// attributes, so re-running the tool on its own output is a fixed point.
void insert_header_doc(rapidxml::xml_document<char>& doc, xml_node<char>* root) {
  if (std::string(kHeaderDoc).find("--") != std::string::npos)
    throw std::runtime_error("Generated XML comment contains forbidden '--'");
  xml_node<char>* c = doc.allocate_node(rapidxml::node_comment, nullptr,
                                        doc.allocate_string(kHeaderDoc));
  root->insert_node(root->first_node(), c);
}

long zakey(int Z, int A) { return (static_cast<long>(Z) << 16) | A; }

std::uint64_t za_isomer_key(int Z, int A, int isomer) {
  return (static_cast<std::uint64_t>(static_cast<unsigned>(Z)) << 32) |
         (static_cast<std::uint64_t>(static_cast<unsigned>(A)) << 16) |
         static_cast<unsigned short>(isomer);
}

struct NuclideAddress {
  int Z = 0;
  int A = 0;
  int isomer = 0;
  double atomic_mass = 0.0;
  double half_life_s = 0.0;
};

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

// Exact identity resolver for transition-level decay feeds. A floating-level
// flag is not an annotation: it distinguishes physically different levels that
// can have identical excitation energies (Pa-234 has both '-' and '+X' at
// 73.92 keV). Therefore no nearest-energy fallback is permitted here.
int exact_level_index(const G4LevelScheme& s, double ex_keV,
                      const std::string& floating_flag,
                      double tol = 0.01) {
  int found = -1;
  for (const G4Level& l : s.levels) {
    if (l.floating_flag != floating_flag ||
        std::fabs(l.excitation_keV - ex_keV) > tol)
      continue;
    if (found >= 0)
      return -1;  // ambiguous exact identity: fail closed
    found = l.index;
  }
  return found;
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

    // Idempotency: strip any prior cascade attributes + glossary/cascade-match
    // comments so re-running the tool on its own output is a fixed point (and no
    // stale fl/tl survive a match that is no longer accepted this run).
    strip_cascade_detail(doc, root);

    // Prepend the human-readable field glossary / cascade-match legend.
    insert_header_doc(doc, root);

    // Symbol -> (Z,A,isomer) from <nuclide> elements. RadioactiveDecay stores
    // a separate P block for each parent state, so dropping isomerNumber here
    // silently assigns ground-state EC/beta-plus data to metastable parents.
    std::map<std::string, NuclideAddress> sym_to_nuclide;
    for (xml_node<char>* n = root->first_node("nuclide"); n;
         n = n->next_sibling("nuclide")) {
      xml_attribute<char>* sa = n->first_attribute("symbol");
      xml_attribute<char>* za = n->first_attribute("atomicNumber");
      xml_attribute<char>* ma = n->first_attribute("massNumber");
      xml_attribute<char>* ia = n->first_attribute("isomerNumber");
      xml_attribute<char>* aa = n->first_attribute("atomicMass");
      xml_attribute<char>* ha = n->first_attribute("halfLife");
      if (sa && za && ma) {
        NuclideAddress a;
        a.Z = std::atoi(za->value());
        a.A = std::atoi(ma->value());
        a.isomer = ia ? std::atoi(ia->value()) : 0;
        a.atomic_mass = aa ? std::atof(aa->value()) : 0.0;
        a.half_life_s = ha ? std::atof(ha->value()) : 0.0;
        sym_to_nuclide[sa->value()] = a;
      }
    }

    std::map<long, double> ground_atomic_mass;
    std::map<long, double> ground_half_life;
    for (const auto& entry : sym_to_nuclide) {
      const NuclideAddress& a = entry.second;
      if (a.isomer == 0) {
        if (std::isfinite(a.atomic_mass) && a.atomic_mass > 0.0)
          ground_atomic_mass[zakey(a.Z, a.A)] = a.atomic_mass;
        if (std::isfinite(a.half_life_s) && a.half_life_s > 0.0)
          ground_half_life[zakey(a.Z, a.A)] = a.half_life_s;
      }
    }

    std::map<long, G4LevelScheme> level_cache;
    auto get_scheme = [&](int Z, int A) -> const G4LevelScheme& {
      const long k = zakey(Z, A);
      auto it = level_cache.find(k);
      if (it == level_cache.end())
        it = level_cache.emplace(k, load_photon_evaporation(Z, A, photoneva)).first;
      return it->second;
    };
    std::map<std::uint64_t, G4ECData> ec_cache;
    auto get_ec = [&](const NuclideAddress& a) -> const G4ECData& {
      const std::uint64_t k = za_isomer_key(a.Z, a.A, a.isomer);
      auto it = ec_cache.find(k);
      if (it == ec_cache.end()) {
        ParentStateHint hint;
        hint.isomer = a.isomer;
        hint.half_life_s = a.half_life_s;
        if (a.isomer == 0)
          hint.excitation_keV = 0.0;
        if (std::isfinite(a.half_life_s) && a.half_life_s > 0.0) {
          for (const auto& sibling_entry : sym_to_nuclide) {
            const NuclideAddress& sibling = sibling_entry.second;
            if (sibling.Z == a.Z && sibling.A == a.A &&
                sibling.isomer != a.isomer && std::isfinite(sibling.half_life_s) &&
                sibling.half_life_s > 0.0 &&
                std::fabs(std::log(sibling.half_life_s / a.half_life_s)) <=
                    std::log(1.05)) {
              hint.half_life_identity_unique = false;
              break;
            }
          }
        }
        const auto gh = ground_half_life.find(zakey(a.Z, a.A));
        if (gh != ground_half_life.end())
          hint.ground_half_life_s = gh->second;
        const auto gm = ground_atomic_mass.find(zakey(a.Z, a.A));
        if (a.isomer > 0 && gm != ground_atomic_mass.end() &&
            std::isfinite(a.atomic_mass) && a.atomic_mass > gm->second) {
          constexpr double kAtomicMassUnitKeV = 931494.10242;
          hint.excitation_keV = (a.atomic_mass - gm->second) * kAtomicMassUnitKeV;
        }
        it = ec_cache.emplace(k, load_ec_shell_fractions(a.Z, a.A, raddecay, hint)).first;
      }
      return it->second;
    };

    // A RadioactiveDecay P block describes the complete feed law for one decay
    // mode. It can be attached directly only when SandiaDecay has exactly one
    // transition of that mode from the parent. Multiple child/isomer
    // transitions require a separately proven terminal-basin partition; this
    // tool deliberately omits them rather than guessing.
    std::map<std::pair<std::string, std::string>, int> same_mode_count;
    for (xml_node<char>* tr = root->first_node("transition"); tr;
         tr = tr->next_sibling("transition")) {
      xml_attribute<char>* pa = tr->first_attribute("parent");
      xml_attribute<char>* ma = tr->first_attribute("mode");
      xml_attribute<char>* ca = tr->first_attribute("child");
      if (!pa || !ma || !ca) continue;
      const std::string mode = ma->value();
      if (mode == "a" || mode == "b-")
        ++same_mode_count[{pa->value(), mode}];
    }

    long n_trans = 0, n_gamma = 0, n_accepted = 0, n_collision_rejected = 0,
         n_struct_changed = 0,
         n_struct_used = 0, n_ambiguous = 0, n_coinc = 0, n_coinc_a2a4 = 0,
         n_icc = 0, n_ec = 0, n_ec_shell = 0, n_ec_level = 0,
         n_ec_summary = 0, n_ec_q_ambiguous = 0,
         n_positron = 0, n_positron_level = 0,
         n_positron_q_ambiguous = 0,
         n_ec_parent_missing = 0, n_positron_parent_missing = 0,
         n_parent_block_missing = 0, n_doc = 0, n_mult = 0,
         n_e0_verified = 0, n_e0_parent_branch = 0, n_n16_e0_verified = 0,
         n_eu152_e0_parent_branch = 0, n_eu152_e0_verified = 0,
         n_eu152_weak_parent_branch = 0, n_eu152_weak_converted = 0,
         n_intensity_fixed = 0, n_u235_alpha_feed_fixed = 0,
         n_u235_parent_branch = 0,
         n_th234_gamma_normalized = 0, n_th234_parent_branch = 0,
         n_u238_alpha_feed_fixed = 0, n_u238_parent_branch = 0,
         n_level_feed_transitions = 0, n_level_feed_records = 0,
         n_level_feed_omitted = 0;
    std::map<std::string, long> level_feed_omit_reasons;

    for (xml_node<char>* tr = root->first_node("transition"); tr;
         tr = tr->next_sibling("transition")) {
      ++n_trans;
      xml_attribute<char>* pa = tr->first_attribute("parent");
      xml_attribute<char>* ca = tr->first_attribute("child");
      if (!pa || !ca) continue;
      xml_attribute<char>* ma = tr->first_attribute("mode");
      bool hbr = false;
      const double branch_ratio = attr_double(tr, "branchRatio", hbr);
      const char* mode = ma ? ma->value() : "";
      const bool n16_o16_branch = std::string(pa->value()) == "N16" &&
          std::string(ca->value()) == "O16";
      if (n16_o16_branch)
        ++n_e0_parent_branch;
      const bool eu152_gd_branch = std::string(pa->value()) == "Eu152" &&
          std::string(ca->value()) == "Gd152" && std::string(mode) == "b-";
      if (eu152_gd_branch)
        ++n_eu152_e0_parent_branch;
      const bool eu152_sm_branch = std::string(pa->value()) == "Eu152" &&
          std::string(ca->value()) == "Sm152" && std::string(mode) == "ec";
      if (eu152_sm_branch)
        ++n_eu152_weak_parent_branch;
      const bool u235_alpha_branch = std::string(pa->value()) == "U235" &&
          std::string(ca->value()) == "Th231" && std::string(mode) == "a";
      if (u235_alpha_branch)
        ++n_u235_parent_branch;
      if (std::string(pa->value()) == "Th234" &&
          std::string(ca->value()) == "Pa234m" && std::string(mode) == "b-")
        ++n_th234_parent_branch;
      const bool u238_alpha_branch = std::string(pa->value()) == "U238" &&
          std::string(ca->value()) == "Th234" && std::string(mode) == "a";
      if (u238_alpha_branch)
        ++n_u238_parent_branch;
      auto cit = sym_to_nuclide.find(ca->value());
      if (cit == sym_to_nuclide.end()) continue;
      const int daughter_Z = cit->second.Z;
      const G4LevelScheme& scheme = get_scheme(cit->second.Z, cit->second.A);

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

      // SandiaDecay's public U-235 alpha products predate the Aug-2022 ENSDF
      // evaluation and are not a coherent daughter-level vector. Replace the
      // entire vector only after both it and the independently parsed RDM feed
      // law have passed exact identity gates. The target is deliberately the
      // G4-topology-backed normalized projection of ENSDF, not the complete
      // ENSDF alpha table; the provenance comment below records the three
      // tentative/unmapped rows that cannot be represented by that feed law.
      std::vector<xml_node<char>*> u235_alpha_nodes;
      xml_node<char>* u235_alpha_insert_before = nullptr;
      if (u235_alpha_branch) {
        if (!hbr || std::fabs(branch_ratio - 1.0) > 1.0e-12)
          throw std::runtime_error(
              "Reviewed U235->Th231 alpha projection saw unexpected branch ratio");
        constexpr double legacy_energy[] = {
            4044.656, 4148.7, 4150.0, 4214.7, 4219.0, 4271.0, 4293.01,
            4295.0, 4353.787, 4361.618, 4366.1, 4397.8, 4400.9, 4403.275,
            4414.0, 4435.0, 4437.819, 4457.302, 4502.0, 4556.0, 4596.4};
        constexpr double legacy_probability[] = {
            0.00001, 0.0002, 0.009, 0.057, 0.009, 0.004, 0.00017,
            0.00009, 0.00033, 0.001, 0.17, 0.55, 0.044, 0.0021,
            0.021, 0.007, 0.0024, 0.0031, 0.017, 0.042, 0.05};
        constexpr double legacy_hindrance[] = {
            83.0, 34.0, 3.6, 2.0, 15.0, 83.0, 696.0, 6400.0, 3000.0,
            430.0, 11.0, 6.1, 20.0, 447.0, 226.0, 1000.0, 733.0,
            806.0, 1340.0, 1370.0, 2340.0};
        constexpr double current_energy[] = {
            3976.0, 4080.6, 4152.5, 4216.0, 4219.0, 4248.0, 4266.2,
            4273.0, 4282.6, 4323.0, 4323.2, 4357.0, 4364.6, 4377.0,
            4395.7, 4415.1, 4438.3, 4502.9, 4556.5, 4597.7};
        constexpr double current_probability[] = {
            0.00006915006016, 0.0001718901495, 0.002933902552,
            0.06025905243, 0.001215101057, 0.0006915006016,
            0.002272101977, 0.0003259902836, 0.001086600945,
            0.001649701435, 0.03536503077, 0.002074501805,
            0.1871001628, 0.003457503008, 0.5708804967,
            0.03052502656, 0.002420202106, 0.012645011,
            0.03773603283, 0.047121041};
        constexpr double current_hindrance[] = {
            12.0, 41.8, 11.5, 1.94, 117.0, 320.0, 150.0, 1200.0,
            410.0, 580.0, 26.0, 870.0, 10.4, 740.0, 6.02, 160.0,
            3100.0, 1850.0, 1570.0, 2550.0};
        constexpr std::size_t n_legacy =
            sizeof(legacy_energy) / sizeof(legacy_energy[0]);
        constexpr std::size_t n_current =
            sizeof(current_energy) / sizeof(current_energy[0]);
        for (xml_node<char>* alpha = tr->first_node("alpha"); alpha;
             alpha = alpha->next_sibling("alpha"))
          u235_alpha_nodes.push_back(alpha);
        if (!u235_alpha_nodes.empty())
          u235_alpha_insert_before = u235_alpha_nodes.back()->next_sibling();
        bool all_legacy = u235_alpha_nodes.size() == n_legacy;
        bool all_current = u235_alpha_nodes.size() == n_current;
        for (std::size_t i = 0; i < u235_alpha_nodes.size(); ++i) {
          bool have_energy = false, have_probability = false,
               have_hindrance = false;
          const double energy =
              attr_double(u235_alpha_nodes[i], "energy", have_energy);
          const double probability =
              attr_double(u235_alpha_nodes[i], "intensity", have_probability);
          const double hindrance =
              attr_double(u235_alpha_nodes[i], "hindrance", have_hindrance);
          if (!have_energy || !have_probability || !have_hindrance)
            throw std::runtime_error(
                "Reviewed U235->Th231 alpha projection saw missing attributes");
          if (all_legacy)
            all_legacy =
                std::fabs(energy - legacy_energy[i]) <= 1.0e-9 &&
                std::fabs(probability - legacy_probability[i]) <= 1.0e-12 &&
                std::fabs(hindrance - legacy_hindrance[i]) <= 1.0e-9;
          if (all_current)
            all_current =
                std::fabs(energy - current_energy[i]) <= 1.0e-9 &&
                std::fabs(probability - current_probability[i]) <= 1.0e-10 &&
                std::fabs(hindrance - current_hindrance[i]) <= 1.0e-9;
        }
        if (!all_legacy && !all_current)
          throw std::runtime_error(
              "Reviewed U235->Th231 alpha projection saw mixed or unexpected product vector");
      }

      // The legacy U-238 alpha vector and the old RDM daughter-feed law come
      // from the same superseded evaluation. Validate the entire particle
      // vector now, but delay mutation until the independent RDM topology/Q/feed
      // vector has also passed its exact identity gate below. This prevents an
      // output containing only half of the reviewed correction.
      std::vector<xml_node<char>*> u238_alpha_nodes;
      if (u238_alpha_branch) {
        if (!hbr || std::fabs(branch_ratio - 1.0) > 1.0e-12)
          throw std::runtime_error(
              "Reviewed U238->Th234 alpha correction saw unexpected branch ratio");
        constexpr double energies[] = {4038.0, 4151.0, 4198.0};
        constexpr double old_intensities[] = {0.00078, 0.209, 0.79};
        constexpr double corrected_intensities[] = {0.00068, 0.2292, 0.77012};
        constexpr double old_hindrances[] = {40.0, 1.4, 1.0};
        constexpr double corrected_hindrances[] = {45.0, 1.4, 1.0};
        bool all_old = true;
        bool all_corrected = true;
        bool all_old_hindrance = true;
        bool all_corrected_hindrance = true;
        std::size_t index = 0;
        for (xml_node<char>* alpha = tr->first_node("alpha"); alpha;
             alpha = alpha->next_sibling("alpha"), ++index) {
          if (index >= 3)
            throw std::runtime_error(
                "Reviewed U238->Th234 alpha correction saw extra alpha products");
          bool have_energy = false, have_intensity = false, have_hindrance = false;
          const double energy = attr_double(alpha, "energy", have_energy);
          const double intensity = attr_double(alpha, "intensity", have_intensity);
          const double hindrance = attr_double(alpha, "hindrance", have_hindrance);
          if (!have_energy || !have_intensity || !have_hindrance ||
              std::fabs(energy - energies[index]) > 1.0e-9)
            throw std::runtime_error(
                "Reviewed U238->Th234 alpha correction saw unexpected energy vector");
          all_old = all_old &&
              std::fabs(intensity - old_intensities[index]) <= 1.0e-12;
          all_corrected = all_corrected &&
              std::fabs(intensity - corrected_intensities[index]) <= 1.0e-12;
          all_old_hindrance = all_old_hindrance &&
              std::fabs(hindrance - old_hindrances[index]) <= 1.0e-12;
          all_corrected_hindrance = all_corrected_hindrance &&
              std::fabs(hindrance - corrected_hindrances[index]) <= 1.0e-12;
          u238_alpha_nodes.push_back(alpha);
        }
        if (index != 3)
          throw std::runtime_error(
              "Reviewed U238->Th234 alpha correction saw missing alpha products");
        if (!all_old && !all_corrected)
          throw std::runtime_error(
              "Reviewed U238->Th234 alpha correction saw mixed or unexpected intensities");
        if (!all_old_hindrance && !all_corrected_hindrance)
          throw std::runtime_error(
              "Reviewed U238->Th234 alpha correction saw mixed or unexpected hindrances");
        if (!(all_old && all_old_hindrance) &&
            !(all_corrected && all_corrected_hindrance))
          throw std::runtime_error(
              "Reviewed U238->Th234 alpha correction saw incoherent intensity/hindrance state");
      }

      bool transition_th234_normalized = false;
      if (std::string(pa->value()) == "Th234" &&
          std::string(ca->value()) == "Pa234m" && std::string(mode) == "b-" &&
          hbr && std::fabs(branch_ratio - 1.0) <= 1.0e-12) {
        // NNDC ENSDF 234Th beta-minus evaluation, NDS 207 (2026), evaluated
        // through Dec-2023: NR=0.0050(4), normalized from
        // %Ig(63.29+62.86)=3.72(7). The Sandia vector uses the superseded
        // NR=0.0066 normalization, uniformly high by 33/25. Require the complete
        // ordered 13-gamma identity before scaling ONLY gamma intensities.
        // https://www.nndc.bnl.gov/nudat3/getdecaydataset.jsp?dsid=234th+b-+decay&nucleus=234PA
        constexpr double energies[] = {
            62.86, 73.92, 83.3, 87.02, 92.38, 92.80, 103.35,
            103.71, 108.0, 112.8, 132.9, 184.8, 63.29};
        constexpr double old_intensities[] = {
            0.0002112, 0.0, 0.000792, 0.0001914, 0.028116,
            0.027715, 0.00004158, 0.0, 0.0001056, 0.002772,
            0.0, 0.000132, 0.0485};
        constexpr std::size_t n_expected =
            sizeof(energies) / sizeof(energies[0]);
        if (ng != n_expected)
          throw std::runtime_error(
              "Reviewed Th234 gamma normalization saw unexpected gamma count");
        bool all_old = true;
        bool all_corrected = true;
        for (std::size_t i = 0; i < n_expected; ++i) {
          if (std::fabs(gs[i].energy - energies[i]) > 1.0e-9)
            throw std::runtime_error(
                "Reviewed Th234 gamma normalization saw unexpected energy vector");
          const double corrected = old_intensities[i] * 25.0 / 33.0;
          if (old_intensities[i] != 0.0) {
            all_old = all_old &&
                std::fabs(gs[i].intensity - old_intensities[i]) <= 1.0e-12;
            all_corrected = all_corrected &&
                std::fabs(gs[i].intensity - corrected) <= 1.0e-10;
          }
        }
        if (!all_old && !all_corrected)
          throw std::runtime_error(
              "Reviewed Th234 gamma normalization saw mixed or unknown intensity vector");
        if (all_old) {
          for (std::size_t i = 0; i < n_expected; ++i) {
            const double corrected = old_intensities[i] * 25.0 / 33.0;
            set_attr_precise(doc, gs[i].node, "intensity", corrected);
            gs[i].intensity = corrected;
          }
        }
        transition_th234_normalized = true;
        ++n_th234_gamma_normalized;
      }

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
      // The conservative structure-aware matcher is trusted as-is. A Jul 2026
      // verification pass over the common field/calibration nuclides + progeny
      // (checked vs the GEANT4 z{Z}.a{A} level files) confirmed every unmatched
      // above-threshold gamma is a genuine SandiaDecay-vs-GEANT4 data difference
      // (line absent from the G4 scheme / ambiguous), not a matcher
      // error, so the ordinary structure matcher has no per-line overrides.
      // The reviewed primary-evaluation identity exceptions below are separate
      // physical-provenance corrections and may deliberately select an exact
      // edge independently of the general matcher result.
      const std::vector<GammaMatch> mm = match_gammas(scheme, sg);

      std::vector<int> from_lvl(ng, -1), to_lvl(ng, -1);
      std::vector<const G4Transition*> gtr(ng, nullptr);
      std::set<int> populated;
      bool transition_intensity_fixed = false;
      for (std::size_t i = 0; i < ng; ++i) {
        if (mm[i].collision_rejected) ++n_collision_rejected;
        if (mm[i].ambiguous) ++n_ambiguous;
        const VerifiedE0Override* e0_override = hbr
            ? verified_e0_identity(pa->value(), ca->value(), mode,
                                   branch_ratio, gs[i].energy)
            : nullptr;
        int matched_index = mm[i].accepted ? mm[i].g4 : -1;
        if (e0_override) {
          matched_index = -1;
          for (std::size_t candidate = 0;
               candidate < scheme.transitions.size(); ++candidate) {
            if (!matches_verified_e0_edge(*e0_override,
                                          scheme.transitions[candidate]))
              continue;
            if (matched_index >= 0)
              throw std::runtime_error(
                  "Reviewed pure-E0 override matched multiple G4 edges");
            matched_index = static_cast<int>(candidate);
          }
          if (matched_index < 0)
            throw std::runtime_error(
                "Reviewed pure-E0 override did not find its exact G4 edge");
        }
        if (matched_index < 0) continue;
        ++n_accepted;
        if (mm[i].structure_changed) ++n_struct_changed;
        if (mm[i].structure_used) ++n_struct_used;
        gtr[i] = &scheme.transitions[matched_index];
        for (std::size_t previous = 0; previous < i; ++previous)
          if (gtr[previous] == gtr[i])
            throw std::runtime_error(
                "Reviewed pure-E0 override collided with another accepted gamma");
        from_lvl[i] = gtr[i]->from_level;
        to_lvl[i] = gtr[i]->to_level;
        populated.insert(gtr[i]->from_level);
        populated.insert(gtr[i]->to_level);
        if (hbr) {
          const GammaIntensityOverride* fix = gamma_intensity_override(
              pa->value(), ca->value(), mode, branch_ratio, gs[i].energy, *gtr[i]);
          if (fix) {
            constexpr double kIntensityIdentityTolerance = 1.0e-10;
            const bool original = std::fabs(gs[i].intensity - fix->original_intensity) <=
                                  kIntensityIdentityTolerance;
            const bool already_fixed =
                std::fabs(gs[i].intensity - fix->corrected_intensity) <=
                kIntensityIdentityTolerance;
            if (!original && !already_fixed)
              throw std::runtime_error(
                  "Reviewed U235 19.59-keV intensity correction saw unexpected input intensity");
            set_attr_precise(doc, gs[i].node, "intensity", fix->corrected_intensity);
            gs[i].intensity = fix->corrected_intensity;
            transition_intensity_fixed = true;
            ++n_intensity_fixed;
          }
        }
        // Internal-conversion coefficients (raw inputs for the vacancy-level
        // x-ray model). Written only when nonzero, on IC-active gammas.
        // K is written for all Z (K x-rays are >10 keV for every EC/IC-relevant
        // progeny). L is written per-subshell (iccL1/iccL2/iccL3) ONLY for
        // progeny whose L x-rays clear the 10 keV file cut (Z >= kLxrayMinZ);
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
        // Copy GEANT4's multipolarity and mixing ratio as advisory metadata for
        // consumers such as angular-correlation calculations. The code is not
        // primary-evaluation provenance. In particular, PhotonEvaporation code 1
        // occurs on zero- and finite-ICC records that SandiaDecay tabulates as
        // measured gamma lines (Ac-232 373.3 keV and Rb-78 189.54 keV are concrete
        // zero-ICC examples). GEANT4 itself selects photon versus conversion from
        // alpha, not this code. A consumer must therefore not treat mult == 1 by
        // itself as a photon veto or as permission to synthesize an electron.
        // Parsed by SandiaDecay as RadParticle::multipole / mixing_ratio.
        if (gtr[i]->multipole_code != 0) {
          set_attr_int(doc, gs[i].node, "mult", gtr[i]->multipole_code);
          if (gtr[i]->mixing_delta != 0.0)
            set_attr(doc, gs[i].node, "delta", gtr[i]->mixing_delta);
          ++n_mult;
        }
        // Unlike `mult`, this flag is authoritative provenance. It comes only
        // from the reviewed primary-evaluation table above and is attached after
        // an accepted structure match has confirmed the expected level topology.
        if (e0_override) {
          set_attr_int(doc, gs[i].node, "e0Verified", 1);
          ++n_e0_verified;
          if (n16_o16_branch) ++n_n16_e0_verified;
          if (eu152_gd_branch) ++n_eu152_e0_verified;
        }
        // Level-scheme topology (G4 level indices) for the level-path cascade
        // realization: which progeny levels this gamma connects.
        set_attr_int(doc, gs[i].node, "fl", gtr[i]->from_level);
        set_attr_int(doc, gs[i].node, "tl", gtr[i]->to_level);
      }

      if (transition_intensity_fixed) {
        const std::string txt =
            " intensity-fix: U235->Th231 19.59-keV gamma intensity "
            "0.61 -> 0.00583; NNDC ENSDF NDS 185 (2022), evaluated Aug-2022, "
            "%IG=0.583(11), alpha=114.7(19); energy unchanged. ";
        tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                          doc.allocate_string(txt.c_str())));
      }
      if (transition_th234_normalized) {
        const std::string txt =
            " intensity-fix: Th234->Pa234m 13-gamma vector scaled by 25/33; "
            "NNDC ENSDF NDS 207 (2026), evaluated through Dec-2023, "
            "NR 0.0050(4) replaces 0.0066; gamma energies, beta intensities and "
            "coincident conditional intensities unchanged. ";
        tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                          doc.allocate_string(txt.c_str())));
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

      // Documentation comment. Three parts:
      //  (1) SD gammas with no accepted G4 match, bucketed by a per-transition
      //      DETECTION THRESHOLD T = max(1e-5/decay, 0.1% of the strongest gamma):
      //        - I >= T   : a real cascade gap, listed explicitly with its reason;
      //        - 0 < I < T: undetectable, collapsed to a count + energy list;
      //        - I <= 0   : not emitted, just counted.
      //  (2) Curated-override provenance (matches / confirmed no-matches / INVALID).
      //  (3) G4 transitions on populated levels that no SD gamma claims (doc-only:
      //      SD simply lacks the line; adding it is out of scope).
      double I_max = 0.0;
      for (std::size_t i = 0; i < ng; ++i) I_max = std::max(I_max, gs[i].intensity);
      const double T = std::max(1e-5, 1e-3 * I_max);

      std::string doc_sd_above;      // I >= T, real gaps
      std::vector<double> below_T;   // 0 < I < T
      int n_zero = 0;                // I <= 0
      int nd_above = 0;
      for (std::size_t i = 0; i < ng; ++i) {
        if (mm[i].accepted) continue;
        const double I = gs[i].intensity;
        if (I <= 0.0) { ++n_zero; continue; }
        if (I < T) { below_T.push_back(gs[i].energy); continue; }
        if (nd_above < 25) {
          const char* reason;
          if (mm[i].collision_rejected) {
            reason = "edge-collision";
          } else if (std::fabs(gs[i].energy - 510.999) < 1.0) {
            reason = "annih511";
          } else {
            double gap = 0.0;
            const G4Transition* nt =
                scheme.loaded ? nearest_g4(scheme, gs[i].energy, gap) : nullptr;
            reason = (nt && std::fabs(gap) <= 3.0) ? "neartol" : "absent";
          }
          char b[96];
          std::snprintf(b, sizeof(b), " %.2f(I=%.3g,%s)", gs[i].energy, I, reason);
          doc_sd_above += b;
        } else if (nd_above == 25) {
          doc_sd_above += " ...";
        }
        ++nd_above;
      }

      std::string doc_g4;
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

      const bool have_note = !doc_sd_above.empty() || !doc_g4.empty() ||
                             !below_T.empty() || n_zero > 0 || !scheme.loaded;
      if (have_note) {
        std::string txt = " cascade-match:";
        if (!scheme.loaded) txt += " (no G4 level file for progeny)";
        if (!doc_sd_above.empty())
          txt += " SD-unmatched-above-threshold[keV]:" + doc_sd_above + ";";
        if (!below_T.empty()) {
          char h[96];
          std::snprintf(h, sizeof(h),
              " %zu SD gammas below detectable threshold T=%.3g/decay:",
              below_T.size(), T);
          txt += h;
          int shown = 0;
          for (double e : below_T) {
            if (shown++ >= 25) { txt += " ..."; break; }
            char b[24]; std::snprintf(b, sizeof(b), " %.1f", e); txt += b;
          }
          txt += ";";
        }
        if (n_zero > 0) {
          char b[72];
          std::snprintf(b, sizeof(b),
              " %d SD gammas with I=0 (not emitted, suppressed);", n_zero);
          txt += b;
        }
        if (!doc_g4.empty()) txt += " G4-listed-not-in-SD[keV]:" + doc_g4 + ";";
        txt += " ";
        if (txt.find("--") != std::string::npos)
          throw std::runtime_error("Generated cascade-match comment contains forbidden '--'");
        tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                          doc.allocate_string(txt.c_str())));
        ++n_doc;
      }

      // Independent alpha/beta-minus daughter-level population metadata. These
      // RadioactiveDecay records are the categorical feed law; Sandia particle
      // endpoints are not used because evaluated XML commonly splits,
      // duplicates or aggregates them. Delivery is deliberately all-or-nothing.
      const std::string mode_string(mode);
      if (mode_string == "a" || mode_string == "b-") {
        struct PendingLevelFeed {
          int level = -1;
          double intensity_percent = 0.0;
          double q_keV = 0.0;
          std::string forbiddenness;
        };
        std::string feed_omit;
        std::vector<PendingLevelFeed> pending_feeds;
        auto parent_it = sym_to_nuclide.find(pa->value());
        const auto count_it = same_mode_count.find({pa->value(), mode_string});
        if (count_it == same_mode_count.end() || count_it->second != 1) {
          // Terminal-basin partitioning of one P-block among multiple Sandia
          // transitions is not yet proven. Omit rather than assigning the same
          // complete feed law to two children.
          feed_omit = "multiple-same-mode-transitions";
        } else if (parent_it == sym_to_nuclide.end()) {
          feed_omit = "missing-parent-nuclide";
        } else if (!scheme.loaded) {
          feed_omit = "missing-progeny-level-file";
        } else {
          const NuclideAddress& parent_address = parent_it->second;
          const NuclideAddress& child_address = cit->second;
          const bool za_ok = mode_string == "a"
              ? (child_address.Z == parent_address.Z - 2 &&
                 child_address.A == parent_address.A - 4)
              : (child_address.Z == parent_address.Z + 1 &&
                 child_address.A == parent_address.A);
          if (!za_ok) {
            feed_omit = "daughter-ZA-mismatch";
          } else {
            const G4ECData& decay = get_ec(parent_address);
            const bool physical_selection = decay.loaded &&
                decay.selection_method != ParentBlockSelectionMethod::None &&
                decay.selection_method != ParentBlockSelectionMethod::OrdinalFallback;
            if (!physical_selection) {
              feed_omit = "parent-block-not-physically-identified";
            } else {
              const DecayModeLevelFeeds& feeds =
                  mode_string == "a" ? decay.alpha : decay.beta_minus;
              if (feeds.summary_count != 1 ||
                  !std::isfinite(feeds.branch_fraction) ||
                  feeds.branch_fraction <= 0.0) {
                feed_omit = "missing-or-duplicate-mode-summary";
              } else if (!hbr || !std::isfinite(branch_ratio) ||
                         std::fabs(branch_ratio - feeds.branch_fraction) > 5.0e-5) {
                feed_omit = "mode-branch-closure";
              } else if (feeds.records.empty()) {
                feed_omit = "no-level-records";
              } else {
                double record_percent_sum = 0.0;
                int positive_records = 0;
                for (const DaughterLevelFeedRecord& record : feeds.records) {
                  if (!std::isfinite(record.daughter_ex_keV) ||
                      !std::isfinite(record.intensity_percent) ||
                      record.intensity_percent < 0.0) {
                    feed_omit = "invalid-level-record";
                    break;
                  }
                  if (record.intensity_percent == 0.0)
                    continue;
                  if (!std::isfinite(record.q_keV) || record.q_keV <= 0.0) {
                    feed_omit = "invalid-level-endpoint";
                    break;
                  }
                  std::string forbiddenness;
                  if (record.forbiddenness.empty()) {
                    forbiddenness.clear(); // untagged RDM record is Allowed
                  } else if (record.forbiddenness == "firstForbidden") {
                    forbiddenness = "1";
                  } else if (record.forbiddenness == "secondForbidden") {
                    forbiddenness = "2";
                  } else if (record.forbiddenness == "uniqueFirstForbidden") {
                    forbiddenness = "1u";
                  } else if (record.forbiddenness == "uniqueSecondForbidden") {
                    forbiddenness = "2u";
                  } else if (record.forbiddenness == "uniqueThirdForbidden") {
                    forbiddenness = "3u";
                  } else {
                    feed_omit = "unknown-beta-shape-token";
                    break;
                  }
                  ++positive_records;
                  record_percent_sum += record.intensity_percent;
                  const int level = exact_level_index(
                      scheme, record.daughter_ex_keV, record.floating_flag);
                  if (level < 0) {
                    feed_omit = "unmapped-or-ambiguous-level-record";
                    break;
                  }
                  auto same = std::find_if(
                      pending_feeds.begin(), pending_feeds.end(),
                      [&](const PendingLevelFeed& p) { return p.level == level; });
                  if (same == pending_feeds.end()) {
                    pending_feeds.push_back(
                        {level, record.intensity_percent, record.q_keV,
                         forbiddenness});
                  } else {
                    // Preserve one metadata outcome per positive RDM record.
                    // A many-record-to-one-level collapse is not proven to be a
                    // duplicate and would erase endpoint/shape identity.
                    feed_omit = "multiple-records-map-to-level";
                    break;
                  }
                }
                if (feed_omit.empty() && positive_records == 0)
                  feed_omit = "no-positive-level-records";
                if (feed_omit.empty() &&
                    std::fabs(record_percent_sum -
                              100.0 * feeds.branch_fraction) > 0.001)
                  feed_omit = "record-sum-closure";
                if (feed_omit.empty()) {
                  // Header/record closure has already proved no material mass is
                  // missing. Normalize by the record sum solely to absorb the
                  // last few rounded decimal places and serialize an exact
                  // transition-conditional probability law.
                  for (PendingLevelFeed& entry : pending_feeds)
                    entry.intensity_percent /= record_percent_sum;
                }
              }
            }
          }
        }

        if (u235_alpha_branch) {
          if (!feed_omit.empty() || pending_feeds.size() != 20 ||
              u235_alpha_nodes.empty())
            throw std::runtime_error(
                "Reviewed U235->Th231 alpha projection could not establish exact feed law");
          std::sort(pending_feeds.begin(), pending_feeds.end(),
                    [](const PendingLevelFeed& a, const PendingLevelFeed& b) {
                      return a.level < b.level;
                    });
          constexpr int levels[] = {
              0, 1, 2, 3, 4, 5, 6, 8, 10, 13,
              14, 16, 17, 18, 21, 22, 26, 30, 37, 49};
          constexpr double q_keV[] = {
              4678.1, 4636.148, 4581.971, 4515.99, 4492.384,
              4472.79, 4456.703, 4441.201, 4437.221, 4402.675,
              4400.5, 4361.019, 4353.187, 4347.5, 4326.589,
              4300.527, 4290.272, 4225.908, 4147.87, 4044.054};
          constexpr double probability[] = {
              0.047121041, 0.03773603283, 0.012645011,
              0.002420202106, 0.03052502656, 0.5708804967,
              0.003457503008, 0.1871001628, 0.002074501805,
              0.001649701435, 0.03536503077, 0.001086600945,
              0.0003259902836, 0.002272101977, 0.0006915006016,
              0.001215101057, 0.06025905243, 0.002933902552,
              0.0001718901495, 0.00006915006016};
          for (std::size_t i = 0; i < 20; ++i) {
            const PendingLevelFeed& feed = pending_feeds[i];
            if (feed.level != levels[i] ||
                std::fabs(feed.q_keV - q_keV[i]) > 1.0e-9 ||
                std::fabs(feed.intensity_percent - probability[i]) > 1.0e-10 ||
                !feed.forbiddenness.empty())
              throw std::runtime_error(
                  "Reviewed U235->Th231 alpha projection saw unexpected level/Q/feed vector");
          }

          constexpr double alpha_energy[] = {
              3976.0, 4080.6, 4152.5, 4216.0, 4219.0, 4248.0, 4266.2,
              4273.0, 4282.6, 4323.0, 4323.2, 4357.0, 4364.6, 4377.0,
              4395.7, 4415.1, 4438.3, 4502.9, 4556.5, 4597.7};
          constexpr double alpha_probability[] = {
              0.00006915006016, 0.0001718901495, 0.002933902552,
              0.06025905243, 0.001215101057, 0.0006915006016,
              0.002272101977, 0.0003259902836, 0.001086600945,
              0.001649701435, 0.03536503077, 0.002074501805,
              0.1871001628, 0.003457503008, 0.5708804967,
              0.03052502656, 0.002420202106, 0.012645011,
              0.03773603283, 0.047121041};
          constexpr double alpha_hindrance[] = {
              12.0, 41.8, 11.5, 1.94, 117.0, 320.0, 150.0, 1200.0,
              410.0, 580.0, 26.0, 870.0, 10.4, 740.0, 6.02, 160.0,
              3100.0, 1850.0, 1570.0, 2550.0};
          constexpr int alpha_level[] = {
              49, 37, 30, 26, 22, 21, 18, 17, 16, 13,
              14, 10, 8, 6, 5, 4, 3, 2, 1, 0};
          for (std::size_t i = 0; i < 20; ++i) {
            const int* level_it = std::find(levels, levels + 20, alpha_level[i]);
            if (level_it == levels + 20 ||
                std::fabs(alpha_probability[i] -
                          probability[level_it - levels]) > 1.0e-15)
              throw std::runtime_error(
                  "Reviewed U235->Th231 alpha projection has inconsistent product/feed constants");
          }

          // No mutation occurs until the complete public-product and feed-law
          // identities above have both passed. Preserve the particle block's
          // position relative to gamma/x-ray products while replacing it whole.
          for (xml_node<char>* alpha : u235_alpha_nodes)
            tr->remove_node(alpha);
          for (std::size_t i = 0; i < 20; ++i) {
            xml_node<char>* alpha = doc.allocate_node(
                rapidxml::node_element, doc.allocate_string("alpha"));
            set_attr_precise(doc, alpha, "energy", alpha_energy[i]);
            set_attr_precise(doc, alpha, "hindrance", alpha_hindrance[i]);
            set_attr_precise(doc, alpha, "intensity", alpha_probability[i]);
            if (u235_alpha_insert_before)
              tr->insert_node(u235_alpha_insert_before, alpha);
            else
              tr->append_node(alpha);
          }
          for (std::size_t i = 0; i < 20; ++i)
            pending_feeds[i].intensity_percent = probability[i];

          const std::string text =
              " intensity-fix: U235->Th231 alpha products and level feeds are "
              "the G4-topology-backed normalized projection of Aug-2022 ENSDF "
              "(NDS 185, 2022); 20 represented probabilities sum to one; "
              "emitted alpha energies and hindrances updated; omitted "
              "unmapped/tentative ENSDF rows 4042 keV (~0.007%, no G4 level), "
              "4212 keV (0.0131%, no RDM feed), and 4296 keV (0.00915%, no "
              "RDM feed); RDM level/Q topology and gamma products unchanged. ";
          tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                            doc.allocate_string(text.c_str())));
          ++n_u235_alpha_feed_fixed;
        }

        if (u238_alpha_branch) {
          if (!feed_omit.empty() || pending_feeds.size() != 3 ||
              u238_alpha_nodes.size() != 3)
            throw std::runtime_error(
                "Reviewed U238->Th234 alpha correction could not establish exact feed law");
          std::sort(pending_feeds.begin(), pending_feeds.end(),
                    [](const PendingLevelFeed& a, const PendingLevelFeed& b) {
                      return a.level < b.level;
                    });
          constexpr int levels[] = {0, 1, 2};
          constexpr double q_keV[] = {4269.9, 4220.35, 4106.85};
          constexpr double old_raw_percent[] = {79.017, 20.905, 0.078017};
          constexpr double corrected_probability[] = {0.77012, 0.2292, 0.00068};
          constexpr double old_percent_sum = 100.000017;
          bool all_old = true;
          bool all_corrected = true;
          for (std::size_t i = 0; i < 3; ++i) {
            const PendingLevelFeed& feed = pending_feeds[i];
            if (feed.level != levels[i] ||
                std::fabs(feed.q_keV - q_keV[i]) > 1.0e-9 ||
                !feed.forbiddenness.empty())
              throw std::runtime_error(
                  "Reviewed U238->Th234 alpha correction saw unexpected level/Q topology");
            all_old = all_old &&
                std::fabs(feed.intensity_percent -
                          old_raw_percent[i] / old_percent_sum) <= 1.0e-10;
            all_corrected = all_corrected &&
                std::fabs(feed.intensity_percent - corrected_probability[i]) <= 1.0e-10;
          }
          if (!all_old && !all_corrected)
            throw std::runtime_error(
                "Reviewed U238->Th234 alpha correction saw mixed or unexpected feed probabilities");

          constexpr double corrected_alpha_by_xml_order[] =
              {0.00068, 0.2292, 0.77012};
          constexpr double corrected_hindrance_by_xml_order[] =
              {45.0, 1.4, 1.0};
          for (std::size_t i = 0; i < 3; ++i) {
            set_attr_precise(doc, u238_alpha_nodes[i], "intensity",
                             corrected_alpha_by_xml_order[i]);
            set_attr_precise(doc, u238_alpha_nodes[i], "hindrance",
                             corrected_hindrance_by_xml_order[i]);
            pending_feeds[i].intensity_percent = corrected_probability[i];
          }
          const std::string text =
              " intensity-fix: U238->Th234 alpha/feed probabilities "
              "[ground,49.55,163.05] -> [0.77012,0.2292,0.00068]; "
              "4038-keV alpha hindrance 40 -> 45; NNDC ENSDF NDS 207 "
              "(2026), cutoff 1-Dec-2023; alpha energies and RDM level/Q "
              "topology unchanged. ";
          tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                            doc.allocate_string(text.c_str())));
          ++n_u238_alpha_feed_fixed;
        }

        if (feed_omit.empty() && !pending_feeds.empty()) {
          std::sort(pending_feeds.begin(), pending_feeds.end(),
                    [](const PendingLevelFeed& a, const PendingLevelFeed& b) {
                      return a.level < b.level;
                    });
          for (const PendingLevelFeed& entry : pending_feeds) {
            xml_node<char>* feed = doc.allocate_node(
                rapidxml::node_element, doc.allocate_string("levelFeed"));
            set_attr_int(doc, feed, "lvl", entry.level);
            set_attr_precise(doc, feed, "probability", entry.intensity_percent);
            set_attr_precise(doc, feed, "q", entry.q_keV);
            if (!entry.forbiddenness.empty())
              feed->append_attribute(doc.allocate_attribute(
                  doc.allocate_string("forbiddenness"),
                  doc.allocate_string(entry.forbiddenness.c_str())));
            tr->append_node(feed);
            ++n_level_feed_records;
          }
          ++n_level_feed_transitions;
        } else {
          if (feed_omit.empty()) feed_omit = "empty-feed-law";
          const std::string text = " cascade-feed: omitted(" + feed_omit + "). ";
          tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                            doc.allocate_string(text.c_str())));
          ++n_level_feed_omitted;
          ++level_feed_omit_reasons[feed_omit];
        }
      }

      // The legacy Eu-152 EC product vector is exceptional: its EC, positron
      // and Sm K-x-ray intensities are absolute fractions of all parent decays and
      // close to the 0.721 transition branch ratio.  SandiaDecay's public
      // contract, and every consumer, instead treats RadParticle::intensity as
      // conditional on this Transition and multiplies by branchRatio.  Convert
      // this reviewed vector once at the data boundary.  Every node's type,
      // energy and original value is identity-gated; the explicit basis marker
      // makes an enriched rerun idempotent without a magnitude heuristic.
      if (eu152_sm_branch) {
        if (!hbr || std::fabs(branch_ratio - 0.721) > 1.0e-12)
          throw std::runtime_error(
              "Reviewed Eu152->Sm152 product-intensity correction saw unexpected branch ratio");
        struct WeakProductIdentity {
          const char* element;
          double energy_keV;
          double parent_absolute_mass;
        };
        constexpr WeakProductIdentity expected[] = {
            {"electronCapture", 105.2,       0.00066},
            {"electronCapture", 117.268,     0.00047},
            {"electronCapture", 144.06,      0.000423},
            {"electronCapture", 224.411,     0.00923},
            {"electronCapture", 261.52,      0.000211},
            {"electronCapture", 294.868,     0.02085},
            {"electronCapture", 344.506,     0.2492},
            {"electronCapture", 502.556,     0.00859},
            {"electronCapture", 581.543,     0.00658},
            {"electronCapture", 640.445,     0.1741},
            {"electronCapture", 788.417,     0.2173},
            {"electronCapture", 833.186,     0.00076},
            {"electronCapture", 851.334,     0.00237},
            {"electronCapture", 1063.847,    0.013},
            {"electronCapture", 1507.8205,   0.009},
            {"electronCapture", 1752.5183,   0.0074},
            {"positron",        485.822664,  0.000028},
            {"positron",        730.520464,  0.00011},
            {"xray",            39.097,      0.0000536},
            {"xray",            39.522,      0.211},
            {"xray",            40.118,      0.383},
            {"xray",            45.294,      0.0374},
            {"xray",            45.414,      0.0724},
            {"xray",            45.741,      0.00107},
            {"xray",            46.578,      0.0239},
            {"xray",            46.705,      0.0091},
        };
        std::vector<xml_node<char>*> weak_nodes;
        for (xml_node<char>* e = tr->first_node("electronCapture"); e;
             e = e->next_sibling("electronCapture"))
          weak_nodes.push_back(e);
        for (xml_node<char>* p = tr->first_node("positron"); p;
             p = p->next_sibling("positron"))
          weak_nodes.push_back(p);
        for (xml_node<char>* x = tr->first_node("xray"); x;
             x = x->next_sibling("xray"))
          weak_nodes.push_back(x);
        if (weak_nodes.size() != sizeof(expected) / sizeof(expected[0]))
          throw std::runtime_error(
              "Reviewed Eu152->Sm152 product-intensity correction saw unexpected product count");

        xml_attribute<char>* basis = tr->first_attribute("productIntensityBasis");
        xml_attribute<char>* legacy_weak_basis =
            tr->first_attribute("weakIntensityBasis");
        if (basis && std::string(basis->value(), basis->value_size()) !=
                         "transitionConditional")
          throw std::runtime_error(
              "Reviewed Eu152->Sm152 product-intensity correction saw unknown basis");
        bool all_parent_absolute = true;
        bool all_transition_conditional = true;
        bool migrated_weak_only = legacy_weak_basis && !basis &&
            std::string(legacy_weak_basis->value(),
                        legacy_weak_basis->value_size()) ==
                "transitionConditional";
        for (std::size_t i = 0; i < weak_nodes.size(); ++i) {
          const WeakProductIdentity& id = expected[i];
          if (std::string(weak_nodes[i]->name(), weak_nodes[i]->name_size()) !=
              id.element)
            throw std::runtime_error(
                "Reviewed Eu152->Sm152 product-intensity correction saw unexpected product type");
          bool has_energy = false, has_intensity = false;
          const double energy = attr_double(weak_nodes[i], "energy", has_energy);
          const double intensity =
              attr_double(weak_nodes[i], "intensity", has_intensity);
          if (!has_energy || !has_intensity ||
              std::fabs(energy - id.energy_keV) > 1.0e-7)
            throw std::runtime_error(
                "Reviewed Eu152->Sm152 product-intensity correction saw unexpected product identity");
          all_parent_absolute = all_parent_absolute &&
              std::fabs(intensity - id.parent_absolute_mass) <= 1.0e-12;
          all_transition_conditional = all_transition_conditional &&
              std::fabs(intensity - id.parent_absolute_mass / branch_ratio) <=
                  5.0e-10;
          if (migrated_weak_only) {
            const bool weak_product = i < 18;
            const double migrated_expected = weak_product
                ? id.parent_absolute_mass / branch_ratio
                : id.parent_absolute_mass;
            migrated_weak_only =
                std::fabs(intensity - migrated_expected) <= 5.0e-10;
          }
        }
        if ((!basis && !all_parent_absolute) ||
            (basis && !all_transition_conditional)) {
          if (!migrated_weak_only)
            throw std::runtime_error(
                "Reviewed Eu152->Sm152 product-intensity correction saw mixed or unexpected values");
        }
        if (legacy_weak_basis && !migrated_weak_only)
          throw std::runtime_error(
              "Reviewed Eu152->Sm152 product-intensity correction saw unexpected legacy basis");
        for (std::size_t i = 0; i < weak_nodes.size(); ++i)
          set_attr_precise(doc, weak_nodes[i], "intensity",
                           expected[i].parent_absolute_mass / branch_ratio);
        del_attr(tr, "weakIntensityBasis");
        set_attr_string(doc, tr, "productIntensityBasis", "transitionConditional");
        const std::string text =
            " intensity-fix: Eu152->Sm152 EC, beta-plus and Sm K-x-ray products converted "
            "from parent-absolute masses to transition-conditional intensities "
            "by dividing by branchRatio=0.721; represented conditional mass "
            "0.999004; the residual Other outcome is retained; product energies, "
            "x-ray assignments, logFT values and RDM level/shell metadata unchanged. ";
        tr->append_node(doc.allocate_node(rapidxml::node_comment, nullptr,
                                          doc.allocate_string(text.c_str())));
        ++n_eu152_weak_converted;
      }

      // electronCapture shell fractions (raw input for the vacancy-level x-ray
      // model), matched by EC Q-value to the parent's RadioactiveDecay data.
      auto pit = sym_to_nuclide.find(pa->value());
      if (pit != sym_to_nuclide.end()) {
        long transition_ec = 0, transition_positron = 0;
        for (xml_node<char>* e = tr->first_node("electronCapture"); e;
             e = e->next_sibling("electronCapture"))
          ++transition_ec;
        for (xml_node<char>* p = tr->first_node("positron"); p;
             p = p->next_sibling("positron"))
          ++transition_positron;
        n_ec += transition_ec;
        n_positron += transition_positron;

        if (!transition_ec && !transition_positron)
          continue;

        const G4ECData& ec = get_ec(pit->second);
        if (ec.loaded) {
          for (xml_node<char>* e = tr->first_node("electronCapture"); e;
               e = e->next_sibling("electronCapture")) {
            bool he; const double q = attr_double(e, "energy", he);
            // Keep level identity separate from the summary fallback. Summary
            // shell fractions are useful for vacancy yields, but their
            // daughter_ex_keV == 0 is not evidence of ground-state feeding.
            bool q_ambiguous = false;
            const ECShellFractions* exact = he
                ? ec.for_q(q, 5.0, &q_ambiguous) : nullptr;
            if (q_ambiguous) ++n_ec_q_ambiguous;
            const ECShellFractions* f = exact ? exact
                                               : (ec.has_summary ? &ec.summary : nullptr);
            if (!f) continue;
            // K vacancy fraction for all Z (K x-rays > 10 keV everywhere this
            // matters). L only for progeny whose L x-rays clear the 10 keV cut;
            // M/N x-rays are sub-10 keV for every element and are never modelled,
            // so ecM/ecN are dropped entirely.
            set_attr(doc, e, "ecK", f->fK);
            if (daughter_Z >= kLxrayMinZ)
              set_attr(doc, e, "ecL", f->fL);
            // Progeny level fed by this EC branch (index in the progeny scheme),
            // for the level-path cascade realization.
            const int lvl = (exact && scheme.loaded)
                ? level_index_for_energy(scheme, exact->daughter_ex_keV) : -1;
            if (lvl >= 0) {
              set_attr_int(doc, e, "lvl", lvl);
              ++n_ec_level;
            }
            if (!exact) ++n_ec_summary;
            ++n_ec_shell;
          }

          // SandiaDecay positron energy is the kinetic endpoint. GEANT4's
          // per-level BetaPlus Q includes the two electron rest masses, so use
          // endpoint + 2*m_e and only a real 5-token BetaPlus record. A summary
          // BetaPlus probability cannot identify a daughter level.
          constexpr double kTwoElectronMassKeV = 2.0 * 510.99895;
          for (xml_node<char>* p = tr->first_node("positron"); p;
               p = p->next_sibling("positron")) {
            bool he; const double endpoint = attr_double(p, "energy", he);
            bool q_ambiguous = false;
            const BetaPlusRecord* b = he
                ? ec.beta_plus_for_q(endpoint + kTwoElectronMassKeV, 5.0,
                                     &q_ambiguous)
                : nullptr;
            if (q_ambiguous) ++n_positron_q_ambiguous;
            const int lvl = (b && scheme.loaded)
                ? level_index_for_energy(scheme, b->daughter_ex_keV) : -1;
            if (lvl >= 0) {
              set_attr_int(doc, p, "lvl", lvl);
              ++n_positron_level;
            }
          }
        } else if (transition_ec || transition_positron) {
          ++n_parent_block_missing;
          n_ec_parent_missing += transition_ec;
          n_positron_parent_missing += transition_positron;
        }
      }
    }

    // If the canonical parent branch was present, data drift must not silently
    // drop or duplicate reviewed provenance. Validate before creating output.
    if (n_e0_parent_branch &&
        (n_e0_parent_branch != 1 || n_n16_e0_verified != 1))
      throw std::runtime_error(
          "Reviewed N16->O16 E0 override did not match exactly once");
    if (n_eu152_e0_parent_branch &&
        (n_eu152_e0_parent_branch != 1 || n_eu152_e0_verified != 1))
      throw std::runtime_error(
          "Reviewed Eu152->Gd152 E0 override did not match exactly once");
    if (n_eu152_weak_parent_branch &&
        (n_eu152_weak_parent_branch != 1 || n_eu152_weak_converted != 1))
      throw std::runtime_error(
          "Reviewed Eu152->Sm152 product-intensity correction did not match exactly once");
    if (n_u235_parent_branch &&
        (n_u235_parent_branch != 1 || n_intensity_fixed != 2 ||
         n_u235_alpha_feed_fixed != 1))
      throw std::runtime_error(
          "Reviewed U235->Th231 gamma (19.59, 64.35) and alpha/feed corrections "
          "did not each match exactly once");
    if (n_th234_parent_branch &&
        (n_th234_parent_branch != 1 || n_th234_gamma_normalized != 1))
      throw std::runtime_error(
          "Reviewed Th234->Pa234m gamma normalization did not match exactly once");
    if (n_u238_parent_branch &&
        (n_u238_parent_branch != 1 || n_u238_alpha_feed_fixed != 1))
      throw std::runtime_error(
          "Reviewed U238->Th234 alpha/feed correction did not match exactly once");

    std::ofstream out(out_path, std::ios::out | std::ios::binary);
    if (!out) { std::cerr << "Cannot open output.\n"; return 1; }
    rapidxml::print<char>(out, doc, 0);

    const double pct = n_gamma ? 100.0 * n_accepted / n_gamma : 0.0;
    long n_sel_ground = 0, n_sel_mass = 0, n_sel_half = 0, n_sel_drift = 0,
         n_sel_ordinal = 0, n_sel_none = 0, n_sel_nonordinal = 0;
    for (const auto& entry : ec_cache) {
      const G4ECData& ec = entry.second;
      switch (ec.selection_method) {
        case ParentBlockSelectionMethod::GroundExcitation: ++n_sel_ground; break;
        case ParentBlockSelectionMethod::MassExcitation: ++n_sel_mass; break;
        case ParentBlockSelectionMethod::HalfLife: ++n_sel_half; break;
        case ParentBlockSelectionMethod::SingleBlockHalfLifeDrift: ++n_sel_drift; break;
        case ParentBlockSelectionMethod::OrdinalFallback: ++n_sel_ordinal; break;
        case ParentBlockSelectionMethod::None: ++n_sel_none; break;
      }
      if (ec.loaded && ec.parent_block != ec.parent_isomer)
        ++n_sel_nonordinal;
    }
    std::printf("Enriched %ld transitions (%ld documented with miss comments):\n",
                n_trans, n_doc);
    std::printf("  gammas: %ld; accepted G4 matches: %ld (%.1f%%)\n",
                n_gamma, n_accepted, pct);
    std::printf("  ambiguous (2+ within tol): %ld; structure changed the pick: %ld; "
                "accepted via structure only: %ld\n",
                n_ambiguous, n_struct_changed, n_struct_used);
    std::printf("  records rejected by accepted-edge collisions: %ld\n",
                n_collision_rejected);
    std::printf("  coincidences scanned: %ld; got a2/a4: %ld\n", n_coinc, n_coinc_a2a4);
    std::printf("  gammas with ICC written: %ld\n", n_icc);
    std::printf("  electronCapture shell fractions: %ld/%ld; exact fed levels: %ld; "
                "summary-only shell fallbacks: %ld; ambiguous Q matches suppressed: %ld\n",
                n_ec_shell, n_ec, n_ec_level, n_ec_summary, n_ec_q_ambiguous);
    std::printf("  positron exact fed levels: %ld/%ld; ambiguous Q matches suppressed: %ld\n",
                n_positron_level, n_positron, n_positron_q_ambiguous);
    std::printf("  relevant requested parent blocks missing: %ld transitions "
                "(%ld EC, %ld positron records)\n",
                n_parent_block_missing, n_ec_parent_missing, n_positron_parent_missing);
    std::printf("  relevant parent-state selection: ground-excitation=%ld "
                "mass-excitation=%ld half-life=%ld single-block-drift=%ld "
                "ordinal-fallback=%ld no-match=%ld; nonordinal blocks=%ld\n",
                n_sel_ground, n_sel_mass, n_sel_half, n_sel_drift, n_sel_ordinal,
                n_sel_none, n_sel_nonordinal);
    std::printf("  gammas with multipolarity written: %ld; verified E0: %ld\n",
                n_mult, n_e0_verified);
    std::printf("  reviewed gamma-intensity corrections: %ld\n", n_intensity_fixed);
    std::printf("  reviewed Eu152 product-intensity conversions: %ld\n",
                n_eu152_weak_converted);
    std::printf("  reviewed U235 alpha/feed vector corrections: %ld\n",
                n_u235_alpha_feed_fixed);
    std::printf("  reviewed Th234 gamma-vector normalizations: %ld\n",
                n_th234_gamma_normalized);
    std::printf("  reviewed U238 alpha/feed vector corrections: %ld\n",
                n_u238_alpha_feed_fixed);
    std::printf("  exact alpha/beta-minus feed laws: %ld transitions, %ld level records; "
                "omitted: %ld\n",
                n_level_feed_transitions, n_level_feed_records,
                n_level_feed_omitted);
    for (const auto& reason : level_feed_omit_reasons)
      std::printf("    feed omitted %-38s %ld\n", reason.first.c_str(), reason.second);
    std::printf("Wrote %s\n", out_path.c_str());
  } catch (std::exception& e) {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return 0;
}
