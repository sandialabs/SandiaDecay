#include "G4LevelData.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace cascade_enrich {

const G4Level* G4LevelScheme::level(int index) const {
  for (const G4Level& l : levels)
    if (l.index == index)
      return &l;
  return nullptr;
}

bool decode_multipolarity(int code, int& L0, int& Lp) {
  L0 = Lp = 0;
  if (code <= 0)
    return false;
  // Map a single multipole digit d (1..9 => E0,E1,M1,E2,M2,E3,M3,E4,M4) to its
  // angular-momentum order L = d/2 (integer): 2->1, 3->1, 4->2, 5->2, ...
  auto order = [](int d) { return d / 2; };
  if (code < 100) {
    L0 = Lp = order(code);
  } else {
    const int nx = code / 100;   // lower multipole
    const int ny = code % 100;   // higher multipole
    L0 = order(nx);
    Lp = order(ny);
  }
  // Guard against a degenerate decode (e.g. an E0 monopole digit).
  if (L0 <= 0 && Lp <= 0)
    return false;
  if (L0 <= 0) L0 = Lp;
  if (Lp <= 0) Lp = L0;
  return true;
}

namespace {

std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> t;
  std::istringstream is(line);
  std::string tok;
  while (is >> tok)
    t.push_back(tok);
  return t;
}

// A level line's 2nd token is the floating-level flag, a non-numeric string
// ('-', '+X', ...). A transition line's 2nd token is the gamma energy (numeric).
bool is_level_line(const std::vector<std::string>& t) {
  if (t.size() < 6)
    return false;
  const std::string& flag = t[1];
  if (flag == "-")
    return true;
  return !flag.empty() && (flag[0] == '+');
}

int two_j_from_jpi(double jpi, bool& known, int& parity) {
  if (std::abs(jpi - 99.0) < 0.5) {
    known = false;
    parity = 0;
    return -1;
  }
  known = true;
  parity = (jpi < 0.0) ? -1 : +1;
  return static_cast<int>(std::lround(2.0 * std::fabs(jpi)));
}

} // namespace

G4LevelScheme load_photon_evaporation(int Z, int A, const std::string& data_dir) {
  G4LevelScheme scheme;
  scheme.Z = Z;
  scheme.A = A;

  std::string path = data_dir;
  if (!path.empty() && path.back() != '/')
    path += '/';
  path += "z" + std::to_string(Z) + ".a" + std::to_string(A);

  std::ifstream in(path);
  if (!in)
    return scheme;  // loaded == false

  int current_level = -1;
  int remaining_gammas = 0;

  std::string line;
  while (std::getline(in, line)) {
    const std::vector<std::string> t = tokenize(line);
    if (t.empty())
      continue;

    const bool want_level = (remaining_gammas <= 0);
    if (want_level && is_level_line(t)) {
      G4Level lvl;
      lvl.index = std::atoi(t[0].c_str());
      lvl.excitation_keV = std::atof(t[2].c_str());
      lvl.halflife_s = std::atof(t[3].c_str());
      const double jpi = std::atof(t[4].c_str());
      lvl.twoJ = two_j_from_jpi(jpi, lvl.j_known, lvl.parity);
      const int n_gammas = std::atoi(t[5].c_str());
      scheme.levels.push_back(lvl);
      current_level = lvl.index;
      remaining_gammas = n_gammas;
      continue;
    }

    // Otherwise this is a transition line of the current level.
    if (t.size() < 6 || current_level < 0)
      continue;

    G4Transition tr;
    tr.from_level = current_level;
    tr.to_level = std::atoi(t[0].c_str());
    tr.gamma_keV = std::atof(t[1].c_str());
    tr.rel_intensity = std::atof(t[2].c_str());
    tr.multipole_code = std::atoi(t[3].c_str());
    tr.mixing_delta = std::atof(t[4].c_str());
    tr.alpha_total = std::atof(t[5].c_str());
    decode_multipolarity(tr.multipole_code, tr.L0, tr.Lp);

    if (tr.alpha_total > 0.0 && t.size() > 6) {
      tr.has_icc = true;
      const std::size_t navail = t.size() - 6;
      for (std::size_t k = 0; k < 10 && k < navail; ++k)
        tr.icc_shell[k] = std::atof(t[6 + k].c_str());
    }

    scheme.transitions.push_back(tr);
    --remaining_gammas;
  }

  scheme.loaded = true;
  return scheme;
}

} // namespace cascade_enrich
