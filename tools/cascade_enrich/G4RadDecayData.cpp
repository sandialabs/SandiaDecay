#include "G4RadDecayData.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

namespace cascade_enrich {

const ECShellFractions* G4ECData::for_level(double ex_keV, double tol) const {
  const ECShellFractions* best = nullptr;
  double best_d = tol;
  for (const ECShellFractions& f : per_level) {
    const double d = std::fabs(f.daughter_ex_keV - ex_keV);
    if (d <= best_d) {
      best_d = d;
      best = &f;
    }
  }
  if (best)
    return best;
  return has_summary ? &summary : nullptr;
}

const ECShellFractions* G4ECData::for_q(double q_keV, double tol) const {
  const ECShellFractions* best = nullptr;
  double best_d = tol;
  for (const ECShellFractions& f : per_level) {
    const double d = std::fabs(f.q_keV - q_keV);
    if (d <= best_d) {
      best_d = d;
      best = &f;
    }
  }
  if (best)
    return best;
  return has_summary ? &summary : nullptr;
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

// Returns 'K'/'L'/'M'/'N' for an EC mode token like "KshellEC", else 0.
char ec_shell(const std::string& mode) {
  if (mode.size() < 7)
    return 0;
  if (mode.compare(1, 6, "shellE") == 0 && mode.back() == 'C') {
    const char s = mode[0];
    if (s == 'K' || s == 'L' || s == 'M' || s == 'N')
      return s;
  }
  return 0;
}

struct Accum { double K = 0, L = 0, M = 0, N = 0; };

void add(Accum& a, char shell, double v) {
  switch (shell) {
    case 'K': a.K += v; break;
    case 'L': a.L += v; break;
    case 'M': a.M += v; break;
    case 'N': a.N += v; break;
    default: break;
  }
}

ECShellFractions normalize(double ex, double q, const Accum& a) {
  ECShellFractions f;
  f.daughter_ex_keV = ex;
  f.q_keV = q;
  const double tot = a.K + a.L + a.M + a.N;
  f.total_intensity = tot;
  if (tot > 0.0) {
    f.fK = a.K / tot;
    f.fL = a.L / tot;
    f.fM = a.M / tot;
    f.fN = a.N / tot;
  }
  return f;
}

} // namespace

G4ECData load_ec_shell_fractions(int Zparent, int Aparent, const std::string& data_dir) {
  G4ECData data;
  data.Z = Zparent;
  data.A = Aparent;

  std::string path = data_dir;
  if (!path.empty() && path.back() != '/')
    path += '/';
  path += "z" + std::to_string(Zparent) + ".a" + std::to_string(Aparent);

  std::ifstream in(path);
  if (!in)
    return data;  // loaded == false

  // Group per-level (5-token) EC records by rounded daughter excitation, and the
  // 3-token summary records separately. Only the first (ground-state) parent
  // block is used.
  std::map<long, Accum> per_level_acc;       // key = round(daughter_ex*1000)
  std::map<long, double> per_level_ex;       // key -> representative ex
  std::map<long, double> per_level_q;        // key -> representative Q
  Accum summary_acc;
  bool seen_ground_parent = false;
  bool past_ground_parent = false;

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    const std::vector<std::string> t = tokenize(line);
    if (t.empty())
      continue;

    if (t[0] == "P") {
      // Parent-state header: P excitation flag halflife
      const double parent_ex = (t.size() > 1) ? std::atof(t[1].c_str()) : 0.0;
      if (!seen_ground_parent && std::fabs(parent_ex) < 1e-3) {
        seen_ground_parent = true;
      } else if (seen_ground_parent) {
        past_ground_parent = true;  // a subsequent P starts an isomer block
      }
      continue;
    }

    if (!seen_ground_parent || past_ground_parent)
      continue;

    const char shell = ec_shell(t[0]);
    if (!shell)
      continue;

    if (t.size() == 3) {
      // Summary record: Mode daughterEx(=0) intensity(fraction)
      add(summary_acc, shell, std::atof(t[2].c_str()));
    } else if (t.size() >= 5) {
      // Per-level record: Mode daughterEx flag intensity(%) Q
      const double ex = std::atof(t[1].c_str());
      const double inten = std::atof(t[3].c_str());
      const double q = std::atof(t[4].c_str());
      const long key = std::lround(ex * 1000.0);
      add(per_level_acc[key], shell, inten);
      per_level_ex[key] = ex;
      per_level_q[key] = q;
    }
  }

  for (const auto& kv : per_level_acc)
    data.per_level.push_back(
        normalize(per_level_ex[kv.first], per_level_q[kv.first], kv.second));

  if (summary_acc.K + summary_acc.L + summary_acc.M + summary_acc.N > 0.0) {
    data.summary = normalize(0.0, 0.0, summary_acc);
    data.has_summary = true;
  }

  data.loaded = true;
  return data;
}

} // namespace cascade_enrich
