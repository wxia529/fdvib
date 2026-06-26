#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace fdvib {

double norm(const Vec3 &v) { return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

void write_dynG(const fs::path &p, const QEInput &q, const std::vector<double> &h) {
    const double alat_ang = norm(q.cell[0]);
    const double alat_bohr = alat_ang / BOHR_TO_ANG;
    std::ostringstream o;
    o << "Dynamical matrix file\nfdvib finite-difference local/Gamma matrix\n";
    o << std::setw(3) << q.ntyp << std::setw(5) << q.nat << std::setw(4) << 0;
    o << std::fixed << std::setprecision(7) << std::setw(12) << alat_bohr;
    for (int i = 0; i < 5; ++i) o << std::setw(12) << 0.0;
    o << "\nBasis vectors\n" << std::setprecision(9);
    for (const auto &v : q.cell) o << "  " << std::setw(15) << v[0]/alat_ang << std::setw(15) << v[1]/alat_ang << std::setw(15) << v[2]/alat_ang << '\n';
    o << std::setprecision(12);
    for (int i = 0; i < q.ntyp; ++i)
        o << std::setw(12) << i+1 << "  '" << std::left << std::setw(7) << q.species[i].symbol << std::right << "' " << std::setw(22) << AMU_RY*q.species[i].mass << '\n';
    o << std::setprecision(10);
    for (int i = 0; i < q.nat; ++i)
        o << std::setw(5) << i+1 << std::setw(5) << q.atoms[i].type << std::setw(18) << q.atoms[i].r[0]/alat_ang
          << std::setw(18) << q.atoms[i].r[1]/alat_ang << std::setw(18) << q.atoms[i].r[2]/alat_ang << '\n';
    o << "\n     Dynamical  Matrix in cartesian axes\n\n     q = (    0.000000000   0.000000000   0.000000000 ) \n\n";
    o << std::fixed << std::setprecision(10);
    const int n3 = 3*q.nat;
    for (int a = 0; a < q.nat; ++a) for (int b = 0; b < q.nat; ++b) {
        o << std::setw(5) << a+1 << std::setw(5) << b+1 << '\n';
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) o << std::setw(14) << h[(3*a+i)*n3 + 3*b+j] << "   " << std::setw(12) << 0.0 << "  ";
            o << '\n';
        }
    }
    write_text(p, o.str());
}

std::vector<Mode> parse_modes(const fs::path &p, int nat) {
    if (nat <= 0) throw std::runtime_error("Invalid atom count while parsing modes: " + p.string());
    std::ifstream in(p);
    if (!in) throw std::runtime_error("Cannot read " + p.string());
    std::regex fr(R"(freq\s*\(\s*\d+\s*\)\s*=\s*[-+0-9.EeDd]+\s*\[THz\]\s*=\s*([-+0-9.EeDd]+)\s*\[cm-1\])", std::regex::icase);
    std::regex pair(R"(\(\s*([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*\))");
    std::vector<Mode> out; std::string line; std::smatch m;
    while (std::getline(in, line)) if (std::regex_search(line, m, fr)) {
        Mode mode; mode.freq = number(m[1]); mode.displacement.reserve(nat);
        for (int i = 0; i < nat; ++i) {
            if (!std::getline(in, line) || !std::regex_search(line, m, pair)) throw std::runtime_error("Bad eigenvector block in " + p.string());
            mode.displacement.push_back({number(m[1]), number(m[3]), number(m[5])});
        }
        out.push_back(std::move(mode));
    }
    if (static_cast<int>(out.size()) != 3 * nat)
        throw std::runtime_error("Expected " + std::to_string(3 * nat) + " modes, found " +
                                 std::to_string(out.size()) + " in " + p.string());
    return out;
}

void write_compact_molden(const fs::path &p, const DynGeometry &g,
                          const std::vector<Mode> &modes, double zero_tol);

fs::path completed_forces(const Settings &s, const std::string &id) {
    const auto marker = s.workdir / "state" / (id + ".complete");
    std::istringstream in(read_text(marker));
    std::string calculation, digest, extra;
    if (!(in >> calculation >> digest) || in >> extra ||
        calculation.rfind(id + '_', 0) != 0 || fs::path(calculation).filename() != calculation)
        throw std::runtime_error("Invalid displacement completion snapshot: " + marker.string());
    return s.workdir / "calculations" / calculation / "forces.dat";
}

void analyze(const Settings &s) {
    const auto q = parse_qe_input(s.workdir / "scf.in.reference");
    auto selected = s.selected;
    if (s.system_type == "gas") { selected.resize(q.nat); std::iota(selected.begin(), selected.end(), 1); }
    const int n3 = 3*q.nat;
    std::vector<double> h(n3*n3, 0.0);
    const double delta_bohr = s.displacement / BOHR_TO_ANG;
    for (int atom1 : selected) for (int axis = 0; axis < 3; ++axis) {
        const auto p = read_forces(completed_forces(s, job_name(atom1,axis,1)), q.nat);
        const auto m = read_forces(completed_forces(s, job_name(atom1,axis,-1)), q.nat);
        const int col = 3*(atom1-1)+axis;
        for (int atomj1 : selected) for (int beta=0; beta<3; ++beta) {
            const int row=3*(atomj1-1)+beta;
            h[row*n3+col]=-(p[atomj1-1][beta]-m[atomj1-1][beta])/(2*delta_bohr);
        }
    }
    double asym=0;
    for (int i=0;i<n3;++i) for(int j=0;j<i;++j) {
        asym=std::max(asym,std::abs(h[i*n3+j]-h[j*n3+i]));
        const double x=0.5*(h[i*n3+j]+h[j*n3+i]); h[i*n3+j]=h[j*n3+i]=x;
    }
    const auto results=s.workdir/"results"; fs::create_directories(results);
    const auto dyn=results/(s.output_prefix+".dynG");
    for (const auto &p : {dyn, results/"dynmat.in", results/"metadata.dat"})
        if (fs::exists(p)) throw std::runtime_error("Refuse to overwrite " + p.string());
    fs::copy_file(s.workdir / "metadata.dat", results / "metadata.dat");
    write_dynG(dyn,q,h);
    std::ostringstream di;
    di << "&INPUT\n  fildyn='" << dyn.filename().string() << "',\n  filout='" << s.output_prefix << ".freq.out',\n"
       << "  asr='no',\n"
       << "  remove_interaction_blocks=" << (s.system_type=="gas"?".false.":".true.") << ",\n/\n";
    write_text(results/"dynmat.in",di.str());
    std::cout << "Hessian max asymmetry: " << std::scientific << asym << " Ry/Bohr^2\n"
              << "Wrote " << dyn << " and " << (results/"dynmat.in") << "\n";
}

DynGeometry read_dyn_geometry(const fs::path &p) {
    std::ifstream in(p); if(!in) throw std::runtime_error("Cannot read "+p.string());
    std::string line;
    for (int i = 0; i < 3; ++i)
        if (!std::getline(in, line)) throw std::runtime_error("Incomplete dynG header: " + p.string());
    std::istringstream hs(line); int ntyp,nat,ibrav; double alat; hs>>ntyp>>nat>>ibrav>>alat;
    if(!hs||ibrav!=0||ntyp<=0||nat<=0||!(alat>0.0)) throw std::runtime_error("Unsupported dynG header");
    if (!std::getline(in,line)) throw std::runtime_error("Incomplete dynG basis section: " + p.string());
    for(int i=0;i<3;++i)
        if (!std::getline(in,line)) throw std::runtime_error("Incomplete dynG basis vectors: " + p.string());
    std::vector<double> tmass(ntyp); std::vector<std::string> tsymbol(ntyp);
    const std::regex quoted(R"('([^']*)')");
    for(int i=0;i<ntyp;++i){
        if (!std::getline(in,line)) throw std::runtime_error("Incomplete dynG species block: " + p.string());
        std::smatch match;
        if (!std::regex_search(line, match, quoted)) throw std::runtime_error("Bad species symbol in dynG");
        tsymbol[i]=trim(match[1]);
        std::istringstream ls(line); std::string token, last;
        while(ls>>token) last=token;
        if(last.empty()) throw std::runtime_error("Bad species mass in dynG");
        tmass[i]=number(last)/AMU_RY;
    }
    DynGeometry g; g.masses.resize(nat); g.r_bohr.resize(nat); g.symbols.resize(nat);
    for(int i=0;i<nat;++i){
        if (!std::getline(in,line)) throw std::runtime_error("Incomplete dynG atom block: " + p.string());
        std::istringstream ls(line); int n,t; Vec3 x; ls>>n>>t>>x[0]>>x[1]>>x[2];
        if(!ls||n!=i+1||t<1||t>ntyp||!std::isfinite(x[0])||!std::isfinite(x[1])||!std::isfinite(x[2]))
            throw std::runtime_error("Bad atom row in dynG: " + p.string());
        g.masses[i]=tmass.at(t-1);g.symbols[i]=tsymbol.at(t-1);for(int k=0;k<3;++k)g.r_bohr[i][k]=x[k]*alat;
    }
    return g;
}

void write_compact_molden(const fs::path &p, const DynGeometry &g,
                          const std::vector<Mode> &modes, double zero_tol) {
    std::vector<const Mode*> keep;
    for (const auto &m : modes) {
        double n = 0;
        for (const auto &d : m.displacement) for (double x : d) n += x*x;
        if (std::abs(m.freq) >= zero_tol && n > 1e-16) keep.push_back(&m);
    }
    std::ostringstream o;
    o << "[Molden Format]\n[FREQ]\n" << std::fixed << std::setprecision(8);
    for (auto m : keep) o << m->freq << '\n';
    o << "[FR-COORD]\n";
    for (std::size_t i=0;i<g.symbols.size();++i) o << std::setw(6) << g.symbols[i] << ' ' << std::setw(15) << g.r_bohr[i][0] << ' '
                                    << std::setw(15) << g.r_bohr[i][1] << ' ' << std::setw(15) << g.r_bohr[i][2] << '\n';
    o << "[FR-NORM-COORD]\n" << std::setprecision(10);
    for (std::size_t i = 0; i < keep.size(); ++i) {
        o << " vibration " << i+1 << '\n';
        for (const auto &d : keep[i]->displacement) o << std::setw(15) << d[0] << ' ' << std::setw(15) << d[1] << ' ' << std::setw(15) << d[2] << '\n';
    }
    write_text(p, o.str());
    std::cout << "Wrote " << keep.size() << " compact Molden modes\n";
}

void modes(const fs::path &results) {
    const auto files = result_files(results, "Normal-mode analysis");
    const auto g=read_dyn_geometry(files.dyn);
    const auto parsed=parse_modes(files.freq,static_cast<int>(g.masses.size()));
    const auto mold=results/(files.dyn.stem().string()+".mold");
    write_compact_molden(mold,g,parsed,1.0e-6);
}

} // namespace fdvib
