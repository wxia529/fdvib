#include "fdvib.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace fdvib {

std::array<double,3> symmetric_eigenvalues(Vec3 diag, Vec3 off) {
    // Stable analytic eigenvalues of a real symmetric 3x3 matrix.
    const double p1=off[0]*off[0]+off[1]*off[1]+off[2]*off[2];
    if(p1==0){ std::array<double,3>a{diag[0],diag[1],diag[2]}; std::sort(a.begin(),a.end()); return a; }
    const double q=(diag[0]+diag[1]+diag[2])/3;
    const double p2=(diag[0]-q)*(diag[0]-q)+(diag[1]-q)*(diag[1]-q)+(diag[2]-q)*(diag[2]-q)+2*p1;
    const double p=std::sqrt(p2/6);
    const double a00=(diag[0]-q)/p,a11=(diag[1]-q)/p,a22=(diag[2]-q)/p;
    const double a01=off[0]/p,a02=off[1]/p,a12=off[2]/p;
    const double det=a00*a11*a22+2*a01*a02*a12-a00*a12*a12-a11*a02*a02-a22*a01*a01;
    const double r=std::clamp(det/2.0,-1.0,1.0), phi=std::acos(r)/3;
    const double e3=q+2*p*std::cos(phi), e1=q+2*p*std::cos(phi+2*PI/3), e2=3*q-e1-e3;
    std::array<double,3>a{e1,e2,e3}; std::sort(a.begin(),a.end()); return a;
}

struct VibThermo { double zpe{},u{},s{},f{}; int imag{},zero{},floored{},used{}; };

VibThermo vibration_thermo(const std::vector<Mode>& modes,double T,const std::string& model,double floor,double zero_tol) {
    VibThermo v;
    for(const auto&m:modes){ double nu=m.freq;
        if(std::abs(nu)<zero_tol){++v.zero;continue;} if(nu<0){++v.imag;continue;}
        if(model=="frequency_floor"&&nu<floor){nu=floor;++v.floored;} ++v.used;
        const double e=nu*CM_TO_EV,x=e/(KB_EV*T);
        v.zpe+=0.5*e; v.u+=0.5*e+e/std::expm1(x);
        v.s+=KB_EV*(x/std::expm1(x)-std::log1p(-std::exp(-x)));
    }
    v.f=v.u-T*v.s; return v;
}

std::vector<Mode> gas_vibrational_modes(const std::vector<Mode> &modes,
                                         const std::string &rotor_type,
                                         double rigid_tolerance) {
    const int rigid_dof = rotor_type == "atom" ? 3 : (rotor_type == "linear" ? 5 : 6);
    if (static_cast<int>(modes.size()) < rigid_dof)
        throw std::runtime_error("Not enough normal modes for gas RRHO rigid-body projection");

    std::vector<std::size_t> order(modes.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return std::abs(modes[a].freq) < std::abs(modes[b].freq);
    });
    std::vector<bool> rigid(modes.size(), false);
    double largest_rigid = 0.0;
    for (int i = 0; i < rigid_dof; ++i) {
        rigid[order[i]] = true;
        largest_rigid = std::max(largest_rigid, std::abs(modes[order[i]].freq));
    }
    if (largest_rigid > rigid_tolerance) {
        std::ostringstream message;
        message << "Expected " << rigid_dof << " rigid-body modes within "
                << rigid_tolerance << " cm^-1, but the largest is "
                << largest_rigid << " cm^-1";
        throw std::runtime_error(message.str());
    }

    std::vector<Mode> vibrations;
    vibrations.reserve(modes.size() - rigid_dof);
    for (std::size_t i = 0; i < modes.size(); ++i) {
        if (rigid[i]) continue;
        if (modes[i].freq <= 0.0) {
            std::ostringstream message;
            message << "Gas RRHO has an imaginary/non-positive vibrational mode at "
                    << modes[i].freq << " cm^-1; optimize the geometry before thermochemistry";
            throw std::runtime_error(message.str());
        }
        vibrations.push_back(modes[i]);
    }
    return vibrations;
}

void thermo(const fs::path &results) {
    if (!fs::is_directory(results)) throw std::runtime_error("Not a result directory: " + results.string());
    const auto c=Config::load(results/"thermo.in");
    c.require_only({"model", "temperature_k", "pressure_bar", "symmetry_number",
                    "electronic_degeneracy", "rotor_type", "low_frequency_model",
                    "frequency_floor_cm1", "zero_tolerance_cm1"}, "thermo.in");
    const std::string model=lower(c.get("model"));
    if(model!="gas_rrho"&&model!="local_harmonic") throw std::runtime_error("THERMO model must be gas_rrho or local_harmonic");
    if(model=="local_harmonic" &&
       (c.has("pressure_bar") || c.has("symmetry_number") ||
        c.has("electronic_degeneracy") || c.has("rotor_type")))
        throw std::runtime_error("local_harmonic must not contain gas-only thermochemistry parameters");
    if (fs::exists(results/"dynmat.in")) {
        const auto d=Config::load(results/"dynmat.in");
        const auto asr=lower(d.get("asr","no"));
        const auto remove=lower(d.get("remove_interaction_blocks",".false."));
        const bool removes_blocks=(remove==".true."||remove=="true"||remove=="t");
        if(model=="local_harmonic"&&!removes_blocks)
            throw std::runtime_error("local_harmonic requires remove_interaction_blocks=.true. in dynmat.in");
        if(model=="gas_rrho"&&(removes_blocks||asr!="zero-dim"))
            throw std::runtime_error("gas_rrho requires asr='zero-dim' and remove_interaction_blocks=.false. in dynmat.in");
    }
    const double T=c.real("temperature_k",-1); if(T<=0) throw std::runtime_error("temperature_k must be > 0");
    const std::string low=lower(c.get("low_frequency_model","harmonic"));
    if(low!="harmonic"&&low!="frequency_floor") throw std::runtime_error("Bad low_frequency_model");
    if(model=="gas_rrho"&&low!="harmonic")
        throw std::runtime_error("gas_rrho requires low_frequency_model='harmonic'; rigid-body modes are excluded by molecular degrees of freedom");
    const double floor=c.real("frequency_floor_cm1",50),zt=c.real("zero_tolerance_cm1",1);
    if (zt < 0.0) throw std::runtime_error("zero_tolerance_cm1 must be non-negative");
    if (low=="frequency_floor" && floor <= 0.0) throw std::runtime_error("frequency_floor_cm1 must be positive");
    std::vector<fs::path> dyns, freqs;
    for(const auto&e:fs::directory_iterator(results)){const auto n=e.path().filename().string();if(n.size()>=5&&n.substr(n.size()-5)==".dynG")dyns.push_back(e.path());if(n.size()>=9&&n.substr(n.size()-9)==".freq.out")freqs.push_back(e.path());}
    if(dyns.size()!=1||freqs.size()!=1) throw std::runtime_error("Result directory must contain exactly one .dynG and one .freq.out");
    const auto &dyn=dyns.front(), &freq=freqs.front();
    const auto g=read_dyn_geometry(dyn);
    const auto modes=parse_modes(freq,static_cast<int>(g.masses.size()));
    std::string type;
    std::array<double,3> I{};
    double mtot=0.0;
    if(model=="gas_rrho") {
        mtot=std::accumulate(g.masses.begin(),g.masses.end(),0.0);
        Vec3 com{0,0,0}; for(size_t i=0;i<g.masses.size();++i)for(int k=0;k<3;++k)com[k]+=g.masses[i]*g.r_bohr[i][k]/mtot;
        Vec3 d{0,0,0},off{0,0,0};
        for(size_t i=0;i<g.masses.size();++i){const double m=g.masses[i]*AMU_KG;double x[3];for(int k=0;k<3;++k)x[k]=(g.r_bohr[i][k]-com[k])*BOHR_TO_ANG*1e-10;
            d[0]+=m*(x[1]*x[1]+x[2]*x[2]);d[1]+=m*(x[0]*x[0]+x[2]*x[2]);d[2]+=m*(x[0]*x[0]+x[1]*x[1]);off[0]-=m*x[0]*x[1];off[1]-=m*x[0]*x[2];off[2]-=m*x[1]*x[2];}
        I=symmetric_eigenvalues(d,off); type=lower(c.get("rotor_type","auto"));
        if(type=="auto") type=(g.masses.size()==1?"atom":(I[0]/std::max(I[2],1e-300)<1e-6?"linear":"nonlinear"));
        if(type!="atom"&&type!="linear"&&type!="nonlinear") throw std::runtime_error("rotor_type must be auto, atom, linear, nonlinear");
        if((type=="atom")!=(g.masses.size()==1)) throw std::runtime_error("rotor_type='atom' is valid only for a monatomic species");
        if(type=="nonlinear"&&g.masses.size()<3) throw std::runtime_error("A molecule with fewer than three atoms cannot be a nonlinear rotor");
    }
    const auto thermo_modes=model=="gas_rrho"?gas_vibrational_modes(modes,type,zt):modes;
    const auto vib=vibration_thermo(thermo_modes,T,low,floor,model=="gas_rrho"?0.0:zt);
    double hcorr=vib.u,stotal=vib.s,gcorr=vib.f,strans=0,srot=0,selec=0,htrans=0,urot=0;
    if(model=="gas_rrho"){
        if(!c.has("pressure_bar")||!c.has("symmetry_number")) throw std::runtime_error("gas_rrho requires explicit pressure_bar and symmetry_number");
        const auto dataset=Config::load(results/"fdvib.in.reference");
        if(lower(dataset.get("system_type"))!="gas") throw std::runtime_error("gas_rrho requires a gas fdvib.in.reference");
        const double pbar=c.real("pressure_bar",-1); const int sigma=c.integer("symmetry_number",0),mult=dataset.integer("multiplicity",0);
        if(pbar<=0||sigma<=0||mult<=0) throw std::runtime_error("Gas pressure, symmetry, and multiplicity must be positive");
        const double qtrans=std::pow(2*PI*(mtot*AMU_KG)*KB_SI*T/(H_SI*H_SI),1.5)*(KB_SI*T/(pbar*BAR_PA));
        strans=KB_EV*(std::log(qtrans)+2.5);htrans=2.5*KB_EV*T;
        if(type=="linear") {
            const double moment=std::max(I[1],I[2]); if(!(moment>0.0)) throw std::runtime_error("Invalid linear-molecule moment of inertia");
            const double qr=8*PI*PI*moment*KB_SI*T/(sigma*H_SI*H_SI);srot=KB_EV*(std::log(qr)+1);urot=KB_EV*T;
        }
        else if(type=="nonlinear") {
            if(!(I[0]>0.0&&I[1]>0.0&&I[2]>0.0)) throw std::runtime_error("Invalid nonlinear-molecule moments of inertia");
            const double qr=std::sqrt(PI)/sigma*std::pow(8*PI*PI*KB_SI*T/(H_SI*H_SI),1.5)*std::sqrt(I[0]*I[1]*I[2]);srot=KB_EV*(std::log(qr)+1.5);urot=1.5*KB_EV*T;
        }
        int degeneracy=mult; const auto ed=lower(c.get("electronic_degeneracy","auto"));
        if(ed!="auto") {
            const double value=number(ed), rounded=std::round(value);
            if(std::abs(value-rounded)>1.0e-10) throw std::runtime_error("electronic_degeneracy must be an integer or auto");
            degeneracy=static_cast<int>(rounded);
        }
        if(degeneracy<=0)throw std::runtime_error("electronic_degeneracy must be positive or auto");
        selec=KB_EV*std::log(static_cast<double>(degeneracy)); hcorr=vib.u+htrans+urot; stotal=vib.s+strans+srot+selec; gcorr=hcorr-T*stotal;
    }
    std::ostringstream o;o<<"# FDVIB thermochemistry\n# model: "<<model<<"\n# low_frequency_model: "<<low<<"\n";
    if(low=="frequency_floor")o<<"# frequency_floor_cm1: "<<floor<<"\n";
    if(model=="gas_rrho")o<<"# rotor_type: "<<type<<"\n# rigid_body_modes_excluded: "<<(type=="atom"?3:(type=="linear"?5:6))<<"\n# expected_vibrational_modes: "<<thermo_modes.size()<<"\n";
    o<<"# imaginary_modes_excluded: "<<vib.imag<<"\n# zero_modes_excluded: "<<vib.zero<<"\n# positive_modes_used: "<<vib.used<<"\n# modes_floored: "<<vib.floored<<"\n";
    o<<"# units: T=K energies=eV entropy=eV/K\n"<<std::fixed;
    if(model=="local_harmonic") {
        o<<"# "<<std::setw(11)<<"T/K"<<std::setw(16)<<"ZPE/eV"<<std::setw(16)<<"U_vib/eV"
         <<std::setw(19)<<"S_vib/eV_K"<<std::setw(16)<<"TS_vib/eV"<<std::setw(16)<<"F_vib/eV"<<'\n';
        o<<"  "<<std::setw(11)<<std::setprecision(3)<<T
         <<std::setw(16)<<std::setprecision(10)<<vib.zpe<<std::setw(16)<<vib.u
         <<std::setw(19)<<std::setprecision(12)<<vib.s
         <<std::setw(16)<<std::setprecision(10)<<T*vib.s<<std::setw(16)<<vib.f<<'\n';
    } else {
        o<<"# "<<std::setw(11)<<"T/K"<<std::setw(15)<<"ZPE/eV"<<std::setw(15)<<"U_vib/eV"
         <<std::setw(15)<<"H_trans/eV"<<std::setw(15)<<"U_rot/eV"
         <<std::setw(18)<<"S_trans/eV_K"<<std::setw(18)<<"S_rot/eV_K"
         <<std::setw(18)<<"S_vib/eV_K"<<std::setw(18)<<"S_elec/eV_K"
         <<std::setw(15)<<"H_corr/eV"<<std::setw(15)<<"G_corr/eV"<<'\n';
        o<<"  "<<std::setw(11)<<std::setprecision(3)<<T
         <<std::setw(15)<<std::setprecision(10)<<vib.zpe<<std::setw(15)<<vib.u
         <<std::setw(15)<<htrans<<std::setw(15)<<urot
         <<std::setw(18)<<std::setprecision(12)<<strans<<std::setw(18)<<srot
         <<std::setw(18)<<vib.s<<std::setw(18)<<selec
         <<std::setw(15)<<std::setprecision(10)<<hcorr<<std::setw(15)<<gcorr<<'\n';
    }
    write_text(results/"thermo.dat",o.str());
    std::cout<<o.str()<<"Written "<<(results/"thermo.dat")<<'\n';
}

} // namespace fdvib
