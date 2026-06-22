#include "fdvib.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

#ifndef FDVIB_VERSION
#define FDVIB_VERSION "unknown"
#endif

namespace {

void usage() {
    std::cerr << "Usage:\n  fdvib prepare fdvib.in [--force]\n  fdvib run fdvib.in\n  fdvib analyze fdvib.in\n  fdvib modes RESULTS_DIR\n  fdvib thermo RESULTS_DIR\n";
}

} // namespace

int main(int argc,char**argv){
    try{
        if(argc==2&&std::string(argv[1])=="--help"){usage();return 0;}
        if(argc==2&&std::string(argv[1])=="--version"){std::cout<<"fdvib "<<FDVIB_VERSION<<'\n';return 0;}
        if(argc<3){usage();return 2;} const std::string cmd=argv[1];
        if(cmd=="prepare") fdvib::prepare(fdvib::settings(argv[2]),argc==4&&std::string(argv[3])=="--force");
        else if(cmd=="run") fdvib::run_jobs(fdvib::settings(argv[2]));
        else if(cmd=="analyze") {
            const auto config_path = fdvib::fs::absolute(argv[2]);
            const auto root = config_path.parent_path();
            const auto current = fdvib::Config::load(config_path);
            const auto workdir = root / current.get("workdir", "fdvib");
            fdvib::analyze(fdvib::settings(workdir / "fdvib.in.reference", root));
        }
        else if(cmd=="modes") fdvib::modes(fdvib::fs::absolute(argv[2]));
        else if(cmd=="thermo") fdvib::thermo(fdvib::fs::absolute(argv[2]));
        else {usage();return 2;}
        return 0;
    }catch(const std::exception&e){std::cerr<<"fdvib: error: "<<e.what()<<'\n';return 1;}
}
