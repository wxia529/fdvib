#include "fdvib.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

#ifndef FDVIB_VERSION
#define FDVIB_VERSION "unknown"
#endif

namespace {

void usage() {
    std::cerr << "Usage:\n  fdvib -inp fdvib.in\n  fdvib modes RESULTS_DIR\n  fdvib thermo RESULTS_DIR -inp thermo.in\n  fdvib shm RESULTS_DIR\n";
}

} // namespace

int main(int argc,char**argv){
    try{
        if(argc==2&&std::string(argv[1])=="--help"){usage();return 0;}
        if(argc==2&&std::string(argv[1])=="--version"){std::cout<<"fdvib "<<FDVIB_VERSION<<'\n';return 0;}
        if(argc<3){usage();return 2;} const std::string cmd=argv[1];
        if(cmd=="-inp" && argc==3) fdvib::calculate(fdvib::settings(argv[2]));
        else if(cmd=="modes" && argc==3) fdvib::modes(fdvib::fs::absolute(argv[2]));
        else if(cmd=="thermo" && argc==5 && std::string(argv[3])=="-inp")
            fdvib::thermo(fdvib::fs::absolute(argv[2]), fdvib::fs::absolute(argv[4]));
        else if(cmd=="shm" && argc==3) fdvib::shm(fdvib::fs::absolute(argv[2]));
        else {usage();return 2;}
        return 0;
    }catch(const std::exception&e){std::cerr<<"fdvib: error: "<<e.what()<<'\n';return 1;}
}
