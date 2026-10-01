#ifndef readGAM_h
#define readGAM_h
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Euka.h"
#include "soibean.h"
#include "assembly.h"
#include "HaploCart.h"
#include "TrailMix.h"
#include "Dup_Remover.h"
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "gam2prof.h" // Mikkel code
#include "version.h"
#include <sys/wait.h>
#include <unistd.h>
#include <cstring>
#include <cerrno>

//#include "gam2prof.cpp" //Mikkel code

#define VERBOSE
#define ADDDATA
//#define DEBUGREADGRAPH
using namespace std;
using namespace vg;

// "tempeh" is a thin subprocess wrapper around the `tempeh` executable built
// from dep/cdx (github.com/JolanBoucher/cdx, built via its own CMake system -
// a separate C++17/CMake project, not part of vgan's Makefile build; "CDX"
// remains the name of the on-disk index format it builds/reads, unrelated to
// the executable's own name). It already fully handles its own argument
// parsing, mode dispatch (build/inspect/coverage, picked from the input
// file's binary signature), and --help text, so this wrapper just locates
// the built binary and execv()s it with the tempeh args forwarded unmodified
// - no FIFO/pipe plumbing needed, since it reads/writes its inputs/outputs
// directly as regular files (unlike, say, the safari subprocess wrapper for
// trailmix, which streams a GAM through a pipe).
static int run_tempeh_subprocess(int argc, char** argv, const string& cwdProg) {
    const string tempeh_bin = getFullPath(cwdProg + "../dep/cdx/build/tempeh");

    vector<char*> child_argv;
    child_argv.reserve(argc + 1);
    child_argv.emplace_back(const_cast<char*>(tempeh_bin.c_str()));
    for (int i = 2; i < argc; ++i) { child_argv.emplace_back(argv[i]); }
    child_argv.emplace_back(nullptr);

    pid_t pid = fork();
    if (pid == -1) {
        cerr << "[tempeh] fork() failed: " << strerror(errno) << endl;
        return 1;
    }
    if (pid == 0) {
        execv(tempeh_bin.c_str(), child_argv.data());
        cerr << "[tempeh] execv() failed on " << tempeh_bin << ": " << strerror(errno) << endl;
        _exit(127);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) { return WEXITSTATUS(status); }
    if (WIFSIGNALED(status)) {
        cerr << "[tempeh] terminated by signal " << WTERMSIG(status) << endl;
        return 128 + WTERMSIG(status);
    }
    return 1;
}

int main(int argc, char *argv[]) {

    /// JUST FOR TESTING WILL BE PLACED ELSEWHERE //////
    /////// Testing MCMC ////////

    // const double alpha = 0.1;
    // vector<long double> current_vec = {0.05, 0.6, 0.1, 0.25};

    // const vector <long double> test = MCMC().generate_proposal(current_vec, alpha);

    // for (auto elemet:test){
    //     cout << elemet << endl;
    // }


    string_view usage=string("\n")+
        "   vgan is a suite of tools for mitochondrial pangenomics. \n"+
        "   We currently support four subcommands: euka (for the classification of eukaryotic taxa) \n"+
	"   and HaploCart (for modern human mtDNA haplogroup classification), keelime (for the \n"+
        "   hybrid assembly for identified species), and soibean (for the classification \n" +
        "   of eukaryotic species). \n" +
        "   The underlying data structure is the VG graph. \n\n" +
        string(argv[0]) +" <command> options\n"+
        "\n"+
        "   Commands:\n"+
        "      euka         Identify eukaryotic taxa "+"\n"+
        "      duprm        Remove PCR duplicates from a GAM file "+"\n"+
        "      gam2prof     Reads a GAM file and produces a deamination profile for the\n"+
        "                   5' and 3' ends (only for euka_db) " +"\n" +
        "      haplocart    Predict human mitochondrial haplogroup  "+"\n"+
        "      soibean      Identify eukaryotic species " +"\n"+
        "      tempeh       Build/inspect a CDX pangenome coordinate index, or compute\n"+
        "                   GAM coverage against one (see 'tempeh --help') "+"\n"+
        "      trailmix     Inference on ancient human mtDNA mixture "+"\n"+
        "      version      Print version                         " +
	"";

    if( argc==1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage "<<usage<<"\n";
        return 1;
    }


    if(string(argv[1]) == "version"){cerr << "vgan "<<VERSION << endl; return 0;}

    else if(string(argv[1]) == "euka"){
        Euka  euka_;

        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<euka_.usage()<<"\n";
            return 1;
        }

	const string cwdProg=getCWD(argv[0]);
        argv++;
        argc--;
        return euka_.run(argc, argv, cwdProg);

    }
    else if(string(argv[1]) == "keelime"){
        assembly  assembly_;

        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<assembly_.usage()<<"\n";
            return 1;
        }

    const string cwdProg=getCWD(argv[0]);
        argv++;
        argc--;
        return assembly_.run(argc, argv, cwdProg);

    }


    // // Mikkel code begins // //
    else if(string(argv[1]) == "gam2prof"){
        

        Gam2prof  gam2prof_;
        
        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<gam2prof_.usage()<<"\n";
            return 1;
        }
        
    const string cwdProg=getCWD(argv[0]);
        argv++;
        argc--;
        return gam2prof_.run(argc, argv, cwdProg);
    
    }
    // // Mikkel code ends // //

    else if(string(argv[1]) == "duprm"){
        Dup_Remover dup_remover;
        if( argc==2 || argc > 3 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<dup_remover.usage()<<"\n";
            return 1;
        }
        const string cwdProg=getCWD(argv[0]);
        const char *gamfile = argv[2];
        dup_remover.remove_duplicates(gamfile);
        return 1;
                                        }

    else if(string(argv[1]) == "soibean"){
        

        soibean soibean_;
        
        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            const string cwdProg=getCWD(argv[0]);
            cerr<<soibean_.usage(cwdProg)<<"\n";
            return 1;
        }
        
        const string cwdProg=getCWD(argv[0]);
            argv++;
            argc--;
            return soibean_.run(argc, argv, cwdProg);
        
        }


else if(string(argv[1]) == "haplocart"){

        Haplocart  haplocart_;

        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<haplocart_.usage()<<"\n";
            return 1;
        }

        const string cwdProg=getCWD(argv[0]);
        argv++;
        argc--;
        return haplocart_.run(argc, argv, cwdProg);

    }

    else if(string(argv[1]) == "trailmix"){

        Trailmix  trailmix_;
        const string cwdProg = getCWD(argv[0]);

        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<trailmix_.usage()<<"\n";
            return 1;
        }

        argv++;
        argc--;
        return trailmix_.run(argc, argv, cwdProg);

    }

    else if(string(argv[1]) == "tempeh"){

        const string cwdProg = getCWD(argv[0]);
        return run_tempeh_subprocess(argc, argv, cwdProg);

    }else{
        cerr<<"invalid command "<<string(argv[1])<<"\n";
        return 1;
	}





   return 0;
}

#endif
