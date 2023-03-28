#ifndef readGAM_h
#define readGAM_h
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Euka.h"
#include "HaploCart.h"
#include "Dup_Remover.h"
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "gam2prof.h" // Mikkel code
#include "version.h"

//#include "gam2prof.cpp" //Mikkel code

#define VERBOSE
#define ADDDATA
//#define DEBUGREADGRAPH
using namespace std;
using namespace vg;



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
        "   We currently support two subcommands: euka (for the classification of eukaryotic taxa) \n"+
	"   and HaploCart (for modern human mtDNA haplogroup classification). \n"+
        "   The underlying data structure is the VG graph. \n\n" +
        string(argv[0]) +" <command> options\n"+
        "\n"+
        "   Commands:\n"+
        "      euka         Identify eukaryotic species "+"\n"+
        "      duprm         Remove PCR duplicates from a GAM file "+"\n"+
        "      haplocart    Predict human mitochondrial haplogroup  "+"\n"+
        //"      gam2prof     Reads a GAM file and produces a deamination profile for the\n"+
        //"                   5' and 3' ends " +"\n" +
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

else{      if(string(argv[1]) == "haplocart"){

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

    }else{
        cerr<<"invalid command "<<string(argv[1])<<"\n";
        return 1;
	}}





   return 0;
}

#endif
