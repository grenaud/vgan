#ifndef readGAM_h
#define readGAM_h
#include "TrailMix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Euka.h"
#include "soibean.h"
#include "readVG.h"
#include "HaploCart.h"
#include "gam2prof.h"
#include "Dup_Remover.h"
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "version.h"

#define VERBOSE
#define ADDDATA
//#define DEBUGREADGRAPH
using namespace std;
using namespace vg;
// using namespace google::protobuf;

void vg_preflight() {
    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    //choose_good_thread_count();
    temp_file::set_system_dir();
    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        exit(1);
                                      }
                    }


int main(int argc, char *argv[]) {

    vg_preflight();

    string_view usage=string("\n")+
        "   vgan is a suite of tools for pangenomics. \n"+
        "   We currently support two subcommands: euka (for the classification of eukaryotic taxa) \n"+
	"   And haplocart (for modern human mtDNA haplogroup classification). \n"+
        "   The underlying data structure is the VG graph. \n\n" +
        string(argv[0]) +" <command> options\n"+
        "\n"+
        "   Commands:\n"+
        "      euka         Classify eukaryotic taxa "+"\n"+
        "      haplocart    Predict human mitochondrial haplogroup  "+"\n"+
        "      soibean      Identify eukaryotic species " +"\n"+
        "      trailmix     Inference on ancient human mtDNA mixture                                        "+"\n"+
	"";

    if( argc==1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage "<<usage<<"\n";
        return 1;
    }


    if(string(argv[1]) == "euka"){
        Euka  euka_;
        const string cwdProg = getCWD(argv[0]);
        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<euka_.usage()<<"\n";
            return 1;
        }

        argv++;
        argc--;
        return euka_.run(argc, argv, cwdProg);

    }else{      if(string(argv[1]) == "haplocart"){
        Haplocart  haplocart_;
        const string cwdProg = getCWD(argv[0]);
        shared_ptr<Trailmix_struct> dta = make_unique<Trailmix_struct>();
        dta->cwdProg = cwdProg;

        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<haplocart_.usage()<<"\n";
            return 1;
        }

        argv++;
        argc--;

        haplocart_.run(argc, argv, dta);
        return 0; // Fix later

    }

        // // Mikkel code begins // //
    else if(string(argv[1]) == "gam2prof"){
        shared_ptr<Trailmix_struct> null_  = nullptr;
        Gam2prof  gam2prof_;
        if( argc==2 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<gam2prof_.usage()<<"\n";
            return 1;
        }
    string cwdProg=getCWD(argv[0]);
        argv++;
        argc--;
        return gam2prof_.run(argc, argv, cwdProg, null_);

    }
    // // Mikkel code ends // //

        else if(string(argv[1]) == "duprm"){
        shared_ptr dta = make_unique<Trailmix_struct>();
        Dup_Remover dup_remover;
        if( argc==2 || argc > 3 ||
            (argc == 3 && (string(argv[2]) == "-h" || string(argv[2]) == "--help") )
            ){
            cerr<<dup_remover.usage()<<"\n";
            return 1;
        }
        const string cwdProg=getCWD(argv[0]);
        const char *gamfile = argv[2];
        dup_remover.remove_duplicates(dta, gamfile);
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


    else{      if(string(argv[1]) == "trailmix"){
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

    }else {if (string(argv[1]) == "--path-supports"){

        // NOTE: This only works if the graph has one connected component.
        // TODO: Have a function check if the graph has one connected component.

        const auto _ = readVG(string(argv[2]));

        return 0;
                                                    }


     else{
        cerr<<"invalid command "<<string(argv[1])<<"\n";
        return 1;
	}}}
   return 0;
}
                                }
#endif
