#include "Euka.h"
#include "index_registry.hpp"
#include <random>


//#define DEBUGGIRAFFE
using namespace vg;

const char get_dummy_qual_score(double &background_error_prob)             {
    // Given a background error probability, return a dummy quality score for the artificial FASTQ reads
    string illumina_encodings = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
    const int Q = -10 * log10(background_error_prob);
    return illumina_encodings[Q];
                                                                           }


void Euka::map_giraffe(string fastq1filename, string fastq2filename, const int n_threads, bool interleaved,
		       const char * fifo_A, const vg::subcommand::Subcommand* sc,
		       const string &tmpdir, const string &cwdProg, const string &prefix, const string &minprefix){
    cerr << "Mapping reads..." << endl;

    int retcode;
    vector<string> arguments;
    arguments.emplace_back("vg");
    arguments.emplace_back("giraffe");

    auto normal_cout = cout.rdbuf();
    ofstream cout(fifo_A);
    std::cout.rdbuf(cout.rdbuf());

    // Map with VG Giraffe to generate a GAM file
    string minimizer_to_use= minprefix + ".min";

    // giraffe's zipcode-based seed clustering (new in this vg version) always
    // needs a zipcode-annotated minimizer index. Whenever the plain `-m`
    // file given above lacks a matching "<prefix>.shortread.zipcodes" file
    // next to it (true for euka_db.min/Ursidae.min/etc, which predate
    // zipcodes), giraffe silently DISCARDS it and rebuilds a fresh
    // minimizer+zipcode pair from the graph itself using
    // IndexingParameters::short_read_minimizer_k/w -- hardcoded to 29/11
    // (vg's generic default) with NO command-line flag to override it. This
    // means whatever k/w the on-disk .min file was actually built with is
    // irrelevant at runtime. The original production euka_db.min (confirmed
    // by reading its on-disk header before the v1.75 migration) was built
    // with k=20/w=10, matching this project's own short-read/ancient-DNA
    // convention used elsewhere (soibean's make_graph_files.sh, HaploCart's
    // k17_w18.min) -- vg's 29/11 default requires reads of at least 39bp to
    // find even one minimizer, silently dropping shorter/more-degraded
    // fragments that are exactly euka/soibean's signal of interest. Setting
    // these here (IndexingParameters is a plain public static struct, safe
    // to set before any in-process giraffe invocation) restores that
    // sensitivity without needing a CLI flag vg doesn't provide.
    IndexingParameters::short_read_minimizer_k = 20;
    IndexingParameters::short_read_minimizer_w = 10;

    if (fastq1filename != "" && fastq2filename != "")
	{

	    arguments.emplace_back("-f");
	    arguments.emplace_back(fastq1filename);
	    arguments.emplace_back("-f");
	    arguments.emplace_back(fastq2filename);
	    arguments.emplace_back("-Z");
	    arguments.emplace_back(prefix + ".giraffe.gbz");
	    arguments.emplace_back("-d");
	    arguments.emplace_back(prefix + ".dist");
	    arguments.emplace_back("-m");
	    arguments.emplace_back(minimizer_to_use);
	    char** argvtopass = new char*[arguments.size()];
	    for (int i=0;i<arguments.size();i++) {
		argvtopass[i] = const_cast<char*>(arguments[i].c_str());
	    }

	    auto* sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);	    
	    auto normal_cerr = cerr.rdbuf();
	    //std::cerr.rdbuf(NULL);
	    (*sc)(arguments.size(), argvtopass);
	    //std::cerr.rdbuf(normal_cerr);

	}

    else if (fastq1filename != "" && fastq2filename == "")
	{
	    arguments.emplace_back("-f");
	    arguments.emplace_back(fastq1filename);
	    arguments.emplace_back("-Z");
	    arguments.emplace_back(prefix + ".giraffe.gbz");
	    arguments.emplace_back("-d");
	    arguments.emplace_back(prefix + ".dist");
	    arguments.emplace_back("-m");
	    arguments.emplace_back(minimizer_to_use);

	    if (interleaved) {
		arguments.emplace_back("-i");
	    }

	    char** argvtopass = new char*[arguments.size()];
	    for (int i=0;i<arguments.size();i++) {
		argvtopass[i] = const_cast<char*>(arguments[i].c_str());
#ifdef DEBUGGIRAFFE		
		cerr<<"argvtopass["<<i<<"] = "<<argvtopass[i] <<endl;
#endif
	    }

	    auto* sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);
	    auto normal_cerr = cerr.rdbuf();
	    //std::cerr.rdbuf(NULL);
	    (*sc)(arguments.size(), argvtopass);
	    //std::cerr.rdbuf(normal_cerr);
	    delete[] argvtopass;

	}

    std::cout.rdbuf(normal_cout);
    std::cerr << "Reads mapped" << endl;
}




