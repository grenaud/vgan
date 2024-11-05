#include "Euka.h"
#include <random>
#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << v[i] << '\t';}cerr << endl;

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
		       const string &tmpdir, const string &cwdProg, const string &prefix, const string &deam3pfreqE, const string &deam5pfreqE,
                       const double posterior_threshold){

    cerr << "Mapping reads..." << endl;

    int retcode;
    vector<string> arguments;
    arguments.emplace_back("vg");
    arguments.emplace_back("safari");

    auto normal_cout = cout.rdbuf();
    ofstream cout(fifo_A);
    std::cout.rdbuf(cout.rdbuf());

    // Map with VG Giraffe to generate a GAM file
    string minimizer_to_use= prefix + ".min";
    string rymer_to_use= prefix + ".ry";

    if (fastq1filename != "" && fastq2filename != "")
	{

	    arguments.emplace_back("-f");
	    arguments.emplace_back(fastq1filename);
	    arguments.emplace_back("-f");
	    arguments.emplace_back(fastq2filename);
	    arguments.emplace_back("-g");
	    arguments.emplace_back(prefix + ".gg");
	    arguments.emplace_back("-d");
	    arguments.emplace_back(prefix + ".dist");
	    arguments.emplace_back("-m");
	    arguments.emplace_back(minimizer_to_use);
	    arguments.emplace_back("-H");
	    arguments.emplace_back(prefix + ".gbwt");
	    arguments.emplace_back("-x");
	    arguments.emplace_back(prefix + ".og");
            arguments.emplace_back("-q");
            arguments.emplace_back(rymer_to_use);
            arguments.emplace_back("--deam-3p");
            arguments.emplace_back(deam3pfreqE);
            arguments.emplace_back("--deam-5p");
            arguments.emplace_back(deam5pfreqE);
            arguments.emplace_back("-j");
            arguments.emplace_back(to_string(posterior_threshold));
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
	    arguments.emplace_back("-g");
	    arguments.emplace_back(prefix + ".gg");
	    arguments.emplace_back("-d");
	    arguments.emplace_back(prefix + ".dist");
	    arguments.emplace_back("-m");
	    arguments.emplace_back(minimizer_to_use);
	    arguments.emplace_back("-H");
	    arguments.emplace_back(prefix + ".gbwt");
	    arguments.emplace_back("-x");
	    arguments.emplace_back(prefix + ".og");
            arguments.emplace_back("-q");
            arguments.emplace_back(rymer_to_use);
            arguments.emplace_back("--deam-3p");
            arguments.emplace_back(deam3pfreqE);
            arguments.emplace_back("--deam-5p");
            arguments.emplace_back(deam5pfreqE);
            arguments.emplace_back("-j");
            arguments.emplace_back(to_string(posterior_threshold));

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
