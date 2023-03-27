
#include "Euka.h"
#include <random>

using namespace vg;

#define PRINTVEC(v) for (const auto &el : v){cerr << el << '\t';}cerr << endl<< endl;   
void Euka::map_giraffe(string fastq1filename, string fastq2filename, const int n_threads, bool interleaved,
                                                    const char * fifo_A, const vg::subcommand::Subcommand* sc,
                                                    const string &tmpdir, const string &cwdProg, const string &prefix){

    int retcode;
    vector<string> arguments;
    arguments.emplace_back("vg");
    arguments.emplace_back("giraffe");

    auto normal_cout = cout.rdbuf();
    ofstream cout(fifo_A);
    std::cout.rdbuf(cout.rdbuf());
    // Map with VG Giraffe to generate a GAM file
    string minimizer_to_use= prefix + ".min";

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
        arguments.emplace_back("-s");
        arguments.emplace_back("100");
        arguments.emplace_back("-D");
        arguments.emplace_back("400");
        arguments.emplace_back("-w");
        arguments.emplace_back("40");
        arguments.emplace_back("-v");
        arguments.emplace_back("2");
        arguments.emplace_back("-r");
        arguments.emplace_back("30");
        arguments.emplace_back("-t");
        arguments.emplace_back(to_string(n_threads));
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
               arguments.emplace_back("-s");
               arguments.emplace_back("100");
               arguments.emplace_back("-D");
               arguments.emplace_back("400");
               arguments.emplace_back("-w");
               arguments.emplace_back("40");
               arguments.emplace_back("-v");
               arguments.emplace_back("2");
               arguments.emplace_back("-r");
               arguments.emplace_back("30");
               arguments.emplace_back("-t");
               arguments.emplace_back(to_string(n_threads));


           if (interleaved) {
               arguments.emplace_back("-i");
                            }

               //PRINTVEC(arguments)

               char** argvtopass = new char*[arguments.size()];
               for (int i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                                    }

               auto* sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);
               auto normal_cerr = cerr.rdbuf();
               //std::cerr.rdbuf(NULL);
               (*sc)(arguments.size(), argvtopass);
               //std::cerr.rdbuf(normal_cerr);
               delete[] argvtopass;

    }

    std::cout.rdbuf(normal_cout);

}
