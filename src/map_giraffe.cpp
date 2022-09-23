#include "HaploCart.h"
#include <random>

//https://stackoverflow.com/questions/71890603/how-do-i-redirect-stderr-to-dev-null-in-c
struct null_streambuf: public std::streambuf
{
  using int_type = std::streambuf::int_type;
  using traits   = std::streambuf::traits_type;

  virtual int_type overflow( int_type value ) override
  {
    return value;
  }
};

const char Haplocart::get_dummy_qual_score(const double background_error_prob) {
    // Given a background error probability, return a dummy quality score for the artificial FASTQ reads
    string illumina_encodings = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
    const int Q = -10 * log10(background_error_prob);
    return illumina_encodings[Q];
                                                                           }


void Haplocart::map_giraffe(string &fastaseq, string &fastq1filename, string &fastq2filename, const int n_threads, bool interleaved,
                                    double background_error_prob, const string & samplename, const char * fifo_A, const vg::subcommand::Subcommand* sc,
                                    const string & tmpdir, const string &graph_dir_path, const bool quiet){



    if (!quiet) {cerr << "Mapping reads..." << endl;}

    int retcode;
    vector<string> arguments;
    arguments.emplace_back("vg");
    arguments.emplace_back("giraffe");

    auto normal_cout = cout.rdbuf();
    ofstream cout(fifo_A);
    std::cout.rdbuf(cout.rdbuf());
    // Map with VG Giraffe to generate a GAM file
    string minimizer_to_use=  graph_dir_path + "k31_w11.min";


    int nonbase_count = 0;
    for (auto base : fastaseq) {
        if (base != 'A' && base != 'C' && base != 'T' && base != 'G' && base != 'a' && base != 'c' && base != 't' && base != 'g') {
            ++nonbase_count;
                                                                      }
                               }


    if (nonbase_count > 7999) {
        if (!quiet) {cerr << "Detecting many ambiguous bases, using alternative minimizer index..." << endl;}
        minimizer_to_use=  graph_dir_path + "k17_w18.min";
                              }


if (fastq1filename != "" && fastq2filename != "")
    {

	//@josh what is this? I commented
        // if(fastq1filename.front() != '/'){fastq1filename = getFullPath(cwdProg + fastq1filename);}
        // if(fastq2filename.front() != '/'){fastq1filename = getFullPath(cwdProg + fastq2filename);}

        arguments.emplace_back("-f");
        arguments.emplace_back(fastq1filename);
        arguments.emplace_back("-f");
        arguments.emplace_back(fastq2filename);
        arguments.emplace_back("-Z");
        arguments.emplace_back(getFullPath(graph_dir_path + "graph.giraffe.gbz"));
        arguments.emplace_back("-d");
        arguments.emplace_back(getFullPath(graph_dir_path + "graph.dist"));
        arguments.emplace_back("-m");
        arguments.emplace_back(minimizer_to_use);
        arguments.emplace_back("-b");
        arguments.emplace_back("fast");
        char** argvtopass = new char*[arguments.size()];
        for (int i=0;i<arguments.size();i++) {
            argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                             }

        auto* sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);
        auto normal_cerr = cerr.rdbuf();
        std::cerr.rdbuf(NULL);
        (*sc)(arguments.size(), argvtopass);
        std::cerr.rdbuf(normal_cerr);

    }

else if (fastq1filename != "" && fastq2filename == "")
    {
		//@josh what is this? I commented

	//if(fastq1filename.front() != '/'){fastq1filename = getFullPath(cwdProg + fastq1filename);}

               arguments.emplace_back("-f");
               arguments.emplace_back(fastq1filename);
               arguments.emplace_back("-Z");
               arguments.emplace_back(getFullPath(graph_dir_path + "graph.giraffe.gbz"));
               arguments.emplace_back("-d");
               arguments.emplace_back(getFullPath(graph_dir_path + "graph.dist"));
               arguments.emplace_back("-m");
               arguments.emplace_back(minimizer_to_use);
               arguments.emplace_back("-b");
               arguments.emplace_back("fast");

           if (interleaved) {
               arguments.emplace_back("-i");
                            }

               char** argvtopass = new char*[arguments.size()];
               for (int i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                                    }

               auto* sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);
               auto normal_cerr = cerr.rdbuf();
               std::cerr.rdbuf(NULL);
               (*sc)(arguments.size(), argvtopass);
               std::cerr.rdbuf(normal_cerr);
               delete[] argvtopass;

    }

else if (fastaseq != ""){
    string fastapath;
    const char dummyq = Haplocart::get_dummy_qual_score(background_error_prob);

    string fasta_cmd;
    const string dummy_fastq_file = fa2fq(fastaseq, dummyq, tmpdir);


    arguments.emplace_back("-f");
    arguments.emplace_back(dummy_fastq_file);
    arguments.emplace_back("-Z");
    arguments.emplace_back(getFullPath(graph_dir_path + "graph.giraffe.gbz"));
    arguments.emplace_back("-d");
    arguments.emplace_back(getFullPath(graph_dir_path + "graph.dist"));
    arguments.emplace_back("-m");
    arguments.emplace_back(minimizer_to_use);
    arguments.emplace_back("-b");
    arguments.emplace_back("fast");

    char** argvtopass = new char*[arguments.size()];
    for (int i=0;i<arguments.size();i++) {
            argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                         }

    if (sc == NULL) {
            sc = vg::subcommand::Subcommand::get(arguments.size(), argvtopass);
                    }

    auto normal_cerr = cerr.rdbuf();
    std::cerr.rdbuf(NULL);
    (*sc)(arguments.size(), argvtopass);
    std::cerr.rdbuf(normal_cerr);

    delete[] argvtopass;
    remove(dummy_fastq_file.c_str());
   }

    std::cout.rdbuf(normal_cout);
    if (!quiet) {std::cerr << "Reads mapped" << endl;}
}

