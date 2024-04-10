#include "TrailMix.h"
#include "Trailmix_struct.h"
//#include "trailmix_functions.h"
//#include "TrailMix.h"
#include "HaploCart.h"
#include <algorithm>
#include <thread>
#include <chrono>
#include <omp.h>
#include <numeric>
#include <variant>
#include <sys/wait.h>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "utils.hpp"
#include "tree.h"
#include "getLCAfromGAM.h"
#include <vg/io/vpkg.hpp>
//#include <gbwtgraph/gbwtgraph.h>
//#include <gbwtgraph/path_cover.h>

using namespace vg;

void run_vg_surject(shared_ptr<Trailmix_struct> &dta)
{
    if (dta->paths_to_surject.empty()){
        cerr << "No paths to surject" << endl;
        return;
                                      }

    string output_file_path = "../src/consensus.bam";

    cerr << "SURJECTING..." << endl;
    vector<string> sjarguments;

    sjarguments.emplace_back("vg");
    sjarguments.emplace_back("surject");
    sjarguments.emplace_back("-x");
    sjarguments.emplace_back(dta->graph_prefix + ".xg");
    sjarguments.emplace_back("-t");
    sjarguments.emplace_back(to_string(dta->n_threads));
    sjarguments.emplace_back("-b");
    //sjarguments.emplace_back(dta->tmpdir + "/to_surject.gam");
    sjarguments.emplace_back("../share/toy_hcfiles/sim.gam");
    sjarguments.emplace_back("-p");
    for (const auto & p : dta->paths_to_surject){
        sjarguments.emplace_back(p);
                                                }

    // Convert the vector of strings to a vector of C strings.
    std::vector<char*> argv;
    for (const auto& arg : sjarguments) {
        argv.push_back(const_cast<char*>(arg.c_str()));
    }

    // execvp expects a null pointer as the last element.
    argv.push_back(nullptr);

    pid_t pid = fork();
    if (pid == -1) {
        throw std::runtime_error("Error in fork");
    }

    if (pid == 0) {
        // This is the child process.

        // Open the output file. O_WRONLY means "open for writing" and O_CREAT means "create file if it does not exist".
        int output_fd = open(output_file_path.c_str(), O_WRONLY | O_CREAT, 0644); // 0644 means "user can read/write, others can read"
        if (output_fd == -1) {
            perror("open");
            exit(EXIT_FAILURE);
        }

        // Redirect stdout to the file.
        if (dup2(output_fd, STDOUT_FILENO) == -1) {
            perror("dup2");
            exit(EXIT_FAILURE);
        }

        // Run vg surject.
        execvp(argv[0], argv.data());

        // If execvp returns, there was an error.
        perror("execvp");
        exit(EXIT_FAILURE);
    } else {
        // This is the parent process. Wait for the child to finish.
        int status;
        waitpid(pid, &status, 0);

        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
            //std::cerr << "vg surject ran successfully" << std::endl;
        } else {
            std::cerr << "vg surject failed" << std::endl;
        }
    }
 exit(0);
}

Trailmix::Trailmix(){

}

Trailmix::~Trailmix(){

}

std::string Trailmix::usage() const {
    std::stringstream usage;

    usage << "Usage: vgan trailmix [options]\n\n"
          << "Deconvolution and phylogenetic placement of hominin sedaDNA mixtures.\n\n"
          << "Options:\n"
          << "Algorithm parameters:\n"
          << "\t-e [FLOAT]\t\tBackground error probability for FASTA input (default: 0.0001)\n"
          << "\t-f [STR]\t\tFASTA consensus input file\n"
          << "\t-fq1 [STR]\t\tFASTQ input file\n"
          << "\t-fq2 [STR]\t\tFASTQ second input file (for paired-end reads)\n"
          << "\t-g [STR]\t\tGAM input file\n"
          << "\t-i\t\t\tInput FASTQ (-fq1) is interleaved\n"
          << "\t-k\t\t\tNumber of distinct contributing haplogroups\n"
          << "\t-o [STR]\t\tOutput file (default: stdout)\n"
          << "\t--auto\t\t\tAutomatically infer the optimal number of distinct sources in the sample\n"
          << "Non-algorithm parameters:\n"
          << "\t-s [STR]\t\tSample name\n"
          << "\t-t\t\t\tNumber of threads\n"
          << "\t-v\t\t\tVerbose mode\n"
          << "\t-z\t\t\tTemporary directory (default: /tmp/)\n";

    return usage.str();
}

const int Trailmix::run(int argc, char *argv[], const string & cwdProg){

    cerr <<
 " ████████ ██████   █████  ██ ██      ███    ███ ██ ██   ██ \n"
 "    ██    ██   ██ ██   ██ ██ ██      ████  ████ ██  ██ ██   \n"
  "    ██    ██████  ███████ ██ ██      ██ ████ ██ ██   ███    \n"
   "    ██    ██   ██ ██   ██ ██ ██      ██  ██  ██ ██  ██ ██   \n"
   "    ██    ██   ██ ██   ██ ██ ███████ ██      ██ ██ ██   ██  \n" << endl;

    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    temp_file::set_system_dir();

    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        exit(1);
                                      }
    bool verbose = true;
    bool graph_dir_specified = false;
    bool webapp = false;
    bool compute_posteriors = true;
    bool output_profs=false;
    bool auto_mode=false;
    bool debug=false;
    bool strand_specific = false;
    string gamfilename, samplename, fastafilename, fastq1filename, fastq2filename, outputfilename, posteriorfilename;
    string mu="125";
    string sigma="0.00001";
    string strand_specific_library_type = "read1_forward";
    string tmpdir = "/tmp/";
    string rng_seed="NONE";
    string graph_dir = "../share/toy_hcfiles/";
    string graph_prefix = "graph";
    string graphfilename = graph_prefix + ".og";
    unsigned int n_threads = 1;
    double background_error_prob;
    int k=1;
    int auto_max=5;
    int chains=4;
    int burnin = 1;
    int iter=3;

    for(int i=1;i<(argc);++i){

    if(string(argv[i]) == "-e"){
            background_error_prob = stod(argv[i+1]);
            if (background_error_prob < 0 || background_error_prob > 1) {
                       throw std::runtime_error("Error, option -e is not a valid probability.");
                                                                        }
            continue;
                                }

    if(string(argv[i]) == "-g"){
            gamfilename = argv[i+1];
            continue;
                                }

    if(string(argv[i]) == "-f"){
            fastafilename = argv[i+1];
            samplename=fastafilename;
            continue;
                                }


    if(string(argv[i]) == "-fq1"){
            fastq1filename = argv[i+1];
            const int idx = fastq1filename.find_last_of("/");
            samplename=fastq1filename.substr(idx + 1);
            continue;
                                }

    if(string(argv[i]) == "-fq2"){
            fastq2filename = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "-sigma"){
            sigma = argv[i+1];
            continue;
                                   }

    if(string(argv[i]) == "-mu"){
            mu = argv[i+1];
            continue;
                                }

    if(string(argv[i]) == "-k"){
            k = stoi(argv[i+1]);
            if (k < 1) {throw std::runtime_error("[TrailMix] Error, k must be a positive integer");}
            continue;
                                 }

    if(string(argv[i]) == "-o"){
            outputfilename = argv[i+1];
            continue;
                                 }

     if(string(argv[i]) == "-r"){
            rng_seed = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "-s"){
            samplename = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "-t"){
            if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw std::runtime_error("[TrailMix] Error, invalid number of threads");}
            if (stoi(argv[i+1]) == -1) {n_threads = std::thread::hardware_concurrency();}
            else if (stoi(argv[i+1]) <= std::thread::hardware_concurrency()) {
                n_threads = stoi(argv[i+1]);
                                                                             }
            else {
                cerr << "[TrailMix] Warning, specified number of threads is greater than the number available. Using " << n_threads << " threads\n";
                n_threads = std::thread::hardware_concurrency();
                 }
            continue;
                                 }


    if(string(argv[i]) == "-v"){
            verbose = true;
            continue;
                               }

    if(string(argv[i]) == "-z"){
            tmpdir = argv[i+1];
            if (tmpdir.back() != '/') {tmpdir += '/';}
            continue;
                               }

    if(string(argv[i]) == "-w"){
            webapp = true;
            outputfilename = tmpdir+"haplogroups.html";
            posteriorfilename = tmpdir+"posteriors.html";
            continue;
                               }

    if(string(argv[i]) == "--output-profs"){
            output_profs=true;
            continue;
                                           }

   if(string(argv[i]) == "--debug"){
            debug=true;
            continue;
                                   }

   if(string(argv[i]) == "--auto"){
            auto_mode=true;
            continue;
                                          }

    if(string(argv[i]) == "--strand-specific"){
            strand_specific=true;
            continue;
                                            }

    if(string(argv[i]) == "--read1-reverse"){
            strand_specific_library_type = "read1_reverse";
            continue;
                                            }

     if(string(argv[i]) == "--auto-max"){
            auto_max = stoi(argv[i+1]);
            continue;
                                        }

     if(string(argv[i]) == "--iter" || string(argv[i]) == "--iterations"){
         iter=stoi(argv[i+1]);
         assert(iter >= 0);
         continue;
     }
     if(string(argv[i]) == "--burnin"){
         burnin=stoi(argv[i+1]);
         assert(burnin >= 0);
         continue;
     }
     if(string(argv[i]) == "--chains"){
         chains=stoi(argv[i+1]);
         assert(burnin >= 0);
         continue;
     }

               }

    if(!graph_dir_specified){
	graph_dir  =  cwdProg +  graph_dir;
    }

    // Load a bunch of stuff
    Haplocart hc;
    shared_ptr dta = make_unique<Trailmix_struct>();
    dta->to_increment.reserve(1000000);
    dta->running_trailmix=true; // Need this up here for loading path supports correctly
    hc.load_pangenome_map(dta);
    hc.precompute_incorrect_mapping_probs(dta);
    hc.load_path_names(dta);
    //hc.load_path_supports(dta);
    const int nbpaths = dta->path_names.size();
    hc.load_mappabilities(dta);
    const vector<double> qscore_vec = get_qscore_vec();
    vector<long double> log_likelihood_vec(nbpaths, 0);
    bool use_background_error_prob = false;
    bool is_consensus_fasta = false;

    // Make sure we did not load an empty file
    assert(qscore_vec.size() > 0);


    //////////////////////////////////////////////// POPULATE OUR STRUCT ///////////////////////////////////////////////////////////


   std::cerr << "Loading GBWT index..." << std::endl << std::flush;
   std::string gbwt_index_file = dta->graph_dir + dta->graph_prefix + ".gbwt";
   //std::string gbwtgraph_index_file = dta->graph_prefix + "graph.gg";
   dta->gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_index_file);
   //dta->gbwtgraph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwtgraph_index_file);
   std::cerr << "GBWT index loaded." << std::endl << std::flush;

/*
   ofstream out("gbwt_threads");
const gbwt::Metadata& metadata = dta->gbwt->metadata;
const std::vector<gbwt::PathName>& path_names = metadata.path_names;
for (std::size_t i = 0; i < metadata.paths(); ++i) {
  if (metadata.hasPathNames()) {
    gbwt::size_type path_length = dta->gbwt->extract(i).size();
    std::string thread_name = "_gbwt_" + metadata.sample(path_names[i].sample) + "_" +
                              metadata.contig(path_names[i].contig) + "_" +
                              std::to_string(path_names[i].phase) + "_" +
                              std::to_string(path_names[i].count);
    out << thread_name << "," << path_length << std::endl;
  }
}
exit(0);
*/

    dta->k = k;
    dta->qscore_vec = qscore_vec;
    dta->background_error_prob = 0.0001;
    dta->cwdProg = cwdProg;
    dta->mu = mu;
    dta->sigma = sigma;
    dta->strand_specific = strand_specific;
    dta->strand_specific_library_type = strand_specific_library_type;
    dta->webapp=webapp;
    dta->debug=debug;
    dta->auto_mode=auto_mode;
    dta->verbose=verbose;
    dta->output_profs=output_profs;
    dta->compute_posteriors=compute_posteriors;
    dta->n_threads=n_threads;
    dta->fastafilename=fastafilename;
    dta->samplename=samplename;
    dta->fastq1filename=fastq1filename;
    dta->fastq2filename=fastq2filename;
    dta->graphfilename=graphfilename;
    dta->tmpdir = tmpdir + "vgan_" + random_string(7) + "/";
    dta->gamfilename=gamfilename;
    dta->outputfilename=outputfilename;
    dta->posteriorfilename=posteriorfilename;
    dta->background_error_prob=background_error_prob;
    dta->nbpaths=nbpaths;
    dta->graph_prefix = graph_prefix;
    dta->graph_dir_specified = graph_dir_specified;
    dta->iter=iter;
    dta->burnin = burnin;
    dta->chains=chains;
    //dta->treePath = dta->graph_prefix + "/iqtree/haps.treefile";
    dta->treePath = dta->graph_dir + "/hominin_mts_prank.best.dnd";
    dta->output_profs=true;
    dta->rng_seed=rng_seed;
    auto tree = make_tree_from_dnd(dta);
    dta->tree = make_unique<spidir::Tree>(tree);
    Trailmix::create_path_node_map(dta);
    hc.readPathHandleGraph(dta);
    assert(dta->nodevector.size() > 0);
    std::vector<Eigen::Matrix4d> sub_vec;
    Eigen::Matrix4d sub_matrix = Eigen::Matrix4d::Zero();
    dta->sub_vec = std::vector<Eigen::Matrix4d>(dta->path_names.size(), sub_matrix);
    dta->graph.for_each_path_handle([&](const handlegraph::path_handle_t &path_handle) {
    std::string path_name = dta->graph.get_path_name(path_handle);
    dta->gbwt_paths.emplace_back(path_name);
    dta->auto_max = auto_max;
    });

/*

    if (dta->auto_mode) {
    dta->node_combos.clear();
    dta->branch_combos.clear();
    dta->branch_combo_lls.clear();
    for (unsigned int run = 0; run < dta->auto_max; ++run) {
        if (dta->debug && run != 1){continue;}
        cerr << "ON RUN NUMBER: " << run + 1 << endl;
        dta->first_of_new_k = true;
        dta->sources = vector<bool>(run + 1, false);
        dta->k = run + 1;
        generate_combinations(dta, run + 1, 0);
        dta->first_of_new_k = false;
        dta->tmpdir = tmpdir + "vgan_" + random_string(7) + "/";
        if (!fs::is_directory(dta->tmpdir) || !fs::exists(dta->tmpdir)) {
            std::filesystem::create_directory(dta->tmpdir);
        }
    }
}
    else {
        dta->current_source=false;
        Trailmix::run_trailmix(dta);
        //std::filesystem::remove_all(dta->tmpdir);
         }

    //Trailmix::write_output(dta);
    //run_vg_surject(dta);
*/

    run_haplocart(dta);
    //dta->gam = analyse_GAM_trailmix(dta);
    //remove(dta->fifo_A);
    Trailmix::run_trailmix(dta);


    return 0;
                                          }




