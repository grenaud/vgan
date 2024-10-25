#include "TrailMix.h"
#include "Trailmix_struct.h"
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
#include "vgan_utils.h"
#include "tree.h"
#include "precompute.h"
#include <vg/io/vpkg.hpp>

using namespace vg;
#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << v[i] << '\t';}cerr << endl;

#define PRINT_MAP_TO_TSV(map, filename) \
do { \
    std::ofstream outFile(filename); \
    if (!outFile) { \
        std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl; \
        break; \
    } \
    for (const auto& pair : map) { \
        outFile << pair.first << '\t' << pair.second << '\n'; \
    } \
    outFile.close(); \
} while (0)


std::vector<bool> parseBoolVector(const std::string &str, char delimiter) {
    std::vector<bool> result;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        if (token == "1") {
            result.push_back(true);
        } else if (token == "0") {
            result.push_back(false);
        } else {
            throw std::runtime_error("Invalid value in source assignment vector. Must be '0' or '1'.");
        }
    }
    return result;
}

void printBoolVector(const std::vector<bool> &vec, const std::string &name) {
    std::cerr << name << ": [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cerr << (vec[i] ? "true" : "false");
        if (i != vec.size() - 1) {
            std::cerr << ", ";
        }
    }
    std::cerr << "]" << std::endl;
}


// Assuming main_surject is the entry point for vg surject similar to main_gamsort
int main_surject(int argc, char** argv);

void printMatrix(const Matrix& matrix) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cerr << matrix[i][j];
            if (j < 3) std::cerr << "\t"; // Tab-delimited
        }
        std::cerr << std::endl; // New line at the end of each row
    }
}

void normalizeMatrix(Matrix& matrix) {
    for (int i = 0; i < 4; ++i) {
        double rowTotal = 0.0;
        // Calculate the sum of each row
        for (int j = 0; j < 4; ++j) {
            rowTotal += matrix[i][j];
        }

        // Normalize each element in the row if the total is greater than 0
        if (rowTotal > 0.0) {
            for (int j = 0; j < 4; ++j) {
                matrix[i][j] /= rowTotal;
            }
        }
    }
}

void writeMatrixToFile(const std::string& filename, const Matrix& matrix) {
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        // Header for the columns
        outfile << "\tA\tC\tG\tT\n"; // Tab-delimited header

        const char nucleotides[] = {'A', 'C', 'G', 'T'};
        for (int i = 0; i < 4; ++i) {
            outfile << nucleotides[i]; // Row label
            for (int j = 0; j < 4; ++j) {
                outfile << "\t" << matrix[i][j]; // Tab-delimited values
            }
            outfile << std::endl; // New line at the end of each row
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}



// function gets the average for per clade damage profils (specific for either c->t or g->a as well as 5' or 3' end). 
vector <double> get_avg(vector<double> dam, int l, int no_clades){

    vector <double> sum(l, 0.0);
    vector<double> av; 
    
    for (int i = 0; i<dam.size(); i++){

        sum[i%l] += dam.at(i); 
    }

    for (int j = 0; j<sum.size(); j++){
        av.emplace_back(sum[j]/no_clades);
    }
    
    return av;
}

void output_combined_damagefile(
    const std::vector<double>& end5ct,
    const std::vector<double>& end3ct,
    const std::vector<double>& end5ga,
    const std::vector<double>& end3ga,
    int lengthToProf,
    const std::vector<std::string>& clade_id_list,
    const std::string&  TM_outputfilename
) {
    // Calculate averages
    std::vector<double> end5ct_av = get_avg(end5ct, lengthToProf, clade_id_list.size());
    std::vector<double> end3ct_av = get_avg(end3ct, lengthToProf, clade_id_list.size());
    std::vector<double> end5ga_av = get_avg(end5ga, lengthToProf, clade_id_list.size());
    std::vector<double> end3ga_av = get_avg(end3ga, lengthToProf, clade_id_list.size());

    // Output for 5' profile
    std::ofstream outdamall(TM_outputfilename + "_5p.prof", std::ios::trunc);
    outdamall << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";
    for (std::size_t pas = 0; pas < end5ct_av.size(); pas++) {
        for (int col = 0; col < 12; col++) {
            if (col == 5) {
                outdamall << end5ct_av[pas] << '\t';
            } else if (col == 11) {
                outdamall << "0\n";
            } else {
                outdamall << "0\t";
            }
        }
    }

    // Output for 3' profile
    std::ofstream outdamall2(TM_outputfilename + "_3p.prof", std::ios::trunc);
    outdamall2 << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";

    std::reverse(end3ga_av.begin(), end3ga_av.end());

    for (std::size_t pas = 0; pas < end3ga_av.size(); pas++) {
        for (int col = 0; col < 12; col++) {
            if (col == 6) {
                outdamall2 << end3ga_av[pas] << '\t';
            } else if (col == 11) {
                outdamall2 << "0\n";
            } else {
                outdamall2 << "0\t";
            }
        }
    }
}

void Trailmix::load_pangenome_map(shared_ptr<Trailmix_struct> &dta){
const string pangenome_map_path = getFullPath(dta->graph_dir+"parsed_pangenome_mapping");
igzstream myfile;
myfile.open(pangenome_map_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
    {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     const string key = tokens[0];
     const int val = stoi(tokens[1])+1;
     dta->pangenome_map.insert(make_pair(key, val));
    }
assert(!dta->pangenome_map.empty());
return;
}

void Trailmix::precompute_incorrect_mapping_probs(shared_ptr<Trailmix_struct> &dta){
for (int Q=0; Q!=100; ++Q) {
    dta->incorrect_mapping_vec.emplace_back(get_p_incorrectly_mapped(Q));
                           }
return;
}

void Trailmix::load_mappabilities(std::shared_ptr<Trailmix_struct> &dta) {
    std::string mappability_path = getFullPath(dta->graph_dir + "mappability.tsv");
    igzstream myfile;
    myfile.open(mappability_path.c_str(), ios::in);
    std::string line;
    while (std::getline(myfile, line)) {
        std::vector<std::string> tokens = allTokensWhiteSpaces(line);
        if (tokens.size() != 4) {
            throw std::runtime_error("Invalid line in mappability file: " + line);
        }
        int start_pos = std::stoi(tokens[1]);
        int end_pos = std::stoi(tokens[2]);
        const double mappability = std::stod(tokens[3]);
        if (start_pos >= end_pos) {
            throw std::runtime_error("Invalid start/end positions in mappability file: " + line);
        }
        dta->mappabilities.resize(end_pos, mappability);
        std::fill(dta->mappabilities.begin() + start_pos, dta->mappabilities.begin() + end_pos, mappability);
        assert(!dta->mappabilities.empty());
    }
    return;
}

void Trailmix::modifyPathNameInPlace(shared_ptr<Trailmix_struct> &dta, string &path_name, bool first){
     string original_path_name = path_name;

    // Special case: a single letter followed by underscore
    std::regex pattern2("^([A-Za-z])_$");
    path_name = std::regex_replace(path_name, pattern2, "$1");

    // Handle the special path format
    std::regex pattern3("\\+([0-9]+)\\+\\(([0-9]+)\\)");
    path_name = std::regex_replace(path_name, pattern3, "_$1__$2_");

    // Now replace special characters
    std::replace(path_name.begin(), path_name.end(), '+', '_');
    std::replace(path_name.begin(), path_name.end(), '\'', '_');
    std::replace(path_name.begin(), path_name.end(), '*', '_');
    std::replace(path_name.begin(), path_name.end(), '@', '_');
    std::replace(path_name.begin(), path_name.end(), '(', '_');
    std::replace(path_name.begin(), path_name.end(), ')', '_');

    if (!path_name.empty() && path_name.back() == '*') {
        path_name.back() = '_';
    }

   if (path_name.size() == 1){path_name = path_name + "_";}

    // Remove .1 or .2 at the end
    std::regex patternEnd("\\.(1|2)$");
    path_name = std::regex_replace(path_name, patternEnd, "");

    if(first){
    dta->originalPathNames[path_name] = original_path_name;
             }

}

void Trailmix::readPHG(shared_ptr<Trailmix_struct> &dta) {

    cerr << "[TrailMix] Deserializing graph: " << dta->graph_dir + dta->graph_prefix + ".og" << endl;
    dta->graph.deserialize(dta->graph_dir + dta->graph_prefix + ".og");
    cerr << "[TrailMix] Done deserializing" << endl;
    dta->minid = dta->graph.min_node_id();
    dta->maxid = dta->graph.max_node_id();

    // create path support file
    size_t N_nodes = dta->graph.get_node_count();

    unsigned int p_index = 0;
    dta->graph.for_each_path_handle([&](const handlegraph::path_handle_t &path_handle) {
    string path_name = dta->graph.get_path_name(path_handle);
    modifyPathNameInPlace(dta, path_name, true);

size_t pos = 0;
        dta->path_names.emplace_back(path_name);
        p_index++;
    });

    dta->nbpaths = dta->path_names.size();

    for(int64 i=dta->minid;i<=dta->maxid;++i){
        NodeInfo * nodetoadd = new NodeInfo(i,dta->nbpaths,0);
        const auto nodehandle = dta->graph.get_handle(i);
        string seqtoadd = dta->graph.get_sequence(nodehandle);
        nodetoadd->seq = seqtoadd;
        dta->nodevector.emplace_back(move(nodetoadd));
    }

if (dta->nodevector.size() == 0) {
    throw std::runtime_error("Error: nodevector is empty. Expected non-empty nodevector.");
}

if (dta->minid == 0) {
    throw std::runtime_error("Error: minid is zero. Expected non-zero minid.");
}

if (dta->maxid == 0) {
    throw std::runtime_error("Error: maxid is zero. Expected non-zero maxid.");
}

if (dta->path_names.size() == 0) {
    throw std::runtime_error("Error: path_names is empty. Expected non-empty path_names.");
}

   return;
}

void run_vg_surject(shared_ptr<Trailmix_struct> &dta, const string &path_to_surject) {
    cout << "[TrailMix] Starting run_vg_surject function with system call" << endl;

    // Extract filename from path_to_surject
    size_t last_slash_idx = path_to_surject.find_last_of("\\/");
    string surject_filename = (last_slash_idx == string::npos) ? path_to_surject : path_to_surject.substr(last_slash_idx + 1);

    // Prepare file paths
    string output_file_path = "./" + surject_filename + ".bam";  // Output in current directory with surject name
    string vg_binary_path = "../dep/vg/bin/vg";

    // Prepare the system command
    stringstream command_stream;
    command_stream << vg_binary_path << " surject"
                   << " -x " << dta->graph_dir << dta->graph_prefix << ".xg"
                   << " -t " << dta->n_threads
                   << " -b"  // BAM output
                   << " -p " << path_to_surject
                   << " -l -P -A " << "./" + surject_filename + "_"+to_string(dta->k)+".gam"
                   << " > " << output_file_path;

    string command = command_stream.str();
    //cout << "Executing command: " << command << endl;

    // Execute the system call
    int result = system(command.c_str());

    // Check result
    if (result == 0) {
        cout << "vg surject ran successfully" << endl;
    } else {
        cout << "vg surject failed with error code " << result << endl;
    }
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
          << "\t-fq1 [STR]\t\t\t     FASTQ input file\n"
          << "\t-fq2 [STR]\t\t\t     FASTQ second input file (for paired-end reads)\n"
          << "\t-g [STR]\t\t\t     GAM input file\n"
          << "\t-i\t\t\t\t     Input FASTQ (-fq1) is interleaved\n"
          << "\t-k\t\t\t\t     Number of distinct contributing haplogroups\n"
          << "\t-o [STR]\t\t\t     Output file (default: stdout)\n"
          << "\t-pt [FLOAT]\t\t\t	Posterior threshold for ancient DNA (aligned with SAFARI) \n"
          << "\t--depth [INT]\t\t\t     Tree depth to whitelist the MCMC search space (-1 for no restriction)\n"
          << "Non-algorithm parameters:\n"
          << "\t--dbprefix <prefix>\t\t     Specify the prefix for the database files\n"
          << "\t--tm-files [STR]\t\t      Specify the TrailMix file directory (default: \"../share/vgan/tmfiles/\") \n"
          << "\t--save-gam [STR]\t\t      Save GAM file to specified location (e.g. for reuse in contamination mode) \n"
          << "\t-t\t\t\t\t     Number of threads\n"
          << "\t-v\t\t\t\t     Verbose mode\n"
          << "\t-z\t\t\t\t     Temporary directory (default: /tmp/)\n"
          << "Markov chain Monte Carlo options:\n"
          << "  \t--chains [INT]\t\t             Define the number of chains for the MCMC (default: 4)\n"
          << "  \t--iter [INT]\t\t             Define the number of iterations for the MCMC (default: 1.000.000)\n"
          << "  \t--randStart [bool]\t\t          Set to get random starting nodes in the tree instead of the signature nodes (default: false)\n"
          << "  \t--burnin [INT]\t\t             Define the burn-in period for the MCMC (default: 100.000)\n"
          << "Initialization options:\n"
          << "  \t--library-type [STR] \t\t     Strand-specific library type (fr: read1 forward, rf: read1 reverse) (default: unstranded)\n"
          << "Contamination Mode Options:\n"
          << "  \t--contamination-mode, --cont-mode    Enable contamination mode, requires exactly two sources\n"
          << "  \t--is-ancient [BOOL_VECTOR]   Assignments of sources, represented as a boolean vector, comma-delimited, where 1 means ancient (e.g. 0,1).\n";

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
    bool debug=false;
    string gamfilename, samplename, fastafilename, fastq1filename, fastq2filename, posteriorfilename;
    string out_gam_filename = "sample.gam";
    string TM_outputfilename = "TM_out";
    string tmpdir = "/tmp/";
    string graph_dir = getFullPath(cwdProg+"../share/vgan/tmfiles/");
    string graph_prefix="graph";
    unsigned int n_threads = 1;
    bool graphdirspecified = false;
    std::vector<bool> isAncient = {true, false};
    string deam5pfreqE  = getFullPath(cwdProg+"../share/vgan/damageProfiles/none.prof");
    string deam3pfreqE  =  getFullPath(cwdProg+"../share/vgan/damageProfiles/none.prof");
    string none_prof = getFullPath(cwdProg+"../share/vgan/damageProfiles/none.prof");
    string posterior_threshold = "0.5";
    bool specifiedDeam=false;
    bool run_mcmc=true;
    bool cont_mode=false;
    bool save_gam=false;
    unsigned int iter=10000;
    unsigned int burnin=100;
    unsigned int chains=1;
    unsigned int k=1;
    string prof_out_file_path = getFullPath(cwdProg+"../");
    int lengthToProf = 5;
    int depth = 10;
    bool randStart=false;
    bool strand_specific = false;
    string rng_seed="NONE";
    string strand_specific_library_type = "";
    bool k_unset=true;

     for(int i=1;i<(argc);++i){

     if(string(argv[i]) == "-g"){
            gamfilename = argv[i+1];
            continue;
                                }

     if(string(argv[i]) == "-f"){
            fastafilename = argv[i+1];
            samplename    = fastafilename;
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

     if(string(argv[i]) == "-l"){
            lengthToProf = stoi(argv[i+1]);
            continue;
        }

     if(string(argv[i]) == "--out_dir"){
            prof_out_file_path = string(argv[i+1]);
            continue;
        }

       if(std::string(argv[i]) == "-k"){
            k_unset=false;
            k = std::stoi(argv[i+1]);
            if (k < 1) {
                throw std::runtime_error("[TrailMix] Error, k must be a positive integer");
            }
            continue;
      }


       if(string(argv[i]) == "-o"){
            TM_outputfilename = argv[i+1];
            continue;
                               }

       if(string(argv[i]) == "--library-type"){
            strand_specific = true;
            strand_specific_library_type = argv[i+1];
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

        if(string(argv[i]) == "-z"){
            tmpdir = getFullPath(argv[i+1]);
            continue;
                               }

    if(string(argv[i]) == "--debug"){
            debug=true;
            continue;
                                   }

        
     if(string(argv[i]) == "--depth"){
            depth=stod(argv[i+1]);
            continue;
                                   }



        if(string(argv[i]) == "--deam5p"  ){
            deam5pfreqE=string(argv[i+1]);
        specifiedDeam=true;
            continue;
        }

        if(string(argv[i]) == "--deam3p"  ){
            deam3pfreqE=string(argv[i+1]);
        specifiedDeam=true;
            continue;
        }

        if(string(argv[i]) == "--randStart"  || string(argv[i]) == "--randstart"){
            randStart=true;
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

        if(string(argv[i]) == "-pt"){
            posterior_threshold = argv[i+1];
            continue;
                                }

         if(string(argv[i]) == "--dbprefix"){
            graph_prefix=argv[i+1];
            continue;
        }


    if(string(argv[i]) == "--tm-files"){
            graph_dir=getFullPath(argv[i+1]);
            continue;
                                      }

     if(string(argv[i]) == "--save-gam"){
            out_gam_filename= argv[i+1];
            save_gam = true;
            continue;
                                        }

       if(string(argv[i]) == "--contamination-mode" || string(argv[i]) == "--cont-mode"){
            if (k != 2){
              throw runtime_error("[TrailMix] Must have two sources to run in contamination mode");
            }
            cont_mode=true;
            continue;
        }

       if (std::string(argv[i]) == "--is-ancient") {
        if(i + 1 < argc) {
            isAncient = parseBoolVector(argv[++i], ','); // Assume parseBoolVector is correct
        } else {
            throw std::runtime_error(std::string("Option ") + argv[i] + " requires an argument.");
        }
        continue;
    }

    }

    shared_ptr dta = make_unique<Trailmix_struct>();
    dta->graph_dir = graph_dir;

    std::cerr << "Loading GBWT index..." << std::endl << std::flush;
    std::string gbwt_index_file = getFullPath(dta->graph_dir + graph_prefix + ".gbwt");
    std::string gbwtgraph_index_file = getFullPath(dta->graph_dir + graph_prefix + ".gg");
    dta->gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_index_file);
    dta->gbwtgraph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwtgraph_index_file);
    std::cerr << "GBWT index loaded." << std::endl << std::flush;

    // Load a bunch of stuff
    //Haplocart hc;

    dta->running_trailmix=true;
    dta->save_gam=save_gam;
    dta->out_gam_filename = out_gam_filename;
    dta->is_ancient_vec = isAncient;
    dta->strand_specific = strand_specific;
    dta->strand_specific_library_type = strand_specific_library_type;
    dta->fastafilename=fastafilename;
    dta->samplename=samplename;
    dta->fastq1filename=fastq1filename;
    dta->fastq2filename=fastq2filename;
    dta->randStart = randStart;
    dta->none_prof=none_prof;
    dta->tmpdir = tmpdir + "vgan_" + random_string(7) + "/";
    dta->rpvg_gamfilename = dta->tmpdir + random_string(7);
    if (!fs::is_directory(dta->tmpdir) || !fs::exists(dta->tmpdir)) {
      std::filesystem::create_directory(dta->tmpdir);
    }
    dta->gamfilename=gamfilename;
    dta->TM_outputfilename=TM_outputfilename;
    dta->n_threads=n_threads;
    dta->k = k;
    dta->graphdirspecified = graphdirspecified;
    dta->treePath = dta->graph_dir + "haps.treefile";
    auto tree = make_tree_from_dnd(dta);

    if (tree.nodes.size() == 0){
        throw runtime_error("The tree is empty");
    }

        for (size_t i = 0; i < tree.nodes.size(); i++) {
        auto & tn = tree.nodes[i];

        if (tn == NULL){cerr << "TREE NODE IS NULL!" << endl; continue;}

    }

    dta->cont_mode=cont_mode;
    dta->tree = &tree;
    dta->qscore_vec = get_qscore_vec();
    dta->running_trailmix=true;
    dta->iter=iter;
    dta->burnin=burnin;
    dta->chains=chains;
    dta->posterior_threshold=posterior_threshold;
    dta->depth=depth;
    dta->graph_prefix=graph_prefix;
    dta->graph_dir = graph_dir;
    readPHG(dta);
    dta->deam5pfreqE = deam5pfreqE;
    dta->deam3pfreqE = deam3pfreqE;
    dta->dmg.initDeamProbabilities(deam5pfreqE,deam3pfreqE);
    dta->dmg_none.initDeamProbabilities(none_prof, none_prof);
    load_mappabilities(dta);
    load_pangenome_map(dta);

for (int Q=0; Q!=100; ++Q) {
    dta->incorrect_mapping_vec.emplace_back(get_p_incorrectly_mapped(Q));
                           }

    Trailmix::run_mcmc(dta);

/*
    // Check the validity of the paths in dta->path_node_map
    for (const auto& path : dta->paths_to_surject) {
        if (dta->path_node_map.find(path) == dta->path_node_map.end()) {
            std::cerr << "Error: Path " << path << " not found in path_node_map." << std::endl;
            throw std::runtime_error("Path not found in path_node_map.");
        } else {
            std::cerr << "Path " << path << " found in path_node_map!" << std::endl;
        }
    }

// Step 1: Map Reads to All Paths
std::map<int, std::set<std::string>> readToPathMap; // Maps read indices to the set of paths they support
for (const auto& path : dta->paths_to_surject) {
    for (size_t index = 0; index < dta->read_probs.second.size(); ++index) {
        const auto& node_vector = dta->read_probs.second[index];
        int node_to_check = dta->path_node_map[path];
        if (std::find(node_vector.begin(), node_vector.end(), node_to_check) != node_vector.end()) {
            readToPathMap[index].insert(path);
        }
    }
}

// Step 2: Surject Reads Uniquely for Each Path
for (const auto& path : dta->paths_to_surject) {
    std::vector<int> reads_to_surject;
    Matrix aggregated_matrix = {}; // Initialize your matrix

    for (const auto& readMapEntry : readToPathMap) {
        int readIndex = readMapEntry.first;
        const auto& pathsSupported = readMapEntry.second;

        // Check if the read supports only the current path
        if (pathsSupported.count(path) > 0) {
            reads_to_surject.push_back(readIndex);

            // Aggregate substitution matrices for the reads
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    aggregated_matrix[i][j] += dta->substitution_matrices[readIndex][i][j];
                }
            }
        }
    }

std::cerr << "Aggregated Matrix (Before Normalization):" << std::endl;
printMatrix(aggregated_matrix); // Implement this function to print the matrix

    // Normalize the aggregated matrix here
    normalizeMatrix(aggregated_matrix);

std::cerr << "Normalized Matrix:" << std::endl;
printMatrix(aggregated_matrix); // Print matrix after normalization

    // Write the normalized matrix to a file
    std::string matrix_filename = dta->originalPathNames[path] + "_substitution_matrix.txt";
    writeMatrixToFile(matrix_filename, aggregated_matrix);

    cerr << "READ PROBS SIZE: " << dta->read_probs.second.size() << endl;
    cerr << "NUMBER OF READS SURJECTED: " << reads_to_surject.size() << endl;

    writeDeconvolvedReads(dta, dta->rpvg_gamfilename, dta->originalPathNames[path]+"_"+to_string(dta->k)+".gam", reads_to_surject);
    run_vg_surject(dta, dta->originalPathNames[path]);
}
*/

    //exit(0);

    return 0;
}
