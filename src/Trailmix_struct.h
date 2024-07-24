#pragma once
#include "forward_declarations.h"
#include "../dep/vg/src/subcommand/subcommand.hpp"
#include "AlignmentInfo.h"
#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "NodeInfo.h"
#include "../dep/spimap/src/Tree.h"
#include "damage.h"
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>
#include <cfloat>

// Define a 4x4 matrix
using Matrix = std::array<std::array<double, 4>, 4>;

struct Trailmix_struct {

    Damage dmg;

    // Boolean variables
    bool include_read_props=true;
    bool reads_already_processed = false;
    bool output_profs=false;
    bool webapp = false;
    bool cont_mode=false;
    bool verbose = false;
    bool quiet = false;
    bool debug=false;
    bool use_background_error_prob = false;
    bool compute_posteriors = false;
    bool is_consensus_fasta = false;
    bool invoked_samplename = false;
    bool running_trailmix = false;
    bool running_single_source = false;
    bool interleaved=false;
    bool dump_json=false;
    bool auto_mode=false;
    bool graphdirspecified = false;
    bool randStart=false;

    // Integer variables
    unsigned int n_threads = 1;
    unsigned int k = 2;
    unsigned int n_samples = 1;
    int iter=100000;
    int burnin = 1000;
    int chains = 4;
    int depth=3;

    // Double variables
    double background_error_prob = 0.0001;

    // String variables
    string graph_dir = "../share/tmfiles/";
    string graph_prefix = "graph";
    string graphfilename = graph_prefix + ".og";
    string treePath;

    // Giraffe parameters
    string cwdProg;
    string fastafilename;
    string fastaseq;
    string samplename;
    string fastq1filename;
    string fastq2filename;
    string gamfilename;
    string posteriorfilename;
    string jsonfilename;
    string outputfilename = "/dev/stdout";
    string TM_outputfilename;
    string tmpdir = "/tmp/";

    // Named pipe variables
    int status, status1, status2, status3, rpvg_status, rpvg_ht_status, pack_status;
    pid_t wpid, wpid2, wpid3, pid1, pid2, pid3, rpvg_pid, rpvg_ht_pid, pack_pid;
    string first_fifo, second_fifo, third_fifo;
    const char* fifo_A;
    const char* fifo_B;
    const char* fifo_C;

    // Alignments
    shared_ptr<vector<AlignmentInfo*>> algnvector = make_shared<vector<AlignmentInfo*>>();
    shared_ptr<vector<AlignmentInfo*>> current_algnvector = make_shared<vector<AlignmentInfo*>>();

        // RPVG parameters
    string rng_seed = "NONE";
    string rpvg_gamfilename;
    bool strand_specific = false;
    string strand_specific_library_type = "rf";

    // Loaded data
    vector<vector<bool>> path_supports;
    string deam3pfreqE;
    string deam5pfreqE;
    std::vector<Matrix> substitution_matrices;
    std::unordered_map<std::string, std::string> originalPathNames;
    set<string> in_pruned_set;
    int n_leaves=0;
    vector<vector<string>> hap_combos;
    vector<pair<string, double>> tpms;
    pair<vector<double>, vector<vector<unsigned int>>> read_probs;
    vector<string> RPVG_hap_names;
    vector<vector<spidir::Node*>> node_combos;
    map<const string, int> pangenome_map;
    vector<NodeInfo*> nodevector;
    vector<double> mappabilities;
    vector<double> incorrect_mapping_vec;
    vector<string> path_names;
    vector<double> qscore_vec;
    vector<unsigned int> quality_scores;
    vector<unsigned int> dup_alignment_counts;
    std::shared_ptr<gbwt::GBWT> gbwt;
    shared_ptr<gbwtgraph::GBWTGraph> gbwtgraph;
    vector<string> gbwt_paths;
    std::map<unsigned int, unsigned int> gbwt_to_path;
    std::map<unsigned int, unsigned int> path_to_gbwt;

    // Inference

    double read_coverage;
    vector<double> log_likelihood_vec;
    unsigned int max_branch_combo_index;
    unsigned int nbpaths;
    unsigned int** baseshift_data_array;
    spidir::Tree* tree;
    vg::Mapping mppg;
    string graph_seq;
    string mapping_seq;
    vector<vector<string>> nodepaths;
    unordered_map<string, int> path_node_map;
    unordered_map<int, string> node_path_map;
    unsigned int minid = 1;
    unsigned int maxid = -1;
    const vg::subcommand::Subcommand* sc = NULL;
    map<string, vector<string>> parents;
    map<string, vector<string>> children;
    bdsg::ODGI graph;
    vector<double> seed;
    std::unordered_map<unsigned int, int> gbwt_increments;
    vector<bool> source_assignments = vector<bool>{true, false}; // True is ancient, false is modern

    // Surjection
    vector<string> paths_to_surject;

                  };

