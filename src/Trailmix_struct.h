#pragma once
#include "forward_declarations.h"
#include "../dep/vg/src/subcommand/subcommand.hpp"
#include "AlignmentInfo.h"
#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "NodeInfo.h"
#include <Eigen/Dense>
#include "../dep/spimap/src/Tree.h"
#include "damage.h"
#include "utils.hpp"
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>
#include <cfloat>
//#include "MCMC.h"

using namespace Eigen;

struct Trailmix_struct {
    // Damage info object
    Damage dmg;

    // Boolean variables
    bool reads_already_processed = false;
    bool output_profs=false;
    bool webapp = false;
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
    bool graph_dir_specified = false;

    // Integer variables
    unsigned int n_threads = 1;
    unsigned int k = 1;
    unsigned int n_samples = 1;
    int auto_max = 5;
    int iter=1;
    int burnin = 0;
    int chains = 4;

    // Double variables
    double background_error_prob = 0.0001;

    // String variables
    string graph_prefix = "graph";
    string graph_dir = "../share/toy_hcfiles/";
    string graphfilename = graph_prefix + "graph.og";
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
    string tmpdir = "/tmp/";

    // RPVG parameters
    string rng_seed = "NONE";
    string rpvg_gamfilename;
    string mu = "125";
    string sigma = "0.000001";
    bool strand_specific = false;
    string strand_specific_library_type = "rf";

    // Named pipe variables
    int status, status1, status2, status3, rpvg_status, rpvg_ht_status, pack_status;
    pid_t wpid, wpid2, wpid3, rpvg_wpid, rpvg_ht_wpid, pid1, pid2, pid3, rpvg_pid, rpvg_ht_pid, pack_pid;
    string first_fifo, second_fifo, third_fifo;
    const char* fifo_A;
    const char* fifo_B;
    const char* fifo_C;

    // Alignments
    shared_ptr<vector<AlignmentInfo*>> algnvector = make_shared<vector<AlignmentInfo*>>();
    shared_ptr<vector<AlignmentInfo*>> current_algnvector = make_shared<vector<AlignmentInfo*>>();

    // Loaded data
    vector<string> MCMC_path_names;
    map<const string, int> pangenome_map;
    vector<NodeInfo*> nodevector;
    vector<double> mappabilities;
    vector<double> incorrect_mapping_vec;
    vector<string> path_names;
    vector<double> qscore_vec;
    vector<unsigned int> quality_scores;
    vector<vector<string>> hap_combos;
    vector<pair<string, double>> tpms;
    vector<unsigned int> dup_alignment_counts;
    pair<vector<double>, vector<vector<unsigned int>>> read_probs;
    vector<string> RPVG_hap_names;
    std::shared_ptr<gbwt::GBWT> gbwt;
    shared_ptr<gbwtgraph::GBWTGraph> gbwtgraph;
    vector<string> gbwt_paths;
    std::map<unsigned int, unsigned int> gbwt_to_path;
    std::map<unsigned int, unsigned int> path_to_gbwt;

    // Inference

    vector<AlignmentInfo*>* gam = NULL;
    long double read_coverage;
    vector<long double> log_likelihood_vec;
    unsigned int max_branch_combo_index;
    vector<vector<long double>> bc_dists;
    vector<long double> branch_combo_lls;
    unsigned int nbpaths;
    unsigned int** baseshift_data_array;
    unordered_map<int, int> path_node_map;
    unordered_map<int, int> node_path_map;
    vector<double> hap_combo_posteriors;
    vector<vector<spidir::Node*>> node_combos;
    vector<vector<unsigned int>> branch_combos;
    vector<unsigned int> best_branch_combo;
    shared_ptr<spidir::Tree> tree;
    vg::Mapping mppg;
    string graph_seq;
    string mapping_seq;
    vector<vector<bool>> path_supports;
    vector<unsigned int> node_counts;
    unsigned int minid = -1;
    unsigned int maxid = -1;
    long double max_ll_all = 42;
    std::vector<Matrix4d> sub_vec;
    vg::subcommand::Subcommand* sc = NULL;
    map<string, vector<string>> parents;
    map<string, vector<string>> children;
    bdsg::ODGI graph;
    map<unsigned int, vector<unsigned int>> mapping_supports;
    vector<tuple<char, char, unsigned int, unsigned int, long double>> to_increment;
    vector<long double> seed;
    vector<long double> optimized_branch_placements;
    vector<long double> optimized_branch_lengths;
    std::unordered_map<unsigned int, int> gbwt_increments;

    // Auto mode
    bool current_source=true;
    vector<bool> sources;
    vector<vector<vector<bool>>> all_sources;
    unsigned int source_pos;
    long double best_model_likelihood = -DBL_MAX;
    long double cur_model_likelihood = -DBL_MAX;
    vector<long double> best_branch_combo_of_best_model;
    vector<unsigned int> max_branch_combo_indices;
    vector<vector<vector<unsigned int>>> best_branch_combos;
    vector<vector<long double>> best_branch_placements;
    vector<vector<long double>> model_lls;
    unsigned int cur_within_k;
    bool first_of_new_k = true;
    vector<tuple<char, char, unsigned int, unsigned int, long double>> to_increment_ancient;
    size_t ancient_counts=0;
    size_t modern_counts=0;

    // Inferred values
    vector<vector<long double>> all_branch_placements;
    vector<vector<long double>> all_branch_lengths;
    vector<vector<double>> optimized_thetas;


    // Surjection
    vector<string> paths_to_surject;

                  };

