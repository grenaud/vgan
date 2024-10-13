#pragma once
#include "Trailmix_struct.h"
#include "sys/wait.h"
#include "libgab.h"
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "readVG.h"
#include "readGAM.h"
#include "MCMC.h"
#include "subcommand/subcommand.hpp"

using namespace std;

int main_pack(int argc, char** argv);
int main_call(int argc, char** argv);
int main_surject(int argc, char** argv);

class Trailmix{
private:

const char get_dummy_qual_score(const double &background_error_prob);
void write_fq_read(auto & dummyFASTQFile, int offset, const int window_size, const string &fastaseq, char dummyqualscore);
const string fa2fq(const string & fastaseq, const char & dummyqualscore, const string & tmpdir);
void map_giraffe(shared_ptr<Trailmix_struct> &dta);
const pair<vector<string>, vector<string>> read_fasta(const string & fastafilename);
void load_mappabilities(std::shared_ptr<Trailmix_struct> &dta);
void load_pangenome_map(shared_ptr<Trailmix_struct> &dta);
vector<string> paths_through_node(const bdsg::ODGI& graph, const bdsg::handle_t& node);
void run_mcmc(shared_ptr<Trailmix_struct>& dta);
void run_gam2prof(shared_ptr<Trailmix_struct>& dta);
const double get_likelihood_of_source(const shared_ptr<Trailmix_struct> &dta, const unsigned int &branch, const unsigned int &branch_parent, \
                                           const double &proposed_branch_position, double &theta);
const double sum_log_likelihoods(vector<double> &log_lik_vec);
void modifyPathNameInPlace(shared_ptr<Trailmix_struct> &dta, string &path_name, bool first=false);
void load_read_probs(shared_ptr<Trailmix_struct> &dta);
void precompute_incorrect_mapping_probs(shared_ptr<Trailmix_struct> &dta);
void readPathHandleGraph (shared_ptr<Trailmix_struct> &dta);
void get_seed_source_estimates(shared_ptr<Trailmix_struct> &dta, const vector<string> &sigpaths);
void compute_read_coverage(shared_ptr<Trailmix_struct> &dta);
void writeDeconvolvedReads(shared_ptr<Trailmix_struct> &dta, const string &input_file, const string &output_file, vector<int> reads_to_include);
void write_output(std::shared_ptr<Trailmix_struct> &dta);
void create_path_node_map(shared_ptr<Trailmix_struct> &dta);
void place_sources(shared_ptr<Trailmix_struct> &dta);
void infer(shared_ptr<Trailmix_struct> &dta);
void pack(shared_ptr<Trailmix_struct> &dta);
void load_hap_combos(shared_ptr<Trailmix_struct> &dta);
void load_tpms(shared_ptr<Trailmix_struct> &dta);
void initializeParams(RunTreeProportionParams &params, shared_ptr<Trailmix_struct>& dta);

const double get_branch_combo_posterior(shared_ptr<Trailmix_struct> &dta, unsigned int branch_combo_idx);

public:
    Trailmix();
    Trailmix(const Trailmix & other);
    ~Trailmix();
    Trailmix & operator= (const Trailmix & other);


    string usage() const;
    const int run(int argc, char *argv[], const string & cwdProg);
    void mcmc_setup(shared_ptr<Trailmix_struct> &dta);
    void run_trailmix(shared_ptr<Trailmix_struct> &dta);
    void readPHG(shared_ptr<Trailmix_struct> &dta);
};

int main_filter(int argc, char** argv);

