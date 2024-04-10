#pragma once
#include "Trailmix_struct.h"
#include "sys/wait.h"
#include "libgab.h"
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "readVG.h"
#include "readGAM.h"
#include "subcommand/subcommand.hpp"

using namespace std;

int main_pack(int argc, char** argv);
int main_call(int argc, char** argv);
int main_surject(int argc, char** argv);

class Trailmix{
private:

void create_node_combos(shared_ptr<Trailmix_struct>& dta);
void run_mcmc(shared_ptr<Trailmix_struct>& dta);
void run_haplocart(shared_ptr<Trailmix_struct>& dta);
void run_gam2prof(shared_ptr<Trailmix_struct>& dta);
const long double get_likelihood_of_source(const shared_ptr<Trailmix_struct> &dta, const unsigned int &branch, const unsigned int &branch_parent, \
                                           const long double &proposed_branch_position, long double &theta);
const long double sum_log_likelihoods(vector<long double> &log_lik_vec);
void load_read_probs(shared_ptr<Trailmix_struct> &dta);
void get_seed_source_estimates(shared_ptr<Trailmix_struct> &dta, const vector<unsigned int> &branch_combo);
void get_likelihood_of_branch_combo(shared_ptr<Trailmix_struct> &dta, const vector<unsigned int> &branch_combo, unsigned int branch_combo_idx, unsigned int N);
void compute_read_coverage(shared_ptr<Trailmix_struct> &dta);
void write_output(std::shared_ptr<Trailmix_struct> &dta);
void create_path_node_map(shared_ptr<Trailmix_struct> &dta);
void place_sources(shared_ptr<Trailmix_struct> &dta);
void infer(shared_ptr<Trailmix_struct> &dta);
void load_hap_combos(shared_ptr<Trailmix_struct> &dta);
void get_posterior(shared_ptr<Trailmix_struct> &dta);
void remove_duplicate_branch_combos(vector<vector<unsigned int>> & branch_combos);
void get_dists_on_branch(const std::vector<unsigned int> &branch_combo, std::shared_ptr<Trailmix_struct> &dta);
void compute_optimal_placement(shared_ptr<Trailmix_struct> &dta);
void enumerate_branch_combos(shared_ptr<Trailmix_struct> &dta);
void branch_combo_cartesian_product( vector<vector<unsigned int> >& v , shared_ptr<Trailmix_struct> &dta);
void generate_combinations(shared_ptr<Trailmix_struct> &dta, int num_sources, int current_pos=0);
void load_tpms(shared_ptr<Trailmix_struct> &dta);
void pack(shared_ptr<Trailmix_struct> &dta);
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
};

int main_filter(int argc, char** argv);
