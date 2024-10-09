#pragma once
#include "forward_declarations.h"
#include "Clade.h"
#include <map>
#include "Trailmix_struct.h"
//#include "tree.h"
#include <vector>

//#define DEBUGHKY
//#define DEBUGPAIN3

using namespace std;

struct RunTreeProportionParams {
    spidir::Node* root = NULL;
    vector<vector<double>> probMatrix;
    spidir::Tree* tr; // Changed to a regular pointer
    vector<unsigned int> sources;
    unsigned int maxIter = 50000;
    unsigned int burn = 5000;
    unsigned int chains = 4;
    bool soibean = true;
    const vector<AlignmentInfo*>* align;
    double logLike;
    std::array<double, 256> freqs;

    // Constructor that takes a pointer to the tree
    RunTreeProportionParams(spidir::Tree* tree) : tr(tree) {}
};


typedef struct PosTree {
    spidir::Node* pos;
    double pos_branch;
    double theta;
    double branch_place_anc;
    double branch_place_der;

} PosTree;


typedef struct MCMCiteration {
    int n_components;
    vector<double> proportions;
    vector<double> max_branch_lens;
    vector<PosTree> positions_tree;
    double logLike;
} MCMCiteration;




class MCMC {

double proposal_sd = 0.01; // proposal SD
double kappa = 1/22;

private:

public:
    void run_trailmix(shared_ptr<Trailmix_struct> &dta);
    const vector<double> generate_proposal(vector<double> &current_vec, const double &alpha, const bool branch_pos);
    const std::vector<double> sample_normal_euka(std::vector<double>& x, double alpha);
    std::vector<double> sample_normal(std::vector<double>& x, double alpha);
    void generateNumbers(double num1);

    template<typename T> const vector<T> softmax(vector<T> &log_vec);
    double calculateLogWeightedAverage(double logValueChild, double weightChild, double logValueParent, double weightParent);
    vector<double> run(int iter, int burnin, double tol, const vector<double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
    const double get_proposal_likelihood(const vector <double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
    std::vector<MCMCiteration> run_tree_proportion(RunTreeProportionParams &params, std::vector<MCMCiteration> &state_t_vec, const bdsg::ODGI& graph, \
                             const vector<vector<string>> &nodepaths, string num, shared_ptr<Trailmix_struct> &dta, bool running_trailmix, int chain);
    void validateInputs(const PosTree &current_position, const double move_distance);
    inline spidir::Node* findLCA(const spidir::Tree* tree, spidir::Node* node1, spidir::Node* node2);
    //bool is_in_pruned(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used);
    //bool is_in_pruned_verbose(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used);
    bool is_in_pruned_set(spidir::Node* p, shared_ptr<Trailmix_struct> &dta);
    void add_nodes_at_depth(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta);
    //const void updatePosition(std::shared_ptr<spidir::Tree> tr, PosTree &current_position, const double move_distance, bool move_forward);
    inline void updatePosition(PosTree &current_position, const double move_distance, bool move_forward);
    const double calculateDistanceToAncestor(spidir::Node* startNode, spidir::Node* ancestor);
    const double getSumPatristicDistances(const PosTree &current_position, const shared_ptr<std::shared_ptr<spidir::Tree>> tr);
    double calculateRhat(const std::vector<double>& means, const std::vector<double>& variances, int chainLength, int numchains);
    inline const std::vector<double> getPatristicDistances(spidir::Tree* tr, spidir::Node* node, int numofLeafs, double posonbranch);
    inline const double calculateEuclideanDistance(const std::vector<double> &vec1, const std::vector<double> &vec2);
    pair<unordered_map<string, vector<vector<double>>>, double> processMCMCiterations(shared_ptr<Trailmix_struct>& dta, const std::vector<MCMCiteration> &MCMCiterations, int k, const string &num, int chain, \
                                                                                        spidir::Tree* tr, int numofleafs);
    double moveForward(PosTree &current_position, double move_distance_abs);
    void checkPositionErrors(const PosTree &current_position, double move_distance);
    inline const double calculateDistanceToLeaf(const std::shared_ptr<spidir::Tree> &tr, const spidir::Node* current, const spidir::Node* leaf, double currentDistance, \
                               std::unordered_set<const spidir::Node*>& visited);
    //void processMCMCiterations(const std::vector<std::vector<MCMCiteration>> &MCMCiterationsVec, int k, string num);


  inline const double computeBaseLogLike(shared_ptr<Trailmix_struct>& dta, const AlignmentInfo* read, RunTreeProportionParams &params, const int basevec, \
                                   const int base, const string &pathName, const double t, double branch_len, bool cont_mode)
    {

        const auto &detail = read->detailMap.at(pathName).at(basevec).at(base);

        const double purinfreq = params.freqs['R'];
        const double pyrinfreq = params.freqs['Y'];
        const double mu = params.freqs['M'];
        const double obaseFreq = params.freqs[detail.readBase];
        const char refb = detail.referenceBase;
        const char readb = detail.readBase;

        if (readb == 'S' || refb == 'S' || readb == 'N' || refb == 'N' || readb == '-' || refb == '-'){return 1e-9;}

#ifdef DEBUGPAIN3  // Debugging block start
        cerr << "Debugging Information START" << endl;

        cerr << "Input Variables:" << endl;
        cerr << "basevec: " << basevec << endl;
        cerr << "base: " << base << endl;
        cerr << "pathName: " << pathName << endl;
        cerr << setprecision(14)<< "t: " << t << endl;
        cerr << setprecision(14)<< "branch_len: " << branch_len << endl;

        cerr << "Intermediate Variables:" << endl;
        cerr<< setprecision(14) << "purinfreq: " << purinfreq << endl;
        cerr<< setprecision(14) << "pyrinfreq: " << pyrinfreq << endl;
        cerr<< setprecision(14) << "mu: " << mu << endl;
        cerr << setprecision(14)<< "obaseFreq: " << obaseFreq << endl;
        cerr << setprecision(14)<< "refb: " << refb << endl;
        cerr << setprecision(14)<< "readb: " << readb << endl;

#endif
        

        double probBaseHKY[4];
        for (int bpd = 0; bpd < 4; bpd++) {
            probBaseHKY[bpd] = 0.0;
        }

        for (int bpo = 0; bpo < 4; bpo++) {
            char rb = "ACGT"[bpo];
            //cerr << "base " << rb << " iteration " << bpo << endl;
            if (rb == refb) {
                // no mutation
                if (rb == 'A' || rb == 'G') {
                    const double A = 1 + purinfreq * (kappa - 1);
                    //cerr << "A "  << A << endl;
                    const double jut1 = params.freqs[rb] + params.freqs[rb] * ((1 / purinfreq) - 1) * exp(-(mu * t));
                    //cerr << "jut1 " << jut1 << endl;
                    const double jut11 = ((purinfreq - params.freqs[rb]) / purinfreq) * exp(-(mu * t * A));
                    //cerr << "jut11 " << jut11 << endl;
                    probBaseHKY[bpo] = jut1 + jut11;
                    if(probBaseHKY[bpo] < 1e-8){
                        probBaseHKY[bpo] = 1e-8;
                    }
                    if (isnan(probBaseHKY[bpo]) || isinf(probBaseHKY[bpo])|| probBaseHKY[bpo] < 1e-8){
                        cerr << "log like A & G match " << probBaseHKY[bpo] << " for base " << refb << endl;
                        throw runtime_error("HKY is invalid for log like A & G match");
                    }

                    

                }
                // case 2: we have a match but the base is a Pyrimidine
                else if (rb == 'C' || rb == 'T') {
                    const double A = 1 + pyrinfreq * (kappa - 1);
                    //cerr << "A "  << A << endl;
                    const double jut1 = params.freqs[rb] + params.freqs[rb] * ((1 / pyrinfreq) - 1) * exp(-(mu * t));
                    //cerr << "jut1 " << jut1 << endl;
                    const double jut11 = ((pyrinfreq - params.freqs[rb]) / pyrinfreq) * exp(-(mu * t * A));
                    //cerr << "jut11 " << jut11 << endl;
                    probBaseHKY[bpo] = jut1 + jut11;
                    if(probBaseHKY[bpo] < 1e-8){
                        probBaseHKY[bpo] = 1e-8;
                    }
                    if (isnan(probBaseHKY[bpo]) || isinf(probBaseHKY[bpo])|| probBaseHKY[bpo] < 1e-8){
                        cerr << "log like C & T match " << probBaseHKY[bpo] << " for base " << refb << endl;
                        throw runtime_error("HKY is invalid for log like C & T match");
                    }
                    

                }
            } else {
                // case 1: we have a mismatch and the base is a Purine
                if ((rb == 'A' && refb == 'G') || (rb == 'G' && refb == 'A')) {
                    const double A = 1 + purinfreq * (kappa - 1);
                    //cerr << "A "  << A << endl;
                    const double jut1 = params.freqs[rb] + params.freqs[rb] * ((1 / purinfreq) - 1) * exp(-(mu * t));
                    //cerr << "jut1 " << jut1 << endl;
                    const double jut11 = (params.freqs[rb] / purinfreq) * exp(-(mu * t * A));
                    //cerr << "jut11 " << jut11 << endl;
                    if(jut1 > jut11){
                        probBaseHKY[bpo] = jut1 - jut11;
                    }else{
                        probBaseHKY[bpo] = jut11 - jut1;
                    }
                    if(probBaseHKY[bpo] < 1e-8){
                        probBaseHKY[bpo] = 1e-8;
                    }
                    
                    if (isnan(probBaseHKY[bpo]) || isinf(probBaseHKY[bpo])|| probBaseHKY[bpo] < 1e-8){
                        cerr << "log like C & T match " << probBaseHKY[bpo] << " for base " << refb << endl;
                        throw runtime_error("HKY is invalid for log like C & T match.");
                    }


                    

                }
                // case 2: we have a mismatch and the base is a Pyrimidine
                else if ((rb == 'C' && refb == 'T') || (rb == 'T' && refb == 'C')) {
                    const double A = 1 + pyrinfreq * (kappa - 1);
                    //cerr << "A "  << A << endl;
                    const double jut1 = params.freqs[rb] + params.freqs[rb] * ((1 / pyrinfreq) - 1) * exp(-(mu * t));
                    //cerr << "jut1 " << jut1 << endl;
                    const double jut11 = (params.freqs[rb] / pyrinfreq) * exp(-(mu * t * A));
                    //cerr << "jut11 " << jut11 << endl;
                    if(jut1 > jut11){
                        probBaseHKY[bpo] = jut1-jut11;
                    }else{
                        probBaseHKY[bpo] = jut11-jut1;
                    }
                    
                    if(probBaseHKY[bpo] < 1e-8){
                        probBaseHKY[bpo] = 1e-8;
                    }
                    if (isnan(probBaseHKY[bpo]) || isinf(probBaseHKY[bpo])|| probBaseHKY[bpo] < 1e-8){
                        cerr << "log like C & T mismatch " << probBaseHKY[bpo] << " for base " << refb << endl;
                        throw runtime_error("HKY is invalid for log like C & T mismatch.");
                    }
                    

                } else {
                    probBaseHKY[bpo] = params.freqs[rb] * (1 - exp(-(mu * t)));
                    if(probBaseHKY[bpo] < 1e-8){
                        probBaseHKY[bpo] = 1e-8;
                    }

                    if (isnan(probBaseHKY[bpo]) || isinf(probBaseHKY[bpo]) || probBaseHKY[bpo] < 1e-8){
                        cerr << "log like the trash " << probBaseHKY[bpo] << " for base " << refb << endl;
                        throw runtime_error("HKY loglikemarg is invalid.");
                    }

                }
            }
        } // end filling probability matrix.

        double log_lik_marg = -std::numeric_limits<double>::infinity();
        for (int bpd = 0; bpd < 4; bpd++) {
            if ("ACGT"[bpd] == readb) {
                log_lik_marg = oplusInitnatl(log_lik_marg, (log(probBaseHKY[bpd])));
                
            } else {
                log_lik_marg = oplusInitnatl(log_lik_marg, (log(probBaseHKY[bpd])));
                
            }
        }
        if (log_lik_marg > 1e-8){
            log_lik_marg = log(0.999999999);
        }


        if (isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
            cerr << "post HKY:" << endl;
            for (int bpd = 0; bpd < 4; bpd++) {
                cerr << setprecision(15) << bpd << "\t" << probBaseHKY[bpd] << endl;
            }
            cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
            cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;
            cerr << "UNUSUAL HKY" << endl;
            throw runtime_error("HKY loglikemarg is invalid.");
                                                                              }

#ifdef DEBUGHKY
        cerr << "post HKY:" << endl;
        for (int bpd = 0; bpd < 4; bpd++) {
            cerr << setprecision(15) << bpd << "\t" << probBaseHKY[bpd] << endl;
        }
        cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
        cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;
        //if (readb != refb){throw std::runtime_error("Mismatch ");}
#endif

        if (isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
            cerr << "cont mode? " << cont_mode << endl;
            cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
            cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;
            cerr << setprecision(15) << "detail map log like no damage " << detail.logLikelihoodNoDamage << endl;

            throw runtime_error("HKY loglikemarg is invalid.");}

//cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;

        return log_lik_marg;
    }



   
    const int factorial(int n){
    if (n < 0) {
        // Invalid input, return an error value or throw an exception
        // depending on your requirements
        // For simplicity, let's return -1.0
        return -1.0;
    } else if (n == 0) {
        // Base case: factorial of 0 is 1
        return 1.0;
    } else {
        // Recursive case: factorial of n is n multiplied by factorial of (n-1)
        return static_cast<double>(n) * factorial(n - 1);
    }
                             };

    unsigned long long nChooseK(int n, int k) {
    if (k > n) {
        return 0;
    }

    // Optimization to reduce the number of multiplications
    if (k > n/2) {
        k = n - k;
    }

    unsigned long long result = 1;
    for (int i = 1; i <= k; ++i) {
        if (result > ULLONG_MAX / n) {
            // Handle overflow, if needed
            throw std::overflow_error("Integer overflow in nChooseK");
        }
        result *= n--;
        result /= i;
    }
    return result;
}

void get_proposal_sd(double& proposal_sd, double& acceptance_rate, int current_iteration, int total_iterations, int burn_in) {
    std::random_device rd;
    std::mt19937 gen(rd());

    // Decide which range to choose
    std::uniform_int_distribution<> range_choice(0,5);
    int choice = range_choice(gen);

    // Choose uniformly within the selected range
    if (choice == 0) {
        std::uniform_real_distribution<> range(1e-4, 1e-3); // Very small range
        proposal_sd = range(gen);
    } else if (choice == 1) {
        std::uniform_real_distribution<> range(1e-2, 1e-1); // Small range
        proposal_sd = range(gen);
    } else if (choice == 2) {
        std::uniform_real_distribution<> range(1e-1, 1.0); // Medium-small range
        proposal_sd = range(gen);
    } else if (choice == 3) {
        std::uniform_real_distribution<> range(1.0, 10.0); // Medium-large range
        proposal_sd = range(gen);
    } else if (choice == 4) {
        std::uniform_real_distribution<> range(10.0, 100.0); // Large range
        proposal_sd = range(gen);
    }
    else if (choice == 5) {
        std::uniform_real_distribution<> range(300.0, 800.0); // Extra Large range
        proposal_sd = range(gen);
    }
}

// Generate thetaVec and max_branch_lens
std::pair<std::vector<double>, std::vector<double>> generateThetaVecAndMaxBranchLens(const std::vector<PosTree>& current_positions) {
    std::vector<double> thetaVec;
    std::vector<double> max_branch_lens;

    for (auto& p : current_positions) {
        thetaVec.emplace_back(max(0.1, p.theta));
        max_branch_lens.emplace_back(p.pos->dist);
        //max_branch_lens.emplace_back(0.0001);
    }

    return make_pair(thetaVec, max_branch_lens);
}

/*
// Generate a vector of random numbers
std::vector<double> generateRandomNumbers(int size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> random_numbers(size);
    double sum = 0.0;

    for (double &num : random_numbers) {
        num = dis(gen);
        sum += num;
    }

    for (double &num : random_numbers) {
        num /= sum;
    }

    return random_numbers;
}
*/

std::vector<double> generateRandomNumbers(int size) {
    std::vector<double> distribution(size, 0.0);

    if (size > 0) {
        double value = 1.0 / size;
        for (double &num : distribution) {
            num = value;
        }
    }

    return distribution;
}



// Initialize the MCMCiteration state
MCMCiteration initializeState(RunTreeProportionParams& params, vector<double> &seed) {
    MCMCiteration state;
    state.n_components = params.sources.size();

    // Generate random numbers
    std::vector<double> random_numbers = generateRandomNumbers(state.n_components);

    // Initialize current_positions
    state.positions_tree = initializePositions(random_numbers, params, seed);

    // Generate thetaVec and max_branch_lens
    auto [thetaVec, max_branch_lens] = generateThetaVecAndMaxBranchLens(state.positions_tree);

    state.proportions = thetaVec;
    state.max_branch_lens = max_branch_lens;

    //std::vector<double> LLvec_init = get_LLvec_init(state, params);
    double likelihood_t = params.logLike;
    state.logLike = likelihood_t;

    return state;
}

// Initialize current_positions
std::vector<PosTree> initializePositions(const std::vector<double>& random_numbers, RunTreeProportionParams& params, vector<double> &seed) {
    std::vector<PosTree> current_positions(random_numbers.size());
    int index = 0;

    for (auto& p : current_positions) {
        p.pos = params.tr->nodes[params.sources[index]];
        p.pos_branch = 0.0001; //1.0 / random_numbers.size();
        double seed_sum = std::accumulate(seed.begin(), seed.end(), 0.0);
        if (seed_sum == 1.0) {
            p.theta = seed[index];
                             }
        else  {
            p.theta = random_numbers[index];
              }
        p.branch_place_anc = 0.5;
        p.branch_place_der = 0.5;
        index++;
    }

    return current_positions;
}


/*
std::vector<PosTree> initializePositions(const std::vector<double>& random_numbers, RunTreeProportionParams& params, vector<double> &seed) {
    std::vector<PosTree> current_positions(random_numbers.size());
    double equal_theta = 1.0 / random_numbers.size();  // Equal division of theta for each source

    for (int index = 0; index < random_numbers.size(); ++index) {
        current_positions[index].pos = params.tr->nodes[params.sources[index]];
        current_positions[index].pos_branch = 1.0 / random_numbers.size();
        current_positions[index].theta = equal_theta;  // Assign equal theta to each source
        current_positions[index].branch_place_anc = 0.5;
        current_positions[index].branch_place_der = 0.5;
    }

    return current_positions;
}
*/

    MCMC();
    MCMC(const MCMC & other);
    ~MCMC();
    MCMC & operator= (const MCMC & other);
};

// computes the softmax of the proposal vector.
// normalises so each proposal vector sums up to 1
template <typename T> const vector<T> MCMC::softmax(vector<T> &log_vec){
        double K = 0.0;

        for (size_t i = 0; i<log_vec.size(); i++){

                K += exp(log_vec.at(i));
#ifdef DEBUGGENERATEVEC
                cerr << "sum to k " << K << endl;
#endif

        }
#ifdef DEBUGGENERATEVEC
        cerr << "K " << K << endl;
#endif

        vector<T> transformed_log_vec;
        for(size_t i = 0; i < log_vec.size(); i++){
                transformed_log_vec.emplace_back(exp(log_vec.at(i))/K);
        }

        return transformed_log_vec;
}


