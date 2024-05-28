#pragma once
#include "forward_declarations.h"
#include "Clade.h"
#include <map>
#include <Eigen/Dense>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
//#include "TrailMix.h"
#include "Trailmix_struct.h"
//#include "tree.h"
#include <vector>



//#define DEBUGHKY
//#define DEBUGPAIN3
//#define DEBUGDIAGNOSTIC
using namespace std;


struct RunTreeProportionParams {
    spidir::Node* root = NULL;
    vector<vector<double>> probMatrix;
    spidir::Tree* tr; // Changed to a regular pointer
    vector<int> sources;
    bool RPVG = false;
    unsigned int maxIter = 1000000;
    unsigned int burn = 100000;
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

double proposal_sd = 1.0; // proposal SD
double kappa = 1/22;
double factor = pow(15.0, 15);


private:


public:
    void run_rpvg_haplotype_transcripts(shared_ptr<Trailmix_struct>& dta);
    void run_trailmix(shared_ptr<Trailmix_struct> &dta);
    const vector<double> generate_proposal(vector<double> &current_vec, const double &alpha, const bool branch_pos);
    double getQuantile(const std::vector<double>& sortedData, double q);
    const double calculateDistanceToLeaf(const spidir::Tree* tr, const spidir::Node* current, const spidir::Node* leaf, double currentDistance, \
                               std::unordered_set<const spidir::Node*>& visited);
    pair<unordered_map<string, vector<vector<double>>>, double> processMCMCiterations(
                             shared_ptr<Trailmix_struct> &dta, const vector<MCMCiteration>& MCMCiterations,int k,const string &num,int chain,\
                             spidir::Tree* tr, int numofleafs, \
                             bool include_read_props);
    const double calculateEuclideanDistance(const std::vector<double> &vec1, const std::vector<double> &vec2);
    double calculateRhat(const std::vector<double>& means, const std::vector<double>& variances, int chainLength, int numChains);
    const double calculateDistanceToAncestor(spidir::Node* startNode, spidir::Node* ancestor);
    bool is_in_pruned(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used);
    bool is_in_pruned_verbose(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used);
    const std::vector<double> sample_normal_euka(std::vector<double>& x, double alpha);
    std::vector<double> sample_normal(std::vector<double>& x, double alpha);
    void generateNumbers(double num1);
    template<typename T> const vector<T> softmax(vector<T> &log_vec);
    void run_rpvg_haplotypes(shared_ptr<Trailmix_struct>& dta);
    const vector<double> run(int iter, int burnin, double tol, const vector<double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
    double get_proposal_likelihood(const vector <double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
    std::vector<MCMCiteration> run_tree_proportion(RunTreeProportionParams params, std::vector<MCMCiteration> state_t_vec, const bdsg::ODGI& graph, \
                                                     vector<vector<string>> nodepaths, string num, int n_threads, int numPaths, int chainindex, double con);
    std::vector<MCMCiteration> run_tree_proportion_TM(RunTreeProportionParams &params, std::vector<MCMCiteration> &state_t_vec, const bdsg::ODGI& graph, \
                               const vector<vector<string>> &nodepaths, string num, shared_ptr<Trailmix_struct> &dta, bool running_trailmix, int chain);

    void updatePosition(PosTree &current_position, double move_distance, bool move_forward);
    inline const std::vector<double> getPatristicDistances(const spidir::Tree* tr, spidir::Node* node, int numofLeafs, double posonbranch);

    double log_diff_exp(double logA, double logB) {
        if (logB >= logA) {
            throw std::runtime_error("logB must be less than logA");
        }
        else if (logA == -INFINITY && logB == -INFINITY) {
            return -INFINITY;
        }
        else {
            return logA + log1p(-exp(logB - logA));
        }
    }


 inline const double computeBaseLogLike(const AlignmentInfo* read, RunTreeProportionParams &params, const int basevec, \
                                   const int base, const string &pathName, const double t)
    {

        const auto &detail = read->detailMap.at(pathName).at(basevec).at(base);

        const double purinfreq = params.freqs['R'];
        const double pyrinfreq = params.freqs['Y'];
        const double mu = params.freqs['M'];
        const double obaseFreq = params.freqs[detail.readBase];
        const char refb = detail.referenceBase;
        const char readb = detail.readBase;

// Check if either refb or readb are not A, C, T, or G
    if (refb != 'A' && refb != 'C' && refb != 'T' && refb != 'G' &&
        readb != 'A' && readb != 'C' && readb != 'T' && readb != 'G') {
        //7return log(0.999);
        cerr << "BASES: " << refb << '\t' << readb << endl;
        throw runtime_error("[TrailMix] INVALID BASE");
    }

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
            cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
            cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;

            throw runtime_error("HKY loglikemarg is invalid.");}

//cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;

        return log_lik_marg;
    }


    double calculateLogWeightedAverage(double logValueChild, double weightChild, double logValueParent, double weightParent) {

        // Log-sum-exp trick to compute log(a + b) from log(a) and log(b)
        double maxLogValue = std::max(logValueChild + std::log(weightChild), logValueParent + std::log(weightParent));
        double logSumExp = maxLogValue + std::log(std::exp(logValueChild + std::log(weightChild) - maxLogValue) + std::exp(logValueParent + std::log(weightParent) - maxLogValue));
        
        double logWeightSum = std::log(weightChild + weightParent);
        
        if (std::isinf(logWeightSum)) {
            std::cerr << "Error: Denominator is zero. Cannot compute log weighted average." << std::endl;
            return -std::numeric_limits<double>::infinity();
        }

        double logWeightedAverage = logSumExp - logWeightSum;

        return logWeightedAverage;
    }


void get_proposal_sd(double& proposal_sd, double& acceptance_rate, int current_iteration, int total_iterations, int burn_in) {
    std::random_device rd;
    std::mt19937 gen(rd());

    // Decide which range to choose
    std::uniform_int_distribution<> range_choice(0,5);
    int choice = range_choice(gen);

    // Choose uniformly within the selected range
    if (choice == 0) {
        std::uniform_real_distribution<> range(1e-3, 1e-2); // Very small range
        proposal_sd = range(gen);
    } else if (choice == 1) {
        std::uniform_real_distribution<> range(1e-2, 1e-1); // Small range
        proposal_sd = range(gen);
    } else if (choice == 2) {
        std::uniform_real_distribution<> range(1e-1, 1.0); // Medium-small range
        proposal_sd = range(gen);
    } else if (choice == 3) {
        std::uniform_real_distribution<> range(1.0, 3.0); // Medium-large range
        proposal_sd = range(gen);
    } else if (choice == 4) {
        std::uniform_real_distribution<> range(3.0, 10.0); // Large range
        proposal_sd = range(gen);
    }
    else if (choice == 5) {
        std::uniform_real_distribution<> range(10.0, 200.0); // Large range
        proposal_sd = range(gen);
    }
}

// Generate thetaVec and max_branch_lens
std::pair<std::vector<double>, std::vector<double>> generateThetaVecAndMaxBranchLens(const std::vector<PosTree>& current_positions) {
    std::vector<double> thetaVec;
    std::vector<double> max_branch_lens;

    for (auto& p : current_positions) {
        thetaVec.emplace_back(max(0.001, p.theta));
        max_branch_lens.emplace_back(p.pos->dist);
    }

    return make_pair(thetaVec, max_branch_lens);
}

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

// Initialize the MCMCiteration state
MCMCiteration initializeState(RunTreeProportionParams& params) {
    MCMCiteration state;
    state.n_components = params.sources.size();

    // Generate random numbers
    std::vector<double> random_numbers = generateRandomNumbers(state.n_components);

    // Initialize current_positions
    state.positions_tree = initializePositions(random_numbers, params);

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
std::vector<PosTree> initializePositions(const std::vector<double>& random_numbers, RunTreeProportionParams& params) {
    std::vector<PosTree> current_positions(random_numbers.size());
    int index = 0;

    for (auto& p : current_positions) {
        p.pos = params.tr->nodes[params.sources[index]];
        p.pos_branch = 0.5;//1.0 / random_numbers.size();
        p.theta = random_numbers[index];
        p.branch_place_anc = 0.5;
        p.branch_place_der = 0.5;
        index++;
    }

    return current_positions;
}


// Function to find the lowest common ancestor (LCA) of two nodes in a given tree
    spidir::Node* findLCA(const spidir::Tree* tree, spidir::Node* node1, spidir::Node* node2) {

     // Define the lambda expression inside the function
    auto isNodeInTree = [&tree](spidir::Node* nodeToCheck) -> bool {
        for (int i = 0; i < tree->nodes.size(); ++i) {
            if (tree->nodes[i] == nodeToCheck) {
                return true;
            }
        }
        return false;
    };

    if (tree == nullptr || node1 == nullptr || node2 == nullptr) {
        throw std::runtime_error("Tree or node pointers are null.");
    }
    if (!isNodeInTree(node1) || !isNodeInTree(node2)) {
        throw std::runtime_error("One or both nodes do not belong to the provided tree.");
    }
    std::unordered_map<spidir::Node*, bool> ancestors;
    // Traverse ancestors of the first node and mark them
    spidir::Node* current = node1;
    while (current != nullptr) {
        ancestors[current] = true;
        current = current->parent;
    }
    // Traverse ancestors of the second node to find the common ancestor
    current = node2;
    while (current != nullptr) {
        if (ancestors.find(current) != ancestors.end()) {
            return current;  // Return the LCA node
        }
        current = current->parent;
    }
    throw std::runtime_error("No common ancestor found.");
}

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




