#pragma once
#include "forward_declarations.h"
#include "Clade.h"
#include <map>
#include <Eigen/Dense>
//#include <nlopt.hpp>
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
    const vector<long double> generate_proposal(vector<long double> &current_vec, const double &alpha, const bool branch_pos);
    const std::vector<long double> sample_normal_euka(std::vector<long double>& x, long double alpha);
    std::vector<double> sample_normal(std::vector<double>& x, double alpha);
    void generateNumbers(double num1);

    template<typename T> const vector<T> softmax(vector<T> &log_vec);
    void run_rpvg_haplotypes(shared_ptr<Trailmix_struct>& dta);
    vector<long double> run(int iter, int burnin, double tol, const vector<long double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
    
    long double get_proposal_likelihood(const vector <long double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);

    std::vector<MCMCiteration> run_tree_proportion(RunTreeProportionParams params, std::vector<MCMCiteration> state_t_vec, const bdsg::ODGI& graph, vector<vector<string>> nodepaths, string num, int n_threads, int numPaths, int chainindex, double con);

    void updatePosition(PosTree &current_position, double move_distance, bool move_forward);
    //double moveBackward(PosTree &current_position, double move_distance_abs);
    //double moveForward(PosTree &current_position, double move_distance_abs);
    //void checkPositionErrors(const PosTree &current_position, double move_distance);
    //const double calculateEuclideanDistance(const std::vector<double> &vec1, const std::vector<double> &vec2);
    //const std::vector<double> getPatristicDistances(const spidir::Tree* tr, const spidir::Node* node);
    //double calculateDistanceToLeaf(const spidir::Tree* tr, const spidir::Node* current, const spidir::Node* leaf, double currentDistance, std::unordered_set<const spidir::Node*>& visited);
    pair<unordered_map<string, vector<vector<double>>>, double> processMCMCiterations(const std::vector<MCMCiteration> MCMCiterationsVec, int k, string num, int chain, const spidir::Tree* tr, int numofleafs);
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



    inline double computeBaseLogLike(const AlignmentInfo* read, RunTreeProportionParams &params, const int basevec, const int base, const string &pathName, const double t, const double branch_len)
    {
        // Store the repeated map lookup in a reference
        const auto detail = read->detailMap.at(pathName).at(basevec).at(base);

        const double purinfreq = params.freqs['R'];
        const double pyrinfreq = params.freqs['Y'];
        const double mu = params.freqs['M'];
        const double obaseFreq = params.freqs[detail.readBase];
        const char refb = detail.referenceBase;
        const char readb = detail.readBase;

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
                        throw runtime_error("HKY loglikemarg is invalid for trash.");
                    }

                }
            }
        } // end filling probability matrix.

        double log_lik_marg = -std::numeric_limits<double>::infinity();
        for (int bpd = 0; bpd < 4; bpd++) {
            if ("ACGT"[bpd] == readb) {
                log_lik_marg = oplusInitnatl(log_lik_marg, (log(probBaseHKY[bpd])  + log((1 -branch_len) ))); //+ log((1- (branch_len*0.1)) )
                
            } else {
                log_lik_marg = oplusInitnatl(log_lik_marg, (log(probBaseHKY[bpd])  + log((branch_len/3)))); //+ log(((branch_len*0.1)/3))
                
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

            throw runtime_error("HKY loglikemarg is invalid.");}

#ifdef DEBUGHKY
        cerr << "post HKY:" << endl;
        for (int bpd = 0; bpd < 4; bpd++) {
            cerr << setprecision(15) << bpd << "\t" << probBaseHKY[bpd] << endl;
        }
        cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
        cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;
        //if (readb != refb){throw std::runtime_error("Mismatch ");}
#endif
        log_lik_marg = log_lik_marg + detail.logLikelihood;
        if (isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
            cerr << setprecision(15) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
            cerr << setprecision(15) << "detail map log like " << detail.logLikelihood << endl;

            throw runtime_error("HKY loglikemarg is invalid.");}
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


   
    const int factorial(int n){
    if (n < 0) {
        // Invalid input, return an error value or throw an exception
        // depending on your requirements
        // For simplicity, let's return -1.0
        return -1.0L;
    } else if (n == 0) {
        // Base case: factorial of 0 is 1
        return 1.0L;
    } else {
        // Recursive case: factorial of n is n multiplied by factorial of (n-1)
        return static_cast<long double>(n) * factorial(n - 1);
    }
                             };

    int getBCIndex(const std::shared_ptr<Trailmix_struct>& dta, const std::vector<unsigned int>& combo) {
    // convert input vector to a multiset for order-agnostic comparison
    std::multiset<unsigned int> comboSet(combo.begin(), combo.end());

    for (size_t i = 0; i < dta->branch_combos.size(); i++) {
        // convert each combo in branch_combos to a multiset
        std::multiset<unsigned int> existingComboSet(dta->branch_combos[i].begin(), dta->branch_combos[i].end());

        // compare sets
        if (existingComboSet == comboSet) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

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


void get_proposal_sd(double acceptanceRate, int currentIteration, int totalIterations) {
    // Define ranges
    double minSD = 1e2;
    double maxSD = 1e2;

    // Increase the base adaptation rate
    double baseAdaptationRate = 0.1;

    // Define cutoff for adaptation: stop adapting after 50% of total iterations
    int cutoffIteration = totalIterations / 2;

    // The total number of iterations, used for adjusting the target acceptance rate
    int maxIterations = totalIterations;

    // Only adapt if we haven't reached the cutoff
    if (currentIteration < cutoffIteration) {

        // Adaptation rate decreases slowly as iterations progress, but doesn't get too small
        double adaptationRate = baseAdaptationRate / (1 + 0.001 * currentIteration);

        // Adjust the target acceptance rate dynamically based on progress in iterations
        double initialTargetAcceptanceRate = 0.50; // initial rate
        double finalTargetAcceptanceRate = 0.23; // final rate
        double targetAcceptanceRate = initialTargetAcceptanceRate -
            ((initialTargetAcceptanceRate - finalTargetAcceptanceRate) * currentIteration / maxIterations);

        // Adjust proposal SD
        if (acceptanceRate > targetAcceptanceRate) {
            proposal_sd = std::max(minSD, proposal_sd / (1.0 + adaptationRate));
        } else if (acceptanceRate < targetAcceptanceRate) {
            proposal_sd = std::min(maxSD, proposal_sd * (1.0 + adaptationRate));
        }
    }
}

long double getQuantile2(const std::vector<long double>& sortedData, double q) {
    if (sortedData.empty()) {
        throw std::runtime_error("Vector is empty");
    }
    
    if (q < 0.0 || q > 1.0) {
        throw std::invalid_argument("Quantile must be between 0 and 1");
    }

    const auto n = sortedData.size();
    const auto index = (n - 1) * q;
    const auto lowerIndex = static_cast<size_t>(floor(index));
    const auto upperIndex = static_cast<size_t>(ceil(index));

    if (lowerIndex == upperIndex) {
        return sortedData[lowerIndex];
    } else {
        const auto frac = index - lowerIndex;
        return (1.0 - frac) * sortedData[lowerIndex] + frac * sortedData[upperIndex];
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

    //std::vector<long double> LLvec_init = get_LLvec_init(state, params);
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

double getQuantile2(const std::vector<double>& sortedData, double q) {
    if (sortedData.empty()) {
        throw std::runtime_error("Vector is empty");
    }
    
    if (q < 0.0 || q > 1.0) {
        throw std::invalid_argument("Quantile must be between 0 and 1");
    }

    const auto n = sortedData.size();
    const auto index = (n - 1) * q;
    const auto lowerIndex = static_cast<size_t>(floor(index));
    const auto upperIndex = static_cast<size_t>(ceil(index));

    if (lowerIndex == upperIndex) {
        return sortedData[lowerIndex];
    } else {
        const auto frac = index - lowerIndex;
        return (1.0 - frac) * sortedData[lowerIndex] + frac * sortedData[upperIndex];
    }
}

// Function to check if a node belongs to a given tree
bool isNodeInTree(const spidir::Tree* tree, spidir::Node* node) {
    for (int i = 0; i < tree->nodes.size(); ++i) {
        if (tree->nodes[i] == node) {
            return true;
        }
    }
    return false;
}
// Function to find the lowest common ancestor (LCA) of two nodes in a given tree
    spidir::Node* findLCA(const spidir::Tree* tree, spidir::Node* node1, spidir::Node* node2) {
    if (tree == nullptr || node1 == nullptr || node2 == nullptr) {
        throw std::runtime_error("Tree or node pointers are null.");
    }
    if (!isNodeInTree(tree, node1) || !isNodeInTree(tree, node2)) {
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


// Function to calculate the distance to an ancestor node
double calculateDistanceToAncestor(spidir::Node* startNode, spidir::Node* ancestor) {
    double distance = 0.0;
    spidir::Node* current = startNode;
    while (current != nullptr && current != ancestor) {
        distance += (current->dist);
        current = current->parent;
    }
    return current == ancestor ? distance : -1.0; // Return -1 if ancestor is not found
}

const std::vector<double> getPatristicDistances(const spidir::Tree* tr, spidir::Node* node, int numofLeafs, double posonbranch) {
    if (!tr || !node) {
        throw std::runtime_error("Tree or node is null.");
    }
    std::vector<double> distances(numofLeafs, std::numeric_limits<double>::max());
    for (size_t i = 0; i < tr->nodes.size(); ++i) {
        const auto &leafNode = tr->nodes[i];
        if (!leafNode || !leafNode->isLeaf()) {
            continue;
        }
        spidir::Node* lca = findLCA(tr, node, leafNode);
        double distanceToLCAFromNode = calculateDistanceToAncestor(node, lca) - posonbranch;
        double distanceToLCAFromLeaf = calculateDistanceToAncestor(leafNode, lca);
        if (distanceToLCAFromNode >= 0.0 && distanceToLCAFromLeaf >= 0.0) {
            distances[i] = distanceToLCAFromNode + distanceToLCAFromLeaf;

        } else {
            //std::cerr << "Error calculating distance for leaf ID " << leafNode->name << std::endl;
        }
    }
    
    return distances;
}


const long double calculateEuclideanDistance(const std::vector<double> &vec1, const std::vector<double> &vec2) {
    if (vec1.empty() || vec2.empty()) {
        throw std::runtime_error("Input vectors are empty.");
    }

    size_t minSize = std::min(vec1.size(), vec2.size());
    if (minSize == 0) {
        throw std::runtime_error("One of the vectors is empty.");
    }

    long double sum = 0.0;
    for (size_t i = 0; i < minSize; ++i) {
        if (vec1[i] == std::numeric_limits<double>::max() || vec2[i] == std::numeric_limits<double>::max()) {
            //std::cerr << "Warning: Skipping comparison at index " << i << " due to invalid value." << std::endl;
            continue;
        }
        double diff = vec1[i] - vec2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}



    MCMC();
    MCMC(const MCMC & other);
    ~MCMC();
    MCMC & operator= (const MCMC & other);
};

// computes the softmax of the proposal vector.
// normalises so each proposal vector sums up to 1
template <typename T> const vector<T> MCMC::softmax(vector<T> &log_vec){
        long double K = 0.0;

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




