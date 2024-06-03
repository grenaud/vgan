#pragma once
#include "MCMC.h"
#include "MCMC_functions.h"
#include "Euka.h"
#include "libgab.h"
#include <algorithm>
#include <gzstream.h>
//#define DEBUGGENERATEVEC
//#define VERBOSE_MCMC
//#define DEBUGMCMCDig
//#define DEBUGMCMCPOS
//#define DEBUGMCMCPOS2
//#define DEBUGPAIN

#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cout << setprecision(10) << v[i] << '\t';}cout << endl << endl;

MCMC::MCMC(){

}

MCMC::~MCMC(){

}

bool MCMC::is_in_pruned_verbose(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used) {
    std::cerr << "Entering is_in_pruned with depth: " << depth << " and node: " << (p ? p->longname : "NULL") << std::endl;

    if (depth == -1) {
        //std::cerr << "Depth is -1, returning true." << std::endl;
        return true;
    }

    // Check if the current node is in the set
    if (dta->in_pruned_set.find(p->longname) != dta->in_pruned_set.end()) {
        //std::cerr << "Node " << p->longname << " is in pruned set." << std::endl;
        // Before returning true, check if the depth is already used
        bool inserted = depths_used.insert(depth).second;
        //std::cerr << "Depth " << depth << (inserted ? " is not " : " is ") << "already used." << std::endl;
        return inserted;
    } else {
        //std::cerr << "Node " << p->longname << " is not in pruned set." << std::endl;
    }

    // If depth is 0, stop the recursion
    if (depth == 0) {
        //std::cerr << "Depth is 0, stopping recursion and returning false." << std::endl;
        return false;
    }

    // Check the parent node
    if (p->parent) {
        //std::cerr << "Checking parent node of " << p->longname << std::endl;
        if (is_in_pruned(p->parent, depth - 1, dta, depths_used)) {
            //std::cerr << "Parent node or its ancestors are in pruned set, returning true." << std::endl;
            return true;
        }
    } else {
        std::cerr << "No parent node for " << p->longname << ", skipping parent check." << std::endl;
    }

    // Check children if it's not a leaf
    if (!p->isLeaf()) {
        std::cerr << "Node " << p->longname << " is not a leaf, checking children." << std::endl;
        for (int i = 0; i < p->nchildren; ++i) {
            spidir::Node* child = p->children[i];
            if (child) {
                std::cerr << "Checking child " << i << " of " << p->longname << std::endl;
                if (is_in_pruned(child, depth - 1, dta, depths_used)) {
                    std::cerr << "Child node or its descendants are in pruned set, returning true." << std::endl;
                    return true;
                }
            } else {
                std::cerr << "Child " << i << " of " << p->longname << " is NULL, skipping." << std::endl;
            }
        }
    } else {
        std::cerr << "Node " << p->longname << " is a leaf, skipping children check." << std::endl;
    }

    std::cerr << "Exiting is_in_pruned for node " << p->longname << " with false." << std::endl;
    return false;
}


bool MCMC::is_in_pruned(spidir::Node* p, const int depth, shared_ptr<Trailmix_struct> &dta, set<int> depths_used) {
    if (depth == -1) {
        return true;
    }
    // Check if the current node is in the set
    if (dta->in_pruned_set.find(p->longname) != dta->in_pruned_set.end()) {
        // Before returning true, check if the depth is already used
        if (!depths_used.insert(depth).second) {
            // Depth already used, return false
            return false;
        }
        return true;
    }
    // If depth is 0, stop the recursion
    if (depth == 0) {
        return false;
    }

    // Check the parent node
    if (p->parent && is_in_pruned(p->parent, depth - 1, dta, depths_used)) {
        return true;
    }

    // Check children if it's not a leaf
    if (!p->isLeaf()) {
        for (int i = 0; i < p->nchildren; ++i) {
            spidir::Node* child = p->children[i];
            if (child && is_in_pruned(child, depth - 1, dta, depths_used)) {
                return true;
            }
        }
    }

    return false;
}


void MCMC::generateNumbers(double num1) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    
    std::normal_distribution<double> dis(num1, 0.1);

    double lower_bound = 0.00001;
    double upper_bound = 0.999999;

    do {
        num1 = dis(gen);
    } while (num1 <= lower_bound || num1 >= upper_bound);

    
}



void MCMC::updatePosition(PosTree &current_position, double move_distance, bool move_forward) {
#ifdef DEBUGMCMCPOS
    //current_position.pos_branch = 0.5;
    std::cerr << "---- Function Called ----" << std::endl;
    
    // Initial validations
    std::cerr<< std::setprecision(14) << "Initial pos_branch: " << current_position.pos_branch << std::endl;
    std::cerr<< std::setprecision(14) << "Initial move_distance: " << move_distance << std::endl;
    std::cerr << std::setprecision(14)<< "Initial move_forward: " << std::boolalpha << move_forward << std::endl;
#endif 
    if (current_position.pos_branch < 0.0 || current_position.pos_branch > 1.0) {
        throw std::runtime_error("Error: Initial pos_branch is out of valid range.");
                                                                                }

    if (current_position.pos == nullptr) {
        throw std::runtime_error("Error: current position pointer is null.");
    }

    if (move_distance < 0.0) {
        throw std::runtime_error("Error: move distance cannot be negative.");
    }

    double move_distance_abs = std::abs(move_distance); // - current_position.pos_branch
#ifdef DEBUGMCMCPOS
    std::cerr<< std::setprecision(14) << "Initial move_distance_abs: " << move_distance_abs << std::endl;
#endif
    while (move_distance_abs > 0.0)
    {
#ifdef DEBUGMCMCPOS
    std::cerr << "Entering while loop with move_distance_abs: " << move_distance_abs << std::endl;
#endif
        if (move_forward)
        {
#ifdef DEBUGMCMCPOS
        std::cerr << "Moving Forward..." << std::endl;

#endif      
            //double updateMoveDist = abs(current_position.pos_branch - move_distance_abs);
            if (current_position.pos_branch + move_distance_abs < 1.0) {
#ifdef DEBUGMCMCPOS
                std::cerr << "Sufficient distance on the current branch to move forward." << std::endl;
#endif
                current_position.pos_branch += move_distance_abs;
#ifdef DEBUGMCMCPOS
                std::cerr << std::setprecision(14) << "New pos_branch: " << current_position.pos_branch << std::endl;
#endif
                 if (current_position.pos_branch < 0.0 || current_position.pos_branch > 1.0) {
                     throw std::runtime_error("Error: pos_branch is out of valid range after increment.");
                                                                                             }
                move_distance_abs = 0.0;
            } else {
#ifdef DEBUGMCMCPOS
                std::cerr << "Insufficient distance on the current branch. Need to move to child node." << std::endl;
#endif
                if (current_position.pos->children == nullptr) {
                    // This is a leaf node, no further movement is possible
                    //current_position.pos_branch = 0.99999; 
                    //move_distance_abs = 0.0;
                    move_forward = false;
                    continue;
                }
                double remaining_distance = move_distance_abs - (1.0 - current_position.pos_branch);
                if (remaining_distance < 0.0) {
                    remaining_distance = 0.0;
                }
#ifdef DEBUGMCMCPOS
                cerr << std::setprecision(14) << "The remaining_distance is "  << remaining_distance << endl;
#endif
                int num_children = current_position.pos->nchildren;
                int random_index = rand() % num_children;
#ifdef DEBUGMCMCPOS
                cerr << "Moving from parent " << current_position.pos->longname << endl;
#endif 
                current_position.pos = current_position.pos->children[random_index];
#ifdef DEBUGMCMCPOS
                cerr << "Moved to child " << current_position.pos->longname << endl;
#endif
                double next_branch_length = current_position.pos->dist;
#ifdef DEBUGMCMCPOS
                cerr << std::setprecision(14)<< "Branch length of the child node is " << next_branch_length << endl;
#endif
                if (next_branch_length < 0.0) {
                    throw std::runtime_error("Error: next branch length cannot be negative.");
                }

                if (remaining_distance > 1.0) {
                    current_position.pos_branch = 1.0;
                    move_distance_abs = remaining_distance - 1.0;
                    
#ifdef DEBUGMCMCPOS
                    cerr << std::setprecision(14)<< "remaining distance exceeded the branch length of " << 1 << endl;
                    cerr << std::setprecision(14)<< "The new move_distance_abs is " << move_distance_abs << endl;
#endif

                } else {
#ifdef DEBUGMCMCPOS
                    cerr << "The remaining_distance is not longer than the next branch length. Calculating next position on the branch. " << endl;
                    cerr << std::setprecision(14)<< "remaining_distance " << remaining_distance << " next_branch_length " << next_branch_length << endl; 
#endif              //current_position.pos_branch = 0.0;
                    current_position.pos_branch = remaining_distance;
                    move_distance_abs = 0.0;
#ifdef DEBUGMCMCPOS
                    cerr << std::setprecision(14)<< "Forward move new branch position 11 " << current_position.pos_branch << endl;
#endif
                }
            }
        } else
        {
#ifdef DEBUGMCMCPOS
             std::cerr << "Moving Backward..." << std::endl;
#endif
            //double current_position_reversed = abs(current_position.pos_branch - move_distance_abs);
#ifdef DEBUGMCMCPOS
            cerr << std::setprecision(14) << "current_position_reversed is " << current_position.pos_branch << endl;
#endif
            if (current_position.pos_branch - move_distance_abs > 0.0) {

                current_position.pos_branch = current_position.pos_branch - move_distance_abs;
#ifdef DEBUGMCMCPOS
                cerr << std::setprecision(14) << "Backwards move on the branch " << current_position.pos_branch << endl;
#endif
                move_distance_abs = 0.0;
            } else
            {

                std::vector<spidir::Node*> possible_nodes;

                if (current_position.pos->parent == nullptr) {
#ifdef DEBUGMCMCPOS
                    cerr << "Root was hit"<< current_position.pos->longname << endl;
                    cerr << "Root children " << current_position.pos->children[0]->longname << " " << current_position.pos->children[1]->longname << " " << current_position.pos->children[2]->longname << endl;
#endif
                    // This is a root node, no further movement is possible
                    move_forward = true; // Change direction if a leaf node is reached
                    int num_children = current_position.pos->nchildren;
                    int random_index = rand() % num_children;

                    current_position.pos = current_position.pos->children[random_index];

                    double next_branch_length = current_position.pos->dist; // come back here: should probably be just 1 because we assume that all the branch lengths are the same in the MCMC

                    if (next_branch_length < 0.0) {
                        throw std::runtime_error("Error: next branch length cannot be negative.");
                    }
                    continue;  
                     
                }else{
                    possible_nodes.push_back(current_position.pos->parent);

                }
                // Add the children except the curent one
                if (current_position.pos->children != nullptr) {
                    for (int i = 0; i < current_position.pos->parent->nchildren; ++i) {
                        if (current_position.pos->parent->children[i] != current_position.pos) {
                            possible_nodes.push_back(current_position.pos->parent->children[i]);
#ifdef DEBUGMCMCPOS
                            cerr << "Current branch " << current_position.pos->longname << endl;
                            cerr << "Possible child node " << current_position.pos->parent->children[i]->longname << endl;
#endif
                        }
                    }
                }
                // If there are no possible nodes to choose from (which is very unlikely but just in case)
                if (possible_nodes.empty()) {
                    throw std::runtime_error("Error: No valid nodes to move to.");
                }

                spidir::Node* chosen_node = possible_nodes[rand() % possible_nodes.size()];

                // Continuous backwards move to the parent node
                if (chosen_node == current_position.pos->parent)
                { //&& chosen_node->parent != nullptr
#ifdef DEBUGMCMCPOS
                    std::cerr << "Current position reversed: " << current_position.pos_branch << std::endl;
                    std::cerr << "Move distance abs: " << move_distance_abs << std::endl;
                    std::cerr << "Calculated remaining distance before check: " << (move_distance_abs - current_position.pos_branch) << std::endl;
#endif
                    double remaining_distance = move_distance_abs - current_position.pos_branch;
                    if (remaining_distance < 0.0) {
                        remaining_distance = 0.0;
#ifdef DEBUGMCMCPOS
                        cerr << "This should never happen" << endl;
#endif
                    }

                    //setting the parent node as the new position 
                    current_position.pos = current_position.pos->parent;
                    double parent_branch_length = current_position.pos->dist;

                    if (parent_branch_length < 0.0) {
                        throw std::runtime_error("Error: parent branch length cannot be negative.");
                    }

                    if (remaining_distance > 1.0) {
                        current_position.pos_branch = 0.0;
                        move_distance_abs = remaining_distance - 1.0;
                        continue;
                        
                    } else {
#ifdef DEBUGMCMCPOS
                        std::cerr << "Remaining distance: " << remaining_distance << std::endl;
                        std::cerr << "Parent branch length: " << parent_branch_length << std::endl;
#endif
                        double new_pos_branch = 1.0 - remaining_distance;
#ifdef DEBUGMCMCPOS
                        std::cerr << "New pos branch: " << new_pos_branch << std::endl;
#endif

                        if (new_pos_branch <= 0.0 || new_pos_branch >= 1.0) {
                            throw std::runtime_error("Error: new position branch is not in the valid range.");
                        }

                        current_position.pos_branch = new_pos_branch;
                        move_distance_abs = 0.0;
                    }
                }else{// Chosen node is the sibling
                    move_forward = true;
                    current_position.pos = chosen_node;
                    double remaining_distance = move_distance_abs - current_position.pos_branch;
                    current_position.pos_branch = 0.0;

#ifdef DEBUGMCMCPOS
                    cerr << "Chose the sibling " << current_position.pos->longname << endl;
                    cerr << "Moving Forward again" << endl;
                    cerr << "The remaining_distance is " << remaining_distance << endl;
#endif

                    if (current_position.pos_branch + remaining_distance < 1.0) {
#ifdef DEBUGMCMCPOS
                        std::cerr << "Sufficient distance on the current branch to move forward." << std::endl;
#endif
                        current_position.pos_branch = remaining_distance;
#ifdef DEBUGMCMCPOS
                        std::cerr << std::setprecision(14) << "New pos_branch: " << current_position.pos_branch << std::endl;
#endif

                        if (current_position.pos_branch < 0.0 || current_position.pos_branch > 1.0) {
                            throw std::runtime_error("Error: pos_branch is out of valid range after increment.");
                        }
                        move_distance_abs = 0.0;
                    } else { // would exceed the leaf 
#ifdef DEBUGMCMCPOS
                        std::cerr << "Insufficient distance on the current branch. Need to move to child node." << std::endl;
#endif
                        if (current_position.pos->children == nullptr) {
                            // This is a leaf node, no further movement is possible
                            //current_position.pos_branch = 0.99999; 
                            //move_distance_abs = 0.0;
                            move_forward = false;
                            continue;
                        }
                        double remaining_distance = move_distance_abs - (1.0 - current_position.pos_branch);
                        if (remaining_distance < 0.0) {
                            remaining_distance = 0.0;
                        }
#ifdef DEBUGMCMCPOS
                        cerr << std::setprecision(14) << "The remaining_distance is "  << remaining_distance << endl;
#endif
                        int num_children = current_position.pos->nchildren;
                        int random_index = rand() % num_children;
#ifdef DEBUGMCMCPOS
                        cerr << "Moving from parent " << current_position.pos->longname << endl;
#endif 
                        current_position.pos = current_position.pos->children[random_index];
#ifdef DEBUGMCMCPOS
                        cerr << "Moved to child " << current_position.pos->longname << endl;
#endif
                        double next_branch_length = current_position.pos->dist;
#ifdef DEBUGMCMCPOS
                        cerr << std::setprecision(14)<< "Branch length of the child node is " << next_branch_length << endl;
#endif
                        if (next_branch_length < 0.0) {
                            throw std::runtime_error("Error: next branch length cannot be negative.");
                        }

                        if (remaining_distance > 1.0) {
                            current_position.pos_branch = 1.0;
                            move_distance_abs = remaining_distance - 1.0;
                            
#ifdef DEBUGMCMCPOS
                            cerr << std::setprecision(14)<< "remaining distance exceeded the branch length of " << 1 << endl;
                            cerr << std::setprecision(14)<< "The new move_distance_abs is " << move_distance_abs << endl;
#endif
                        
                        }
                    }
                
                }
            }
        }
#ifdef DEBUGMCMCPOS
        std::cerr << std::setprecision(14) << "Updated move_distance_abs: " << move_distance_abs << std::endl;
        std::cerr << std::setprecision(14)<< "Updated pos_branch: " << current_position.pos_branch << std::endl;
#endif
        
    }
    if (current_position.pos_branch < 0.0 || current_position.pos_branch > 1.0) {
            throw std::runtime_error("Error: pos_branch is out of valid range after movement.");
        }
    return;
}


const std::vector<double> MCMC::sample_normal_euka(std::vector<double>& x, double alpha) {
    std::vector<double> result;
    std::vector<std::normal_distribution<double>> dists;
    static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count()); // initialize with a random seed
    for (size_t i = 0; i < x.size(); ++i) {
        double stdev = alpha; // define alpha as a fraction of max_branch_lens
        std::normal_distribution<double> dist(x[i], stdev);
        double sample = dist(generator);
        //cerr << "SAMPLE: " << sample << endl;
        result.emplace_back(std::max(0.0000001, std::min(sample, alpha)));
    }
    return result;
}


// samples the proportion of the sources
std::vector<double> MCMC::sample_normal(std::vector<double>& x, double alpha) {
    constexpr double epsilon = 1e-9;
    std::vector<double> result;
    std::vector<std::normal_distribution<double>> dists;
    static std::mt19937 generator(std::random_device{}());

    if (x.empty()) {
        throw std::invalid_argument("vector can't be empty");
    }

    double sum = 0.0L; // Variable to store the sum of sampled values

    for (size_t i = 0; i < x.size(); ++i) {
        std::normal_distribution<double> dist(x[i], 0.1);
        double sample;

        do {
            sample = dist(generator);
        } while (sample < 0.0L || sample > 1.0L);

        result.emplace_back(sample);
        sum += sample;
    }

    // Normalize the sampled values to ensure they sum up to 1
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] /= sum;
    }

    return result;
}

std::vector<MCMCiteration> MCMC::run_tree_proportion_TM(RunTreeProportionParams &params, std::vector<MCMCiteration> &state_t_vec, const bdsg::ODGI& graph, \
                               const vector<vector<string>> &nodepaths, string num, shared_ptr<Trailmix_struct> &dta, bool running_trailmix, int chain) {


    const unsigned int n_sources = params.sources.size();
    cerr << "number of sources " << n_sources << endl;
    if (n_sources > 10){throw runtime_error("We cannot handle this many sources");}

    unsigned int total_proposals = 0;
    unsigned int n_accept = 0;
    double acceptance_rate = 0.5;
    MCMCiteration state_t_1;
    double likelihood_t_1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<PosTree> current_positions(n_sources);
    std::vector<double> random_numbers(n_sources);
    MCMCiteration state_t = initializeState(params, dta->seed);
    auto initial_state = state_t;
    state_t_vec.emplace_back(initial_state);

    ofstream mcmcout(num+"Result"+to_string(n_sources)+"_chain"+to_string(chain)+".mcmc");
    std::ios_base::sync_with_stdio(false);
    mcmcout.tie(nullptr);

    for (int sou = 1; sou< n_sources + 1; ++sou){
        mcmcout << "Source_" << sou << '\t' << "Node_" << sou << '\t' << "Loglikelihood" << '\t' << "proportion_" << sou << '\t' << "branch_position_derived_" << sou;
        if (sou != n_sources){mcmcout << '\t';}
    }
    mcmcout << endl;

    ofstream mcmcdetail(num+"Trace"+to_string(n_sources)+".detail.mcmc");
    for (int sou = 1; sou< n_sources + 1; ++sou){
        mcmcdetail << "Source_" << sou << '\t' << "Node_" << sou << '\t' << "Likelihood" << '\t' << "proportion_" << sou << '\t' << \
                      "branch_position_derived_" << sou << '\t' << "Move";
        if (sou != n_sources){mcmcdetail << '\t';}
    }
    mcmcdetail << endl;

    for (unsigned int iteration = 0; iteration <= params.maxIter; iteration++) {

        if (iteration % 1000 == 0){
            cerr << "ITERATION: " << iteration << endl;
                                  }

        state_t_1 = state_t;

 get_proposal_sd(proposal_sd, acceptance_rate, iteration, params.maxIter, params.burn);

 srand(time(0)); // Initialize random seed

// Choose a random component index
unsigned int componentIndex = std::rand() % state_t_1.n_components;

float randomFloat = 0.01 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/(0.99 - 0.01)));

            if (state_t_1.positions_tree[componentIndex].pos = dta->tree->root){
                int num_children = state_t_1.positions_tree[componentIndex].pos->nchildren;
                int random_index = rand() % num_children;
                state_t_1.positions_tree[componentIndex].pos = state_t_1.positions_tree[componentIndex].pos->children[random_index];
                                                                               }

            // Determine the mean based on the iteration
            double proposal_mean = state_t_1.positions_tree[componentIndex].pos_branch;
            std::normal_distribution<double> distribution_bl(proposal_mean, proposal_sd);
            double proposed_position = distribution_bl(gen);
            if (proposed_position < 0.0) {
                if (state_t_1.positions_tree[componentIndex].pos->parent == nullptr) {

                    const auto rootChildren = state_t_1.positions_tree[componentIndex].pos->children;
                    std::uniform_int_distribution<int> distribution(0.0, state_t_1.positions_tree[componentIndex].pos->nchildren);
                    const unsigned int selectedChildIndex = distribution(gen);
                    state_t_1.positions_tree[componentIndex].pos = rootChildren[selectedChildIndex];
                    double proposed_position_abs = abs(proposed_position);
                    state_t_1.positions_tree[componentIndex].pos_branch = proposed_position_abs;
                } else {
                        updatePosition(state_t_1.positions_tree[componentIndex], -proposed_position, false);
                }
            } else {
                if (proposed_position > state_t_1.positions_tree[componentIndex].pos_branch) {
                    if (state_t_1.positions_tree[componentIndex].pos->isLeaf() && state_t_1.positions_tree[componentIndex].pos_branch == 0.9999999) {
                        continue;
                    } else {
                        updatePosition(state_t_1.positions_tree[componentIndex], proposed_position, true);
                    }
                } else {
                    state_t_1.positions_tree[componentIndex].pos_branch = proposed_position;
                       }
            }

        std::vector<double> tmp_theta;
        for (auto& p : state_t_1.positions_tree) {
            tmp_theta.emplace_back(p.theta);
        }

        tmp_theta = sample_normal(tmp_theta, 0.1);
        for (int idx = 0; idx < current_positions.size(); ++idx) {
            state_t_1.positions_tree[idx].theta = tmp_theta[idx];
        }


        state_t_1.proportions = tmp_theta;
        vector<string> pathNames;
        vector<string> parentpathNames;

        if (state_t_1.positions_tree.empty()){throw runtime_error("TREE POSITIONS ARE EMPTY");}

         for (auto & p : state_t_1.positions_tree){
                pathNames.emplace_back(p.pos->longname);

            if (p.pos->parent != nullptr) {
                parentpathNames.emplace_back(p.pos->parent->longname);

            }else{
                //cerr << "placing back parent: " << p.pos->longname << endl;
                 if (p.pos->longname == ""){
                cerr << "Node number: " << p.pos->name << endl;
                throw runtime_error("Longname is empty 2");
                                           }
                parentpathNames.emplace_back(p.pos->longname);
            }
        }

        unsigned int read_counter=0;
        double logLike = 0.0;

        #pragma omp parallel for num_threads(dta->n_threads) reduction(+:logLike)
        for (auto read : *(params.align))
        {
            read_counter++;
            bool safe=true;

            double readLogLike = 0.0;
            double readLogLikeP = 0.0;

////////////////////////////////////////////////////////////// SINGLE SOURCE //////////////////////////////////////////////////////////

            if (pathNames.size() == 1)
            {
                const double t = state_t_1.positions_tree[0].pos->dist;
                const double t1 = state_t_1.positions_tree[0].pos_branch * t;
                const double t2 = t - t1;

                bool loopentered=false;

if (read->detailMap.find(pathNames[0]) == read->detailMap.end()) {
    cerr << "Missing key in detailMap for path: " << pathNames[0] << endl;
    throw runtime_error("Missing key");
} else {
    for (unsigned int basevec = 0; basevec < read->detailMap[pathNames[0]].size(); ++basevec) {
        for (int base = 0; base < read->detailMap[pathNames[0]][basevec].size(); ++base) {
            loopentered = true;
            if (read->detailMap[pathNames[0]][basevec][base].pathSupport) {
                 readLogLike += computeBaseLogLike(read, params, basevec, base, pathNames[0], t2, t);
            } else {
                 readLogLike += read->detailMap[pathNames[0]][basevec][base].logLikelihood;
            }
        }
    }
}

if (read->detailMap.find(parentpathNames[0]) == read->detailMap.end()) {
    cerr << "Missing key in detailMap for parent path: " << parentpathNames[0] << endl;
    throw runtime_error("Missing key");
} else {
    for (unsigned int basevec = 0; basevec < read->detailMap[parentpathNames[0]].size(); ++basevec) {
        for (int base = 0; base < read->detailMap[parentpathNames[0]][basevec].size(); ++base) {
            if (read->detailMap[parentpathNames[0]][basevec][base].pathSupport) {
                readLogLikeP += computeBaseLogLike(read, params, basevec, base, parentpathNames[0], t1, t);
            } else {
                 readLogLikeP += read->detailMap[parentpathNames[0]][basevec][base].logLikelihood;
            }
        }
    }
}

// Apply bounds to readLogLikeP
const double minLogValue = log(0.0000001); // Lower bound to prevent extremely small values
const double maxLogValue = log(0.9999999); // Upper bound

if (readLogLikeP > maxLogValue) {
    readLogLikeP = maxLogValue;
} else if (readLogLikeP < minLogValue) {
    readLogLikeP = minLogValue;
}


        // Compute intermediate values as before
        double interc2 = log(state_t_1.positions_tree[0].pos_branch) + readLogLike;
        double interp2 = log(1 - state_t_1.positions_tree[0].pos_branch) + readLogLikeP;
        double max_val = std::max(interc2, interp2);
        double inter = max_val + log(exp(interc2 - max_val) + exp(interp2 - max_val));

        logLike += inter;

//if (!pruned){logLike = -std::numeric_limits<double>::max();}

            }

/////////////////////////////////////////////////////////////// END SINGLE SOURCE //////////////////////////////////////////////////////

            else
            {
////////////////////////////////////////////////////////////// BEGIN MULTI SOURCE ///////////////////////////////////////////////////////

                double inter = -std::numeric_limits<double>::infinity();

//////////////////////////////////////////////// CONT MODE ///////////////////////////////////////

if (dta->cont_mode && pathNames.size() == 2) {
//throw runtime_error("NOT DOING CONT MODE YET");
    // Handling the specific case when there are exactly two sources and cont_mode is true
    for (unsigned int y = 0; y < pathNames.size(); ++y) {
        double readLogLike = 0.0;
        double readLogLikeP = 0.0;
        double t = state_t_1.positions_tree[y].pos->dist;
        // Handle special cases for t
        if (t == 0.0) { t = 0.00000001; }
        double t1 = state_t_1.positions_tree[y].pos_branch * t;
        double t2 = t - t1;
        if (state_t_1.positions_tree[y].pos_branch == 1.0) {
            state_t_1.positions_tree[y].pos_branch = 0.999999;
        }

        for (unsigned int basevec = 0; basevec < read->detailMap[pathNames[y]].size(); ++basevec) {
            for (unsigned int base = 0; base < read->detailMap[pathNames[y]][basevec].size(); ++base) {

                //cerr << "paths: " << pathNames[y] << '\t' << parentpathNames[y] << endl;

                // Use .logLikelihood or .logLikelihoodNoDamage based on ancient/modern source
                double logLikelihoodValue = dta->source_assignments[y] ?
                                            read->detailMap[pathNames[y]][basevec][base].logLikelihood :
                                            read->detailMap[pathNames[y]][basevec][base].logLikelihoodNoDamage;

                if (read->detailMap[pathNames[y]][basevec][base].pathSupport) {
                    readLogLike += computeBaseLogLike(read, params, basevec, base, pathNames[y], t2, t);
                } else {
                    readLogLike += logLikelihoodValue;
                }

                // Similar logic for parent paths and readLogLikeP
        double parentLogLikelihoodValue = dta->source_assignments[y] ?
                                          read->detailMap[parentpathNames[y]][basevec][base].logLikelihood :
                                          read->detailMap[parentpathNames[y]][basevec][base].logLikelihoodNoDamage;

        if (read->detailMap[parentpathNames[y]][basevec][base].pathSupport) {
            readLogLikeP += computeBaseLogLike(read, params, basevec, base, parentpathNames[y], t1, t);
                                                                            }
         else {
            readLogLikeP += parentLogLikelihoodValue;
               }

        }
      }

// Apply bounds to readLogLikeP
const double minLogValue = log(0.0000001); // Lower bound to prevent extremely small values
const double maxLogValue = log(0.9999999); // Upper bound

if (readLogLikeP > maxLogValue) {
    readLogLikeP = maxLogValue;
} else if (readLogLikeP < minLogValue) {
    readLogLikeP = minLogValue;
}

        // Compute intermediate values as before
        double interc2 = log(state_t_1.positions_tree[y].pos_branch) + readLogLike;
        double interp2 = log(1 - state_t_1.positions_tree[y].pos_branch) + readLogLikeP;
        double inter2 = oplusnatl(interc2, interp2);
        auto intersum = inter2 + log(state_t_1.proportions[y]);
        inter = oplusInitnatl(inter, intersum);
    }
    logLike += inter;
}
                 //////////////////////////////////////// END CONT MODE /////////////////////////////////////////////////////////
else {
                for (unsigned int y = 0; y < pathNames.size(); ++y)
                {

                    double readLogLike = 0.0;
                    double readLogLikeP = 0.0;

                    double t = state_t_1.positions_tree[y].pos->dist;
                    if (t == 0.0){ //this will only happen at the root node. This node does not have a branch length because we have an unrooted tree.
                        t = 0.00000001;
                    }
                    const double t1 = state_t_1.positions_tree[y].pos_branch * t;
                    const double t2 = t - t1;
                    if (state_t_1.positions_tree[y].pos_branch == 1.0){
                        state_t_1.positions_tree[y].pos_branch == 0.999999;
                    }

                    // RIGHT HERE

for (unsigned int basevec = 0; basevec < read->detailMap[pathNames[y]].size(); ++basevec) {
    // Pre-check to ensure the key exists before accessing its value
    auto itPath = read->detailMap.find(pathNames[y]);
    if (itPath == read->detailMap.end()) {
        cerr << "Missing key in detailMap: " << pathNames[y] << endl;
        break; // Exit the loop if the key is missing
    }

    for (unsigned int base = 0; base < itPath->second[basevec].size(); ++base) {
        if (itPath->second[basevec][base].pathSupport) {
            readLogLike += computeBaseLogLike(read, params, basevec, base, pathNames[y], t2, t);
        } else {
            readLogLike += itPath->second[basevec][base].logLikelihood;
        }
    }
}


auto itParentPath = read->detailMap.find(parentpathNames[y]);

if (itParentPath == read->detailMap.end()) {
    cerr << "Missing key in detailMap: " << parentpathNames[y] << endl;
} else {
    for (unsigned int basevec = 0; basevec < itParentPath->second.size(); ++basevec) {
        for (unsigned int base = 0; base < itParentPath->second[basevec].size(); ++base) {
            if (itParentPath->second[basevec][base].pathSupport) {
                readLogLike += computeBaseLogLike(read, params, basevec, base, parentpathNames[y], t1, t);
            } else {
                readLogLikeP += itParentPath->second[basevec][base].logLikelihood;
            }

         }
      }
}

// Apply bounds to readLogLikeP
const double minLogValue = log(0.0000001); // Lower bound to prevent extremely small values
const double maxLogValue = log(0.9999999); // Upper bound

if (readLogLikeP > maxLogValue) {
    readLogLikeP = maxLogValue;
} else if (readLogLikeP < minLogValue) {
    readLogLikeP = minLogValue;
}

        // Compute intermediate values as before
        double interc2 = log(state_t_1.positions_tree[y].pos_branch) + readLogLike;
        double interp2 = log(1 - state_t_1.positions_tree[y].pos_branch) + readLogLikeP;
        double inter2 = oplusnatl(interc2, interp2);
        auto intersum = inter2 + log(state_t_1.proportions[y]);
        inter = oplusInitnatl(inter, intersum);

//if (!pruned){
//    inter += log(0.0000000000001);logLike += log(0.0000000000001);
//            }

      }
logLike += inter;

   }

}
//////////////////////////////////////////////////////////  END MULTI SOURCE ///////////////////////////////////////////////////////////////////////

            }
//////////////////////////////////////////////////////////  END LOOP OVER READS ///////////////////////////////////////////////////////////////////////

#ifdef VERBOSE_MCMC
         //if (logLike >= 0.0){throw runtime_error("SINGLE SOURCE SHOULDNT HAPPEN");}
         if (iteration % 100 == 1){
            cerr << endl;
            cerr << setprecision(16) << "proposal SD: " << proposal_sd << "  acceptance rate: " << acceptance_rate << endl;
            cerr << "proposed log lik: " << setprecision(16) << logLike << "      current:  " << state_t_1.logLike << endl;
            for (auto & p : state_t_1.positions_tree){
                cerr << "BEST SOURCE: " << dta->originalPathNames[p.pos->longname] << " proportion: " << p.theta << "  pos on branch: " << p.pos_branch << endl;
                                                   }
                                  }
#endif

        likelihood_t_1 = logLike;
        state_t_1.logLike = likelihood_t_1;
        double acceptance_prob = (state_t_1.logLike - state_t.logLike > 0) ? 1.0 : exp((state_t_1.logLike - state_t.logLike));
        //acceptance_prob = max(acceptance_prob, 0.0001);
        //acceptance_prob = min(acceptance_prob, 0.9999);
        //cerr << "AP: " << acceptance_prob << endl;
        //if (    (state_t_1.logLike - state_t.logLike > 0) && ( (state_t_1.logLike - state_t.logLike) > 200.0)     ){acceptance_prob = 0.0000000001;}
        const double u = dis(gen);

        if (u <= acceptance_prob || iteration == 0) {

            if (iteration > params.burn){
                int source_counter=0;
                for (auto & p : state_t_1.positions_tree){

                    mcmcout << setprecision(14) << dta->originalPathNames[p.pos->longname] << '\t' << p.pos->name  << '\t' << state_t.logLike << '\t' << p.theta << '\t' << p.pos_branch;
                    mcmcdetail << std::setprecision(14) << dta->originalPathNames[p.pos->longname]  << '\t' << p.pos->name  << '\t' << state_t_1.logLike << '\t' << p.theta << \
                                                   '\t' <<  p.pos_branch << '\t' << "accepted";

                    source_counter++;

                    if (source_counter < state_t.positions_tree.size()){mcmcdetail << '\t'; mcmcout << '\t';}

                }
                mcmcout << endl;
                mcmcdetail << endl;
            }

            n_accept++;
            total_proposals++;
            state_t = state_t_1;
            if (iteration >= params.burn) {
            state_t_vec.emplace_back(state_t_1);
                                          }
        } else {

            if (iteration > params.burn) {
    int total_positions = state_t.positions_tree.size();

    // Output for current state's positions
    for (auto &p : state_t.positions_tree) {
        mcmcout << std::setprecision(14) << dta->originalPathNames[p.pos->longname] << '\t' << p.pos->name << '\t'
                << state_t.logLike << '\t' << p.theta << '\t' << p.pos_branch;

        if (&p != &state_t.positions_tree.back()) mcmcout << '\t'; // No tab after last element
    }
    mcmcout << std::endl;

    // Output for previous state's positions
    for (auto &p : state_t_1.positions_tree) {
        mcmcdetail << std::setprecision(14) << dta->originalPathNames[p.pos->longname] << '\t'
                   << p.pos->name << '\t' << state_t.logLike << '\t'
                   << p.theta << '\t' << p.pos_branch << '\t' << "rejected";

        if (&p != &state_t_1.positions_tree.back()) mcmcdetail << '\t'; // No tab after last element
    }
    mcmcdetail << std::endl;
}

            total_proposals++;
        }

        acceptance_rate = static_cast<double>(n_accept) / total_proposals;
    }

//for (auto & p : state_t.positions_tree){
//    if (p.pos->longname != ""){
//        dta->paths_to_surject.emplace_back(p.pos->longname);
//                              }
//}

if (state_t_vec.empty()){
cerr << "State_t_vec is empty!" << endl;
//throw runtime_error("State_t_vec is empty!");
}

return state_t_vec;

}


std::vector<MCMCiteration> MCMC::run_tree_proportion(RunTreeProportionParams params, std::vector<MCMCiteration> state_t_vec, const bdsg::ODGI& graph, \
                                                     vector<vector<string>> nodepaths, string num, int n_threads, int numPaths, int chainindex, double con) {

    const unsigned int n_sources = params.sources.size();
    cerr << "number of sources " << n_sources << endl;
    //cerr << "con " << con << endl;
    
    int total_proposals = 0;
    int n_accept = 0;
    double acceptance_rate = 0.5;
    MCMCiteration state_t_1;
    double likelihood_t_1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<PosTree> current_positions(n_sources);
    std::vector<double> random_numbers(n_sources);
    MCMCiteration state_t = initializeState(params);
    int accept_count = 0;

    double proposal_sd; 
    double initSD;

    if(numPaths <= 30.0){
        initSD = 3.0;
    }else{
        initSD = numPaths * (3.0/30.0);
    }

    double step = (initSD - 0.1)/ (params.burn - 1);
    double step2 = (0.1 - 1e-5)/((params.maxIter - params.burn)- 1);

    ogzstream mcmcout;
    mcmcout.open((num+"Result"+to_string(n_sources)+to_string(chainindex)+".mcmc").c_str(), ios::out);
    for (int sou = 1; sou< n_sources + 1; ++sou){
        mcmcout << "Source_" << sou << '\t' << "Log-likelihood" << '\t' << "proportion" << '\t' << "branch_position_derived" << '\t'; 
    }
    mcmcout << endl;

    ogzstream mcmcdetail;

    mcmcdetail.open((num+"Trace"+to_string(n_sources)+to_string(chainindex)+".detail.mcmc").c_str(), ios::out);
    for (int sou = 1; sou< n_sources + 1; ++sou){
        mcmcdetail << "Source_" << sou << '\t' << "Log-likelihood" << '\t' << "proportion_" << sou << '\t' << "branch_position_derived_" << sou << '\t' << "Move" << '\t'; 
    }
    mcmcdetail << endl;

    for (unsigned int iteration = 0; iteration <= params.maxIter; iteration++)
    {
        // Ensure that 'params.burn' is less than 'params.maxIter'
        if (params.burn >= params.maxIter) {
            throw runtime_error("Number of brun in iteration exceedes the number of total iterations. Exiting. ");
        }

        double step = (initSD - 0.1) / std::max(static_cast<unsigned int>(1), params.burn - 1);

        double step2 = (0.1 - 1e-5) / std::max(static_cast<unsigned int>(1), (params.maxIter - params.burn) - 1);


        state_t_1 = state_t;

        if (iteration < params.burn){

             proposal_sd = std::max(1e-5, initSD - iteration * step);
            //cerr << proposal_sd << endl;
        }else{
            if (iteration % 100000 == 0){
                proposal_sd = 1;
            }else{
                proposal_sd = std::max(1e-5, 0.1 - (iteration - params.burn) * step2);
            }
            
            //proposal_sd = 0.1;
        }

        //cerr <<std::setprecision(14)<< proposal_sd << endl;

        if (iteration != 0){
            
        


            for (unsigned int i = 0; i < state_t_1.n_components; i++) {
                
                std::normal_distribution<double> distribution_bl(0, proposal_sd); //state_t_1.positions_tree[i].pos_branch,
                double proposed_position = distribution_bl(gen);
#ifdef DEBUGMCMCPOS2
                cerr << std::setprecision(14) << "sd " << proposal_sd << endl;
                cerr << std::setprecision(14)<<"proposed position " << proposed_position << endl;
                cerr << std::setprecision(14)<<"state t position " << state_t.positions_tree[i].pos_branch << endl;
                cerr << std::setprecision(14)<<"state_t_1 position " << state_t_1.positions_tree[i].pos_branch << endl;
#endif
                if (proposed_position < 0.0) {
                    //if(proposed_position > 0.0){
                        //double newP = proposed_position - state_t_1.positions_tree[i].pos_branch;
                        //proposed_position = newP;
#ifdef DEBUGMCMCPOS2
                        //cerr<< std::setprecision(14) << "Moving Backwards on the same branch with positive move " << proposed_position << endl;
#endif 
                    //}
    //                 if (state_t_1.positions_tree[i].pos->parent == nullptr && state_t_1.positions_tree[i].pos_branch <= 0.0000001) {

    //                     const auto rootChildren = state_t_1.positions_tree[i].pos->children;
    //                     std::uniform_int_distribution<int> distribution(0.0, state_t_1.positions_tree[i].pos->nchildren);
    //                     const unsigned int selectedChildIndex = distribution(gen);
    //                     state_t_1.positions_tree[i].pos = rootChildren[selectedChildIndex];
    //                     double proposed_position_abs = abs(proposed_position);
    //                     //updatePosition(state_t_1.positions_tree[i], proposed_position, true);
    //                     state_t_1.positions_tree[i].pos_branch = proposed_position_abs;
    // #ifdef DEBUGMCMCPOS2
    //                     cerr<< std::setprecision(14) << "The position is at root nullptr; new pos branch is " << state_t_1.positions_tree[i].pos_branch << endl;
    // #endif                    
    //                 } else {
#ifdef DEBUGMCMCPOS2
                        cerr << "Backwards move!" << endl;
#endif
                        updatePosition(state_t_1.positions_tree[i], -proposed_position, false);
#ifdef DEBUGMCMCPOS2
                        cerr<< std::setprecision(14) << "new position on the tree (function1) " << state_t_1.positions_tree[i].pos_branch << endl;
#endif
                }else{
                        /*    
                    }
                } else {
                    if (proposed_position > state_t_1.positions_tree[i].pos_branch) {
                        if (state_t_1.positions_tree[i].pos->isLeaf() && state_t_1.positions_tree[i].pos_branch >= 0.99999999999) {
    #ifdef DEBUGMCMCPOS2
                        cerr << "Hit the leaf" << endl;
                        cerr << "branch position after hitting the leaf " << state_t_1.positions_tree[i].pos_branch << endl;
    #endif                        
                            //state_t_1.positions_tree[i].pos_branch = 0.999999;
                            continue;
                        } else {*/
    #ifdef DEBUGMCMCPOS2
                        cerr << "Forward move!" << endl;
    #endif
                            updatePosition(state_t_1.positions_tree[i], proposed_position, true);
    #ifdef DEBUGMCMCPOS2
                        cerr<< std::setprecision(14) << "new position on the tree (function2) " << state_t_1.positions_tree[i].pos_branch << endl;
    #endif

                        }
    //                 } else {
    //                     state_t_1.positions_tree[i].pos_branch = proposed_position;
    //                     //updatePosition(state_t_1.positions_tree[i], proposed_position, false);
    // #ifdef DEBUGMCMCPOS2
    //                     cerr<< std::setprecision(14) << "new position on the tree " << state_t_1.positions_tree[i].pos_branch << endl;
    // #endif
    //                 }
    //             }
                if (state_t_1.positions_tree[i].pos_branch == 1.0){
                    state_t_1.positions_tree[i].pos_branch == 0.99999999;
                    //cerr << "The tree position is shit" << endl;
                }
                //cerr <<"child in loop "<< state_t_1.positions_tree[i].pos->longname << endl;
                //cerr <<"parent in loop "<< state_t_1.positions_tree[i].pos->parent->longname << endl;
            }
        }
        
        std::vector<double> tmp_theta;
        for (auto& p : state_t_1.positions_tree) {
            tmp_theta.emplace_back(p.theta);
            
        }


        
     
        tmp_theta = sample_normal(tmp_theta, 0.01);
        for (int idx = 0; idx < current_positions.size(); ++idx) {
            state_t_1.positions_tree[idx].theta = tmp_theta[idx];
            //generateNumbers(state_t_1.positions_tree[idx].pos_branch);
        }



        state_t_1.proportions = tmp_theta;
        vector<string> pathNames;
        vector<string> parentpathNames;
        for (auto & p : state_t_1.positions_tree){

            pathNames.emplace_back(p.pos->longname);
            //cerr << "first pathNames " << p.pos->longname << endl; 
            
            if (p.pos->parent != nullptr) {
                parentpathNames.emplace_back(p.pos->parent->longname);

            }else{
                parentpathNames.emplace_back(p.pos->longname);
            }
        }
        
#ifdef DEBUGMCMC
        cerr << "state_t_1.proportions " << state_t_1.proportions[0] << " " << state_t_1.proportions[1] << endl;
        cerr << "pathNames " << pathNames[0] << " " << pathNames[1] << endl;
#endif

#ifdef DEBUGPAIN

        //state_t_1.positions_tree[0].pos_branch = 0.9999;
        //state_t_1.positions_tree[0].pos_branch = 1 - 0.9999;

        //pathNames[0] = "N3Ursidae";
        //pathNames[0] = "N4Ursidae";
        //pathNames[0] = "N10Ursidae";

        
        //cerr << pathNames[0] << endl;
        int mismatch = 0;
        int mismatchP = 0;
        int unsup = 0;
        int unsupP = 0;

        cerr << "Calculating log likelihood for child " << pathNames[0] << " at branch position " << state_t_1.positions_tree[0].pos_branch << endl;
        cerr << "Calculating log likelihood for parent" << parentpathNames[0] << " at branch position " << 1- state_t_1.positions_tree[0].pos_branch << endl;
#endif

        double logLike = 0.0;
        #pragma omp parallel for num_threads(n_threads) reduction(+:logLike)
        for (auto read : *(params.align))
        {   
            
            double readLogLike = 0.0;
            double readLogLikeP = 0.0;
            // single source 
            if (pathNames.size() == 1)
            {
                int counter = 0;
                if (state_t_1.positions_tree[0].pos_branch == 1.0){
                    state_t_1.positions_tree[0].pos_branch == 0.9999999999;
                }
                // getting the proportion of t:
                double t = 0.0;
                t = state_t_1.positions_tree[0].pos->dist;
                if (t == 0.0){ //this will only happen at the root node. This node does not have a branch length because we have an unrooted tree.
                    t = 0.00001;
                }
                double t1 = state_t_1.positions_tree[0].pos_branch * t;
                double t2 = t - t1;
                
#ifdef DEBUGPAIN
                cerr << "Sequence " << read->name << " " << read->seq << endl;
                cerr << "pathName " << pathNames[0] << endl;
                cerr << "Size of vector at key '" << pathNames[0] << "': " << read->detailMap[pathNames[0]].size() << endl;
                if (read->detailMap.find(pathNames[0]) != read->detailMap.end()){
                    cerr << "Found it " << endl;
                }else{
                    cerr << "can't find it " << endl;
                }
#endif
                for(int basevec = 0; basevec < read->detailMap[pathNames[0]].size(); ++basevec)
                {
                    //cerr << "I enter the loop 1" << endl;
                    for (int base = 0; base < read->detailMap[pathNames[0]][basevec].size(); ++base)
                    {
#ifdef DEBUGPAIN
                        cerr << "I enter the loop 2" << endl;
                         cerr << "read "<< read->detailMap[pathNames[0]][basevec][base].readBase << " reference " << read->detailMap[pathNames[0]][basevec][base].referenceBase << endl;
#endif
                        if (read->detailMap[pathNames[0]][basevec][base].pathSupport)
                        {
                            counter++;
                            readLogLike += computeBaseLogLike(read, params, basevec, base, pathNames[0], t2, t);

                            
                            
#ifdef DEBUGPAIN
                            
                            if (read->detailMap[pathNames[0]][basevec][base].readBase != read->detailMap[pathNames[0]][basevec][base].referenceBase){
                                mismatch++;
                                cerr << "Mismatch ? " << read->detailMap[pathNames[0]][basevec][base].readBase << " reference " << read->detailMap[pathNames[0]][basevec][base].referenceBase << endl;
                            }
                            if (isnan(readLogLike) || isnan(readLogLikeP)){
                                cerr << "broke in computation for path" << endl;
                            }
#endif
                        }
                        else
                        {
                            counter++;
                            readLogLike += read->detailMap[pathNames[0]][basevec][base].logLikelihood;
                            
#ifdef DEBUGPAIN
                            
                            if (read->detailMap[pathNames[0]][basevec][base].readBase != read->detailMap[pathNames[0]][basevec][base].referenceBase){unsup++;}
                            if (isnan(readLogLike) || isnan(readLogLikeP)){
                                cerr << "broke in accessing pre calc path " << endl;
                            }
#endif
                        }
                        if(read->detailMap[parentpathNames[0]][basevec][base].pathSupport)
                        {
                            counter++;
                            readLogLikeP += computeBaseLogLike(read, params,basevec, base, parentpathNames[0], t1, t);


#ifdef DEBUGPAIN
                            
                            if (read->detailMap[parentpathNames[0]][basevec][base].readBase != read->detailMap[parentpathNames[0]][basevec][base].referenceBase)
                                {mismatchP++;}
                            if (isnan(readLogLike) || isnan(readLogLikeP)){
                                cerr << "broke in computation for parent path" << endl;
                            }
#endif
                        }
                        else
                        {
                            counter++;
                            readLogLikeP += read->detailMap[parentpathNames[0]][basevec][base].logLikelihood;
                            
                           
#ifdef DEBUGPAIN
                            if (read->detailMap[pathNames[0]][basevec][base].readBase != read->detailMap[pathNames[0]][basevec][base].referenceBase){unsupP++;}
                            
                            if (isnan(readLogLike) || isnan(readLogLikeP)){
                                cerr << "broke in accessing pre calc parent path" << endl;
                            }
#endif
                        }

#ifdef DEBUGPAIN

                        {
                        #pragma omp critical
                            //cerr << pathNames[0] << endl;
                        if (readLogLike > 0 || std::isnan(readLogLike) || std::isinf(readLogLike) || readLogLikeP > 0 || std::isnan(readLogLikeP) || std::isinf(readLogLikeP)){


                            cerr << std::setprecision(14)<< "readLogLike " << readLogLike << " readLogLikeP " << readLogLikeP << endl;
                            //cerr << "path " << pathNames << " parent " << parentpathNames << endl;
                            cerr << std::setprecision(14)<< "t " << t << " t1 " << t1 << " t2 " << t2 << endl;
                            cerr << std::setprecision(14)<< "pre calc log like " << read->detailMap[pathNames[0]][basevec][base].logLikelihood << endl;
                            cerr << std::setprecision(14)<< "parent pre calc log like " << read->detailMap[parentpathNames[0]][basevec][base].logLikelihood << endl;
                            double test = computeBaseLogLike(read, params,basevec, base, pathNames[0], t2, t);
                            double what = computeBaseLogLike(read, params,basevec, base, parentpathNames[0], t1, t);
                            cerr << std::setprecision(14)<< "path calc " << test << endl;
                            cerr << std::setprecision(14)<< "parent path calc " << what << endl;
                            throw std::runtime_error("Error: Log-likelihood is nan.");
                        }}

#endif

                        
                    }
                }
                

                double interProp = calculateLogWeightedAverage(readLogLike, state_t_1.positions_tree[0].pos_branch, readLogLikeP, (1 - state_t_1.positions_tree[0].pos_branch));

#ifdef DEBUGPAIN
                
                //cerr << std::setprecision(14)<< "prop for child " <<  readLogLikeW << " at position " << state_t_1.positions_tree[0].pos_branch <<  " prop for parent " <<  readLogLikePW << " at position " << 1 - state_t_1.positions_tree[0].pos_branch << endl;
                cerr << std::setprecision(14)<< "weighted average log " << interProp << endl;
                cerr << std::setprecision(14)<< "log for child " <<  readLogLike << " at position " << state_t_1.positions_tree[0].pos_branch <<  " log for parent " <<  readLogLikeP << " at position " << 1 - state_t_1.positions_tree[0].pos_branch << endl;
                
#endif

                logLike += interProp;
                // if (counter != (read->seq.size() * 2)){
                //     cerr << "Counter counted " << counter << " the read size is " << read->seq.size() * 2 << endl;
                //     throw runtime_error("The MCMC did not go throug all necessary bases");
                // }
                //logLike += (readLogLike + log(state_t_1.positions_tree[0].pos_branch));
               
                
                    
                    
#ifdef DEBUGMCMC
                cerr << "intermediate logLike " << logLike << endl;
                cerr << "read->pathMap[parentpathNames[0]] " << read->pathMap[parentpathNames[0]] << endl;
                cerr << "Position on Branch " << state_t_1.positions_tree[0].pos_branch << endl;
#endif

            }
            else
            {
                int counter = 0;
                double inter = -std::numeric_limits<double>::infinity(); 
                for (int y = 0; y < pathNames.size(); ++y)
                {
                    double readLogLike = 0.0;
                    double readLogLikeP = 0.0;
                    double t = 0.0;
                    t = state_t_1.positions_tree[y].pos->dist;
                    if (t == 0.0){ //this will only happen at the root node. This node does not have a branch length because we have an unrooted tree.
                        t = 0.00001;
                    }
                    double t1 = state_t_1.positions_tree[y].pos_branch * t;
                    double t2 = t - t1;
                    if (state_t_1.positions_tree[y].pos_branch == 1.0){
                        state_t_1.positions_tree[y].pos_branch == 0.99999;
                    }
                    //#pragma omp parallel for num_threads(15) private(read)
                    for (int basevec = 0; basevec < read->detailMap[pathNames[y]].size(); ++basevec){

                        for (int base = 0; base < read->detailMap[pathNames[y]][basevec].size(); ++base)
                        {              
                            
                            if (read->detailMap[pathNames[y]][basevec][base].pathSupport)
                            {
                                counter++;
                                readLogLike += computeBaseLogLike(read, params,basevec, base, pathNames[y], t2, t);
                                
                                
                            }
                            else
                            {
                                counter++;
                                readLogLike += read->detailMap[pathNames[y]][basevec][base].logLikelihood;    
                            }
                            if(read->detailMap[parentpathNames[y]][basevec][base].pathSupport)
                            {
                                counter++;
                                readLogLikeP += computeBaseLogLike(read, params,basevec, base, parentpathNames[y], t1, t);
                                
                                
                            }
                            else
                            {
                                counter++;
                                readLogLikeP += read->detailMap[parentpathNames[y]][basevec][base].logLikelihood;
                            }

//#ifdef DEBUGPAIN

                            if (readLogLike > 0.0 || std::isnan(readLogLike) || std::isinf(readLogLike) || readLogLikeP > 0.0 || std::isnan(readLogLikeP) || std::isinf(readLogLikeP)){
                            //cerr << pathNames[y] << endl;
                            //if (parentpathNames[y] == "NC_062361.1_Hippotragus_niger_roosevelti_voucher_ZMUC_H.R.Siegismund_1646_haplogroup_Eastern_1_mitocho" || pathNames[y] == "NC_062361.1_Hippotragus_niger_roosevelti_voucher_ZMUC_H.R.Siegismund_1646_haplogroup_Eastern_1_mitocho"){
                                cerr << std::setprecision(14)<<"state_t_1.proportions[y]" << state_t_1.proportions[y] << endl;
                                cerr << std::setprecision(14)<< "readLogLike " << readLogLike << " readLogLikeP " << readLogLikeP << endl;
                                //cerr << "path " << pathNames << " parent " << parentpathNames << endl;
                                cerr << std::setprecision(14)<< "t " << t << " t1 " << t1 << " t2 " << t2 << endl;
                                cerr << std::setprecision(14)<< "pre calc log like " << read->detailMap[pathNames[y]][basevec][base].logLikelihood << endl;
                                cerr << std::setprecision(14)<< "parent pre calc log like " << read->detailMap[parentpathNames[y]][basevec][base].logLikelihood << endl;
                                double test = computeBaseLogLike(read, params,basevec, base, pathNames[y], t2, t);
                                double what = computeBaseLogLike(read, params,basevec, base, parentpathNames[y], t1, t);
                                cerr << std::setprecision(14)<< "path calc " << test << endl;
                                cerr << std::setprecision(14)<< "parent path calc " << what << endl;

                                throw runtime_error("Problem in the likelihood compuation! Intermediate log likelihood is -nan, -inf or positive.");}
//#endif

                        }
                            
                            
                    }
                    double interc2 = log(state_t_1.positions_tree[y].pos_branch) + readLogLike;
                    //cerr << "interc2 " << interc2 << endl;
                    double interp2 = log((1 - state_t_1.positions_tree[y].pos_branch)) + readLogLikeP;
                    //cerr << "interp2 " << interp2 << endl;
                    //double inter2 = calculateLogWeightedAverage(readLogLike, state_t_1.positions_tree[y].pos_branch, readLogLikeP,(1 - state_t_1.positions_tree[y].pos_branch));
                    //cerr << "inter2 " << inter2 << endl;
                    double inter2 = oplusnatl(interc2, interp2);
                    inter = oplusInitnatl(inter , (inter2 + log(state_t_1.proportions[y])));
                    //cerr << "inter " << inter << endl;

                }
                
                
                logLike += inter;
                //cerr << "Increasing logLike " << logLike << endl;
                //cerr << "Total number of bases " << counter << endl;
                // if (counter != (read->seq.size() * 2)){
                //     cerr << "Counter counted " << counter << " the read size is " << read->seq.size()*2 << endl;
                //     throw runtime_error("The MCMC did not go throug all necessary bases");
                // }


            }


            
        }

#ifdef DEBUGPAIN
        cerr << std::setprecision(14) << "Log likelihood " << logLike << endl;

        cerr << "Path name child " << pathNames[0] << " Pathname parent " << parentpathNames[0] << endl;
        cerr << "Path name child 2 " << pathNames[1] << " Pathname parent 2" <<parentpathNames[1] << endl;


#endif
        
        bool initialized=false;
        vector<double> initialPatristicDistances;

        // This should be done once before the loop starts, at the place where initialPatristicDistances is first accessible.
        if (!initialized) {
            initialPatristicDistances = vector<double>(15, 1.0);
            initialized = true;
        }
        likelihood_t_1 = logLike;
#ifdef DEBUGMCMC
        cerr << "logLike " << likelihood_t_1 << endl;
#endif
        state_t_1.logLike = likelihood_t_1;

        double acceptance_prob = (state_t_1.logLike - state_t.logLike > 0) ? 1.0 : exp(state_t_1.logLike - state_t.logLike);
        double u = dis(gen);
        

        if (u <= acceptance_prob || iteration == 0) {

            for (auto p : state_t_1.positions_tree){
                mcmcdetail << std::setprecision(14) << p.pos->longname  << "\t" << state_t_1.logLike << '\t' << p.theta << '\t' <<  p.pos_branch << '\t' << "accepted" << '\t';

            }
            mcmcdetail << endl;
            
             
            if (iteration > params.burn){
                accept_count++;
                for (auto & p : state_t.positions_tree){
                    
                    mcmcout << setprecision(14) << p.pos->longname << '\t' << state_t.logLike << '\t' << p.theta << '\t' <<  p.pos_branch << '\t';
                }
                state_t_vec.emplace_back(state_t);
                mcmcout << endl;
                
#ifdef DEBUGMCMC        
                cerr << std::setprecision(16) << "\t" << state_t.logLike << '\t' << "proposal sd" <<  '\t' << proposal_sd << '\t' << acceptance_rate << endl;
#endif
            }
                                                   

            n_accept++;
            total_proposals++;
            state_t = state_t_1;
            
        } else {

            for (auto p : state_t_1.positions_tree){
                mcmcdetail << std::setprecision(14) << p.pos->longname  << "\t" << state_t_1.logLike << '\t' << p.theta << '\t' <<  p.pos_branch << '\t'<< "rejected" << '\t';

            }
            mcmcdetail << endl;
            if (iteration > params.burn)
            {
#ifdef DEBUGMCMC        
            cerr << "Rejected Log likelihood " << state_t_1.logLike << "  Current: " << state_t.logLike << "  Proposal SD: " << proposal_sd << endl;
#endif
           for (auto & p : state_t.positions_tree){
                mcmcout << setprecision(14) << p.pos->longname << '\t' << state_t.logLike << '\t' << p.theta << '\t' <<  p.pos_branch << '\t';
#ifdef DEBUGMCMC
                cerr << "proposal sd: " << proposal_sd << endl;
                cerr << "Current sources: " << std::setprecision(16) << (p.pos->longname.size() < 65 ? p.pos->longname : p.pos->longname.substr(0, 65)) \
                                         << p.pos->name << "  accept rate: " << acceptance_rate << " proportions " << state_t_1.proportions[0] << '\t' << state_t_1.proportions[1] <<endl;
#endif
                }
                mcmcout << endl;

                state_t_vec.emplace_back(state_t);
             }
                

                                          

            total_proposals++;
        }

        //get_proposal_sd(acceptance_rate, iteration, params.maxIter);
        acceptance_rate = static_cast<double>(n_accept) / total_proposals;
        //cerr << "acception rate " << acceptance_rate << endl;
    }

    return state_t_vec;

}




//////////////// START EUKA MCMC ///////////////////////////////////
// function to generate the proposal vector for the MCMC runs
const vector<double> MCMC::generate_proposal(vector<double> &current_vec, const double &alpha, const bool branch_pos){
	// checking that vector sums up to 1

	double check = 0.0;
	for (size_t p = 0; p<current_vec.size(); p++){
#ifdef DEBUGGENERATEVEC
		cerr << current_vec.at(p) << endl;
#endif
		check += current_vec.at(p);
	}
#ifdef DEBUGGENERATEVEC
	cerr << "check "<< check << endl;
#endif
	double lower = 0.99;
	double upper = 1.01;


#ifdef DEBUGGENERATEVEC
	for (const auto &element:current_vec){
        cerr << "Current vec before log: "<< element << endl;
    }
#endif

        //if(!branch_pos){assert(check > lower && check < upper);}

	vector <double> current_vec_log;
	// transform vector into log space.
	for (size_t q = 0; q<current_vec.size(); q++){
		current_vec_log.emplace_back(log(current_vec.at(q)));

	}
#ifdef DEBUGGENERATEVEC
	for (const auto &element:current_vec_log){
        cerr << "Current vec after log: "<< element << endl;

    }
#endif
    // generate random number generator
	std::random_device rd;
	std::mt19937 g(rd());
	// looping through the vector of log transformed distributions.
	//For each element we will sample from a normal distribution to get a new draw.
        vector<double> projected_vec;

        if (!branch_pos) {

        vector<double> proposal_vec;
	for (double element:current_vec_log){
		std::normal_distribution<double> d{element,alpha};

#ifdef DEBUGGENERATEVEC
		cerr << "exp proposal_vec " << exp(d(g)) << endl;
#endif
		proposal_vec.emplace_back(d(g));
	}

	// transforming the proposal_vec with the softmax function
	projected_vec = MCMC::softmax(proposal_vec);
                        }

        else{
        projected_vec = MCMC::sample_normal_euka(current_vec, alpha);
            }

#ifdef DEBUGGENERATEVEC
	for (auto element:projected_vec){
        cerr << "exp projected_vec: " << element << endl;
    }
#endif
	check = accumulate(projected_vec.begin(), projected_vec.end(), 0.0);
	lower = 0.99;
	upper = 1.01;
	// Check that we have successfully projected back onto the unit simplex

	if(!branch_pos){assert(check > lower && check < upper);}

	return projected_vec;
}


double MCMC::get_proposal_likelihood(const vector <double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id){

	double proposal_log_likelihood = 0.0;
	vector<double> total_f_list;
    for (int i=0; i<clade_list_id.size(); i++){

    	total_f_list.emplace_back(clade_list_id.at(i));
        total_f_list.emplace_back(proposal_vec.at(i));
    }

    for (int j = 0; j<total_f_list.size(); j+=2){
        // Make sure we have an even number of elements
    	assert(total_f_list.size() % 2 == 0);

    	int vec_len = total_f_list.size()/2;
        double frac = total_f_list.at(j+1);
        double total_frac = 0.0;
    	for (int k = 1; k<clade_vec->at(total_f_list.at(j)*6+1)->clade_like.size(); k++){
		total_frac += log(((frac * clade_vec->at(total_f_list.at(j)*6+1)->clade_like.at(k)) + (clade_vec->at(total_f_list.at(j)*6+1)->clade_not_like.at(k)  * (1/334))));
    		//total_frac += log((frac * clade_vec->at(total_f_list.at(j)*6+1)->clade_like.at(k)) + (clade_vec->at(total_f_list.at(j)*6+1)->clade_not_like.at(k)/ vec_len));
#ifdef DEBUGGENERATEVEC
    		cerr << "total_frac " << total_frac << endl;
#endif

    	}

    	proposal_log_likelihood += total_frac;

    	}

#ifdef DEBUGGENERATEVEC
    	cerr << "proposal_log_likelihood: " << proposal_log_likelihood << endl;

#endif

    return proposal_log_likelihood;

	
}



const vector<double > MCMC::run(int iter, int burnin, double tol, const vector<double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id){

	MCMC mcmc;

	vector <double> current_best = init_vec;
	//vector <double> current_best = {0.0909090909,0.0909090909, 0.0909090909, 0.0909090909, 0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909};
	vector <double> proposal_vec;

    double current_log_likelihood = -9999999;
    double acceptance_prob = 0;
    double u=0;
    std::random_device rd;
    std::mt19937 gen(rd());

    struct mcmc_moves{
    	double log_like;
    	vector <double> abund_vec;
    };
#ifdef DEBUGGENERATEVEC
    ofstream outputFile("proposal_log_likelihood.txt", ios::trunc);
    outputFile << "-1" << "\t" << current_log_likelihood << '\t' << current_best[0] << '\t' << current_best[1] << endl;
#endif
    mcmc_moves move[iter];
    int no = 0;

    cerr << "Computing MCMC:" << endl;
	for (int iteration = 0; iteration<iter; iteration++){

		printprogressBarCerr( float(iteration + 1)/float(iter) );

		proposal_vec = mcmc.generate_proposal(current_best, 0.1, false);
        double proposal_log_likelihood = mcmc.get_proposal_likelihood(proposal_vec, clade_vec, clade_list_id);

       // we are not adding proposal vectors to the struct before after the burin in period
		if (iteration > burnin){

        move[iteration].log_like = proposal_log_likelihood;
        move[iteration].abund_vec = proposal_vec;
         
    	}
        else{

        	continue; 
        }

        

#ifdef DEBUGGENERATEVEC
        outputFile << iteration << '\t' << proposal_log_likelihood << '\t' << proposal_vec[0] << '\t' << proposal_vec[1] <<'\t' << proposal_vec[2] << '\t' << proposal_vec[3] <<'\t' << '\t' << proposal_vec[4] << '\t' << proposal_vec[5] <<'\t' << '\t' << proposal_vec[6] << '\t' << proposal_vec[7] <<'\t' << '\t' << proposal_vec[8] << '\t' << proposal_vec[9] <<'\t' << '\t' << proposal_vec[10] << '\t';

        cerr << "proposal_log_likelihood " << proposal_log_likelihood << endl; 
#endif
        acceptance_prob = min((double)(1.0), exp(proposal_log_likelihood-current_log_likelihood));
#ifdef DEBUGGENERATEVEC
        cerr << "acceptance_prob " << acceptance_prob << endl; 
#endif
        std::uniform_real_distribution<> dis(0, 1);
        u = dis(gen);
#ifdef DEBUGGENERATEVEC
        cerr << u << endl;
        //proposal_log_likelihood > current_log_likelihood
#endif
        if (u <= acceptance_prob || iteration == 0) {
#ifdef DEBUGGENERATEVEC
        	outputFile << "accept" << endl;
#endif 
#ifdef VERBOSE_MCMC
            cerr << "ACCEPTING proposal." << endl;
            cerr << "Proposal log likelihood: " << proposal_log_likelihood << endl;
#endif
            current_log_likelihood = proposal_log_likelihood;
            current_best = proposal_vec;
                                                    }
        else {
#ifdef DEBUGGENERATEVEC
        	outputFile << "reject" << endl; 
#endif
#ifdef VERBOSE_MCMC
            cerr << "REJECTING proposal" << endl;
#endif
            current_best = current_best;
             }
#ifdef VERBOSE_MCMC
        cerr << "\n\nCurrent log likelihood: " << current_log_likelihood << endl;
#endif
	}
	cerr<<endl;

#ifdef VERBOSE_MCMC
	cerr << "MCMC completed. The final log likelihood is: " << current_log_likelihood << endl; 
#endif


	
	int per85 = 0.85 * (iter - burnin); 
	int per95 = 0.95 * (iter - burnin);
	double sums[proposal_vec.size()] = {0.0};
	vector<double> posterior_estimate;
	vector<double> sorted_clade; 

	// The posterior mean as well as the confidence intervalls will be calculated per clade (each clade is independent from each other)
	//We are lọoping through the abundance vector with j
	// 
	for (int j=0; j<proposal_vec.size(); j++){

		for (int i = burnin+1; i<iter; i++){
			
			sums[j] += move[i].abund_vec[j];
#ifdef DEBUGOUTPUT
			cerr << "before " <<move[i].abund_vec[j] <<endl;  
#endif		
			sorted_clade.emplace_back(move[i].abund_vec[j]);

		}
		// sorting the struct for highest posterior density intervall (HDI)

		sort(sorted_clade.begin(), sorted_clade.end());

		const int m = sorted_clade.size()/2;
		
		double p = sorted_clade[m];


		// posterior point estimate for each clade (each fraction of the abdundance vector)
		//double p = sums[j]/(iter-burnin);
#ifdef DEBUGOUTPUT
		cerr << "posterior_estimate " << p << endl; 
#endif 
		posterior_estimate.emplace_back(p);

		double low_end_85 = quant(sorted_clade, 0.15);
		double high_end_85 = quant(sorted_clade, 0.85);
		double low_end_95 = quant(sorted_clade, 0.05);
		double high_end_95 = quant(sorted_clade, 0.95);

		posterior_estimate.emplace_back(low_end_85);
		posterior_estimate.emplace_back(high_end_85);
		posterior_estimate.emplace_back(low_end_95); 
		posterior_estimate.emplace_back(high_end_95);

		sorted_clade.clear();


		
	}
	
	cerr<<".. done"<<endl;
	return posterior_estimate;

}
