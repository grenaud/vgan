#include "TrailMix.h"
#include "HaploCart.h"
#include "gam2prof.h"
#include "time.h"
#include <algorithm>
#include <functional>
#include "rpvg_main.hpp"
#include "damage.h"
#include "sys/wait.h"
#include "MCMC.h"
#include "trailmix_functions.h"
#include "precompute.h"
#include "writeDeconvolvedReads.h"
#include "Dup_Remover.h"

#define PRINT_KEYS(map) do { \
    for (const auto& pair : map) { \
        std::cerr << pair.first << "\n"; \
    } \
} while(0)

#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << setprecision(10) << v[i] << '\t';}cerr << endl << endl;
#define RPVG

#define PRINT_SET(SET) \
    do { \
        std::cerr << "Set contents:\n"; \
        for (const auto& elem : SET) { \
            std::cerr << " - " << elem << '\n'; \
        } \
    } while (0)

void Trailmix::run_mcmc(shared_ptr<Trailmix_struct>& dta) {

    string prefix=dta->TM_outputfilename;

    vector<string> rpvgarguments;
    RunTreeProportionParams params(dta->tree);
    params.chains = dta->chains;
    std::vector<std::vector<MCMCiteration>> MCMCiterationsVec(params.chains);
    MCMC mcmc;

    for (int i = dta->minid; i < dta->minid + dta->nodevector.size(); ++i) {
        vector<string> paths;
        uint64_t node = i;
        bdsg::handle_t node_handle = dta->graph.get_handle(node);
        paths = Trailmix::paths_through_node(dta->graph, node_handle);
        dta->nodepaths.emplace_back(paths);
    }

   initializeParams(params, dta);

        string first_fifo = dta->tmpdir + random_string(9);
        dta->fifo_A = first_fifo.c_str();

dta->pid1 = fork();

if (dta->pid1 == -1) {
    throw runtime_error("Error in fork");
}

if (dta->pid1 == 0) {
    if (dta->gamfilename == "") {
        Trailmix::map_giraffe("", dta);
    }
    exit(0);
}
int status;
waitpid(dta->pid1, &status, 0);

dta->algnvector = readGAM(dta);

/*
if (dta->quiet == false && dta->fastafilename=="") {cerr << "Removing PCR duplicates ..." << '\n';}
shared_ptr<vector<bool>> thing = Dup_Remover().remove_duplicates_internal(dta->algnvector, dta->n_threads, dta->quiet);

auto& vec = *(dta->algnvector); // Access the vector object
auto it = vec.begin();
for (size_t i = 0; it != vec.end();) { // Use the vector's iterator and bounds
    if (!(*thing)[i]) {
        // Erase returns the iterator to the next element, so no increment in this branch
        it = vec.erase(it);
    } else {
        ++it;
        ++i; // Increment the index only when not erasing
    }
}
*/

#ifdef RPVG
dta->rpvg_pid = fork();

if (dta->rpvg_pid == -1) {
    throw runtime_error("Error in RPVG fork");
}
#endif

if (dta->rpvg_pid == 0) {
    rpvgarguments = buildRpvgArguments(dta);
     char** rpvgargvtopass = new char*[rpvgarguments.size()];
        for (size_t i = 0; i < rpvgarguments.size(); i++) {
            rpvgargvtopass[i] = const_cast<char*>(rpvgarguments[i].c_str());
        }
   #ifdef RPVG
    int retcode1 = rpvg_main(rpvgarguments.size(), rpvgargvtopass);
    rpvgarguments.clear();
    exit(0);
   #endif
}

waitpid(dta->rpvg_pid, &status, 0);

#ifdef RPVG
dta->rpvg_ht_pid = fork();

if (dta->rpvg_ht_pid == -1) {
    throw runtime_error("Error in RPVG fork");
}
#endif

if (dta->rpvg_ht_pid == 0) {
    rpvgarguments = buildRpvgArgumentsHaplotypeTranscripts(dta);
    char** rpvgargvtopass = new char*[rpvgarguments.size()];
        for (size_t i = 0; i < rpvgarguments.size(); i++) {
            rpvgargvtopass[i] = const_cast<char*>(rpvgarguments[i].c_str());
        }
   #ifdef RPVG
    int retcode1 = rpvg_main(rpvgarguments.size(), rpvgargvtopass);
    rpvgarguments.clear();
    exit(0);
   #endif
}

waitpid(dta->rpvg_ht_pid, &status, 0);

//cerr << "TEMPDIR: " << dta->tmpdir << endl;

           #ifdef RPVG
           while ((dta->rpvg_pid = wait(&dta->rpvg_status)) > 0);
           while ((dta->rpvg_ht_pid = wait(&dta->rpvg_ht_status)) > 0);
           #endif

            dta->rpvg_status = 0;
            dta->rpvg_ht_status = 0;

            if (dta->debug) {
                cerr << "Loading RPVG output" << endl;
            }

           #ifdef RPVG
            load_hap_combos(dta);
            load_tpms(dta);
            assert(!dta->tpms.empty());
            load_read_probs(dta);
            Trailmix::create_path_node_map(dta);

// Convert dta->tpms to a hash table for faster lookups
std::unordered_map<decltype(dta->tpms)::value_type::first_type, decltype(dta->tpms)::value_type::second_type> tpms_map;
for (const auto& p : dta->tpms) {
    tpms_map[p.first] = p.second;
}

bool in_pruned = false;

for (int m = 0; m < dta->path_names.size(); ++m) {
    in_pruned = false; // Reset in_pruned at the beginning of each iteration

    auto it = tpms_map.find(dta->path_names[m]);
    if (it != tpms_map.end()) {
        if (!isnan(it->second) && it->second > 1.0) {
            in_pruned = true;
        }
    }

bool in_pruned=false;

for (int m = 0; m < dta->path_names.size(); ++m) {
    in_pruned = false; // Reset in_pruned at the beginning of each iteration

    auto it = tpms_map.find(dta->path_names[m]);
    if (it != tpms_map.end()) {
        if (!isnan(it->second) && it->second > 1.0) {
            in_pruned = true;
            //std::cerr << "TPM for '" << dta->path_names[m] << "' is " << it->second << " and is valid and greater than 1." << std::endl;
        }
    }

    if (in_pruned) {
        string toadd = dta->path_names[m];
        dta->in_pruned_set.insert(toadd);
        //std::cerr << "Inserting '" << dta->path_names[m] << "' into in_pruned_set as its TPM is valid and greater than 1." << std::endl;
    }
}

}
           #endif


//std::cerr << "PRUNED SET SIZE: " << dta->in_pruned_set.size() << std::endl;
//cerr << "PRUNED SET CONTENTS: " << endl;
//PRINT_SET(dta->in_pruned_set);

            if (dta->debug) {
                cerr << "Number of hap combos: " << dta->hap_combos.size() << endl;
            }

//cerr << "Precomputing..." << endl;

                auto gam = precompute_GAM_trailmix(dta);

//cerr << "Done with precomputations" << endl;


double max_posterior = -1.0; // start with a very low value; assuming posterior probabilities are non-negative
vector<int> sigNodes;

vector<spidir::Node*> max_single_combo; // This will store the combo with the highest posterior

for (auto& combo : dta->hap_combos) {
    vector<spidir::Node*> single_combo;
    double current_posterior = -1.0;
    bool invalidCombo = false;

std::unordered_set<std::string> seen_elements;

for (size_t i = 0; i < combo.size(); ++i) {

// Check for duplicates
        if (seen_elements.find(combo[i]) != seen_elements.end()) {
            invalidCombo = true;
            break; // Break out of loop if duplicate found
        }

    if (i < combo.size() - 2) {
        modifyPathNameInPlace(dta, combo[i]);
        auto found_itr = find(dta->path_names.begin(), dta->path_names.end(), combo[i]);

        if (found_itr == dta->path_names.end()) {
            std::cerr << "Path not found: " << combo[i] << std::endl;
            invalidCombo = true;
            break;
        }

        else{
        single_combo.emplace_back(dta->tree->nodes[dta->path_node_map[combo[i]]]);
            }

    } else if (i == combo.size() - 1) {
        current_posterior = stod(combo[i]);
    }
}

    if (invalidCombo) {
        //cerr << "Invalid combo, continuing" << endl;
        continue;  // Skip the current combo and move on to the next one
    }

    //cerr << "CURRENT POSTERIOR: " << current_posterior << endl;
    // If current_posterior is the new max, store the current combo
    if (current_posterior > max_posterior) {
        max_posterior = current_posterior;
        max_single_combo = single_combo;
    }
}


for (auto& node : max_single_combo) {
    cerr << "SIGNATURE NODE: " << dta->originalPathNames[node->longname] << endl;
    sigNodes.emplace_back(node->name);
}

if (!sigNodes.empty()){
    vector<string> topass;
    for (auto& node : max_single_combo) {
      topass.emplace_back(node->longname);
    }

    Trailmix::get_seed_source_estimates(dta, topass);
    //dta->seed = std::vector<double>(dta->k, 1.0 / dta->k);
                      }

    if (!sigNodes.empty()){
    params.sources = sigNodes;
    //std::random_device rd;  // Random number generator
    //std::mt19937 gen(rd()); // Seed the generator
    //std::uniform_int_distribution distrib(0, dta->tree->nodes.size() - 1);

    //params.sources.resize(dta->k); // Resize vector to hold k elements
    //for (int i = 0; i<dta->k; ++i) {
    //    params.sources[i] = distrib(gen); // Assign random number within the range
   // }
                          }

else {
    std::random_device rd;  // Random number generator
    std::mt19937 gen(rd()); // Seed the generator
    std::uniform_int_distribution distrib(0, dta->tree->nodes.size() - 1);

    params.sources.resize(dta->k); // Resize vector to hold k elements
    for (int i = 0; i < dta->k; ++i) {
        params.sources[i] = distrib(gen); // Assign random number within the range
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




            // This map will store the vectors of statistics for all chains, with branch names as keys
            unordered_map<string, vector<vector<vector<double>>>> branchStatsMap;
            vector<double> chainLogLikes;
              // Open diagnostics file for this chain
                ofstream diagnostics;
                diagnostics.open(prefix + "Diagnostics.txt");
                diagnostics << "Source\tMax Log-likelihood\tFor Chain\tProportion Rhat\tBranch Placement Rhat" << endl;

            for (unsigned int chain = 0; chain < params.chains; ++chain) {
                initializeParams(params, dta);
                cerr << "On chain: " << chain + 1 << "  out of " << params.chains << endl;
                params.probMatrix = convertMapsToVector(gam);
                params.align = gam;
                std::vector<MCMCiteration> chainiter = move(mcmc.run_tree_proportion_TM(params, MCMCiterationsVec[chain], dta->graph, dta->nodepaths, \
                                                                                dta->TM_outputfilename, dta, true, chain));


                  // Process the MCMC iterations to get the statistics map for this chain
                pair<unordered_map<string, vector<vector<double>>>, double> intermStatsMapPair = \
                    move(mcmc.processMCMCiterations(dta, chainiter, dta->k, prefix, chain, dta->tree, dta->n_leaves, dta->include_read_props));


                auto &intermStatsMap = intermStatsMapPair.first;

                chainLogLikes.emplace_back(intermStatsMapPair.second);
                // For each branch name, add the current chain's statistics vector to the main map
                for (const auto& branchStat : intermStatsMap) {
                    const auto& branchName = branchStat.first;
                    const vector<vector<double>>& statsForCurrentChain = branchStat.second;

                    // Ensure that we have a vector to store the stats for this branch
                    if (branchStatsMap.find(branchName) == branchStatsMap.end()) {
                        branchStatsMap[branchName] = vector<vector<vector<double>>>();
                    }

                    // Add the stats for the current chain to the corresponding branch in the map
                    branchStatsMap[branchName].emplace_back(statsForCurrentChain);
                }

                cerr << "Finished running chain: " << chain+1 << endl;
            }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Now, calculate R-hat statistics for each branch across all chains
            int numChains = params.chains;
            int chainLength = dta->iter - dta->burnin; // Assume all chains have the same length


            for (const auto& branchStat : branchStatsMap) {
                const auto& branchName = branchStat.first;
                const auto& allChainStats = branchStat.second;

                if(allChainStats.empty()){throw runtime_error("all Chains vector is empty");}
                std::vector<double> Propmeans(numChains);
                std::vector<double> Propvariances(numChains);
                std::vector<double> Posmeans(numChains);
                std::vector<double> Posvariances(numChains);

                // Collect the statistics for each chain for this branch the branch can have only one name for now

                for (int chain = 0; chain < numChains; ++chain) {

cerr << "Processing chain: " << chain << endl;

if (allChainStats.size() <= chain) {
    continue;
}
if (allChainStats[chain].empty()) {
    continue;
}
if (allChainStats[chain][0].size() < 4) {
    continue;
}
    // Now it's safe to access the elements
    Propmeans[chain] = allChainStats[chain][0][0];
    Propvariances[chain] = allChainStats[chain][0][1];
    Posmeans[chain] = allChainStats[chain][0][2];
    Posvariances[chain] = allChainStats[chain][0][3];

}

                double maxLogLike = chainLogLikes[0];

                int maxIndex = 0;
                for (int h = 0; h <chainLogLikes.size(); ++h){
                    if(chainLogLikes[h] > maxLogLike){
                        maxLogLike = chainLogLikes[h];
                        maxIndex = h;
                    }
                }

                // Calculate R-hat for proportions

                double PropRhat = mcmc.calculateRhat(Propmeans, Propvariances, chainLength, params.chains);
                if(isnan(PropRhat)){
                    cerr << "Warning: R-hat for the proportion estimate is not computed because we are handling a single source." << endl;
                }
                else if (PropRhat == -1){
                    cerr << "Warning: The R-hat for proportion cannot be computed for a single chain." << endl;
                }
                else if(PropRhat > 1.05){
                    std::cerr << "Warning: R-hat for proportion of branch " << branchName << " is above 1.05, indicating that the chains have not converged for the parameter." << std::endl;
                }

                // Calculate R-hat for positions

                double PosRhat = mcmc.calculateRhat(Posmeans, Posvariances, chainLength, params.chains);
                if(isnan(PosRhat)){
                    cerr << "Warning: R-hat cannot be computed. The value for the parameter is identical in every iteration." << endl;
                }
                else if (PosRhat == -1){
                    cerr << "Warning: The R-hat for position cannot be computed for a single chain." << endl;
                }
                else if(PosRhat > 1.05){
                    std::cerr << "Warning: R-hat for position of branch " << branchName << " is above 1.05, indicating that the chains have not converged for the parameter." << std::endl;
                }

                ofstream diagnostics2;
                diagnostics2.open(prefix + "Diagnostics.txt", std::ios_base::app);
                // Write the R-hat statistics for this branch to the diagnostics file
                diagnostics2 << branchName << '\t' << maxLogLike << '\t' << maxIndex << '\t' << PropRhat << '\t' << PosRhat << endl;

            }

}
