#pragma once
#include "MCMC.h"
#include "Euka.h"
#include "libgab.h"
#include <algorithm>
#include <gzstream.h>

pair<unordered_map<string, vector<vector<double>>>, double> MCMC::processMCMCiterations_sb(const std::vector<MCMCiteration> &MCMCiterations, \
                                                                                           int k, string num, int chain, spidir::Tree* tr, \
                                                                                           int numofleafs) {

    unordered_map<string, vector<vector<double>>> branchStatisticsMap;
    vector<vector<double>> chainStatistic;
    std::ofstream estimatesFile, branchestimateFile;
    unsigned int chainIndex = 0;
    estimatesFile.open(num+"ProportionEstimates"+to_string(k)+".txt", std::ios::app | std::ios::out);
    branchestimateFile.open(num+"BranchEstimate"+to_string(k)+".txt", std::ios::app | std::ios::out);

    // Check if the files are open
    if (!estimatesFile.is_open() || !branchestimateFile.is_open()) {
        std::cerr << "Failed to open the files." << std::endl;
        // Handle error appropriately
    }
    // Write column names to both files.
    estimatesFile << "Source\tChain\tMean Proportion Estimate\t5% CI\tMedian Proportion Estimate\t95% CI\tEffective Sample Size\tAutocorrelation\tVariance\n";
    branchestimateFile << "Source\tChain\tMean Branch Position\t5% CI\tMedian Branch Position\t95% CI\tEffective Sample Size\tAutocorrelation\tVariance\tEffective Sample Size for the source estimation\n";
    double chainloglike = MCMCiterations[0].logLike;
    for (int source = 0; source < k; ++source) {

        vector<double> sourceStatistic = {};
        vector<double> proportionVec = {};
        vector<double> positionVec = {};
        string branchName;
        vector<double> initialPatristicDistances;
        vector<double> euc_distances;
        initialPatristicDistances = vector<double>(numofleafs, 1.0);
        size_t totalIterations = MCMCiterations.size();
        size_t startIdx = totalIterations > 500 ? totalIterations - 500 : 0;

        for (size_t idx = 0; idx < totalIterations; ++idx) {
            const auto& iteration = MCMCiterations[idx];

            if(iteration.logLike > chainloglike){
                chainloglike = iteration.logLike;
            }
            branchName = iteration.positions_tree[source].pos->longname;
            if (branchStatisticsMap.find(branchName) == branchStatisticsMap.end()) {
                // If not, initialize it with an empty vector of vectors
                branchStatisticsMap[branchName] = vector<vector<double>>();
            }   
            
           
            proportionVec.emplace_back(iteration.proportions[source]);
            positionVec.emplace_back(iteration.positions_tree[source].pos_branch);
            double t1 = iteration.positions_tree[source].pos->dist * iteration.positions_tree[source].pos_branch;
            double posonbranch = iteration.positions_tree[source].pos->dist - t1;

            
            const vector<double> patristic_distances = getPatristicDistances(tr,iteration.positions_tree[source].pos, numofleafs, posonbranch);
            //PRINTVEC(patristic_distances)
            double euc_dist = calculateEuclideanDistance(patristic_distances, initialPatristicDistances);
            euc_distances.emplace_back(euc_dist);

            

        }

        
        //unsorted vec
        double meanTheta = mean(proportionVec);
        double meanPos = mean(positionVec);
        double Theta_autoc = autocorrelation(proportionVec, 1);
        double Theta_ess = effectiveSampleSize(proportionVec);
        if(Theta_ess < 100){
            cerr << "Warning: The effective sample size for the proportion estimation of chain "<< chain << " is below 100. The estimation of the proportion for the branch " << branchName << " can not be ensured. A rerun using a higher number of iterations is recommended." << endl;
        }
        double Theat_var = variance(proportionVec, meanTheta);
        double Pos_autoc = autocorrelation(positionVec, 1);
        double Pos_ess = effectiveSampleSize(positionVec);
        //cout << "Pos vect size " << positionVec.size() << endl;
        if(Pos_ess < 100){
            cerr << "Warning: The effective sample size for the estimation of the branch position for chain "<< chain << " is below 100. The estimation of the position for the branch " << branchName << " can not be ensured. A rerun using a higher number of iterations is recommended." << endl;
        }
        double dist_ess = effectiveSampleSize(euc_distances);
        //cout << "euc_distances size " << euc_distances.size() << endl;
        if(dist_ess < 100){
            cerr << "Warning: The effective sample size for the estimation of the branch for chain "<< chain << " is below 100. The estimation of the branch " << branchName << " as a source can not be ensured." << endl;
        }
        double Pos_var = variance(positionVec, meanPos);

        //now we need to sort for quantiles
        sort(positionVec.begin(), positionVec.end());
        sort(proportionVec.begin(), proportionVec.end());
        double Theta_fq = getQuantile2(proportionVec, 0.05);
        double Theta_tq = getQuantile2(proportionVec, 0.95);
        
        double Theta_median = getQuantile2(proportionVec, 0.5);
        double Pos_fq = getQuantile2(positionVec, 0.05);
        double Pos_median = getQuantile2(positionVec, 0.5);
        double Pos_tq = getQuantile2(positionVec, 0.95);

        //estimation of autocorrelation and effective sample size. 
        //rule of thump: higher ESS and lower correlation is good, indicating better mixing between sample. 
        

        
        estimatesFile << branchName << '\t' << chain <<'\t' << meanTheta << '\t' << Theta_fq << '\t' << Theta_median << '\t' << Theta_tq << '\t' << Theta_ess << '\t' << Theta_autoc <<'\t' << Theat_var <<'\n';
        branchestimateFile << branchName << '\t' << chain <<'\t' << meanPos << '\t' << Pos_fq << '\t' << Pos_median << '\t' << Pos_tq <<'\t' << Pos_ess << '\t' << Pos_autoc <<'\t' << Pos_var <<'\t'<<dist_ess<<'\n';
        


        sourceStatistic.emplace_back(meanTheta);
        sourceStatistic.emplace_back(Theat_var);
        sourceStatistic.emplace_back(meanPos);
        sourceStatistic.emplace_back(Pos_var);


        branchStatisticsMap[branchName].emplace_back(sourceStatistic);
    }

    return make_pair(branchStatisticsMap, chainloglike);


}


inline double MCMC::computeBaseLogLike(const AlignmentInfo* read, RunTreeProportionParams &params, const int basevec, const int base, \
                                       const string &pathName, const double t, const double branch_len)
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


std::vector<MCMCiteration> MCMC::run_tree_proportion_sb(RunTreeProportionParams params, std::vector<MCMCiteration> state_t_vec, \
const bdsg::ODGI& graph, vector<vector<string>> nodepaths, string num, int n_threads, int numPaths, int chainindex, double con) {

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
                            readLogLike += computeBaseLogLike(read, params, basevec, base, pathNames[0], t2, con);

                            
                            
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
                            readLogLikeP += computeBaseLogLike(read, params,basevec, base, parentpathNames[0], t1, con);


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
                            cerr << "path " << pathNames << " parent " << parentpathNames << endl;
                            cerr << std::setprecision(14)<< "t " << t << " t1 " << t1 << " t2 " << t2 << endl;
                            cerr << std::setprecision(14)<< "pre calc log like " << read->detailMap[pathNames[0]][basevec][base].logLikelihood << endl;
                            cerr << std::setprecision(14)<< "parent pre calc log like " << read->detailMap[parentpathNames[0]][basevec][base].logLikelihood << endl;
                            double test = computeBaseLogLike(read, params,basevec, base, pathNames[0], t2, t2); 
                            double what = computeBaseLogLike(read, params,basevec, base, parentpathNames[0], t1, t1);
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
                                readLogLike += computeBaseLogLike(read, params,basevec, base, pathNames[y], t2, con);//;
                                
                                
                            }
                            else
                            {
                                counter++;
                                readLogLike += read->detailMap[pathNames[y]][basevec][base].logLikelihood;    
                            }
                            if(read->detailMap[parentpathNames[y]][basevec][base].pathSupport)
                            {
                                counter++;
                                readLogLikeP += computeBaseLogLike(read, params,basevec, base, parentpathNames[y], t1, con); //
                                
                                
                            }
                            else
                            {
                                counter++;
                                readLogLikeP += read->detailMap[parentpathNames[y]][basevec][base].logLikelihood;
                            }

#ifdef DEBUGPAIN

                            if (readLogLike > 0.0 || std::isnan(readLogLike) || std::isinf(readLogLike) || readLogLikeP > 0.0 || std::isnan(readLogLikeP) || std::isinf(readLogLikeP)){
                            //cerr << pathNames[y] << endl;
                            //if (parentpathNames[y] == "NC_062361.1_Hippotragus_niger_roosevelti_voucher_ZMUC_H.R.Siegismund_1646_haplogroup_Eastern_1_mitocho" || pathNames[y] == "NC_062361.1_Hippotragus_niger_roosevelti_voucher_ZMUC_H.R.Siegismund_1646_haplogroup_Eastern_1_mitocho"){
                                cerr << std::setprecision(14)<<"state_t_1.proportions[y]" << state_t_1.proportions[y] << endl;
                                cerr << std::setprecision(14)<< "readLogLike " << readLogLike << " readLogLikeP " << readLogLikeP << endl;
                                cerr << "path " << pathNames << " parent " << parentpathNames << endl;
                                cerr << std::setprecision(14)<< "t " << t << " t1 " << t1 << " t2 " << t2 << endl;
                                cerr << std::setprecision(14)<< "pre calc log like " << read->detailMap[pathNames[y]][basevec][base].logLikelihood << endl;
                                cerr << std::setprecision(14)<< "parent pre calc log like " << read->detailMap[parentpathNames[y]][basevec][base].logLikelihood << endl;
                                double test = computeBaseLogLike(read, params,basevec, base, pathNames[y], t2, con); 
                                double what = computeBaseLogLike(read, params,basevec, base, parentpathNames[y], t1, con);
                                cerr << std::setprecision(14)<< "path calc " << test << endl;
                                cerr << std::setprecision(14)<< "parent path calc " << what << endl;

                                throw runtime_error("Problem in the likelihood compuation! Intermediate log likelihood is -nan, -inf or positive.");}
#endif

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
