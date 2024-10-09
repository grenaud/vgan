#pragma once
#include "soibean.h"
#include "bdsg/odgi.hpp"
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
#include "vg/io/vpkg.hpp"
#include "vg/io/protobuf_iterator.hpp"
#include "utility.hpp"
#include "alignment.hpp"
#include "AlignmentInfo.h"
#include "libgab.h"
#include "vgan_utils.h"
#include "MCMC.h"
//#define DEBUG_LLS
//#define DEBUGANALYSEGAM
//#define DEBUGDAMAGE
//#define DEBUGLLM
//#define LOG_INTERMEDIATE_VALUES

#ifdef LOG_INTERMEDIATE_VALUES
    #define LOG_VALUE(var, value, name, read_base, ref_base, bpd) \
        std::cerr << "Variable: " << #var << " = " << std::setprecision(20) << value \
                  << ", Name: " << name \
                  << ", Read Base: " << read_base \
                  << ", Ref Base: " << ref_base \
                  << ", BPD: [" << bpd[0] << ", " << bpd[1] << ", " << bpd[2] << ", " << bpd[3] << "]" \
                  << std::endl;
#else
    #define LOG_VALUE(var, value, name, read_base, ref_base, bpd)
#endif


using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

static vector<AlignmentInfo*>* precompute_GAM(
    const bdsg::ODGI& graph,
    string gamfilename,
    vector<Clade*>* clade_vec,
    const vector<vector<string>>& nodepaths,
    vector<string>& pathNames,
    const vector<double>& qscore_vec,
    bool trailmix,
    int minid,
    bool singlesource,
    vector<vector<diNucleotideProb>>& subDeamDiNuc,
    shared_ptr<Trailmix_struct> dta = NULL
)
{
    int n_reads = 0;
    bool print_dm = true;

    if (gamfilename.empty()) {
        throw std::runtime_error(gamfilename + " is empty. Aborting.");
    }
    if (nodepaths.empty()) {
        throw std::runtime_error("Node paths vector is empty. Aborting.");
    }
    if (pathNames.empty()) {
        throw std::runtime_error("Path names vector is empty. Aborting.");
    }
    if (qscore_vec.empty()) {
        throw std::runtime_error("Quality score vector is empty. Aborting.");
    }

    vector<AlignmentInfo*>* read_vec = new vector<AlignmentInfo*>();
    std::ifstream gam_file(gamfilename, std::ios::binary);
    vg::io::ProtobufIterator<vg::Alignment> iter(gam_file);

    std::unordered_map<std::string, double> results_map;

    // Define nucleotide to index mapping
    std::unordered_map<char, int> nucleotide_index = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    for (; iter.has_current(); iter.advance()) {
        ++n_reads;
        if (n_reads % 10 == 0) {
            cerr << "On read: " << (int)(n_reads / 2) << endl;
        }
        const vg::Alignment& a = *iter;

        if (a.identity() == 0) {
            continue;
        }

        double p_correctly_mapped = 1 - dta->incorrect_mapping_vec[a.mapping_quality()];
        p_correctly_mapped = max(p_correctly_mapped, 0.001);

        if (a.identity() != 0) {
            using BaseInfo = AlignmentInfo::BaseInfo;
            auto&& [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());
            int baseIX = a.path().mapping()[0].position().is_reverse() ? a.sequence().size() - 1 : 0;
            int baseOnRead = baseIX;
            size_t Lseq = a.sequence().size();

            unordered_map<string, vector<vector<BaseInfo>>> detailMap;

            for (size_t c = 0; c < pathNames.size(); ++c) {
                detailMap[pathNames[c]] = vector<vector<BaseInfo>>();
            }

#ifdef DEBUGANALYSEGAM
            unordered_map<string, tuple<string, int, int, int, int, int, int, int, double>> sancheck;
            for (int c = 0; c < pathNames.size(); ++c) {
                sancheck[pathNames[c]] = make_tuple(a.name(), 0, 0, 0, 0, 0, 0, 0, 0.0L);
            }
#endif

            for (int i = 0; i < mppg_sizes.size(); ++i) {

#ifdef DEBUGANALYSEGAMMAP
                cerr << "mapping size " << mppg_sizes.size() << endl;
                cerr << "iter " << i << endl;
#endif

                std::random_device rd;
                std::mt19937 gen(rd());
                // Decide which range to choose
                std::uniform_int_distribution<> range_choice(0, 4);

                int nID = 0;
                int pangenome_base = 100;
                double mappability = 1.0;
                std::unordered_map<std::string, int> pathNameToID;
                std::vector<std::string> idToPathName;

                for (int i = 0; i < dta->path_names.size(); ++i) {
                    pathNameToID[dta->path_names[i]] = i;
                    idToPathName.push_back(dta->path_names[i]);
                }
                std::unordered_set<int> probPaths;

                if (mppg_sizes.size() != a.path().mapping().size()) {
                    if (i >= mppg_sizes.size() - (mppg_sizes.size() - a.path().mapping().size())) {
                        probPaths.insert(-1);  // Insert an invalid ID to represent an empty path
                    } else {
                        nID = a.path().mapping()[i].position().node_id();
                        pangenome_base = dta->pangenome_map.at(std::to_string(nID));
                        mappability = max(1.0, dta->mappabilities[pangenome_base]);
                        for (const auto& path : nodepaths.at(nID - minid)) {
                            probPaths.insert(pathNameToID[path]);
                        }
                    }
                } else {
                    nID = a.path().mapping()[i].position().node_id();
                    pangenome_base = dta->pangenome_map.at(std::to_string(nID));
                    mappability = max(1.0, dta->mappabilities[pangenome_base]);
                    for (const auto& path : nodepaths.at(nID - minid)) {
                        probPaths.insert(pathNameToID[path]);
                    }
                }

                string nodeSeq;
                string partReadSeq;
                int startIndex = 0;
                if (a.path().mapping()[0].position().is_reverse()) {
                    startIndex = (baseIX - mppg_sizes.at(i) - 1 >= 0) ? (baseIX - mppg_sizes.at(i) - 1) : 0;
                    nodeSeq = graph_seq.substr(startIndex, mppg_sizes.at(i));
                    partReadSeq = read_seq.substr(startIndex, mppg_sizes.at(i));
                } else {
                    nodeSeq = graph_seq.substr(baseIX, mppg_sizes.at(i));
                    partReadSeq = read_seq.substr(baseIX, mppg_sizes.at(i));
                }

                std::vector<BaseInfo> readInfo;

                auto modifyPathNameInPrecompute = [&](
                    std::shared_ptr<Trailmix_struct>& dta,
                    std::string& path_name,
                    bool first
                ) {
                    std::string original_path_name = path_name;

                    // Special case: a single letter followed by underscore
                    std::regex pattern2("^([A-Za-z])_$");
                    path_name = std::regex_replace(path_name, pattern2, "$1");

                    // Handle the special path format
                    std::regex pattern3("\\+([0-9]+)\\+\\(([0-9]+)\\)");
                    path_name = std::regex_replace(path_name, pattern3, "_$1__$2_");

                    // Replace special characters
                    std::replace(path_name.begin(), path_name.end(), '+', '_');
                    std::replace(path_name.begin(), path_name.end(), '\'', '_');
                    std::replace(path_name.begin(), path_name.end(), '*', '_');
                    std::replace(path_name.begin(), path_name.end(), '@', '_');
                    std::replace(path_name.begin(), path_name.end(), '(', '_');
                    std::replace(path_name.begin(), path_name.end(), ')', '_');

                    if (!path_name.empty() && path_name.back() == '*') {
                        path_name.back() = '_';
                    }

                    if (path_name.size() == 1) {
                        path_name += "_";
                    }

                    // Remove .1 or .2 at the end
                    std::regex patternEnd("\\.(1|2)$");
                    path_name = std::regex_replace(path_name, patternEnd, "");

                    if (first) {
                        dta->originalPathNames[path_name] = original_path_name;
                    }
                };

                MCMC mcmc;

#pragma omp parallel for num_threads(dta->n_threads) private(readInfo)

                for (size_t m = 0; m < pathNames.size(); ++m) {
                    modifyPathNameInPrecompute(dta, pathNames[m], false);

                    set<int> depths_used;
                    auto p = dta->tree->nodes[dta->path_node_map[pathNames[m]]];
                    // bool in_pruned = mcmc.is_in_pruned(p, dta->depth, dta, depths_used);
                    bool in_pruned = mcmc.is_in_pruned_set(p, dta);
                    bool in_pruned_parent = true;
                    if (p->parent) {
                        in_pruned_parent = mcmc.is_in_pruned_set(p->parent, dta);
                    }

                    if (!in_pruned && !in_pruned_parent) {
                        continue;
                    }

                    // The path is supported
                    bool path_supported = probPaths.find(pathNameToID[pathNames[m]]) != probPaths.end();

                    if (path_supported) {
                        // cerr << "PATH SUPPORTED" << endl;

                        for (int s = 0; s < nodeSeq.size(); ++s) {
                            int base_quality = a.quality()[s];
                            if (nodeSeq[s] == 'N' || partReadSeq[s] == 'N') {

#ifdef DEBUGANALYSEGAM
                                get<4>(sancheck[pathNames[m]]) += 1;
#endif

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log((1-qscore_vec[base_quality]) / 3) + log(p_correctly_mapped) ;
                                info.logLikelihoodNoDamage = log((1-qscore_vec[base_quality]) / 3) + log(p_correctly_mapped) ;

                                if (info.logLikelihood >= 0.0) {
                                    throw std::runtime_error("info.loglikelihood is greater than 0 N");
                                }
                                readInfo.emplace_back(info);

                            } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S') {

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log((1 - qscore_vec[base_quality]) / 3) + log(p_correctly_mapped) ;
                                info.logLikelihoodNoDamage = log((1 - qscore_vec[base_quality]) / 3) + log(p_correctly_mapped) ;


                                if (info.logLikelihood >= 0.0) {
                                    throw std::runtime_error("info.loglikelihood is greater than 0 SSS");
                                }

                                readInfo.emplace_back(info);

#ifdef DEBUGANALYSEGAM
                                get<5>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);  // log((qscore_vec[base_quality])/3);
#endif

                            } else if (nodeSeq[s] == '-' || partReadSeq[s] == '-') {

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log(0.02);
                                info.logLikelihoodNoDamage = log(0.02);
                                readInfo.emplace_back(info);

#ifdef DEBUGANALYSEGAM
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);  // log((qscore_vec[base_quality])/3);
#endif
                            } else {
                                double probBasePreDamage[4];
                                for (int bpd = 0; bpd < 4; bpd++) {
                                    probBasePreDamage[bpd] = 1e-2;
                                }

                                // LEVEL 1 //
                                // Filling a substitution matrix with the probabilities of observing a base pre-damage.
                                for (int bpo = 0; bpo < 4; bpo++) {
                                    if ("ACGT"[bpo] == nodeSeq[s]) {  // No mutation
                                        probBasePreDamage[bpo] = 1 - qscore_vec[base_quality];
                                    } else {  // Mutation
                                        probBasePreDamage[bpo] = qscore_vec[base_quality] / 3;
                                    }
                                }


                                // LEVEL 2 //
                                // Initializing a post-damage substitution matrix. First filled with 0.
                                double probBasePostDamage[4];
                                double probBasePostDamage_none[4];

                                for (int bpd = 0; bpd < 4; bpd++) {
                                    probBasePostDamage[bpd] = 1e-6;
                                    probBasePostDamage_none[bpd] = 1e-6;
                                }

                                // Filling post-damage substitution matrix with the new probabilities of observing a base.
                                // The probabilities are based on the damage profile for --deam5p and --deam3p.
                                // If no damage profile is provided, the pre-damage substitution matrix is multiplied by 0 and stays the same.

                                for (int bpd = 0; bpd < 4; bpd++) {
                                    // Post-damage = prob pre-damage * damage rate from bpo to bpd
                                    for (int bpo = 0; bpo < 4; bpo++) {
                                        double probability = subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                        double probability_contaminant = dta->dmg_none.subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                        probBasePostDamage[bpd] += probBasePreDamage[bpo] * probability;
                                        probBasePostDamage_none[bpd] += probBasePreDamage[bpo] * probability_contaminant;
                                    }

                                   probBasePostDamage_none[bpd] = max(0.0001, probBasePostDamage_none[bpd]);
                                   probBasePostDamage_none[bpd] = min(1.0, probBasePostDamage_none[bpd]);
                                   probBasePostDamage[bpd] = max(0.0001, probBasePostDamage[bpd]);
                                   probBasePostDamage[bpd] = min(1.0, probBasePostDamage[bpd]);
                                }


#ifdef DEBUGDAMAGE
                                {
#pragma omp critical

if (partReadSeq[s] != nodeSeq[s]){
                                    cerr << "POST damage matrix " << endl;
                                    for (int bpd = 0; bpd < 4; bpd++) {
                                        cerr << setprecision(14) << bpd << "\t" << probBasePostDamage[bpd] << endl;
                                    }

                                    cerr << endl;

                                    cerr << "POST damage matrix NONE" << endl;
                                    for (int bpd = 0; bpd < 4; bpd++) {
                                        cerr << setprecision(14) << bpd << "\t" << probBasePostDamage_none[bpd] << endl;
                                    }

                                    cerr << endl << endl << endl;
                                }
                                }
#endif


auto nucleotideToIndex = [](char nucleotide) {
    switch (toupper(static_cast<unsigned char>(nucleotide))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 0;
    }
};
                                double log_lik_marg = -std::numeric_limits<double>::infinity();
                                double log_lik_marg_no_damage = -std::numeric_limits<double>::infinity();

                                double log_prob_damage = log(probBasePostDamage[nucleotideToIndex(partReadSeq[s])]);
                                double log_prob_no_damage = log(probBasePostDamage_none[nucleotideToIndex(partReadSeq[s])]);

auto logSumExp = [](double log_a, double log_b) {
    if (log_a == -std::numeric_limits<double>::infinity()) return log_b;
    if (log_b == -std::numeric_limits<double>::infinity()) return log_a;
    double max_val = std::max(log_a, log_b);
    return max_val + log(exp(log_a - max_val) + exp(log_b - max_val));
};


log_lik_marg_no_damage = log_prob_no_damage; //logSumExp(log_lik_marg_no_damage, log_prob_no_damage);
//if (partReadSeq[s] != nodeSeq[s]) {
    LOG_VALUE(log_lik_marg_no_damage, log_lik_marg_no_damage, a.name(), partReadSeq[s], nodeSeq[s], probBasePostDamage_none);
//}

log_lik_marg = log_prob_damage; //logSumExp(log_lik_marg, log_prob_damage);
//if (partReadSeq[s] != nodeSeq[s]) {
    LOG_VALUE(log_lik_marg, log_lik_marg, a.name(), partReadSeq[s], nodeSeq[s], probBasePostDamage);
//}


#ifdef DEBUGLLM
/*
    if (log_lik_marg != log_lik_marg_no_damage) {
        cerr << endl << endl;
        cerr << "LSEQ: " << Lseq << "  baseIX: " << baseIX << endl;
        cerr << std::setprecision(20) << "log_lik_marg:" << log_lik_marg << " p=" << exp(log_lik_marg) << endl;
        cerr << std::setprecision(20) << "log_lik_marg no damage " << log_lik_marg_no_damage << " p=" << exp(log_lik_marg_no_damage) << endl;

        // Log readBase and refBase
        cerr << "readBase: " << partReadSeq[s] << " refBase: " << nodeSeq[s] << endl;

        cerr << endl << endl;
    }
*/

#endif

                                // if (isnan(log_lik_marg)){throw runtime_error("LOG_LIK_MARG IS NAN");}

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;

                                info.logLikelihood = log_lik_marg + log(p_correctly_mapped) + log(mappability);
                                info.logLikelihoodNoDamage = log_lik_marg_no_damage + log(p_correctly_mapped) + log(mappability);

                                if (info.logLikelihood >= 1e-2) {
                                    cerr << "info.logLikelihood: " << info.logLikelihood << endl;
                                    throw std::runtime_error("info.loglikelihood is greater than 0 D");
                                }

                                readInfo.emplace_back(info);

                                if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                    baseOnRead--;
                                } else if (!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                    baseOnRead++;
                                }
                            }
                        }

                        // The path is not supported
                    } else if (!path_supported) {
                        unsigned int counter = 0;
                        for (size_t s = 0; s < nodeSeq.size(); ++s) {
                            int base_quality = a.quality()[s];

                            // int choice = range_choice(gen);

                            if (s % 4 != 0) {
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log(1 - (qscore_vec[base_quality])) + log(p_correctly_mapped);
                                if (isinf(info.logLikelihood)) {
                                    cerr << "mappability: " << mappability << endl;
                                    cerr << "p_correctly_mapped: " << p_correctly_mapped << endl;
                                    throw runtime_error("THIS IS INF");
                                }

                                info.logLikelihoodNoDamage = log(1 - (qscore_vec[base_quality])) + log(p_correctly_mapped);

                                if (isnan(info.logLikelihood) || info.logLikelihood > 1e-5) {
                                    throw runtime_error("UNSUPPORTED SHOULD NOT HAPPEN");
                                }

                                {
#pragma omp critical
                                    readInfo.emplace_back(move(info));
                                }

#ifdef DEBUGANALYSEGAM
                                get<6>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(1 - (qscore_vec[base_quality]));
#endif
                            } else {
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[90]));
                                info.logLikelihoodNoDamage = log((qscore_vec[90]));
                                {
#pragma omp critical
                                    readInfo.emplace_back(move(info));
                                }

#ifdef DEBUGANALYSEGAM
                                get<7>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);  // log((qscore_vec[base_quality])/3);
#endif
                            }

                            if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                baseOnRead--;
                            } else if (!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                baseOnRead++;
                            }
                        }

                    } else {
                        cerr << "Something is wrong!" << endl;
                    }

                    {
#pragma omp critical
                        // cerr << "ADDING INFO FOR: " << pathNames[m] << endl;
                        detailMap[pathNames[m]].emplace_back(move(readInfo));
                    }

                    if (a.path().mapping()[0].position().is_reverse()) {
                        baseOnRead = baseIX;
                    } else {
                        baseOnRead = baseIX;
                    }
                }

                if (a.path().mapping()[0].position().is_reverse()) {
                    baseIX = startIndex;
                    baseOnRead = baseIX;
                } else {
                    baseIX += mppg_sizes.at(i);
                    baseOnRead = baseIX;
                }
            }

#ifdef DEBUGANALYSEGAMMAP

            cerr << "Number of mappings " << a.path().mapping().size() << endl;
            cerr << "Reported mappings " << mppg_sizes.size() << endl;

            cerr << "Read Sequence " << a.sequence() << endl;
            cerr << "Sequence size " << a.sequence().size() << endl;

#endif

            vector<string> keysWithHighestValue;

#ifdef DEBUGANALYSEGAM
            std::map<std::string, std::pair<double, int>> sumsMap;

            for (const auto& kv : sancheck) {
                double logLikelihoodSum = std::get<8>(kv.second);  // Log Likelihood for this read
                int gapSum = std::get<3>(kv.second);               // Gap for this read

                if (sumsMap.find(kv.first) == sumsMap.end()) {
                    sumsMap[kv.first] = std::make_pair(logLikelihoodSum, gapSum);
                } else {
                    sumsMap[kv.first].first += logLikelihoodSum;
                    sumsMap[kv.first].second += gapSum;
                }
            }

            std::ofstream outfile("output2.txt", std::ios::app);

            // Writing Log Likelihood Sums
            outfile << "Log Likelihood Sums:\n";
            std::vector<std::pair<std::string, double>> likelihoodSums;
            for (const auto& kv : sumsMap) {
                likelihoodSums.emplace_back(kv.first, kv.second.first);
            }
            std::sort(likelihoodSums.begin(), likelihoodSums.end(), [](const auto& a, const auto& b) {
                return a.second > b.second;  // Sort in descending order
            });
            for (const auto& entry : likelihoodSums) {
                outfile << "Path: " << entry.first << "\tSum of Log Likelihood: " << entry.second << '\n';
            }

            // Writing Gap Sums
            outfile << "\nGap Sums:\n";
            std::vector<std::pair<std::string, int>> gapSums;
            for (const auto& kv : sumsMap) {
                gapSums.emplace_back(kv.first, kv.second.second);
            }
            std::sort(gapSums.begin(), gapSums.end(), [](const auto& a, const auto& b) {
                return a.second > b.second;  // Sort in descending order
            });
            for (const auto& entry : gapSums) {
                outfile << "Path: " << entry.first << "\tSum of Gaps: " << entry.second << '\n';
            }

            outfile.close();
#endif

            AlignmentInfo* ai = new AlignmentInfo();
            ai->seq = a.sequence();
            ai->path = a.path();
            ai->mapping_quality = a.mapping_quality();
            ai->quality_scores = a.quality();
            ai->is_paired = a.read_paired();
            ai->mostProbPath = keysWithHighestValue;
            print_dm = false;
            ai->detailMap = move(detailMap);

            read_vec->emplace_back(move(ai));
        }  // end of if identity statement
    }      // end of iteration through reads in gam file

    gam_file.close();

    // std::ofstream outFile("pruned_set.txt");
    // for (const auto& elem : dta->in_pruned_set) {
    //         outFile << elem << std::endl;
    // }
    // outFile.close();

    return read_vec;
}  // end of static function

static vector<AlignmentInfo*>* precompute_GAM_trailmix(shared_ptr<Trailmix_struct>& dta) {
    assert(!dta->path_names.empty());
    assert(!dta->qscore_vec.empty());
    vector<Clade*>* null_clade_vec = NULL;

    return precompute_GAM(dta->graph,
                          dta->fifo_A,
                          null_clade_vec,
                          dta->nodepaths,
                          dta->path_names,
                          dta->qscore_vec,
                          true,
                          dta->minid,
                          true,
                          dta->dmg.subDeamDiNuc,
                          dta);
}

