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

using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

static vector<AlignmentInfo*>* precompute_GAM(
    const bdsg::ODGI& graph,
    string gamfilename,
    vector<Clade*>* clade_vec,
    const vector<NodeInfo*> &nodevector,
    const vector<vector<string>>& nodepaths,
    vector<string> &pathNames,
    const vector<double> &qscore_vec,
    bool trailmix,
    int minid,
    bool singlesource,
    vector <vector<diNucleotideProb> > &subDeamDiNuc,
    shared_ptr<Trailmix_struct> dta = NULL
                                          )
{
    int n_reads = 0;
    bool print_dm=true;

    if (gamfilename.empty()) {
        throw std::runtime_error(gamfilename + " is empty. Aborting.");
    }
    if (nodevector.empty()) {
        throw std::runtime_error("Node vector is empty. Aborting.");
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

    for (; iter.has_current(); iter.advance())
    {
        ++n_reads;
           if (n_reads % 10 == 0){
            cerr << "On read: " << (int)(n_reads/2) << endl;
                                 }
        const vg::Alignment& a = *iter;
        if(a.identity() == 0){continue;}

        double p_correctly_mapped = 1-dta->incorrect_mapping_vec[a.mapping_quality()];
        p_correctly_mapped = max(p_correctly_mapped, 1e-9);

        if (a.identity() != 0)
        {
            using BaseInfo = AlignmentInfo::BaseInfo;
            auto&& [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());
            //PRINTVEC(mppg_sizes)
            //string graph_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
            //string read_seq =  "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
            //vector<int> mppg_sizes(8, 1);

            int baseIX = a.path().mapping()[0].position().is_reverse() ? a.sequence().size() - 1 : 0;
            int baseOnRead = baseIX;
            size_t Lseq = graph_seq.size();

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

            for (int i = 0; i < mppg_sizes.size(); ++i)
            {

#ifdef DEBUGANALYSEGAMMAP
                cerr << "mapping size " << mppg_sizes.size() << endl;
                cerr << "iter " << i << endl;
#endif


std::unordered_set<std::string> probPaths;
int nID = 0;
int pangenome_base = 100;
double mappability = 1.0;

if (mppg_sizes.size() != a.path().mapping().size()) {
    if (i >= mppg_sizes.size() - (mppg_sizes.size() - a.path().mapping().size())) {
        probPaths.insert("");
    } else {
        nID = a.path().mapping()[i].position().node_id();
        pangenome_base = dta->pangenome_map.at(to_string(nID));
        mappability = dta->mappabilities[pangenome_base];
        for (const auto& path : nodepaths.at(nID - minid)) {
            probPaths.insert(path);
        }
    }
} else {
    nID = a.path().mapping()[i].position().node_id();
      pangenome_base = dta->pangenome_map.at(to_string(nID));
      mappability = dta->mappabilities[pangenome_base];
    for (const auto& path : nodepaths.at(nID - minid)) {
        probPaths.insert(path);
    }
}


if (mappability < 0.1){mappability=0.1;}

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
    std::shared_ptr<Trailmix_struct> &dta, 
    std::string &path_name, 
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
bool in_pruned = mcmc.is_in_pruned(p, dta->depth, dta, depths_used);

if(!in_pruned){
continue;
}

                    // the path is supported
                    bool path_supported = probPaths.find(pathNames[m]) != probPaths.end();

                    if (path_supported)
                    {

                        //cerr << "PATH SUPPORTED" << endl;
                        
                        for (int s = 0; s < nodeSeq.size(); ++s)
                        {
                            int base_quality = a.quality()[s];
                            if (nodeSeq[s] == 'N' || partReadSeq[s] == 'N')
                            {

#ifdef DEBUGANALYSEGAM          
                                get<4>(sancheck[pathNames[m]]) += 1;
#endif

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);}
                                 if (info.logLikelihood >= 0.0) {
                                        throw std::runtime_error("info.loglikelihood is greater than 0 N");
                                                                }
                                readInfo.emplace_back(info);

                            } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S')
                            {

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);}
                                 if (info.logLikelihood >= 0.0) {
                                        throw std::runtime_error("info.loglikelihood is greater than 0 SSS");
                                                                }

                                readInfo.emplace_back(info);
#ifdef DEBUGANALYSEGAM
                                get<5>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);//log((qscore_vec[base_quality])/3);
#endif

                            }else if (nodeSeq[s] == '-' || partReadSeq[s] == '-')
                            {

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(qscore_vec[base_quality]/3) + log(p_correctly_mapped) + log(mappability);}
                                readInfo.emplace_back(info);

#ifdef DEBUGANALYSEGAM
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);//log((qscore_vec[base_quality])/3);
#endif
                            } else
                            {
                                double probBasePreDamage [4];
                                for(int bpd=0;bpd<4;bpd++){
                                    probBasePreDamage[bpd]=1e-6;
                                }

                                // LEVEL 1 //
                                // filling a substitution matrix with the propbabilities of observing a base pre damage.
                                for(int bpo=0;bpo<4;bpo++){
                                    if("ACGT"[bpo] == nodeSeq[s]){// no mutation
                                        probBasePreDamage[bpo]= 1 - qscore_vec[base_quality];
                                    }else{ // mutation
                                        probBasePreDamage[bpo]= qscore_vec[base_quality]/3;
                                    }

                                }
#ifdef DEBUGDAMAGE 
                                cerr << "Pre damage matrix" << endl;
                                for(int bpo=0;bpo<4;bpo++){
                                    cerr<<setprecision(14)<<bpo<<"\t"<<probBasePreDamage[bpo]<<endl;
                                }
#endif

                                // LEVEL 2 //
                                //initializing a post damage substitution matrix. First filled with 0.
                                double probBasePostDamage [4];

                                for(int bpd=0;bpd<4;bpd++){
                                    probBasePostDamage[bpd]=1e-10;
                                }

                                // filling post-damage substitution matrix with the new probabilities of observing a base. The probabilities are based on the damage profile for --deam5p and --deam3p
                                // If no damage profile is provided the pre-damage substitution matrix is multiplied by 0 and stays the same.

                                for(int bpd=0;bpd<4;bpd++){
                                    //post damage = prob pre damage * damage rate from bpo to bpd
                                    for(int bpo=0;bpo<4;bpo++){
                                        double probability = max(0.0001, subDeamDiNuc[Lseq][baseIX].p[bpo][bpd]);
                                        probability = min(0.99999, probability);
                                        if (probability < -1e-8 || probability > 1.0 + 1e-8) {
                                            cerr << "Lseq: " << Lseq << endl;
                                            cerr << "baseIX: " << baseIX << endl;
                                            cerr << "bpo: " << bpo << endl;
                                            cerr << "bpd: " << bpd << endl;

                                            throw std::runtime_error("Invalid probability value detected: " + std::to_string(probability) +
                                                                     ". Probability must be between 0 and 1.");
                                                                                             }
probBasePostDamage[bpd] += probBasePreDamage[bpo] * probability;


                                    }
                                }



#ifdef DEBUGDAMAGE
{
#pragma omp critical
                                cerr << "After damage matrix " << endl;
                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<setprecision(14)<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }

                                cerr << "After damage matrix no damage" << endl;
                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<setprecision(14)<<bpd<<"\t"<<probBaseNoDamage[bpd]<<endl;
                                }

}
#endif


                                double log_lik_marg = -std::numeric_limits<double>::infinity();
                                double log_lik_marg_no_damage = -std::numeric_limits<double>::infinity();

                                for(int bpd=0; bpd<4; bpd++){

                                    log_lik_marg = oplusInitnatl(log_lik_marg, log(max(0.000001,probBasePostDamage[bpd])));
                                    log_lik_marg_no_damage = oplusInitnatl(log_lik_marg_no_damage, max(0.000001, log(probBasePreDamage[bpd])));
}


if (log_lik_marg > log(0.999999)){log_lik_marg = log(0.999999);}
if (log_lik_marg_no_damage > log(0.999999)){log_lik_marg_no_damage = log(0.999999);}
if (log_lik_marg < log(0.000001)){log_lik_marg = log(0.000001);}
if (log_lik_marg_no_damage < log(0.000001)){log_lik_marg_no_damage = log(0.000001);}

#ifdef DEBUGDAMAGE                                
                                cerr<<"postdamage:"<<endl;

                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }

                                cerr<<"log_lik_marg:"<<log_lik_marg<<" p="<<exp(log_lik_marg)<<endl;
                                cerr <<"log_lik_marg no damage " <<  log_lik_marg_no_damage << " p="<<exp(log_lik_marg_no_damage) << endl;
                                

#endif
                                //if (isnan(log_lik_marg)){throw runtime_error("LOG_LIK_MARG IS NAN");}


                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;

                                info.logLikelihood = log_lik_marg + log(p_correctly_mapped) + log(mappability);

                                if (dta->cont_mode){info.logLikelihoodNoDamage = log_lik_marg_no_damage  + log(p_correctly_mapped) + log(mappability);}
                                if (info.logLikelihood >= 1e-8) {
                                        cerr << "info.logLikelihood: " << info.logLikelihood << endl;
                                        throw std::runtime_error("info.loglikelihood is greater than 0 D");
                                                                  }

                                readInfo.emplace_back(info);

                                if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-')
                                {
                                    baseOnRead--;
                                } else if(!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                    baseOnRead++;
                                }
                            }
                        }

                        // the path is not supported
                    } else if (probPaths.find(pathNames[m]) == probPaths.end()) {
                        unsigned int counter = 0;
                        for (size_t s = 0; s < nodeSeq.size(); ++s)
                        {
                            int base_quality = a.quality()[s];

                           // ELEPHANT
                           if (s % 2 != 0)
                            {
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log(1 - (qscore_vec[base_quality])) + log(p_correctly_mapped) + log(mappability);
                                if(isinf(info.logLikelihood)){
                                    cerr << "mappability: " << mappability << endl;
                                    cerr << "p_correctly_mapped: " << p_correctly_mapped << endl;
                                    throw runtime_error("THIS IS INF");
                                                             }
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(1 - (qscore_vec[base_quality])) + log(p_correctly_mapped) + log(mappability);}
                                if (isnan(info.logLikelihood) || info.logLikelihood > 1e-5){
                                   throw runtime_error("UNSUPPORTED SHOULD NOT HAPPEN");
                                }

{
                                #pragma omp critical
                               readInfo.emplace_back(move(info));
}

#ifdef DEBUGANALYSEGAM
                                get<6>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(1 - (qscore_vec[base_quality])) + log(p_correctly_mapped) + log(mappability);
#endif
                            } else
                            {
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[base_quality])) + log(p_correctly_mapped) + log(mappability);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log((qscore_vec[base_quality])) + log(p_correctly_mapped) + log(mappability);}
{
                                #pragma omp critical
                                readInfo.emplace_back(move(info));
}

#ifdef DEBUGANALYSEGAM
                                get<7>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5); //log((qscore_vec[base_quality])/3);
#endif
                            }

                            if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                baseOnRead--;
                            } else if(!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-'){
                                baseOnRead++;
                            }


                        }


                    } else {
                        cerr << "Something is wrong!" << endl;

                    }

                   {
                    #pragma omp critical
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
    double logLikelihoodSum = std::get<8>(kv.second); // Log Likelihood for this read
    int gapSum = std::get<3>(kv.second); // Gap for this read

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
    return a.second > b.second; // Sort in descending order
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
    return a.second > b.second; // Sort in descending order
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
        } // end of if identity statement
    } // end of iteration through reads in gam file

    gam_file.close();

    return read_vec;
} // end of static function

 static vector<AlignmentInfo*>* precompute_GAM_trailmix(shared_ptr<Trailmix_struct> &dta){

     assert(!dta->path_names.empty());
     assert(!dta->nodevector.empty());
     assert(!dta->qscore_vec.empty());
     vector<Clade *> * null_clade_vec = NULL;

     return precompute_GAM(dta->graph, \
                        dta->fifo_A, \
                        null_clade_vec, \
                        dta->nodevector, \
                        dta->nodepaths, \
                        dta->path_names, \
                        dta->qscore_vec, \
                        true, \
                        dta->minid, \
                        true, \
                        dta->dmg.subDeamDiNuc,
                        dta
                        );

 }
