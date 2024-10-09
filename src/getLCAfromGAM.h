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
#include "vg/io/alignment_emitter.hpp"
#include "utility.hpp"
#include "alignment.hpp"
#include "AlignmentInfo.h"
#include "libgab.h"
#include "vgan_utils.h"

//#define DEBUGANALYSEGAM

//#define DEBUGANC
//#define DEBUGANALYSEGAM
//#define DEBUGANALYSEGAMMAP
//#define DEBUGDAMAGE

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

//#ifdef DEBUGANALYSEGAM
    cerr<<"In analyse_GAM"<<endl;
//#endif
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
            cerr << "On read: " << n_reads << endl;
                                 }
        const vg::Alignment& a = *iter;

        double p_correctly_mapped = (1-dta->incorrect_mapping_vec[a.mapping_quality()]);


        if (a.identity() != 0)
        {
            using BaseInfo = AlignmentInfo::BaseInfo;

            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());

auto calculateSubstitutionMatrix = [](const std::string& read_seq, const std::string& graph_seq, int N) -> Matrix {
    std::unordered_map<char, int> nucleotide_index = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    Matrix matrix = {}; // Initialize all elements to 0

    // Helper function to update the matrix
    auto updateMatrix = [&](size_t i) {
        char read_char = read_seq[i];
        char graph_char = graph_seq[i];

        if (nucleotide_index.find(read_char) != nucleotide_index.end() &&
            nucleotide_index.find(graph_char) != nucleotide_index.end()) {
            int read_index = nucleotide_index[read_char];
            int graph_index = nucleotide_index[graph_char];
            matrix[graph_index][read_index] += 1.0;
        }
    };

    // Process all bases if N is -1
    if (N == -1) {
        for (size_t i = 0; i < read_seq.length() && i < graph_seq.length(); ++i) {
            updateMatrix(i);
        }
    }
    // Process first N bases
    else {
        for (size_t i = 0; i < N && i < read_seq.length() && i < graph_seq.length(); ++i) {
            updateMatrix(i);
        }

        // Process last N bases
        size_t len = read_seq.length();
        for (size_t i = 0; i < N && i < len && (len - i - 1) < graph_seq.length(); ++i) {
            updateMatrix(len - i - 1);
        }
    }

    return matrix;
};

            auto matrix = calculateSubstitutionMatrix(read_seq, graph_seq, -1);
            dta->substitution_matrices.emplace_back(matrix);

            int baseIX = a.path().mapping()[0].position().is_reverse() ? a.sequence().size() - 1 : 0;
            int baseOnRead = baseIX;
            size_t Lseq = graph_seq.size();

            unordered_map<string, double> pathMap;
            unordered_map<string, vector<vector<BaseInfo>>> detailMap;
                    // Initialize detailMap for each path in pathNames with uniform size vectors

#ifdef DEBUGANALYSEGAM
            unordered_map<string, tuple<string, int, int, int, int, int, int, int, double>> sancheck;
            for (int c = 0; c < pathNames.size(); ++c) {
                sancheck[pathNames[c]] = make_tuple(a.name(), 0, 0, 0, 0, 0, 0, 0, 0.0L);
            }
#endif

            for (size_t c = 0; c < pathNames.size(); ++c) {
                pathMap[pathNames[c]] = 0.0;
                detailMap[pathNames[c]] = vector<vector<BaseInfo>>();
            }

//cerr << "DETAIL MAP SIZE: " << detailMap.size() << endl;

if (detailMap.empty()){throw runtime_error("[MCMC] Detail map is empty");}
            
#ifdef DEBUGANALYSEGAMMAP
            //cerr << "This is the mapping report of read " << a.name() << endl;
            //cerr << "This is the sequence of it " << a.sequence() << endl;
            //cerr << "This is the size of it " << a.sequence().size() << endl;
            cerr << "Reconstructed read sequence " << read_seq << " size " << read_seq.size() << endl;
            cerr << "Reconstructed graph sequence " << graph_seq << " size " << graph_seq.size() <<endl;
            //cerr << "The size of the mappings is " << a.path().mapping().size() << endl;
            //cerr << "The size of the mappings generated " << mppg_sizes.size() << endl;


            //for (int j= 0; j < mppg_sizes.size(); ++j){
            //    cerr << mppg_sizes.at(j) << endl;
           // }
#endif

            for (int i = 0; i < mppg_sizes.size(); ++i)
            {

#ifdef DEBUGANALYSEGAMMAP
                cerr << "mapping size " << mppg_sizes.size() << endl;
                cerr << "iter " << i << endl;
#endif


                                std::unordered_set<std::string> probPaths;
int nID = 0;

if (mppg_sizes.size() != a.path().mapping().size()) {
    if (i >= mppg_sizes.size() - (mppg_sizes.size() - a.path().mapping().size())) {
        probPaths.insert("");
    } else {
        nID = a.path().mapping()[i].position().node_id();
        for (const auto& path : nodepaths.at(nID - minid)) {
            probPaths.insert(path);
        }
    }
} else {
    nID = a.path().mapping()[i].position().node_id();
    for (const auto& path : nodepaths.at(nID - minid)) {
        probPaths.insert(path);
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
                //#pragma omp parallel for num_threads(5) private(readInfo)

for (size_t m = 0; m < pathNames.size(); ++m) {
                    readInfo.clear();

                    // the path is supported
                    bool path_supported = probPaths.find(pathNames[m]) != probPaths.end();
                    //const bool path_supported = (find(probPaths.begin(), probPaths.end(), pathNames[m]) != probPaths.end());

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
                                pathMap[pathNames[m]] += log((qscore_vec[base_quality]/3));
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[base_quality]/3) * p_correctly_mapped);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log((qscore_vec[base_quality]/3));}
                                 if (info.logLikelihood >= 0.0) {
                                        throw std::runtime_error("info.loglikelihood is greater than 0 N");
                                                                }
                                readInfo.emplace_back(info);

                            } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S')
                            {

                                pathMap[pathNames[m]] += log((qscore_vec[base_quality])/3);
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[base_quality]/3) * p_correctly_mapped);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log((qscore_vec[base_quality]/3));}
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
                                pathMap[pathNames[m]] += log((qscore_vec[base_quality])/3);
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood =  log(p_correctly_mapped * (qscore_vec[base_quality])/3);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log((qscore_vec[base_quality])/3);}
                                readInfo.emplace_back(info);

#ifdef DEBUGANALYSEGAM
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);//log((qscore_vec[base_quality])/3);
#endif
                            } else
                            {
                                double probBasePreDamage [4];
                                for(int bpd=0;bpd<4;bpd++){
                                    probBasePreDamage[bpd]=0.0;
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
                                cerr << "After damage matrix " << endl;
                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<setprecision(14)<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }
#endif


                                double log_lik_marg = -std::numeric_limits<double>::infinity(); //sum of the likelihoods but in log space, this will be added to log_lik
                                double log_lik_marg_no_damage = -std::numeric_limits<double>::infinity();
                                
                                for(int bpd=0; bpd<4; bpd++){
    
    if("ACGT"[bpd] == partReadSeq[s]){ // no sequencing error
        
        log_lik_marg = oplusInitnatl(log_lik_marg, log(probBasePostDamage[bpd]));
        log_lik_marg_no_damage = oplusInitnatl(log_lik_marg, max(0.0001, log(probBasePreDamage[bpd])));

    } else { // sequencing error
        
        log_lik_marg = oplusInitnatl(log_lik_marg, log(probBasePostDamage[bpd]));
        log_lik_marg_no_damage = oplusInitnatl(log_lik_marg, max(0.0001, log(probBasePreDamage[bpd])));
    }
}

    //if(log_lik_marg > log(0.9999999)){
    //      log_lik_marg = log(0.9999999);
    //                                 }

    //if(isnan(log_lik_marg) || isinf(log_lik_marg) || isnan(log_lik_marg_no_damage) || isinf(log_lik_marg_no_damage) || log_lik_marg > 1e-8 || log_lik_marg_no_damage > 1e-8){
    //      throw runtime_error("calculated log like is nan");
    //                                                                     }
                                


if (log_lik_marg > log(0.999)){log_lik_marg = log(0.999);}


#ifdef DEBUGDAMAGE                                
                                cerr<<"postdamage:"<<endl;

                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }

                                cerr<<"log_lik_marg:"<<log_lik_marg<<" p="<<exp(log_lik_marg)<<endl;
                                //cout <<"log_lik_marg " <<  log_lik_marg << endl;
                                
                                cerr <<  pathMap[pathNames[m]] << endl;
#endif
                                //if (isnan(log_lik_marg)){throw runtime_error("LOG_LIK_MARG IS NAN");}


                                pathMap[pathNames[m]] += log_lik_marg + log(p_correctly_mapped);

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log_lik_marg + log(p_correctly_mapped);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log_lik_marg_no_damage;}
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
                      //} else if (find(probPaths.begin(), probPaths.end(), pathNames[m]) == probPaths.end()) {
                    } else if (probPaths.find(pathNames[m]) == probPaths.end()) {
                        //cerr << "PATH UNSUPPORTED" << endl;
                        unsigned int counter = 0;
                        for (size_t s = 0; s < nodeSeq.size(); ++s)
                        {
                            int base_quality = a.quality()[s];

                            if (abs(baseOnRead) % 6 == 0)
                            {   
                                pathMap[pathNames[m]] += log((1-(qscore_vec[base_quality]/3))*p_correctly_mapped);
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((1-(qscore_vec[base_quality]/3)) * p_correctly_mapped);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(1-(qscore_vec[base_quality]/3));}
                                if (isnan(info.logLikelihood) || info.logLikelihood > 1e-5){
                                   throw runtime_error("UNSUPPORTED SHOULD NOT HAPPEN");
                                }

                                readInfo.emplace_back(info);

#ifdef DEBUGANALYSEGAM
                                get<6>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(0.5);//log(1 - (qscore_vec[base_quality]));
#endif
                            } else 
                            {
                                pathMap[pathNames[m]] += log(((qscore_vec[base_quality]/3))*p_correctly_mapped);
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log(((qscore_vec[base_quality]/3))*p_correctly_mapped);
                                if (dta->cont_mode){info.logLikelihoodNoDamage = log(0.1);}
                                readInfo.emplace_back(info);

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

                    detailMap[pathNames[m]].emplace_back(readInfo);

                    if (a.path().mapping()[0].position().is_reverse()) {
                        baseOnRead = baseIX;
                    } else {
                        baseOnRead = baseIX;
                    }

                    
                }
                //cout << readInfo.size() << endl;

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
if (!trailmix){
            double highestValue = numeric_limits<double>::lowest();
            for (const auto& pair : pathMap) {
                double currentValue = pair.second;

                highestValue = max(highestValue, currentValue);
                //cout << highestValue << endl;
            }

            for (const auto& pair : pathMap) {

                if (pair.second == highestValue) {
                    keysWithHighestValue.emplace_back(pair.first);
                    //cout << pair.first << endl;
                }
            }
              }


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
            ai->pathMap = move(pathMap);
            if (print_dm){
                //PRINT_DETAIL_MAP_SIZES(detailMap, "detailMap.txt");
                         }
            print_dm = false;
            ai->detailMap = move(detailMap);

            read_vec->emplace_back(ai);
        } // end of if identity statement

    } // end of iteration through reads in gam file

/*(
    // Step 1: Create a map to accumulate the sums
    std::map<std::string, double> accumulatedMap;

    // Step 2: Iterate over each read_vec entry and accumulate the values
    for (const auto& read : *read_vec) {
        for (const auto& [key, value] : read->pathMap) {
            accumulatedMap[key] += value;
        }
    }

    // Step 3: Convert the accumulated map to a vector of pairs
    std::vector<std::pair<std::string, double>> vec(accumulatedMap.begin(), accumulatedMap.end());

    // Step 4: Sort the vector based on the values
    std::sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    // Step 5: Write the sorted values to the debug.txt file
    std::ofstream outFile(random_string(7) + "_debug.txt");
    for (const auto& [key, value] : vec) {
        outFile << setprecision(16) << key << ": " << value << std::endl;
    }
*/
// Preparing to write detailMap information to detailmap.txt
std::vector<std::pair<std::string, double>> detailMapVec;

// Extracting logLikelihood values from detailMap
for (const auto& read : *read_vec) {
    for (const auto& [path, detailVector] : read->detailMap) {
        if (!detailVector.empty() && !detailVector[0].empty()) {
            detailMapVec.emplace_back(path, detailVector[0][0].logLikelihood);
        }
    }
}

// Sorting detailMapVec based on logLikelihood in descending order
std::sort(detailMapVec.begin(), detailMapVec.end(), [](const auto& a, const auto& b) {
    return a.second > b.second;
});

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
