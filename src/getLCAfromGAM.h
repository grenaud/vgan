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

//#define DEBUGANC
//#define DEBUGANALYSEGAM
//#define DEBUGANALYSEGAMMAP
//#define DEBUGDAMAGE
//#define DEBUGPAIN2


using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

static vector<AlignmentInfo*>* analyse_GAM(
    const bdsg::ODGI& graph,
    string gamfilename,
    vector<Clade*>* clade_vec,
    const vector<NodeInfo*> nodevector,
    vector<vector<string>>& nodepaths,
    vector<string> pathNames,
    const vector<double> qscore_vec,
    bool trailmix, 
    int minid,
    int PATHTHRES,
    int PENALTY,
    bool getDetail,
    string outputname,
    bool singlesource, vector <vector<diNucleotideProb> > &subDeamDiNuc)
{
    int n_reads = 0;

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

    cerr << "GAM FILE NAME: " << gamfilename << endl;

#ifdef DEBUGANALYSEGAM
    cerr<<"analyse_GAM"<<endl;
#endif
    vector<AlignmentInfo*>* read_vec = new vector<AlignmentInfo*>();
    std::ifstream gam_file(gamfilename, std::ios::binary);
    vg::io::ProtobufIterator<vg::Alignment> iter(gam_file);
    std::ofstream detailtsv;
    if(getDetail){
        detailtsv.open(outputname +"_MatchInfo.tsv");
        detailtsv << "Path name\tRead name\tNode sequence\tIndex on node\tNode Base\tRead Sequence\tIndex on read\n";
    }
    // cutting path names string to be 100. 
    std::unordered_map<std::string, long double> results_map;
    int maxLength = 101;
    vector<string> pathNamesS;
    for (int g = 0; g<pathNames.size(); ++g){
        if (pathNames[g].length() > maxLength){
            pathNames[g].resize(maxLength);
            pathNamesS.push_back(pathNames[g]);

        }else{pathNamesS.push_back(pathNames[g]);}
    }
    if (pathNames.size() != pathNamesS.size()){throw runtime_error("Path computation went wrong!");}


    for (; iter.has_current(); iter.advance())
    {
        ++n_reads;
        const vg::Alignment& a = *iter;



        if (a.identity() != 0)
        {
         
            using BaseInfo = AlignmentInfo::BaseInfo;
            

            unordered_map<string, int> pathCounts;
            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());
            int baseIX = a.path().mapping()[0].position().is_reverse() ? a.sequence().size() - 1 : 0;
            int baseOnRead = baseIX;
            int Lseq = graph_seq.size();

            unordered_map<string, double> pathMap;
            unordered_map<string, bool> supportMap;
            unordered_map<string, vector<vector<BaseInfo>>> detailMap;

#ifdef DEBUGANALYSEGAM
            
            unordered_map<string, tuple<string, int, int, int, int, int, int, int, long double>> sancheck;
            for (int c = 0; c < pathNames.size(); ++c) {
                sancheck[pathNames[c]] = make_tuple(a.name(), 0, 0, 0, 0, 0, 0, 0, 0.0L);
            }
#endif

            for (int c = 0; c < pathNames.size(); ++c) {
                pathMap[pathNames[c]] = 0.0;
                supportMap[pathNames[c]] = false;
                detailMap[pathNames[c]] = vector<vector<BaseInfo>>();
            }
            
#ifdef DEBUGANALYSEGAMMAP
            cerr << "This is the mapping report of read " << a.name() << endl;
            cerr << "This is the sequence of it " << a.sequence() << endl;
            cerr << "This is the size of it " << a.sequence().size() << endl;
            cerr << "Reconstructed read sequence " << read_seq << " size " << read_seq.size() << endl;
            cerr << "Reconstructed graph sequence " << graph_seq << " size " << graph_seq.size() <<endl;
            cerr << "The size of the mappings is " << a.path().mapping().size() << endl;
            cerr << "The size of the mappings generated " << mppg_sizes.size() << endl;


            for (int j= 0; j < mppg_sizes.size(); ++j){
                cerr << mppg_sizes.at(j) << endl;
            }
#endif

            for (int i = 0; i < mppg_sizes.size(); ++i)
            {

#ifdef DEBUGANALYSEGAMMAP
                cerr << "Mappig size " << mppg_sizes.size() << endl;
                cerr << "iter " << i << endl;
#endif

                // extract all possible paths for a current node ID. There are discrptancies between the generated mapping lengths and the provided mapping length.
                //nodepaths has all path names for each node IDs listed.
                vector<string> probPaths;
                int nID = 0;
                if (mppg_sizes.size() != a.path().mapping().size())
                {
                    if (i >= mppg_sizes.size() - (mppg_sizes.size() - a.path().mapping().size())){
                        probPaths = {"No_support"};
                    }else{
                        nID = a.path().mapping()[i].position().node_id();
                        for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                            probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                        }
                    }
                    

                }else{
                    nID = a.path().mapping()[i].position().node_id();
                    for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                        probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                    }
                }

                

                string nodeSeq, partReadSeq;
                int startIndex = 0;
                if (a.path().mapping()[0].position().is_reverse()) {
                    startIndex = (baseIX - mppg_sizes.at(i) - 1 >= 0) ? (baseIX - mppg_sizes.at(i) - 1) : 0;
                    nodeSeq = graph_seq.substr(startIndex, mppg_sizes.at(i));
                    partReadSeq = read_seq.substr(startIndex, mppg_sizes.at(i));
                } else {
                    nodeSeq = graph_seq.substr(baseIX, mppg_sizes.at(i));
                    partReadSeq = read_seq.substr(baseIX, mppg_sizes.at(i));
                }
                int mismatch = 0;
                int unsup = 0;
                for (int m = 0; m < pathNames.size(); ++m)
                {
                    //cerr <<"path name " << pathNames[m] << endl;
                    std::vector<BaseInfo> readInfo;
                    // the path is supported
                    if (find(probPaths.begin(), probPaths.end(), pathNames[m]) != probPaths.end())
                    {
                        supportMap[pathNames[m]] = true;


                        for (int s = 0; s < nodeSeq.size(); ++s)
                        {
			    if(getDetail){
                                //cerr << pathNames.size() << " " << PATHTHRES << endl; 
			        if(nodeSeq[s] == partReadSeq[s] && probPaths.size() <= PATHTHRES){
                                    detailtsv 
 			                << pathNames[m] << '\t' 
        		                << a.name() << '\t' 
                                        << nodeSeq << '\t' 
        			        << s << '\t' 
        			        << nodeSeq[s] << '\t' 
         			        << a.sequence() << '\t' 
        			        << baseIX + s << '\n'; 

			        }
                            }

                            if (nodeSeq[s] != partReadSeq[s]){mismatch++;
#ifdef DEBUGANALYSEGAM
                                //get<2>(sancheck[pathNames[m]]) += 1;
                                //get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif
		            }
#ifdef DEBUGDAMAGE
                            cerr << "Base Index " << baseOnRead << endl;
                            cerr << "Graph Base " << nodeSeq[s] << endl;
                            cerr << "Read Base " << partReadSeq[s] << endl;
#endif
                            int base_quality = a.quality()[s];
                            if (nodeSeq[s] == 'N' || partReadSeq[s] == 'N')
                            {
                                pathMap[pathNames[m]] += log(0.25);

#ifdef DEBUGANALYSEGAM          
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                //info.logLikelihood = log((qscore_vec[base_quality]/3));
                                info.logLikelihood = log(0.25);
                                 if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("N info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "N info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);

                            } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S')
                            {

                                pathMap[pathNames[m]] += log((qscore_vec[base_quality])/3);
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[base_quality]/3));
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("Softclip info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "Softclip info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);
#ifdef DEBUGANALYSEGAM
                                get<5>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif

                            }else if (nodeSeq[s] == '-' || partReadSeq[s] == '-')
                            {

                                pathMap[pathNames[m]] += log(0.02);
                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                //info.logLikelihood = log((qscore_vec[base_quality])/3);
                                info.logLikelihood = log(0.02);
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("GAP info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "GAP info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
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
                                        probBasePreDamage[bpo]= (qscore_vec[base_quality]/3);
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
                                    probBasePostDamage[bpd]=0.0;
                                }


                                // filling post-damage substitution matrix with the new probabilities of observing a base. The probabilities are based on the damage profile for --deam5p and --deam3p

                                // If no damage profile is provided the pre-damage substitution matrix is multiplied by 0 and stays the same.
                                //cerr << "base index " << baseIX << " read size " << a.sequence().size() << endl;
                                for(int bpd=0;bpd<4;bpd++){
                                    //post damage = prob pre damage * damage rate from bpo to bpd
                                    for(int bpo=0;bpo<4;bpo++){
                                        probBasePostDamage[bpd]+=probBasePreDamage[bpo]*subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                        //cerr << setprecision(14) <<bpd<<"\t"<< subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                    }
                                }
#ifdef DEBUGDAMAGE 
                                cerr << "After damage matrix " << endl;
                                for(int bpd=0;bpd<4;bpd++){

                                    cerr<<setprecision(14)<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }
#endif
                                double log_lik_marg = -std::numeric_limits<double>::infinity(); //sum of the likelihoods but in log space, this will be added to log_lik
                                
                                for(int bpd=0;bpd<4;bpd++){
                                    if( "ACGT"[bpd]== partReadSeq[s]){// no sequencing error
                                    //                                                    Prob of bpd                no seq error
                                        log_lik_marg = oplusInitnatl( log_lik_marg , (log(probBasePostDamage[bpd]) ) );

                                    }else{ // we are already marginalising over all possible bases, if the base is not matching it must be a seq error.
                                    //                                                    Prob of bpd                seq error
                                        log_lik_marg = oplusInitnatl( log_lik_marg , (log(    probBasePostDamage[bpd] ))) ;
                                    }
                                }

                                if(log_lik_marg > log(0.9999999)){
                                    log_lik_marg = log(0.9999999);
                                }

                                if(isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
                                    throw runtime_error("calculated log like is nan");
                                }
                                

#ifdef DEBUGDAMAGE                                
                                cerr<<"postdamage:"<<endl;

                                for(int bpd=0;bpd<4;bpd++){
                                    cerr<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                }

                                cerr<<"log_lik_marg:"<<log_lik_marg<<" p="<<exp(log_lik_marg)<<endl;
                                //cout <<"log_lik_marg " <<  log_lik_marg << endl;
                                
                                cerr <<  pathMap[pathNames[m]] << endl;
#endif
                                pathMap[pathNames[m]] += log_lik_marg;

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log_lik_marg;
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("damage info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr<< std::setprecision(25) << "damage info.logLikelihood " << info.logLikelihood << " "<< log_lik_marg << endl;
#endif
                                readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                get<1>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(1 -(qscore_vec[base_quality])/3);
#endif



                                
                            
                                if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-' )
                                {
                                    baseOnRead--;
                                } else if(!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                    baseOnRead++;
                                }
                            }
                        }

                        // the path is not supported
                        
                    } else if (find(probPaths.begin(), probPaths.end(), pathNames[m]) == probPaths.end()) {
                        int counter = 0;
                        for (int s = 0; s < nodeSeq.size(); ++s)
                        {
                            unsup++;
                            int base_quality = a.quality()[s];
                            if (nodeSeq[s] == 'N' || partReadSeq[s] == 'N')
                            {

                                pathMap[pathNames[m]] += log(0.25);

#ifdef DEBUGANALYSEGAM          
                                get<4>(sancheck[pathNames[m]]) += 1;
#endif
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                //info.logLikelihood = log((qscore_vec[base_quality]/3));
                                info.logLikelihood = log(0.25);
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("N info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "N info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);

                            } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S')
                            {

                                pathMap[pathNames[m]] += log((qscore_vec[base_quality])/3);
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                info.logLikelihood = log((qscore_vec[base_quality]/3));
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("Softclip info.logLikelihood is nan ");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "Softclip info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);
#ifdef DEBUGANALYSEGAM
                                get<5>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif

                            }else if (nodeSeq[s] == '-' || partReadSeq[s] == '-')
                            {

                                pathMap[pathNames[m]] += log(0.02);
                                BaseInfo info;
                                info.readBase = '-';
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = false;
                                //info.logLikelihood = log((qscore_vec[base_quality])/3);
                                info.logLikelihood = log(0.02);
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("GAP info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr << std::setprecision(14)<< "GAP info.logLikelihood " << info.logLikelihood << endl;
#endif
                                readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                get<3>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif
                            
                            }else
                            {
                                if (abs(baseOnRead) % PENALTY  == 0)
                                {   
                                    pathMap[pathNames[m]] += log(1 - (qscore_vec[base_quality]));
                                    BaseInfo info;
                                    info.readBase = '-';
                                    info.referenceBase = nodeSeq[s];
                                    info.pathSupport = false;
                                    info.logLikelihood = log(1 - (qscore_vec[base_quality]));
                                    if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("no sup info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                    cerr << std::setprecision(14)<< "no sup info.logLikelihood " << info.logLikelihood << endl;
#endif
                                    readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                    get<6>(sancheck[pathNames[m]]) += 1;
                                    get<8>(sancheck[pathNames[m]]) += log(1 - (qscore_vec[base_quality]));
#endif
                                } else 
                                {

                                    pathMap[pathNames[m]] += log((qscore_vec[base_quality]/3));

                                    
                                    BaseInfo info;
                                    info.readBase = '-';
                                    info.referenceBase = nodeSeq[s];
                                    info.pathSupport = false;
                                    info.logLikelihood = log((qscore_vec[base_quality])/3);
                                    if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("no sup info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                    cerr << std::setprecision(14)<< "no sup info.logLikelihood " << info.logLikelihood << endl;
#endif
                                    readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                    get<7>(sancheck[pathNames[m]]) += 1;
                                    get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif
                                }
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
                    
                    detailMap[pathNamesS[m]].push_back(readInfo);

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



            long double highestValue = numeric_limits<long double>::lowest();
            for (const auto& pair : pathMap) {
                long double currentValue = pair.second;

                highestValue = max(highestValue, currentValue);
                //cout << highestValue << endl;
            }

            vector<string> keysWithHighestValue;
            for (const auto& pair : pathMap) {

                if (pair.second == highestValue) {
                    keysWithHighestValue.push_back(pair.first);
                    //cout << pair.first << endl;
                }
            }


            
            if(singlesource){

                for (const auto& pair : pathMap) {
                

                    results_map[pair.first] += pair.second;
                }

            }else{
                long double freq = 0.5;

                // Iterate over all distinct pairs in the map
                for (auto it1 = pathMap.begin(); it1 != pathMap.end(); ++it1) {
                    auto it2 = it1;
                    ++it2;
                    for (; it2 != pathMap.end(); ++it2) {
                        std::string new_key = it1->first + "-" + it2->first; // forming a new key
                        long double result = oplusnatl((log(freq) + it1->second), (log(1-freq) + it2->second));
                        results_map[new_key] += result; // storing the result in the new map
                    }
                }
            }




#ifdef DEBUGANALYSEGAM
            std::ofstream outfile("output2.txt", std::ios::app);

            for(const auto& kv : sancheck) {
                outfile << "Path: " << kv.first << '\t';
                outfile << "Read_name: " << std::get<0>(kv.second) << '\t';
                outfile << "Matches: " << std::get<1>(kv.second) << '\t';
                outfile << "Mismatches: " << std::get<2>(kv.second) << '\t';
                outfile << "Gap: " << std::get<3>(kv.second) << '\t';
                outfile << "N: " << std::get<4>(kv.second) << '\t';
                outfile << "S: " << std::get<5>(kv.second) << '\t';
                outfile << "Matches_missing_node: " << std::get<6>(kv.second) << '\t';
                outfile << "Mismatches_missing_node: " << std::get<7>(kv.second) << '\t';
                outfile << "Sum_of_bases: "<< std::get<1>(kv.second)+std::get<2>(kv.second)+std::get<3>(kv.second)+std::get<4>(kv.second)+std::get<5>(kv.second)+std::get<6>(kv.second)+std::get<7>(kv.second)<< '\t';
                outfile << "Log_likelihood: " << std::get<8>(kv.second) << '\n';
                
            }

            outfile.close();
#endif


            AlignmentInfo* ai = new AlignmentInfo();
            ai->seq = a.sequence();
            ai->name = a.name();
            ai->path = a.path();
            ai->mapping_quality = a.mapping_quality();
            ai->quality_scores = a.quality();
            ai->is_paired = a.read_paired();
            ai->mostProbPath = keysWithHighestValue;
            ai->pathMap = pathMap;
            ai->supportMap = supportMap;
            ai->detailMap = detailMap;

            read_vec->push_back(ai);


            // std::string specificKey = "NC_062361.1_Hippotragus_niger_roosevelti_voucher_ZMUC_H.R.Siegismund_1646_haplogroup_Eastern_1_mitochondrion__complete_genome";
            // auto detailMapIt = detailMap.find(specificKey);

            // if (detailMapIt != detailMap.end()) {
            //     const auto& detailVectors = detailMapIt->second; // The vector of vectors of BaseInfo

            //     std::cout << "Sizes of vectors for key '" << specificKey << "':" << std::endl;

            //     // Iterate through the vector of vectors to get their sizes
            //     for (size_t i = 0; i < detailVectors.size(); ++i) {
            //         std::cout << "Size of vector " << i << ": " << detailVectors[i].size() << std::endl;
            //     }
            // } else {
            //     std::cout << "Key '" << specificKey << "' not found in detailMap." << std::endl;
            // }



           // Assuming SOME_TOLERANCE is defined
            const double SOME_TOLERANCE = 1e-6;

            // Iterate through pathMap and compare values with those in detailMap
            for (const auto& pathPair : pathMap) {
                const auto& pathName = pathPair.first;
                double pathMapValue = pathPair.second; // The accumulated log-likelihood from pathMap
                //cerr << "pathName " << pathName << " value " << pathMapValue << endl;
 
                // Check if this path exists in detailMap
                auto detailMapIt = detailMap.find(pathName);
                if (detailMapIt != detailMap.end()) {
                    const auto& detailVectors = detailMapIt->second; // The vector of vectors of BaseInfo
                    double detailMapValue = 0.0; // Initialize the sum for this path

                    // Iterate through the vector of vectors of BaseInfo to sum logLikelihoods
                    for (const auto& baseInfoVec : detailVectors) {
                        for (const auto& baseInfo : baseInfoVec) {
                            detailMapValue += baseInfo.logLikelihood; // Sum logLikelihoods
                        }
                    }

                    // Now compare pathMapValue and detailMapValue
                    if (std::abs(pathMapValue - detailMapValue) > SOME_TOLERANCE) {
                        std::cerr << "Difference detected in path: " << pathName << std::endl;
                        std::cerr << "pathMap value: " << pathMapValue << ", detailMap sum value: " << detailMapValue << std::endl;
                        throw runtime_error("Something is going wrong in the precomputation the pathMap values and the detailMap Value are not the same.");
                    }
                } else {
                    std::cerr << "Path: " << pathName << " is not found in detailMap" << std::endl;
                }
            }


        } // end of if identity statement

    } // end of iteration through reads in gam file
    detailtsv.close();
    gam_file.close();

     for (const auto& pair : results_map) {
            //std::cout << pair.first << " " << pair.second << '\n';
        }

    auto max_it = max_element(results_map.begin(), results_map.end(),
                              [](const auto& a, const auto& b) { return a.second < b.second; });

    // Now to find the "index"
    int idx = 0;
    for (auto it = results_map.begin(); it != results_map.end(); ++it, ++idx) {
        if (it == max_it) {
            break;
        }
    }
#ifdef DEBUGANALYSEGAM
    cerr << "Path with maximum cumulative value: " << max_it->first 
         << ", Value: " << max_it->second 
         << ", Index: " << idx << endl;
#endif








    return read_vec;
} // end of static function 

// static vector<AlignmentInfo*>* analyse_GAM_trailmix(shared_ptr<Trailmix_struct> &dta){

//     vector<Clade*>* null_clade_vec = new vector<Clade*>;
//     vector<vector<string>> null_node_paths;
//     assert(!dta->path_names.empty());

//     return analyse_GAM(dta->graph, \
//                        dta->fifo_A, \
//                        null_clade_vec, \
//                        dta->nodevector, \
//                        null_node_paths, \
//                        dta->path_supports, \
//                        dta->MCMC_path_names, \
//                        dta->qscore_vec, \
//                        1e-8, \
//                        false, \
//                        true, 
//                        2
//                        );

// }
