#pragma once

#include "soibean.h"
#include "bdsg/hash_graph.hpp"
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

// =====================================================================================================================
// SOIBEAN FIXES -- ported from the TrailMix fixes (precompute.h / MCMC.h)
//
// soibean reads alignments with its own copy of the TrailMix logic (this function), and that copy had the same bugs.
// Every change is marked [FIX n]. The old lines are kept as comments ("// OLD: ...") right above the corrected code.
//
// [FIX 1] Which graph node does each piece of the alignment belong to?
//     reconstruct_graph_sequence() returns one entry in mppg_sizes per EDIT, not per node mapping. A mismatch or an
//     indel splits a node match into several edits, so mppg_sizes usually has MORE entries than the read has mappings.
//     The old code used a.path().mapping()[i] for entry i. That is the wrong node as soon as a read has a mismatch, and
//     for all trailing entries it gave up ("No_support": every path counted as unsupported). Now entry_to_mapping[i]
//     gives the mapping that owns the i-th edit, so the node (and therefore the supporting paths) is right.
//
// [FIX 2] Slicing graph_seq / read_seq.
//     Both strings are laid out in alignment-path order, which is read order (5'->3') for forward AND reverse
//     alignments. The old code started reverse reads at the 3' end and walked backwards (startIndex = baseIX - size - 1),
//     which cuts the wrong columns. Now each chunk is simply the next mppg_sizes[i] columns (colStart).
//
// [FIX 3] Per-base indices.
//     The base quality was read at a.quality()[s], where s is the position inside the current chunk, not in the read.
//     The damage profile was read at subDeamDiNuc[graph_seq.size()][baseIX]: baseIX is the start of the chunk (the same
//     value for every base in it) and graph_seq.size() also counts gap columns. Now both use the real read position of
//     the base (readPosOf, which skips deletion columns) and the fragment length is the read length.
//
// [FIX 4] The damage model. This was the main problem.
//     The correct model is a chain of three events:
//         reference base --mutation (HKY, depends on branch time t)--> b_s --deamination--> b_d --sequencing error--> read base
//     b_s and b_d are never observed, so the likelihood of the base we see is the sum over all of their values:
//         P(read | ref, t) = SUM_bs  P_HKY(bs | ref, t) * P(read | bs)
//         P(read | bs)     = SUM_bd  P_deam(bd | bs, position in fragment) * P_seqerr(read | bd, quality)
//     The old code (a) applied sequencing error BEFORE damage, (b) then "summed over the read base" with if/else branches
//     that did the same thing, so it summed the four outcome probabilities of a distribution, which is always exactly 1,
//     i.e. log(1) = 0 (clamped to log(0.9999999)): the damage information was computed and then thrown away, and
//     (c) never summed over the post-mutation base b_s.
//     Now this function stores P(read | bs) for the four possible bs in BaseInfo::damageSeqErrProb. It does not depend on t,
//     so it is computed once per base here; MCMC::computeBaseLogLike (MCMC.h) multiplies it by the HKY term and sums.
//     The scalar BaseInfo::logLikelihood (used for pathMap, i.e. the starting estimates, and for the consistency check
//     at the end of this function) is now log P(read | ref, t = 0) = log damageSeqErrProb[ref], so mismatches and damage
//     now influence it. Unlike TrailMix there is no per-source ancient/modern switch: soibean has one damage profile
//     (--deam5p/--deam3p; the default is "no damage") and always uses it.
// =====================================================================================================================
static vector<AlignmentInfo*>* analyse_GAM(
    const bdsg::HashGraph& graph,
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


    // base -> index (A,C,G,T = 0..3; anything else maps to 0, as in TrailMix's precompute_GAM)
    auto nucToIX = [](char c) {
        switch (toupper(static_cast<unsigned char>(c))) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default:  return 0;
        }
    };

    for (; iter.has_current(); iter.advance())
    {
        ++n_reads;
        const vg::Alignment& a = *iter;



        if (a.identity() != 0)
        {
         
            using BaseInfo = AlignmentInfo::BaseInfo;
            

            unordered_map<string, int> pathCounts;
            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());
            // [FIX 1,2,3] see the explanation above analyse_GAM
            // OLD:             int baseIX = a.path().mapping()[0].position().is_reverse() ? a.sequence().size() - 1 : 0;
            // OLD:             int baseOnRead = baseIX;
            // OLD:             int Lseq = graph_seq.size();
            int Lseq = a.sequence().size();   // fragment length = read length (graph_seq.size() also counts gap columns)

            // reconstruct_graph_sequence() pushes one mppg_sizes entry per EDIT, not per mapping, and every mismatch/indel splits a
            // node match into several edits. Entry i therefore belongs to the mapping (graph node) that owns the i-th edit.
            vector<int> entry_to_mapping;
            for (int mi = 0; mi < a.path().mapping_size(); ++mi) {
                for (int ei = 0; ei < a.path().mapping(mi).edit_size(); ++ei) {
                    entry_to_mapping.push_back(mi);
                }
            }
            if (entry_to_mapping.size() != mppg_sizes.size()) {
                throw std::runtime_error("analyse_GAM: edit count does not match the reconstructed segment count for read " + a.name());
            }
            int colStart = 0;          // first column of the current chunk in graph_seq/read_seq
            int readPosAtChunk = 0;    // 0-based read position of the chunk's first base

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
                // [FIX 1] entry i of mppg_sizes belongs to the mapping that owns the i-th EDIT (not to mapping i)
                // OLD:                 vector<string> probPaths;
                // OLD:                 int nID = 0;
                // OLD:                 if (mppg_sizes.size() != a.path().mapping().size())
                // OLD:                 {
                // OLD:                     if (i >= mppg_sizes.size() - (mppg_sizes.size() - a.path().mapping().size())){
                // OLD:                         probPaths = {"No_support"};
                // OLD:                     }else{
                // OLD:                         nID = a.path().mapping()[i].position().node_id();
                // OLD:                         for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                // OLD:                             probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                // OLD:                         }
                // OLD:                     }
                // OLD:
                // OLD:
                // OLD:                 }else{
                // OLD:                     nID = a.path().mapping()[i].position().node_id();
                // OLD:                     for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                // OLD:                         probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                // OLD:                     }
                // OLD:                 }
                vector<string> probPaths;
                int nID = a.path().mapping()[entry_to_mapping[i]].position().node_id();
                for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                    probPaths.emplace_back(nodepaths.at(nID - minid).at(j));
                }

                

                // [FIX 2,3] slice the next chunk (read order on both strands) and work out the read position of each column
                // OLD:                 string nodeSeq, partReadSeq;
                // OLD:                 int startIndex = 0;
                // OLD:                 if (a.path().mapping()[0].position().is_reverse()) {
                // OLD:                     startIndex = (baseIX - mppg_sizes.at(i) - 1 >= 0) ? (baseIX - mppg_sizes.at(i) - 1) : 0;
                // OLD:                     nodeSeq = graph_seq.substr(startIndex, mppg_sizes.at(i));
                // OLD:                     partReadSeq = read_seq.substr(startIndex, mppg_sizes.at(i));
                // OLD:                 } else {
                // OLD:                     nodeSeq = graph_seq.substr(baseIX, mppg_sizes.at(i));
                // OLD:                     partReadSeq = read_seq.substr(baseIX, mppg_sizes.at(i));
                // OLD:                 }
                string nodeSeq = graph_seq.substr(colStart, mppg_sizes.at(i));
                string partReadSeq = read_seq.substr(colStart, mppg_sizes.at(i));

                // 0-based read position of every column of this chunk ('-' in the read = deletion = no read base); used for the
                // per-base quality score and the position-dependent deamination probabilities
                std::vector<int> readPosOf;
                readPosOf.reserve(partReadSeq.size());
                int readPosEnd = readPosAtChunk;
                for (size_t s2 = 0; s2 < partReadSeq.size(); ++s2) {
                    readPosOf.push_back(readPosEnd);
                    if (partReadSeq[s2] != '-') {
                        ++readPosEnd;
                    }
                }
                auto qualityIndex = [&](size_t s) { return std::min<int>(readPosOf[s], (int)a.quality().size() - 1); };
                auto damageIndex  = [&](size_t s) { return std::min<int>(readPosOf[s], Lseq - 1); };
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
        			        << /* OLD: baseIX + s (start of the chunk, not the read index) */ readPosOf[s] << '\n'; 

			        }
                            }

                            if (nodeSeq[s] != partReadSeq[s]){mismatch++;
#ifdef DEBUGANALYSEGAM
                                //get<2>(sancheck[pathNames[m]]) += 1;
                                //get<8>(sancheck[pathNames[m]]) += log((qscore_vec[base_quality])/3);
#endif
		            }
#ifdef DEBUGDAMAGE
                            // OLD: cerr << "Base Index " << baseOnRead << endl;   (baseOnRead no longer exists, see [FIX 3])
                            cerr << "Base Index " << readPosOf[s] << endl;
                            cerr << "Graph Base " << nodeSeq[s] << endl;
                            cerr << "Read Base " << partReadSeq[s] << endl;
#endif
                            // [FIX 3] s is the position inside the chunk, not in the read
                            // OLD:                             int base_quality = a.quality()[s];
                            int base_quality = a.quality()[qualityIndex(s)];
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
                                // [FIX 4] Damage + sequencing error, marginalized over the (unobserved) post-mutation base. See the explanation
                                // above analyse_GAM. The old block below built the sequencing-error distribution BEFORE damage, then summed the four
                                // outcomes with identical if/else branches (= exactly 1, log = 0), so damage never influenced anything.
                                // OLD:                                 double probBasePreDamage [4];
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:                                     probBasePreDamage[bpd]=0.0;
                                // OLD:                                 }
                                // OLD:
                                // OLD:                                 // LEVEL 1 //
                                // OLD:                                 // filling a substitution matrix with the propbabilities of observing a base pre damage.
                                // OLD:                                 for(int bpo=0;bpo<4;bpo++){
                                // OLD:                                     if("ACGT"[bpo] == nodeSeq[s]){// no mutation
                                // OLD:                                         probBasePreDamage[bpo]= 1 - qscore_vec[base_quality];
                                // OLD:                                     }else{ // mutation
                                // OLD:                                         probBasePreDamage[bpo]= (qscore_vec[base_quality]/3);
                                // OLD:                                     }
                                // OLD:
                                // OLD:                                 }
                                // OLD: #ifdef DEBUGDAMAGE 
                                // OLD:                                 cerr << "Pre damage matrix" << endl;
                                // OLD:                                 for(int bpo=0;bpo<4;bpo++){
                                // OLD:                                     cerr<<setprecision(14)<<bpo<<"\t"<<probBasePreDamage[bpo]<<endl;
                                // OLD:                                 }
                                // OLD: #endif
                                // OLD:
                                // OLD:                                 // LEVEL 2 //
                                // OLD:                                 //initializing a post damage substitution matrix. First filled with 0.
                                // OLD:                                 double probBasePostDamage [4];
                                // OLD:
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:                                     probBasePostDamage[bpd]=0.0;
                                // OLD:                                 }
                                // OLD:
                                // OLD:
                                // OLD:                                 // filling post-damage substitution matrix with the new probabilities of observing a base. The probabilities are based on the damage profile for --deam5p and --deam3p
                                // OLD:
                                // OLD:                                 // If no damage profile is provided the pre-damage substitution matrix is multiplied by 0 and stays the same.
                                // OLD:                                 //cerr << "base index " << baseIX << " read size " << a.sequence().size() << endl;
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:                                     //post damage = prob pre damage * damage rate from bpo to bpd
                                // OLD:                                     for(int bpo=0;bpo<4;bpo++){
                                // OLD:                                         probBasePostDamage[bpd]+=probBasePreDamage[bpo]*subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                // OLD:                                         //cerr << setprecision(14) <<bpd<<"\t"<< subDeamDiNuc[Lseq][baseIX].p[bpo][bpd];
                                // OLD:                                     }
                                // OLD:                                 }
                                // OLD: #ifdef DEBUGDAMAGE 
                                // OLD:                                 cerr << "After damage matrix " << endl;
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:
                                // OLD:                                     cerr<<setprecision(14)<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                // OLD:                                 }
                                // OLD: #endif
                                // OLD:                                 double log_lik_marg = -std::numeric_limits<double>::infinity(); //sum of the likelihoods but in log space, this will be added to log_lik
                                // OLD:
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:                                     if( "ACGT"[bpd]== partReadSeq[s]){// no sequencing error
                                // OLD:                                     //                                                    Prob of bpd                no seq error
                                // OLD:                                         log_lik_marg = oplusInitnatl( log_lik_marg , (log(probBasePostDamage[bpd]) ) );
                                // OLD:
                                // OLD:                                     }else{ // we are already marginalising over all possible bases, if the base is not matching it must be a seq error.
                                // OLD:                                     //                                                    Prob of bpd                seq error
                                // OLD:                                         log_lik_marg = oplusInitnatl( log_lik_marg , (log(    probBasePostDamage[bpd] ))) ;
                                // OLD:                                     }
                                // OLD:                                 }
                                // OLD:
                                // OLD:                                 if(log_lik_marg > log(0.9999999)){
                                // OLD:                                     log_lik_marg = log(0.9999999);
                                // OLD:                                 }
                                // OLD:
                                // OLD:                                 if(isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
                                // OLD:                                     throw runtime_error("calculated log like is nan");
                                // OLD:                                 }
                                // OLD:
                                // OLD:
                                // OLD: #ifdef DEBUGDAMAGE                                
                                // OLD:                                 cerr<<"postdamage:"<<endl;
                                // OLD:
                                // OLD:                                 for(int bpd=0;bpd<4;bpd++){
                                // OLD:                                     cerr<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                                // OLD:                                 }
                                // OLD:
                                // OLD:                                 cerr<<"log_lik_marg:"<<log_lik_marg<<" p="<<exp(log_lik_marg)<<endl;
                                // OLD:                                 //cout <<"log_lik_marg " <<  log_lik_marg << endl;
                                // OLD:
                                // OLD:                                 cerr <<  pathMap[pathNames[m]] << endl;
                                // OLD: #endif
                                // OLD:                                 pathMap[pathNames[m]] += log_lik_marg;
                                // OLD:
                                // OLD:                                 BaseInfo info;
                                // OLD:                                 info.readBase = partReadSeq[s];
                                // OLD:                                 info.referenceBase = nodeSeq[s];
                                // OLD:                                 info.pathSupport = true;
                                // OLD:                                 info.logLikelihood = log_lik_marg;
                                // OLD:                                 if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("damage info.logLikelihood is nan");}
                                // OLD: #ifdef DEBUGPAIN2
                                // OLD:                                 cerr<< std::setprecision(25) << "damage info.logLikelihood " << info.logLikelihood << " "<< log_lik_marg << endl;
                                // OLD: #endif
                                // OLD:                                 readInfo.push_back(info);
                                int readBaseIX = nucToIX(partReadSeq[s]);
                                int refBaseIX  = nucToIX(nodeSeq[s]);
                                // damageSeqErr[bs] = P(read base | post-mutation base bs) = SUM_bd P_deam(bd | bs, position) * P_seqerr(read | bd, quality)
                                double damageSeqErr[4];
                                for (int bs = 0; bs < 4; bs++) {
                                    double acc = 0.0;
                                    for (int bd = 0; bd < 4; bd++) {
                                        double pSeqErr = (bd == readBaseIX) ? (1 - qscore_vec[base_quality]) : (qscore_vec[base_quality] / 3);
                                        acc += subDeamDiNuc[Lseq][damageIndex(s)].p[bs][bd] * pSeqErr;
                                    }
                                    damageSeqErr[bs] = max(0.0001, min(1.0, acc));
                                }
                                // scalar used for pathMap (starting estimates) and the consistency check below: no mutation (t = 0), so the
                                // post-mutation base is the reference base of this path. MCMC::computeBaseLogLike uses the full 4-vector instead.
                                double log_lik_marg = log(damageSeqErr[refBaseIX]);
                                if(isnan(log_lik_marg) || isinf(log_lik_marg) || log_lik_marg > 1e-8){
                                    throw runtime_error("calculated log like is nan");
                                }
                                pathMap[pathNames[m]] += log_lik_marg;

                                BaseInfo info;
                                info.readBase = partReadSeq[s];
                                info.referenceBase = nodeSeq[s];
                                info.pathSupport = true;
                                info.logLikelihood = log_lik_marg;
                                for (int bs = 0; bs < 4; bs++) {
                                    info.damageSeqErrProb[bs] = damageSeqErr[bs];
                                }
                                if(isnan(info.logLikelihood) || isinf(info.logLikelihood) || info.logLikelihood > 1e-8){throw runtime_error("damage info.logLikelihood is nan");}
#ifdef DEBUGPAIN2
                                cerr<< std::setprecision(25) << "damage info.logLikelihood " << info.logLikelihood << " "<< log_lik_marg << endl;
#endif
                                readInfo.push_back(info);

#ifdef DEBUGANALYSEGAM
                                get<1>(sancheck[pathNames[m]]) += 1;
                                get<8>(sancheck[pathNames[m]]) += log(1 -(qscore_vec[base_quality])/3);
#endif



                                
                            
                                // [FIX 3] no running base counter any more (positions come from readPosOf)
                                // OLD:                                 if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-' )
                                // OLD:                                 {
                                // OLD:                                     baseOnRead--;
                                // OLD:                                 } else if(!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                                // OLD:                                     baseOnRead++;
                                // OLD:                                 }
                            }
                        }

                        // the path is not supported
                        
                    } else if (find(probPaths.begin(), probPaths.end(), pathNames[m]) == probPaths.end()) {
                        int counter = 0;
                        for (int s = 0; s < nodeSeq.size(); ++s)
                        {
                            unsup++;
                            // [FIX 3] s is the position inside the chunk, not in the read
                            // OLD:                             int base_quality = a.quality()[s];
                            int base_quality = a.quality()[qualityIndex(s)];
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
                                // [FIX 3] was the running counter baseOnRead (counted from the wrong end on reverse reads); now the real read position
                                // OLD:                                 if (abs(baseOnRead) % PENALTY  == 0)
                                if (abs(readPosOf[s]) % PENALTY  == 0)
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

                            // [FIX 3] no running base counter any more (positions come from readPosOf)
                            // OLD:                             if (a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-') {
                            // OLD:                                 baseOnRead--;
                            // OLD:                             } else if(!a.path().mapping()[0].position().is_reverse() && partReadSeq[s] != '-'){
                            // OLD:                                 baseOnRead++;
                            // OLD:                             }


                        }


                    } else {
                        cerr << "Something is wrong!" << endl;

                    }
                    
                    detailMap[pathNamesS[m]].push_back(readInfo);

                    // [FIX 3] no running base counter any more (positions come from readPosOf)
                    // OLD:                     if (a.path().mapping()[0].position().is_reverse()) {
                    // OLD:                         baseOnRead = baseIX;
                    // OLD:                     } else {
                    // OLD:                         baseOnRead = baseIX;
                    // OLD:                     }

                    
                }
                //cout << readInfo.size() << endl;

                // [FIX 2,3] advance to the next chunk of graph_seq/read_seq and of the read
                // OLD:                 if (a.path().mapping()[0].position().is_reverse()) {
                // OLD:                     baseIX = startIndex;
                // OLD:                     baseOnRead = baseIX;
                // OLD:                 } else {
                // OLD:                     baseIX += mppg_sizes.at(i);
                // OLD:                     baseOnRead = baseIX;
                // OLD:                 }
                colStart += mppg_sizes.at(i);
                readPosAtChunk = readPosEnd;

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
