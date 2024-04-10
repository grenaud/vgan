#pragma once

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



using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

long double calcPathLogLike(
    const bdsg::ODGI& graph,
    const vector<AlignmentInfo*>* align,
    vector<vector<string>>& nodepaths,
    vector<string> pathNames,
    const vector<double> qscore_vec,
    double givenMu,
    bool entire_graph,
    bool trailmix, 
    int minid, vector<long double> distr)
{
    int n_reads = 0;

    if (align->empty()) {
        throw std::runtime_error("Alignment vector is empty. Aborting.");
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
    
    // declare the cumulativePathMap outside the loop
    map<string, long double> cumulativePathMap;
    long double logLike = 0.0L;
    // initialize it with all the paths and zero values
    for (const auto& path : pathNames) {
        cumulativePathMap[path] = 0.0;
    }
    for (AlignmentInfo* a : *align)
    {
        ++n_reads;
        
            
        unordered_map<string, int> pathCounts;
        auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a->path, a->seq);

        int baseIX = a->path.mapping()[0].position().is_reverse() ? a->seq.size() - 1 : 0;
        int baseOnRead = baseIX;

        unordered_map<string, long double> pathMap;
        unordered_map<string, bool> supportMap;



        for (int c = 0; c < pathNames.size(); ++c)
        {
            pathMap[pathNames[c]] = 0.0;
            supportMap[pathNames[c]] = false;
        }


        for (int i = 0; i < mppg_sizes.size(); ++i)
        {
            vector<string> probPaths;
            int nID = 0;
            if (mppg_sizes.size() != a->path.mapping().size())
            {
                if (i == mppg_sizes.size() -1){
                    probPaths = {"No_support"};
                }else{
                    nID = a->path.mapping()[i].position().node_id();
                    for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                        probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                    }
                }
                

            }else{
                nID = a->path.mapping()[i].position().node_id();
                for (int j = 0; j < nodepaths.at(nID - minid).size(); ++j) {
                    probPaths.emplace_back(nodepaths.at(nID - minid).at(j)); 
                }
            }

            

            string nodeSeq, partReadSeq;
            int startIndex = 0;
            if (a->path.mapping()[0].position().is_reverse()) {
                startIndex = (baseIX - mppg_sizes.at(i) - 1 >= 0) ? (baseIX - mppg_sizes.at(i) - 1) : 0;
                nodeSeq = graph_seq.substr(startIndex, mppg_sizes.at(i));

                partReadSeq = read_seq.substr(startIndex, mppg_sizes.at(i));

            } else {
                nodeSeq = graph_seq.substr(baseIX, mppg_sizes.at(i));
       
                partReadSeq = read_seq.substr(baseIX, mppg_sizes.at(i));

            }
            
            
            assert(distr.size() == pathNames.size());
            for (int m = 0; m < pathNames.size(); ++m)
            {

                

                // the path is supported
                if (find(probPaths.begin(), probPaths.end(), pathNames[m]) != probPaths.end())
                {
                    supportMap[pathNames[m]] = true;
                    
                    for (int s = 0; s < nodeSeq.size(); ++s)
                    {
                        int base_quality = a->quality_scores[s];
                        if (nodeSeq[s] == 'N' || partReadSeq[s] == 'N')
                        {

                            continue;

                        } else if (nodeSeq[s] == 'S' || partReadSeq[s] == 'S')
                        {

                            pathMap[pathNames[m]] += log(qscore_vec[base_quality]);

                        }else if (nodeSeq[s] == '-' || partReadSeq[s] == '-')
                        {

                            pathMap[pathNames[m]] += log(qscore_vec[base_quality]);

                        } else if (nodeSeq[s] == partReadSeq[s])
                        {

                            pathMap[pathNames[m]] += log((1 - qscore_vec[base_quality]));
                            
                        } else if (nodeSeq[s] != partReadSeq[s])
                        {

                            pathMap[pathNames[m]] += log(qscore_vec[base_quality] / 3);

                        }else{
                            cerr << "Something is wrong!" << endl;
                        } 

                        if (a->path.mapping()[0].position().is_reverse())
                        {
                            baseOnRead--;
                        } else {
                            baseOnRead++;
                        }
                    }

                    // the path is not supported
                    
                } else if (find(probPaths.begin(), probPaths.end(), pathNames[m]) == probPaths.end())
                {
                    int counter = 0;
                    for (int s = 0; s < nodeSeq.size(); ++s)
                    {
                        int base_quality = a->quality_scores[s];

                        
                        if (abs(baseOnRead) % 4 == 4)
                        {   
                            pathMap[pathNames[m]] += log(1 - qscore_vec[base_quality]);

                        } else 
                        {

                            pathMap[pathNames[m]] += log(qscore_vec[base_quality]);
                        }



                        if (a->path.mapping()[0].position().is_reverse()) {
                            baseOnRead--;
                        } else {
                            baseOnRead++;
                        }


                    }


                } else {
                    cerr << "Something is wrong!" << endl;

                }

                if (a->path.mapping()[0].position().is_reverse()) {
                    baseOnRead = baseIX;
                } else {
                    baseOnRead = baseIX;
                }


            }


            if (a->path.mapping()[0].position().is_reverse()) {
                baseIX = startIndex;
                baseOnRead = baseIX;
            } else {
                baseIX += mppg_sizes.at(i);
                baseOnRead = baseIX;
            }

        }
        
        if (pathNames.size() == 1)
        {
            long double inter = pathMap[pathNames[0]];
            logLike += inter;
        }
        else
        {
            // long double inter = 0.0L;
            // for (int y = 0; y < pathNames.size(); ++y)
            // {
            //     inter = oplusInitnatl(inter, pathMap[pathNames[y]]);
            // }
            long double inter = oplusnatl(log((distr[0]) + pathMap[pathNames[0]]), (log(distr[1]) + pathMap[pathNames[1]]));
            logLike += inter;

        }
        
        
         
        
    } // end of iteration through gam file 

    return logLike;
} // end of function 


