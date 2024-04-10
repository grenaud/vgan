#pragma once
#include <gbwtgraph/gbwtgraph.h>
#include <gbwt/metadata.h>
#include <gbwtgraph/path_cover.h>
#include <vg/io/vpkg.hpp>
#include <vg/vg.pb.h>
#include <vg.hpp>
#include <handlegraph/handle_graph.hpp>
#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
#include "handle.hpp"
#include "io/save_handle_graph.hpp"
#include "NodeInfo.h"
#include "HaploCart.h"
#include "Euka.h"
#include "soibean.h"
#include "progress_bar.hpp"
#include <iostream>
#include <thread>
#include <chrono>
using namespace std::chrono;
using namespace vg;
using namespace vg::algorithms;
using namespace google::protobuf;


tuple<vector<NodeInfo *>, int, bdsg::ODGI,vector<vector<bool>>, vector<string>> Euka::readPathHandleGraph (string & ogfilename, int n_threads, string &gbtwfilename, string &db_prefix,  vector<Clade *> * & clade_vec) {

    bdsg::ODGI graph;
    graph.deserialize(ogfilename);
    const int minid = graph.min_node_id();
    const int maxid = graph.max_node_id();
    // load gbwt graph file
    auto gbwt_graph = vg::io::VPKG::load_one<gbwt::GBWT>(gbtwfilename);
    // create path support file 
    size_t N_nodes = graph.get_node_count();

    vector<string> pathNames;
    unsigned int p_index = 0;
    graph.for_each_path_handle([&](const handlegraph::path_handle_t &path_handle) {
        string path_name = graph.get_path_name(path_handle);
        //cerr << p_index << '\t' << path_name << endl;
        pathNames.emplace_back(path_name);
        p_index++;
    });

    size_t N_paths = pathNames.size();




    vector<vector<bool>> node_path_matrix(N_nodes, vector<bool>(N_paths, false));

    for (size_t path_id = 0; path_id < N_paths; ++path_id) {
        // Get the nodes in the path
        gbwt::vector_type path_nodes = gbwt_graph->extract(path_id);
        for (const auto& node_id : path_nodes) {
            //cout << node_id << " "; 
            // Convert the GBWT node ID to the node handle in the graph
            bdsg::handle_t handle = graph.get_handle(node_id);
            // Determine the index in the node_path_matrix
            int64_t index = graph.get_id(handle) - 1;
            //cout << index << endl;
            if (index >= 0 && index < N_nodes) {
                node_path_matrix[index][path_id] = true;
                //cout << "I get positive too " << endl; 
            }
            
        }
        //cout << endl;
    
    }

    
    vector<NodeInfo *> nodevector;
    const int cladeid=0;
    int nbpaths;
    if (db_prefix == "soibean_db"){nbpaths = pathNames.size();}
    else if (db_prefix == "euka_db"){nbpaths = pathNames.size();}
    else{
        for (int i = 1; i<clade_vec->size(); i+=6){
            if(clade_vec->at(i)->name == db_prefix){
                nbpaths = clade_vec->at(i)->noPaths;
            }
        }
        

    }

    ProgressBar p(maxid-minid);
    p.SetFrequencyUpdate(200);

    auto start = high_resolution_clock::now();

    for(int64 i=minid;i<=maxid;++i){

        NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
        const auto nodehandle = graph.get_handle(i);
        nodetoadd->seq = graph.get_sequence(nodehandle);
        long unsigned int j;
        for (j = 0; j < nbpaths; ++j) {
            //nodetoadd->pathsgo[j] = node_path_matrix[i - 1][j];
        }
        nodevector.emplace_back(move(nodetoadd));

        p.Progressed(i-minid);

        //std::this_thread::sleep_for (std::chrono::milliseconds(10));
    }

    auto stop = high_resolution_clock::now();

    cerr << "Reading graph takes " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;

    assert(nodevector.size() != 0);
    assert(minid != 0);
    // Check if at least one true value exists
    bool hasTrueValue = false;
    for (const auto& row : node_path_matrix) {
        for (bool value : row) {
            if (value) {
                hasTrueValue = true;
                break;
            }
        }
        if (hasTrueValue) {
            break;
        }
    }

    //assert(hasTrueValue);
    assert(pathNames.size() != 0);

    


   return make_tuple(nodevector, minid, graph, node_path_matrix, pathNames);
}
