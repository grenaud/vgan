#pragma once
#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
#include "handle.hpp"
#include "io/save_handle_graph.hpp"
#include "NodeInfo.h"
#include "HaploCart.h"
#include "Euka.h"
#include "progress_bar.hpp"
#include <iostream>
#include <thread>
#include <chrono>
using namespace std::chrono;
using namespace vg;
using namespace vg::algorithms;
using namespace google::protobuf;

tuple<vector<NodeInfo *>, int, bdsg::ODGI> Euka::readPathHandleGraph (string ogfilename, int n_threads, string pathsupport) {

    bdsg::ODGI graph;
    graph.deserialize(ogfilename);
    const int minid = graph.min_node_id();
    const int maxid = graph.max_node_id();
    vector<NodeInfo *> nodevector;
    const int cladeid=0;
    const int nbpaths=graph.get_path_count();
    const vector<vector<bool>> path_supports = Euka().load_path_supports_Euka(pathsupport);

    ProgressBar p(maxid-minid);
    p.SetFrequencyUpdate(200);

    auto start = high_resolution_clock::now();

    for(int64 i=minid;i<=maxid;++i){

                                 NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
                                 const auto nodehandle = graph.get_handle(i);
                                 nodetoadd->seq = graph.get_sequence(nodehandle);
                                 nodevector.emplace_back(nodetoadd);

                                 p.Progressed(i-minid);

                                 //std::this_thread::sleep_for (std::chrono::milliseconds(10));
    }

   auto stop = high_resolution_clock::now();

   cerr << "Reading graph takes " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;

   return make_tuple(nodevector, minid, graph);
}
