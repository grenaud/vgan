#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
#include "handle.hpp"
#include "io/save_handle_graph.hpp"
#include "NodeInfo.h"
#include "HaploCart.h"

using namespace vg;
using namespace vg::algorithms;
using namespace google::protobuf;

const tuple<vector<NodeInfo *>, const int, const bdsg::ODGI> Haplocart::readPathHandleGraph (const string & ogfilename, const int n_threads, const string &hcfiledir) {
    bdsg::ODGI graph;
    graph.deserialize(ogfilename);
    const int minid = graph.min_node_id();
    const int maxid = graph.max_node_id();
    vector<NodeInfo *> nodevector;
    const int cladeid=0;
    const int nbpaths=graph.get_path_count();
    const vector<vector<bool>> path_supports = Haplocart::load_path_supports(hcfiledir);
    for(int64 i=minid;i<=maxid;++i){
                                 NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
                                 const auto nodehandle = graph.get_handle(i);
                                 nodetoadd->seq = graph.get_sequence(nodehandle);
                                 long unsigned int j;
                                 #pragma omp parallel for private(j) num_threads(n_threads)
                                 for (j=0; j<nbpaths; ++j)
                                     {
                                        {nodetoadd->pathsgo[j] = path_supports[i][j];}
                                     }
                                 nodevector.emplace_back(nodetoadd);
                                   }

   return make_tuple(nodevector, minid, graph);
}


