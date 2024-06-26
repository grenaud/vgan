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

void Haplocart::readPathHandleGraph (shared_ptr<Trailmix_struct> &dta) {
    dta->graph.deserialize(dta->graph_dir+dta->graphfilename);
    dta->minid = dta->graph.min_node_id();
    dta->maxid = dta->graph.max_node_id();
    vector<NodeInfo *> nodevector;
    const int nbpaths=dta->graph.get_path_count();
    load_path_supports(dta);
    for(int64 i=dta->minid;i<=dta->maxid;++i){
                                 NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,0);
                                 const auto nodehandle = dta->graph.get_handle(i);
                                 nodetoadd->seq = dta->graph.get_sequence(nodehandle);
                                 long unsigned int j;
                                 //#pragma omp parallel for private(j) num_threads(dta->n_threads)
                                 for (j=0; j<nbpaths; ++j)
                                     {
                                        {nodetoadd->pathsgo[j] = dta->path_supports[i][j];}
                                     }
                                 dta->nodevector.emplace_back(nodetoadd);
                                   }
}


