#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/odgi.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
#include "handle.hpp"
#include "io/save_handle_graph.hpp"
#include "NodeInfo.h"
#include "HaploCart.h"
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>

using namespace vg;
using namespace vg::algorithms;
using namespace google::protobuf;

void Haplocart::readPathHandleGraph(shared_ptr<Trailmix_struct> &dta) {

     if (!dta->running_trailmix){
        std::string gbwt_index_file = dta->hcfiledir + "graph.gbwt";
        std::string gbwtgraph_index_file = dta->hcfiledir + "graph.gg";
        dta->gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_index_file);
        dta->gbwtgraph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwtgraph_index_file);
     }


      std::string filePath = dta->graph_dir + dta->graphfilename;

    if (!std::filesystem::exists(filePath)) {
        throw std::runtime_error("File does not exist: " + filePath);
    }

    dta->graph.deserialize(filePath);
    dta->minid = dta->graph.min_node_id();
    dta->maxid = dta->graph.max_node_id();

    // Create the path support matrix
    size_t N_nodes = dta->graph.get_node_count();
    size_t N_paths = dta->path_names.size();

    if (N_paths == 0){throw runtime_error("NO PATHS");}

    vector<vector<bool>> node_path_matrix(N_nodes, vector<bool>(N_paths, false));
    for (size_t path_id = 0; path_id < N_paths; ++path_id) {
        // Get the nodes in the path
        gbwt::vector_type path_nodes = dta->gbwt->extract(path_id);
        for (const auto& node_id : path_nodes) {
            // Convert the GBWT node ID to the node handle in the graph
            bdsg::handle_t handle = dta->graph.get_handle(node_id);
            // Determine the index in the node_path_matrix
            int64_t index = dta->graph.get_id(handle) - 1;
            if (index >= 0 && index < N_nodes) {
                node_path_matrix[index][path_id] = true;
            }
        }
    }

    vector<NodeInfo *> nodevector;
    const int cladeid = 0;
    const int nbpaths = dta->graph.get_path_count();

    for (int64 i = dta->minid; i <= dta->maxid; ++i) {
        NodeInfo* nodetoadd = new NodeInfo(i, nbpaths, cladeid);
        const auto nodehandle = dta->graph.get_handle(i);
        nodetoadd->seq = dta->graph.get_sequence(nodehandle);
        long unsigned int j;
        for (j = 0; j < nbpaths; ++j) {
            nodetoadd->pathsgo[j] = node_path_matrix[i-1][j];
        }
        dta->nodevector.emplace_back(std::move(nodetoadd));
    }

    return;
}

