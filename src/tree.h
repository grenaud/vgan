#pragma once
#include "forward_declarations.h"
#include "Trailmix_struct.h"
#include "../dep/spimap/src/Tree.h"
#include "../dep/spimap/src/newick.h"
#include "vgan_utils.h"

spidir::Tree make_tree_from_dnd(shared_ptr<Trailmix_struct> &dta) {
    spidir::Tree tree = NULL;
    spidir::readNewickTree(dta->treePath.c_str(), &tree);
    for (int i = 0; i < tree.nodes.size(); i++) {
         if (tree.nodes[i]->isLeaf()){dta->n_leaves++;}

     }
    return tree;
}

std::unordered_map<int, int> inverse_map(const std::unordered_map<int, int>& original) {
    std::unordered_map<int, int> inverse;
    for (const auto& [key, value] : original) {
        inverse.insert(std::make_pair(value, key));
    }
    return inverse;
}

void writeDebugFile(const unordered_map<string, int>& path_node_map, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing.\n";
        return;
    }

    for (const auto& pair : path_node_map) {
        file << pair.first << ": " << pair.second << '\n';
    }

    file.close();
}

void Trailmix::create_path_node_map(shared_ptr<Trailmix_struct> &dta) {
    unordered_map<string, int> path_node_map;
    unordered_map<int, string> node_path_map; // Map from node index to path name

    unordered_set<string> path_names_set(dta->path_names.begin(), dta->path_names.end());

    for (size_t i = 0; i < dta->tree->nodes.size(); i++) {
        if (!dta->tree->nodes[i]) continue; // Safety check

        string node_longname = dta->tree->nodes[i]->longname;
        modifyPathNameInPlace(dta, node_longname);

        // Efficient search with unordered_set
        if (path_names_set.find(node_longname) != path_names_set.end()) {
            path_node_map[node_longname] = i;
            node_path_map[i] = node_longname; // Populate the inverse map
        }
    }

    dta->path_node_map = move(path_node_map);
    dta->node_path_map = move(node_path_map);
}


