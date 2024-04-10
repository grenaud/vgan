#pragma once
#include "forward_declarations.h"
#include "Trailmix_struct.h"
#include "../dep/spimap/src/Tree.h"
#include "../dep/spimap/src/newick.h"

spidir::Tree make_tree_from_dnd(shared_ptr<Trailmix_struct> &dta)                 {

spidir::Tree tree = NULL;
spidir::readNewickTree(dta->treePath.c_str(), &tree);

//    // iterate over all nodes in the tree
//    for (int i = 0; i < tree.nodes.size(); i++) {
//        cerr << "Name: " << tree.nodes[i]->name << ", Long Name: " << tree.nodes[i]->longname << endl;
//    }

return tree;
                                                                                  }


