#include "soibean.h"
#include <iostream>
#include <handlegraph/handle_graph.hpp>

using namespace std;

// Function to find the paths that go through a given node
vector<string> soibean::paths_through_node(const bdsg::ODGI& graph, const bdsg::handle_t& node) {
    vector<string> paths;
    graph.for_each_step_on_handle(node, [&](const bdsg::step_handle_t& step) {
        // Get the path associated with this step
        bdsg::path_handle_t path_handle = graph.get_path_handle_of_step(step);
        // Get the path name and add it to the result vector
        paths.push_back(graph.get_path_name(path_handle));
    });
    return paths;
}

unordered_map<int, int> soibean::makePathToNode(spidir::Tree &tree, vector<string>& path_names)
{
    unordered_map<int, int> path_node_map;
    auto it = path_names.begin();  // Initialize the iterator

    if (tree.nodes.size() == 0) {
        cout << "The tree is empty" << endl;
        return path_node_map;
    }

    for (size_t i = 0; i < tree.nodes.size(); i++)
    {
        if (it != path_names.end())
        {
            const int found_idx = std::distance(path_names.begin(), it);
            path_node_map[found_idx] = i;
            ++it;  // Increment the iterator
        }
        else
        {
            break;  // Exit the loop if all path_names have been processed
        }
    }

    return path_node_map;
}
