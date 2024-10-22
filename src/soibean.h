#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gzstream.h>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include "libgab.h"
#include "Euka.h"
#include "HaploCart.h"
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Clade.h"
#include "getLCAfromGAM.h"
#include "subcommand/subcommand.hpp"
#include "bdsg/odgi.hpp"
#include "../dep/spimap/src/Tree.h"
#include "../dep/spimap/src/newick.h"

using namespace std;

class soibean{

private:

	vector<string> paths_through_node(const bdsg::ODGI& graph, const bdsg::handle_t& node);
	std::vector<std::vector<double>> convertMapsToVector(vector<AlignmentInfo*>* & gam);
	unordered_map<int, int> makePathToNode(spidir::Tree &tree, vector<string>& path_names);
	std::vector<unsigned int> generateRandomNumbers(const int maxNum, const int k);
        double calculateRhat(const std::vector<double>& means, const std::vector<double>& variances, int chainLength);



public:
	soibean();
	soibean(const soibean & other);
	~soibean();

	soibean & operator= (const soibean & other);

    const string usage(const std::string& cwdProg) const;
    const int run(int argc, char *argv[] , const string & cwdProg);

};


