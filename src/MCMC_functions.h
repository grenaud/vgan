#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << setprecision(10) << v[i] << '\t';}cerr << endl << endl;
#define PRINT_KEYS(map) do { \
    for (const auto& pair : map) { \
        std::cerr << pair.first << "\n"; \
    } \
} while(0)


void modifyPathNameInPlace(string &path_name){
     string original_path_name = path_name;

    // Special case: a single letter followed by underscore
    std::regex pattern2("^([A-Za-z])_$");
    path_name = std::regex_replace(path_name, pattern2, "$1");

    // Handle the special path format
    std::regex pattern3("\\+([0-9]+)\\+\\(([0-9]+)\\)");
    path_name = std::regex_replace(path_name, pattern3, "_$1__$2_");

    // Now replace special characters
    std::replace(path_name.begin(), path_name.end(), '+', '_');
    std::replace(path_name.begin(), path_name.end(), '\'', '_');
    std::replace(path_name.begin(), path_name.end(), '*', '_');
    std::replace(path_name.begin(), path_name.end(), '@', '_');
    std::replace(path_name.begin(), path_name.end(), '(', '_');
    std::replace(path_name.begin(), path_name.end(), ')', '_');

    if (!path_name.empty() && path_name.back() == '*') {
        path_name.back() = '_';
    }

   if (path_name.size() == 1){path_name = path_name + "_";}

    // Remove .1 or .2 at the end
    std::regex patternEnd("\\.(1|2)$");
    path_name = std::regex_replace(path_name, patternEnd, "");

}


// Function to calculate the distance to an ancestor node
const double MCMC::calculateDistanceToAncestor(spidir::Node* startNode, spidir::Node* ancestor) {

    double distance = 0.0;
    spidir::Node* current = startNode;

    //std::cerr << "Calculating distance from node ID " << startNode->name
    //          << " to ancestor ID " << ancestor->name << std::endl;

    while (current != nullptr && current != ancestor) {
        //std::cerr << "Current node ID: " << current->name << ", distance: " << current->dist << std::endl;
        distance += (current->dist);
        current = current->parent;
    }

    if (distance < 0.0){throw runtime_error("NEGATIVE DISTANCE!");}

    if (distance == 0.0 && startNode != ancestor){throw runtime_error("ZERO DISTANCE!");}

    if (current == ancestor) {
        //std::cerr << "Total distance: " << distance << std::endl;
        return distance;
    } else {
        std::cerr << "Ancestor not found in path. Returning -1." << std::endl;
        return -1.0; // Return -1 if ancestor is not found
    }
}

double getQuantile2(const std::vector<double>& sortedData, double q) {
    if (sortedData.empty()) {
        throw std::runtime_error("Vector is empty");
    }

    if (q < 0.0 || q > 1.0) {
        throw std::invalid_argument("Quantile must be between 0 and 1");
    }

    const auto n = sortedData.size();
    const auto index = (n - 1) * q;
    const auto lowerIndex = static_cast<size_t>(floor(index));
    const auto upperIndex = static_cast<size_t>(ceil(index));

    if (lowerIndex == upperIndex) {
        return sortedData[lowerIndex];
    } else {
        const auto frac = index - lowerIndex;
        return (1.0 - frac) * sortedData[lowerIndex] + frac * sortedData[upperIndex];
    }
}

double MCMC::calculateRhat(const std::vector<double>& means, const std::vector<double>& variances, int chainLength, int numChains) {
    if (numChains < 2) {
        std::cerr << "Error: At least two chains are required to compute R-hat." << std::endl;
        return -1;
    }
    if (means.size() != variances.size()) {
        throw std::runtime_error("Error: Mismatched sizes of means and variances vectors.");
    }
    if (chainLength <= 0) {
        std::cerr << "Error: Chain length must be positive." << std::endl;
        return -1;
    }

    // Compute within-chain variance W
    double W = 0.0;
    for (double variance : variances) {
        if (variance != 0) {
            W += variance;
        }
    }
    W /= numChains;

    // Compute between-chain variance B
    double grandMean = 0.0;
    int validMeansCount = 0;
    for (double mean : means) {
        if (mean != 0) {
            grandMean += mean;
            validMeansCount++;
        }
    }
    grandMean /= validMeansCount;

    double B = 0.0;
    for (double mean : means) {
        if (mean != 0) {
            B += std::pow(mean - grandMean, 2);
        }
    }
    B *= chainLength / (numChains - 1);

    // Compute R-hat
    double varEstimate = ((chainLength - 1.0) * W + B) / chainLength;
    double Rhat = std::sqrt(varEstimate / W);

    return Rhat;
}


// A helper function to recursively calculate the distance from any node to a leaf.
const double MCMC::calculateDistanceToLeaf(const spidir::Tree* tr, const spidir::Node* current, const spidir::Node* leaf, double currentDistance, \
                               std::unordered_set<const spidir::Node*>& visited) {
    // Base case: if the current node is the leaf, return the current distance
    if (current == leaf) {
        return currentDistance;
    }

    // Mark the current node as visited
    visited.insert(current);

    // Explore all children using index-based access
    for (size_t i = 0; i < current->nchildren; ++i) {
        const spidir::Node* child = current->children[i];
        if (visited.find(child) == visited.end()) {
            double distance = calculateDistanceToLeaf(tr, child, leaf, currentDistance + child->dist, visited);
            if (distance >= 0) {
                // If a valid distance was found, return it
                return distance;
            }
        }
    }

    // If the current node has a parent and it hasn't been visited, explore it
    if (current->parent && visited.find(current->parent) == visited.end()) {
        return calculateDistanceToLeaf(tr, current->parent, leaf, currentDistance + current->dist, visited);
    }
}


inline const std::vector<double> MCMC::getPatristicDistances(const spidir::Tree* tr, spidir::Node* node, int numofLeafs, double posonbranch) {
    if (!tr || !node) {
        throw std::runtime_error("Tree or node is null.");
    }
    std::vector<double> distances(numofLeafs, std::numeric_limits<double>::max());
    for (size_t i = 0; i < tr->nodes.size(); ++i) {
        const auto &leafNode = tr->nodes[i];
        if (!leafNode || !leafNode->isLeaf()) {
            continue;
        }
        spidir::Node* lca = findLCA(tr, node, leafNode);
        double distanceToLCAFromNode = calculateDistanceToAncestor(node, lca) + (posonbranch * node->dist);
        double distanceToLCAFromLeaf = calculateDistanceToAncestor(leafNode, lca);
        if (distanceToLCAFromNode >= 0.0 && distanceToLCAFromLeaf >= 0.0) {
            distances[i] = distanceToLCAFromNode + distanceToLCAFromLeaf;

        } else {
            std::cerr << "Error calculating distance for leaf ID " << leafNode->name << std::endl;
            cerr << "distances: " << distanceToLCAFromNode << "  " << distanceToLCAFromLeaf << endl;
        }
    }
    return distances;
}

const double MCMC::calculateEuclideanDistance(const std::vector<double> &vec1, const std::vector<double> &vec2) {
    if (vec1.empty() || vec2.empty()) {
        throw std::runtime_error("Input vectors are empty.");
    }

    size_t minSize = std::min(vec1.size(), vec2.size());
    if (minSize == 0) {
        throw std::runtime_error("One of the vectors is empty.");
    }

    double sum = 0.0;
    for (size_t i = 0; i < minSize; ++i) {
        if (vec1[i] == std::numeric_limits<double>::max() || vec2[i] == std::numeric_limits<double>::max()) {
            std::cerr << "Warning: Skipping comparison at index " << i << " due to invalid value." << std::endl;
            continue;
        }
        double diff = vec1[i] - vec2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

void printUnorderedMap(const std::unordered_map<std::string, std::vector<std::vector<double>>>& myMap) {
    for (const auto& pair : myMap) {
        std::cerr << "Key: " << pair.first << "\n";
        std::cerr << "Values:\n";
        for (const auto& innerVec : pair.second) {
            std::cerr << "  [ ";
            for (const auto& value : innerVec) {
                std::cerr << value << " ";
            }
            std::cerr << "]\n";
        }
        std::cerr << "\n";
    }
}


pair<unordered_map<string, vector<vector<double>>>, double> MCMC::processMCMCiterations(
    shared_ptr<Trailmix_struct> &dta,
    const vector<MCMCiteration>& MCMCiterations,
    int k,
    const string &num,
    int chain,
    spidir::Tree* tr,
    int numofleafs) {

    unordered_map<string, vector<vector<double>>> branchStatisticsMap;
    double chainloglike = numeric_limits<double>::lowest();
    ofstream readEstimatesFile, baseEstimatesFile, branchestimateFile;

    // Open files for writing statistics
    baseEstimatesFile.open(num + "BaseProportionEstimates.txt", ios::app | ios::out);
    branchestimateFile.open(num + "BranchEstimate.txt", ios::app | ios::out);

    if (!baseEstimatesFile.is_open() || !branchestimateFile.is_open()){
        cerr << "Failed to open files for writing." << endl;
        cerr << "num: " << num << endl;
        throw runtime_error("File opening failed.");
    }

    baseEstimatesFile << "Source\tChain\tMean Prop Est\t5% CI\tMedian Prop Est\t95% CI\tESS\tAutocorrelation\tVariance\n";
    branchestimateFile << "Source\tChain\tMean Branch Position\t5% CI\tMedian Branch Position\t95% CI\tESS\tAutocorrelation\tVariance\tESS for the source estimation\n";

    if (MCMCiterations.empty()) {
        cerr << "MCMC Iterations Vector is Empty!" << endl;
        throw runtime_error("MCMC Iterations Vector is Empty!");
    }

    cerr << "Starting processing MCMC iterations..." << endl;
    //cerr << "NUMBER OF MCMCITERATIONS: " << MCMCiterations.size() << endl;

    for (int source = 0; source < k; ++source) {
        cerr << "Processing source: " << source + 1 << endl;

        vector<double> baseProportionVec, readProportionVec, position_vec, euc_distances;
        string branchName;
        vector<double> initialPatristicDistances = vector<double>(numofleafs, 1.0);
        vector<string> branch_names;

        int iterationCounter=0;
        for (const auto& iteration : MCMCiterations) {
            iterationCounter++;
            if (iterationCounter % 500 == 0){cerr << "Processing MCMC iter: " << iterationCounter << endl;}
            chainloglike = max(chainloglike, iteration.logLike);

            if (dta){
            branchName = dta->originalPathNames[iteration.positions_tree[source].pos->longname];
                    }
            else{
            branchName = iteration.positions_tree[source].pos->longname;
                }
            //modifyPathNameInPlace(branchName);
            branch_names.emplace_back(branchName);

            auto baseprops = iteration.proportions;
            baseProportionVec.emplace_back(baseprops[source]);
            position_vec.emplace_back(iteration.positions_tree[source].pos_branch);

            double posonbranch = iteration.positions_tree[source].pos->dist * iteration.positions_tree[source].pos_branch;
            const vector<double> patristic_distances = getPatristicDistances(tr, iteration.positions_tree[source].pos, numofleafs, posonbranch);
            euc_distances.emplace_back(calculateEuclideanDistance(patristic_distances, initialPatristicDistances));

            if (branchStatisticsMap.find(branchName) == branchStatisticsMap.end()) {
                branchStatisticsMap[branchName] = vector<vector<double>>();
            }
        }


        double meanBaseTheta = mean(baseProportionVec);
        double meanPos = mean(position_vec);
        cerr << "autocorrelations..." << endl;
        double Theta_autoc = autocorrelation(baseProportionVec, 1);
        double Pos_autoc = autocorrelation(position_vec, 1);
        cerr << "...Done. ESS: " << endl;
        double Theta_ess = effectiveSampleSize(baseProportionVec);
        double Pos_ess = effectiveSampleSize(position_vec);
        double dist_ess = effectiveSampleSize(euc_distances);
        cerr << "...Done" << endl;
        double theta_var = variance(baseProportionVec, meanBaseTheta);
        double Pos_var = variance(position_vec, meanPos);

        sort(baseProportionVec.begin(), baseProportionVec.end());
        sort(position_vec.begin(), position_vec.end());

        double base_theta_median = getQuantile2(baseProportionVec, 0.5);
        double Pos_median = getQuantile2(position_vec, 0.5);
        double base_theta_fq = getQuantile2(baseProportionVec, 0.05);
        double base_theta_tq = getQuantile2(baseProportionVec, 0.95);
        double Pos_fq = getQuantile2(position_vec, 0.05);
        double Pos_tq = getQuantile2(position_vec, 0.95);

        baseEstimatesFile << branchName << "\t" << chain << "\t" << meanBaseTheta << "\t" << base_theta_fq << "\t" << base_theta_median << "\t" << \
                             base_theta_tq << "\t" << Theta_ess << "\t" << Theta_autoc << "\t" << theta_var << "\n";

        branchestimateFile << branchName << "\t" << chain << "\t" << meanPos << "\t" << Pos_fq << "\t" << Pos_median << "\t" << Pos_tq << \
                         "\t" << Pos_ess << "\t" << Pos_autoc << "\t" << Pos_var << "\t" << dist_ess << "\n";

        vector<double> sourceStatistic = {meanBaseTheta, theta_var, meanPos, Pos_var};
        branchStatisticsMap[branchName].emplace_back(sourceStatistic);
    }

    return {branchStatisticsMap, chainloglike};
}
