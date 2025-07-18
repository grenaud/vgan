#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gzstream.h>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <utility>
#include <queue>
#include <optional>
#include "libgab.h"
#include "Euka.h"
#include "HaploCart.h"
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Clade.h"
#include "getLCAfromGAM.h"
#include "subcommand/subcommand.hpp"
#include "bdsg/odgi.hpp"


using namespace std;

struct frags {
    string name;
    string sequence;
    string graphSeq;
    vector<int> nodeIDs;
    vector<int> offsets;
    vector<int> baseQ;
    vector<int> coverage;
    std::vector<std::array<double, 5>> probReadPostDamage;
	vector<array<int,5>> countsRead;
    pair<bool, int> cutbool;
    string rymer_sequence;


    void display() const {
        std::cerr << "Name: " << name << "\nSequence: " << sequence << "\nNode IDs: ";
        for (int id : nodeIDs) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
    }
};

// struct GraphData {
//     std::unordered_map<uint64_t, uint64_t> node_depths;

//     // Get depth by node ID
//     uint64_t get_depth_by_node_id(uint64_t node_id) {
//         auto it = node_depths.find(node_id);
//         return it != node_depths.end() ? it->second : 0; // Return depth if exists, 0 if not
//     }
// };
// Define the GraphData struct
struct GraphData {
    std::unordered_map<uint64_t, uint64_t> node_depths;

    // Get depth by node ID
    uint64_t get_depth_by_node_id(uint64_t node_id) const {
        auto it = node_depths.find(node_id);
        return it != node_depths.end() ? it->second : 0; // Return depth if exists, 0 if not
    }
};

struct pair_hash {
    template <class T1, class T2>
    size_t operator () (const std::pair<T1,T2>& pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ hash2;
    }
};

struct GraphAss {
    std::vector<std::map<int, double>> adj_list;  // stores outgoing edges
    std::vector<int> in_degree;
    std::unordered_set<std::pair<int, int>, pair_hash> used_edges;

    GraphAss(int n) : adj_list(n), in_degree(n, 0) {}

    void add_edge(int u, int v, double s) {
        adj_list[u].insert({v, s});
        in_degree[v]++;
    }

    bool is_edge_used(int u, int v) const {
        return used_edges.find({u, v}) != used_edges.end();
    }

    void mark_edge_used(int u, int v) {
        used_edges.insert({u, v});
    }
    void unmark_edge_used(int u, int v) {
            used_edges.erase({u, v});
        }
    
    bool edgeExists(int node1, int node2) {
        if (node1 < 0 || node1 >= adj_list.size()) {
            std::cerr << "Invalid node index: " << node1 << std::endl;
            return false;
        }
        // Using map's find method to search for node2 in the adjacency list of node1
        return adj_list[node1].find(node2) != adj_list[node1].end();
    }
    int numVertices() const {
        return adj_list.size();
    }

    void print_adjacency_list() {
        for (size_t i = 0; i < adj_list.size(); ++i) {
            std::cout << "Adjacency list for vertex " << i << ": ";
            for (const auto& pair : adj_list[i]) {
                int neighbor = pair.first;
                double edge_weight = pair.second;
                std::cout << "(" << neighbor << ", " << edge_weight << ") ";
            }
            std::cout << std::endl;
        }

    }

};

struct OverlapResult {
        int length;   // Length of the overlap
        double score; // Score of the overlap
        string sequence;

    };

struct MergeResult {
    bool canMerge;
    OverlapResult best_overlap;
    bool isOverlapAtEnd;
    int numN;
};

class assembly{

private:
    void reindex_odgi_graph(bdsg::ODGI& graph, GraphData& new_graph_data, int startNode, int endNode);
    void writeContigsToFasta(const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>>& contigs, const std::string& filename);
    void saveToFastaGz(const std::string& fasta, const std::string& filename);
	string reverse_complement(const std::string& dna);
	bool basesMatch(char a, char b, double& mismatch_penalty);
    bool isRYMatch(char a, char b);
    bool basesMatchWithDamage(char a, char b);
	void convertToRYmerSpace(std::vector<frags>& reads);
    int calculateMinOverlapLength(int numberOfReads,int averNumber);
    double calculate_match_score(char a, char b);
    pair<unordered_map<int,vector<int>>, int> initial_overlap(const std::vector<frags>& reads,  const std::unordered_map<uint64_t, uint64_t>& node_depths);
	void find_overlaps(const std::vector<frags>& reads, int min_overlap_length, GraphAss& g, std::map<std::pair<int, int>, OverlapResult>& overlap_map, std::unordered_map<int, std::vector<int>> density_map, int min_score, int averNumber, int lenMin, bool specifiedDeam);
	std::string merge_sequences(const std::string& a, const std::string& b, int overlap_length);

    std::vector<int> sortOverlapsByScore(int vertex, const GraphAss& graph, const std::map<std::pair<int, int>, OverlapResult>& overlap_map);
    void updateScoringMatrix(std::vector<std::unordered_map<char, double>>& scoringMatrix,
    const std::string& contig,
    const std::string& newSeq,
    size_t startIdx,
    const std::vector<std::array<double, 5>>& probReadPostDamage);
    void updateCountMatrix(std::vector<std::unordered_map<char, int>>& countMatrix,
                                 const std::string& contig,
                                 const std::string& newSeq,
                                 size_t startIdx, const std::vector<std::array<int, 5>>& countReads,  int best_overlap_length, bool contigfirst);
    std::vector<std::unordered_map<char, double>> createScoringMatrix(const std::string& contig,
    const std::vector<std::array<double, 5>>& probReadPostDamage);
    std::vector<std::unordered_map<char, int>> createCountMatrix(const std::string& contig, const std::vector<std::array<int, 5>>& countReads);
    void adjustFinalCut(pair<bool, int>& finalcut, const pair<bool, int>& neighborCut);
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> mergeAllPaths(const GraphAss& graph, const std::map<std::pair<int, int>, OverlapResult>& overlap_map, const std::vector<frags>& reads, int min_overlap_length, bool specifiedDeam, bool minimeta, GraphData new_graph_data, int lenMin); 
    //std::vector<std::pair<int, int>> findPartialOverlappingContigs(const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, string>>& contigs);
    bool isSubset(const std::string& a, const std::vector<int>& aNodeIds,
                        const std::string& b, const std::vector<int>& bNodeIds,
                        int allowedMismatches = 1);
    size_t findAlignmentPosition(const std::string& contig, const std::string& read, int allowedMismatches);
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> removeSubsetContigs(const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>>& contigs);
    
    std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>> mergeContigs(const std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>& contig1,
const std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>& contig2,
int overlapLength, bool isOverlapAtEnd, string mode
    );
    MergeResult tryMergeContigs(
	    const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>>& contigs,
	    int i,
	    int j,
	    int min_overlap_length,
	    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap,
	    std::unordered_set<int>& mergedIndices, bool specifiedDeam, int lenMin
    );

    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> mergeOverlappingContigs(
		std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int> >> contigs,
		    int min_overlap_length,
		    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap, bool specifiedDeam, string mode, int lenMin
    );
    void printScoringMatrix(const std::vector<std::unordered_map<char, double>>& scoringMatrix, string &outfile, size_t startAt = 0);
    void printCountMatrix(const std::vector<std::unordered_map<char, int>>& countMatrix, std::string &outfile, size_t startAt = 0);

    std::optional<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> checkAndMergeContigs(
	    std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>& contig_i,
	    std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>& contig_j,
	    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap, int lenMin, int scoreMin, bool specifiedDeam, string mode);
    
    

public:
    set<int> deadEnds;
    std::map<std::string, int> pathScores;
    OverlapResult get_overlap_length_and_score(const std::string& a, const std::string& b, int min_overlap_length, int min_score);
	static bool compareByFirstNodeID(const frags& a, const frags& b, const GraphData& graphData);

	assembly();
	assembly(const assembly & other);
	~assembly();

	assembly & operator= (const assembly & other);

    string usage() const;
    const int run(int argc, char *argv[] , const string & cwdProg);

};
