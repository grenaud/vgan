#include "assembly.h"
#include "soibean.h"
#include <omp.h>
#include <sys/wait.h>
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/vg.pb.h>
#include <vg.hpp>
#include <boost/algorithm/string.hpp>
#include <handlegraph/handle_graph.hpp>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "miscfunc.h"
#include "vgan_utils.h"
#include "Clade.h"
#include "MCMC.h"
#include "damage.h"

//#define DEBUGASS
//#define DEBUGDAM
//#define DEBUGASS2
//#define DEBUGASS3
//#define DEBUGASS4
//#define DEBUGMERGEALL
#define INDELERRORPROB 1.0e-5 // https://link.springer.com/article/10.1186/gb-2011-12-11-r112
#define SIMILARITY 0.98


using namespace vg;


assembly::assembly(){

}

assembly::~assembly(){

}

string assembly::usage() const {
    return string("") +

    "\n"+
    "\n"+


       "            ____         ____                                        ____     \n"+   
       " |    ..'' |            |            |        |       .'. .`.       |         \n"+     
       " |..''     |______      |______      |        |     .'   `   `.     |______   \n"+     
       " |``..     |            |            |        |   .'           `.   |          \n"+    
       " |    ``.. |___________ |___________ |_______ | .'               `. |___________ \n"+  

    "\n"+
    "\n"+
    "Usage: keelime [options]\n" +
        "\n"+
        "Example:\n"+
        "\tvgan keelime -fq1 [input.fq.gz] --dbprefix [taxon name] -o [output prefix] -p [path name] --mode normal\n"+
        "\n"+
		"Example with damage profiles:\n"+
		"\tvgan keelime -fq1 [input.fq.gz] --dbprefix [taxon name] -o [output prefix] -p [path name] --mode normal --deam5p ../share/damageProfiles/dhigh5.prof --deam3p ../share/damageProfiles/dhigh3.prof\n"+
        
        "Options:\n" +
        //"  --print                     Prints available references paths names\n" +
        "  --keelime_dir <directory>   Specify the directory containing the soibean files\n" +
        "  --dbprefix <prefix>         Specify the prefix for the database files\n" +
        "  -o [STR]                    Output file prefix (default: keeOut) "+"\n"+
        "  -fq1 <filename>             Specify the input FASTQ file (single-end or first pair)\n" +
        "  -fq2 <filename>             Specify the input FASTQ file (second pair)\n" +
        "  -g <filename>               Specify the input GAM file (SAFARI output)\n" +
        "  -M [STR]                    Alternative minimizer prefix input\n"+
        "  -t <threads>                Specify the number of threads to use (default: 1)\n" +
        "  -z <directory>              Specify the temporary directory for intermediate files (default: /tmp/)\n" +
        "  -p [STR]                    Specify the name of the reference path\n" +
        
        "Assembly parameters:"+"\n"+
        "  --mode [STR]                normal: positions with contradicting base counts will be resolved by a 65% majority rule (failing positions will be masked) (default mode).\n"+
        "                              reckless: positions with contradicting base counts will be resolved by majority rule (equal base counts for one position will be masked).\n"+
        "                              strict: positions with contradicting base counts will be resolved by a 90% majority rule (failing positions will be masked).\n"+
		"\n"+
		"  -c  [INT]                   Minimum number of bases covering a position in the consensus genome (default: 1).\n"+
        "  -mL [INT]                   Minimum length for overlap between contigs (default: 10)\n"+
        "  -mS [INT]                   Minimum score for the overlap between contigs (default: 15)\n"+
        "  -uR [BOOL]                  Boolean for unknown reference sequence (no N bridges over uncovered bases/ unknown node IDs) (default false)\n"+
        "  -uC [BOOL]                  Boolean to add all assembled contigs to the consensus fasta without a nodeID assignment (default: false).\n"+
       
        "Damage options:"+"\n"+
        "  --deam5p [.prof]            5p deamination frequency for eukaryotic species (default: no damage)"+"\n"+
        "  --deam3p [.prof]            3p deamination frequency for eukaryotic species (default: no damage)"+"\n"+

        
        "\n"+
        "\n";
        
}
void assemblySetup() {

 preflight_check();
  configure_memory_allocator();
  enable_crash_handling();
  temp_file::set_system_dir();

  if (!vg::io::register_libvg_io()) {
      cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
      exit(1);
                                    }

    std::cerr << 
    
       "            ____         ____                                        ____     \n"   
       " |    ..'' |            |            |        |       .'. .`.       |         \n"     
       " |..''     |______      |______      |        |     .'   `   `.     |______   \n"     
       " |``..     |            |            |        |   .'           `.   |          \n"    
       " |    ``.. |___________ |___________ |_______ | .'               `. |___________ \n"  

              << std::endl;
}

////////////// Assembly function storage - temporary ///////////////////////////

void assembly::reindex_odgi_graph(bdsg::ODGI& graph, GraphData& new_graph_data, int startNode, int endNode) {
    std::queue<std::pair<bdsg::handle_t, uint64_t>> bfs_queue;  // Queue of pairs of handle and current depth
    std::unordered_map<bdsg::handle_t, uint64_t> visited_depths;  // Maps handles to their max depth encountered

    bdsg::handle_t start_handle = graph.get_handle(startNode);
    bdsg::handle_t end_handle = graph.get_handle(endNode);

    bfs_queue.push({start_handle, 1});  // Start with depth 1 for the starting node
    visited_depths[start_handle] = 1;

    while (!bfs_queue.empty()) {
        auto [current_handle, current_depth] = bfs_queue.front();
        bfs_queue.pop();

        uint64_t current_id = graph.get_id(current_handle);
        
        // Update the node_depths map with the maximum depth found so far
        if (new_graph_data.node_depths.find(current_id) == new_graph_data.node_depths.end() ||
            new_graph_data.node_depths[current_id] < current_depth) {
            new_graph_data.node_depths[current_id] = current_depth;
        }

        // Check if the current node is the end node
        if (current_handle == end_handle) {
            std::cerr << "End node " << endNode << " reached at depth " << current_depth << ". Terminating BFS." << std::endl;
            break;
        }

        // Process all edges leading from the current node
        graph.follow_edges(current_handle, false, [&](const bdsg::handle_t& next_handle) {
            if (visited_depths.find(next_handle) == visited_depths.end() || visited_depths[next_handle] < current_depth + 1) {
                bfs_queue.push({next_handle, current_depth + 1});
                visited_depths[next_handle] = current_depth + 1;  // Update the depth for this handle
            }
        });
    }

    std::cerr << "Depth assignment complete. Total nodes processed: " << new_graph_data.node_depths.size() << std::endl;
}


void assembly::writeContigsToFasta(const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>>& contigs, const std::string& filename) {
    gzFile fastaFile = gzopen((filename + ".gz").c_str(), "wb");
    
    if (!fastaFile) {
        std::cerr << "Error opening file " << filename << ".gz" << std::endl;
        perror("Additional error info");  // Print system error message
        return;
    }

    for (size_t i = 0; i < contigs.size(); ++i) {
        const auto& contig = contigs[i];
        const std::string& sequence = std::get<0>(contig);

        // Debugging output for header writing attempt
        //std::cerr << "Attempting to write header for contig " << i + 1 << std::endl;
        int headerResult = gzprintf(fastaFile, ">contig_%zu\n", i + 1);
        if (headerResult <= 0) {
            int err;
            const char* headerErrorMsg = gzerror(fastaFile, &err);
            std::cerr << "Failed to write header for contig " << i + 1 << ": " << headerErrorMsg << std::endl;
            continue; // Skip writing sequence if header fails
        }

        // Writing the sequence in chunks
        const size_t CHUNK_SIZE = 1024;  // Adjust as needed
        for (size_t pos = 0; pos < sequence.length(); pos += CHUNK_SIZE) {
            size_t end = std::min(pos + CHUNK_SIZE, sequence.length());
            int sequenceResult = gzprintf(fastaFile, "%s", sequence.substr(pos, end - pos).c_str());
            if (sequenceResult <= 0) {
                int err;
                const char* sequenceErrorMsg = gzerror(fastaFile, &err);
                std::cerr << "Failed to write sequence chunk starting at " << pos << " for contig " << i + 1 << ": " << sequenceErrorMsg << std::endl;
                break;  // Exit the loop on failure
            }
        }
        gzprintf(fastaFile, "\n");  // Write a newline after the whole sequence

        // Flush after each contig to ensure data is written to the file system
        if (gzflush(fastaFile, Z_SYNC_FLUSH) != Z_OK) {
            std::cerr << "Failed to flush gzip buffer after writing contig " << i + 1 << std::endl;
        }
    }
    
    // Close the file
    if (gzclose(fastaFile) != Z_OK) {
        std::cerr << "Failed to close file properly " << filename << ".gz" << std::endl;
    }
}

void assembly::saveToFastaGz(const std::string& fasta, const std::string& filename) {
    gzFile gzFastaFile = gzopen((filename + ".gz").c_str(), "wb");
    if (!gzFastaFile) {
        std::cerr << "Error opening file " << filename << ".gz" << std::endl;
        return;
    }
    // Adding a header line for the FASTA sequence
    std::string header = ">Consensus\n"; // Ensure there is a newline after the header
    if (gzputs(gzFastaFile, header.c_str()) < 0) {
        std::cerr << "Error writing header to file " << filename << ".gz" << std::endl;
    }
    
    if (gzputs(gzFastaFile, fasta.c_str()) < 0) {
        std::cerr << "Error writing to file " << filename << ".gz" << std::endl;
    }
    if (gzputs(gzFastaFile, "\n") < 0) {
        std::cerr << "Error writing newline to file " << filename << ".gz" << std::endl;
    }

    gzclose(gzFastaFile);
    //std::cout << "FASTA written to " << filename << ".gz" << std::endl;
}

bool assembly::compareByFirstNodeID(const frags& a, const frags& b, const GraphData& graphData) {
    // Helper function to get depth of a fragment's first node ID
    auto get_depth = [&](const frags& f) {
        if (f.nodeIDs.empty()) {
            return static_cast<uint64_t>(0); // Return 0 depth if nodeIDs is empty
        }
        return graphData.get_depth_by_node_id(f.nodeIDs.front());
    };

    uint64_t depthA = get_depth(a);
    uint64_t depthB = get_depth(b);

    // Primary comparison by depth
    if (depthA != depthB) {
        return depthA < depthB;
    }

    // Check if either nodeIDs vector is empty and sort such that non-empty comes first
    if (a.nodeIDs.empty() || b.nodeIDs.empty()) {
        return b.nodeIDs.empty() && !a.nodeIDs.empty();
    }

    // Secondary comparison by the front of nodeIDs
    if (a.nodeIDs.front() != b.nodeIDs.front()) {
        return a.nodeIDs.front() < b.nodeIDs.front();
    }

    // Tertiary comparison by offsets, considering potential emptiness
    if (a.offsets.empty() || b.offsets.empty()) {
        return b.offsets.empty() && !a.offsets.empty();
    }
    if (a.offsets.front() != b.offsets.front()) {
        return a.offsets.front() < b.offsets.front();
    }

    // Quaternary comparison by cutbool
    if (a.cutbool.first != b.cutbool.first) {
        return a.cutbool.first > b.cutbool.first; // true (1) comes before false (0)
    }

    // Detailed element-by-element comparison of nodeIDs vectors
    auto size = std::min(a.nodeIDs.size(), b.nodeIDs.size());
    for (size_t i = 0; i < size; ++i) {
        if (a.nodeIDs[i] != b.nodeIDs[i]) {
            return a.nodeIDs[i] < b.nodeIDs[i];
        }
    }

    // If one nodeIDs list is shorter and all previous elements are equal, the shorter list is considered less
    if (a.nodeIDs.size() != b.nodeIDs.size()) {
        return a.nodeIDs.size() < b.nodeIDs.size();
    }

    // Quinary comparison by read length
    if (a.sequence.size() != b.sequence.size()) {
        return a.sequence.size() > b.sequence.size();
    }

    // All criteria equal; preserve original order (stable sort)
    return false;
}





std::string assembly::reverse_complement(const std::string& dna) {
    std::string result;
    result.reserve(dna.size());  // Reserve space to improve performance

    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        switch (*it) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            case 'S': result += 'S'; break;
            case '-': result += '-'; break;
            case 'N': result += 'N'; break;
            default:  // In case of an invalid character
                std::cerr << "Invalid DNA base encountered: " << *it << std::endl;
                return "";
        }
    }
    return result;
}

bool assembly::basesMatch(char a, char b, double& mismatch_penalty) {
    if (a == b) {
        mismatch_penalty = 0.0;
        return true; // Perfect match
    }

    // Gap penalty
    if (a == '-' || b == '-') {
        mismatch_penalty = 0.0;
        return false;
    }
	
    // Gap penalty
    if (a == 'N' || b == 'N') {
        mismatch_penalty = 0.0;
        return false;
    }

    // Define RY space matching rules
    switch (a) {
        case 'R':
        case 'r':
            if (b == 'A' || b == 'a' || b == 'G' || b == 'g' || b == 'R' || b == 'r') {
                mismatch_penalty = 0.0;
                return true;
            }
            break;
        case 'Y':
        case 'y':
            if (b == 'C' || b == 'c' || b == 'T' || b == 't' || b == 'Y' || b == 'y') {
                mismatch_penalty = 0.0;
                return true;
            }
            break;
        case 'A':
        case 'a':
        case 'G':
        case 'g':
            if (b == 'R' || b == 'r') {
                mismatch_penalty = 0.0;
                return true;
            }
            break;
        case 'C':
        case 'c':
        case 'T':
        case 't':
            if (b == 'Y' || b == 'y') {
                mismatch_penalty = 0.0;
                return true;
            }
            break;
    }

    // Ancient damage mismatches
    if ((a == 'C' && b == 'T') || (a == 'c' && b == 't') || (a == 'G' && b == 'A') || (a == 'g' && b == 'a') || (a == 'T' && b == 'C') || (a == 't' && b == 'c') || (a == 'A' && b == 'G') || (a == 'a' && b == 'g')) {
        mismatch_penalty = 0.5;
        return false;
    }

    // All other mismatches
    mismatch_penalty = 3.0;
    return false;
}



bool assembly::isRYMatch(char a, char b) {
    if (a == b) return true;

    // Define matching groups
    std::string Rs = "AGag";
    std::string Ys = "CTct";

    if ((Rs.find(a) != std::string::npos && Rs.find(b) != std::string::npos) ||
        (Ys.find(a) != std::string::npos && Ys.find(b) != std::string::npos)) {
        return true;
    }

    // Allow matching if either character is a '-'
    if (a == '-' || b == '-') {
        return true;
    }

    return false;
}


void assembly::convertToRYmerSpace(std::vector<frags>& reads) {
    for (auto& read : reads) {
        const std::string& original = read.sequence;
        
        // If the sequence is shorter than 10 bases, skip conversion
        if (original.length() < 10) {
            read.rymer_sequence = original; // Assign original sequence without modification
            continue; // Skip to the next read
        }

        std::string rymer_sequence;

        // Convert the first 5 bases
        size_t prefix_length = std::min(size_t(5), original.length());
        for (size_t i = 0; i < prefix_length; ++i) {
            char base = original[i];
            if (base == 'a' || base == 'A' || base == 'g' || base == 'G') {
                rymer_sequence += 'R';
            } else if (base == 'c' || base == 'C' || base == 't' || base == 'T') {
                rymer_sequence += 'Y';
            } else {
                rymer_sequence += base;  // Non-standard bases are added as is
            }
        }

        // Copy the middle section unchanged
        rymer_sequence += original.substr(5, original.length() - 10);

        // Convert the last 5 bases
        for (size_t i = original.length() - 5; i < original.length(); ++i) {
            char base = original[i];
            if (base == 'a' || base == 'A' || base == 'g' || base == 'G') {
                rymer_sequence += 'R';
            } else if (base == 'c' || base == 'C' || base == 't' || base == 'T') {
                rymer_sequence += 'Y';
            } else {
                rymer_sequence += base;  // Non-standard bases are added as is
            }
        }

        read.rymer_sequence = rymer_sequence;
    }
}

double assembly::calculate_match_score(char a, char b) {
    // Both are gaps
    if (a == '-' && b == '-') {
        return 1.0; // Treat as a match
    }

    // One of them is a gap
    if (a == '-' || b == '-') {
        return 0.0; // Penalize if only one is a gap
    }

    // One or both are N (unknown)
    if (a == 'N' || b == 'N') {
        return 0.0; // Mildly positive or neutral
    }

    // Exact match between A, C, G, and T
    if ((a == 'A' || a == 'C' || a == 'G' || a == 'T') && a == b) {
        return 3.0;
    } else if ((a == 'R' && (b == 'A' || b == 'G')) || (b == 'R' && (a == 'A' || a == 'G'))) {
        return 2.0; // R matches A or G
    } else if ((a == 'Y' && (b == 'C' || b == 'T')) || (b == 'Y' && (a == 'C' || a == 'T'))) {
        return 2.0; // Y matches C or T
    } else if ((a == 'R' && b == 'R') || (a == 'Y' && b == 'Y')) {
        return 1.0; // R matches R or Y matches Y
    } else if ((a == 'A' && b == 'G') || (a == 'G' && b == 'A') ||
               (a == 'C' && b == 'T') || (a == 'T' && b == 'C')) {
        return -1.0; // Less penalty for mismatch between A and G or C and T
    }

    return -3.0; // Mismatch between any other bases
}



int assembly::calculateMinOverlapLength(int numberOfReads, int averNumber) {
    int minOverlap = 2;
    int maxOverlap = 8;

    if (numberOfReads <= 1) {
        //std::cerr << "Number of reads is 1 or less. Returning minOverlap: " << minOverlap << std::endl;
        return minOverlap;
    }

    // Base overlap is the minimum overlap
    int baseOverlap = minOverlap;
    

    int increment = averNumber > 0 ? (maxOverlap - minOverlap) / averNumber : 0;
    if (increment < 1){
        increment = 1;
    }
    int proposedOverlap = baseOverlap + (numberOfReads - 1) * increment;
    int finalOverlap = std::max(minOverlap, std::min(proposedOverlap, maxOverlap));

    //std::cerr << "Returning finalOverlap: " << finalOverlap << std::endl;

    return finalOverlap;
}


std::string assembly::merge_sequences(const std::string& a, const std::string& b, int overlap_length) {
        return a + b.substr(overlap_length);
    return a;  
}



OverlapResult assembly::get_overlap_length_and_score(const std::string& a, const std::string& b, int min_overlap_length, int min_score) {
    int a_length = a.size();
    int b_length = b.size();
    OverlapResult best_overlap = {0, 0.0, ""}; // To store the best overlap found
#ifdef DEBUGASS3
    std::cerr << "Finding overlap between sequences:\n";
    std::cerr << "Sequence A: " << a << "\n";
    std::cerr << "Sequence B: " << b << "\n";
    std::cerr << "Minimum overlap length: " << min_overlap_length << "\n";
#endif

    // Iterate over possible starting points in sequence A
    for (int start = 0; start <= a_length - min_overlap_length; ++start) {
        int current_overlap_length = std::min(b_length, a_length - start);
        double current_score = 0.0;
        bool match = true;  // Assume a match unless proven otherwise
        int mismatch_count = 0;  // Counter for mismatches
#ifdef DEBUGASS3
        std::cerr << "Checking overlap starting at position " << start << " with overlap length " << current_overlap_length << "\n";
#endif
        // Check if the entire current segment matches, allowing for up to three mismatches
        for (int i = 0; i < current_overlap_length; ++i) {
            char base_a = a[start + i];
            char base_b = b[i];
            double mismatch_penalty;

            if (!basesMatch(base_a, base_b, mismatch_penalty)) {
                mismatch_count += mismatch_penalty;
                current_score -= mismatch_penalty;
#ifdef DEBUGASS3
                std::cerr << "Mismatch at A[" << start + i << "] = " << base_a << " and B[" << i << "] = " << base_b << " (Mismatch count: " << mismatch_count << ")\n";
#endif
                //if (mismatch_count > 3.0) {
                //     match = false;
                //     std::cerr << "Too many mismatches (" << mismatch_count << ") at position " << start + i << "\n";
                //     break;
                // }
            } else {
#ifdef DEBUGASS3
                std::cerr << "Match at A[" << start + i << "] = " << base_a << " and B[" << i << "] = " << base_b << "\n";
#endif
            }
        }

        // Only calculate score if a valid overlap with at most three mismatches was found
        if (match) {
#ifdef DEBUGASS3
            std::cerr << "Valid overlap found with " << mismatch_count << " mismatches. Calculating score...\n";
#endif
            for (int i = 0; i < current_overlap_length; ++i) {
                char base_a = a[start + i];
                char base_b = b[i];
                double score = calculate_match_score(base_a, base_b);
                current_score += score;
#ifdef DEBUGASS3
                std::cerr << "Match score for A[" << start + i << "] = " << base_a << " and B[" << i << "] = " << base_b << " is " << score << "\n";
#endif
            }
#ifdef DEBUGASS3
            std::cerr << "Current overlap score: " << current_score << "\n";
#endif
            // Store the best overlap result if it meets the minimum score requirement
            if (current_score >= min_score && current_overlap_length >= min_overlap_length) {
                OverlapResult result = {current_overlap_length, current_score, b};

                // Update the best overlap found
                if (result.length > best_overlap.length || (result.length == best_overlap.length && result.score > best_overlap.score)) {
                    best_overlap = result;
#ifdef DEBUGASS3
                    std::cerr << "Best overlap updated: Length=" << best_overlap.length << ", Score=" << best_overlap.score << "\n";
                    std::cerr << "Best overlap sequences: " << a.substr(start, current_overlap_length) << " " << b.substr(0, current_overlap_length) << "\n";
                    cerr << "Mismatch count " << mismatch_count << endl;
#endif
                }
            }
        } else {
#ifdef DEBUGASS3
            std::cerr << "No valid overlap found starting at position " << start << "\n";
#endif
        }
    }
#ifdef DEBUGASS3
    // Print the best valid overlap
    std::cerr << "Best overlap found: Length=" << best_overlap.length << ", Score=" << best_overlap.score << "\n";
#endif
    // Return the best overlap found or zero overlap and score if no suitable overlap is found
    return best_overlap;
}





// Assume node_depths is a map from node ID to its depth from the start node
std::pair<std::unordered_map<int, std::vector<int>>, int> assembly::initial_overlap(const std::vector<frags>& reads, const std::unordered_map<uint64_t, uint64_t>& node_depths) {
    int n = reads.size();
    std::vector<int> nodeDepths;
    std::unordered_map<int, std::vector<int>> density_map;

    std::cerr << "Processing " << n << " reads." << std::endl;

    // Convert node IDs to depths
    for (int i = 0; i < n; ++i) {
        int nodeID = reads[i].nodeIDs.front();  // Get the first node ID from each read
        if (node_depths.find(nodeID) != node_depths.end()) {
            nodeDepths.emplace_back(node_depths.at(nodeID));
        } else {
            nodeDepths.emplace_back(-1);  // Handle case where node ID is not in the depth map
        }
    }

    int inc = 0; // Initialize inc to avoid potential undefined behavior
    for (int i = 0; i < n; i++) {
        if (nodeDepths[i] == -1) continue;  // Skip if depth is not known

        int current_depth = nodeDepths[i];
        int max_depth = current_depth + 10;
        std::vector<int> potOverlaps;

        for (int j = 0; j < n; j++) {
            if (nodeDepths[j] == -1 || j == i) continue;  // Skip if depth is not known or it's the same read

            int target_depth = nodeDepths[j];
            if (target_depth >= current_depth && target_depth <= max_depth) {
                potOverlaps.emplace_back(j);
                inc++;
            }
        }

        density_map[i] = potOverlaps;
        //std::cerr << "Read " << i << " has " << potOverlaps.size() << " potential overlaps." << std::endl;
    }

    if (inc == 0) {
        //std::cerr << "Warning: Increment value 'inc' is zero, potential division by zero." << std::endl;
        inc = 1;  // Prevent division by zero
    }

    int avNoReads = inc / n;
    std::cerr << "Average number of reads calculated: " << avNoReads << std::endl;

    return {density_map, avNoReads};
}



void assembly::find_overlaps(const std::vector<frags>& reads, int min_overlap_length, GraphAss& g, std::map<std::pair<int, int>, OverlapResult>& overlap_map, std::unordered_map<int, std::vector<int>> density_map, int min_score, int averNumber, int lenMin, bool specifiedDeam) {
    
    for (const auto& pair : density_map) {
        for (int j : pair.second) {
            //min_overlap_length = calculateMinOverlapLength(pair.second.size(), averNumber);
			
			

            OverlapResult results1 = get_overlap_length_and_score(reads[pair.first].rymer_sequence, reads[j].rymer_sequence, min_overlap_length, lenMin);
            OverlapResult results2 = get_overlap_length_and_score(reads[j].rymer_sequence, reads[pair.first].rymer_sequence, min_overlap_length, lenMin);
			
			if(!specifiedDeam){
                results1 = get_overlap_length_and_score(reads[pair.first].sequence, reads[j].sequence, min_overlap_length, lenMin);
                results2 = get_overlap_length_and_score(reads[j].sequence, reads[pair.first].sequence, min_overlap_length, lenMin);
			}

            bool endOfFirstInSecond = std::find(reads[j].nodeIDs.begin(), reads[j].nodeIDs.end(), reads[pair.first].nodeIDs.back()) != reads[j].nodeIDs.end();
            bool endOfSecondInFirst = std::find(reads[pair.first].nodeIDs.begin(), reads[pair.first].nodeIDs.end(), reads[j].nodeIDs.back()) != reads[pair.first].nodeIDs.end();

            if (endOfFirstInSecond && results1.length >= min_overlap_length && results1.score >= min_score) {
                g.add_edge(pair.first, j, results1.score);
                overlap_map[{pair.first, j}] = results1;
            }

            if (endOfSecondInFirst && results2.length >= min_overlap_length && results2.score >= min_score) {
                g.add_edge(j, pair.first, results2.score);
                overlap_map[{j, pair.first}] = results2;
            }
        }
    }

    std::cerr << "Completed finding overlaps." << std::endl;
}

// Function to sort overlaps by score

std::vector<int> assembly::sortOverlapsByScore(int vertex, const GraphAss& graph, const std::map<std::pair<int, int>, OverlapResult>& overlap_map) {
    std::vector<std::pair<int, OverlapResult>> overlaps;

    for (const auto& neighbor : graph.adj_list[vertex]) {
        int target = neighbor.first;
        OverlapResult overlap = overlap_map.at({vertex, target});
        overlaps.push_back({target, overlap});
    }

    // Sort overlaps based on score in descending order
    std::sort(overlaps.begin(), overlaps.end(), [](const std::pair<int, OverlapResult>& a, const std::pair<int, OverlapResult>& b) {
        return a.second.score > b.second.score;
    });

    std::vector<int> sortedNeighbors;
    for (const auto& overlap : overlaps) {
        sortedNeighbors.push_back(overlap.first);
    }

    return sortedNeighbors;
}

std::vector<std::unordered_map<char, double>> assembly::createScoringMatrix(
    const std::string& contig,
    const std::vector<std::array<double, 5>>& probReadPostDamage
) {
    std::vector<std::unordered_map<char, double>> scoringMatrix(contig.size());

    for (size_t i = 0; i < contig.size(); ++i) {
        scoringMatrix[i]['A'] = std::log(probReadPostDamage[i][0]);
        scoringMatrix[i]['C'] = std::log(probReadPostDamage[i][1]);
        scoringMatrix[i]['G'] = std::log(probReadPostDamage[i][2]);
        scoringMatrix[i]['T'] = std::log(probReadPostDamage[i][3]);
        scoringMatrix[i]['-'] = log(probReadPostDamage[i][4]);
    }

    // Check for invalid values in the scoring matrix
    for (const auto& map : scoringMatrix) {
        for (const auto& [base, score] : map) {
            if (std::isinf(score) || std::isnan(score)) {
                std::cerr << "Error: Invalid value in scoring matrix during creation: " << score << "\n";
            }
        }
    }
#ifdef DEBUGDAM
    // Print the initial scoring matrix with the bases from the contig
    std::cerr << "Initial scoring matrix with bases from the contig:" << std::endl;
    for (size_t i = 0; i < contig.size(); ++i) {
        std::cerr << "Position " << i << " (" << contig[i] << "): A=" << std::exp(scoringMatrix[i]['A'])
                  << ", C=" << std::exp(scoringMatrix[i]['C'])
                  << ", G=" << std::exp(scoringMatrix[i]['G'])
                  << ", T=" << std::exp(scoringMatrix[i]['T']) << std::endl;
    }
#endif
    return scoringMatrix;
}

// Function to create a scoring matrix for each base of the contig
std::vector<std::unordered_map<char, int>> assembly::createCountMatrix(const std::string& contig, const std::vector<std::array<int, 5>>& countsRead) {
    std::vector<std::unordered_map<char, int>> countMatrix(contig.size());

    for (size_t i = 0; i < contig.size(); ++i) {
        countMatrix[i]['A'] = countsRead[i][0];
        countMatrix[i]['C'] = countsRead[i][1];
        countMatrix[i]['G'] = countsRead[i][2];
        countMatrix[i]['T'] = countsRead[i][3];
        countMatrix[i]['-'] = countsRead[i][4];
        
    }
    // Check for invalid values in the scoring matrix
    for (const auto& map : countMatrix) {
        for (const auto& [base, score] : map) {
            if (std::isinf(score) || std::isnan(score)) {
                std::cerr << "Error: Invalid value in count matrix during creation: " << score << "\n";
            }
        }
    }

    return countMatrix;
}


void assembly::updateScoringMatrix(
    std::vector<std::unordered_map<char, double>>& scoringMatrix,
    const std::string& contig,
    const std::string& newSeq,
    size_t startIdx,
    const std::vector<std::array<double, 5>>& probReadPostDamage
) {
    size_t originalSize = scoringMatrix.size();
    //std::cerr << "original Size: " << originalSize << std::endl;

    // Ensure the scoring matrix can accommodate the extended contig
    if (scoringMatrix.size() < contig.size()) {
        scoringMatrix.resize(contig.size());
        for (size_t i = originalSize; i < contig.size(); ++i) {
            scoringMatrix[i]['A'] = 0.0;
            scoringMatrix[i]['C'] = 0.0;
            scoringMatrix[i]['G'] = 0.0;
            scoringMatrix[i]['T'] = 0.0;
            scoringMatrix[i]['-'] = 0.0;
        }
    }
#ifdef DEBUGDAM
    // Print the initial scoring matrix with the bases
    std::cerr << "Initial scoring matrix with bases from the contig:" << std::endl;
    for (size_t i = 0; i < contig.size(); ++i) {
        std::cerr << "Position " << i << " (" << contig[i] << "): A=" << std::exp(scoringMatrix[i]['A'])
                  << ", C=" << std::exp(scoringMatrix[i]['C'])
                  << ", G=" << std::exp(scoringMatrix[i]['G'])
                  << ", T=" << std::exp(scoringMatrix[i]['T']) << std::endl;
    }
#endif
    // Update the scoring matrix with the new sequence and log probabilities
    for (size_t i = 0; i < newSeq.size(); ++i) {
        size_t matrixIndex = startIdx + i;
#ifdef DEBUGDAM
        std::cerr << "startIdx: " << startIdx << ", matrixIndex: " << matrixIndex << std::endl;
        std::cerr << "contig: " << contig << ", size: " << contig.size() << std::endl;
        std::cerr << "newSeq: " << newSeq << ", size: " << newSeq.size() << std::endl;
        std::cerr << "i: " << i << std::endl;
        std::cerr << "start2 (startIdx - newSeq.size()): " << probReadPostDamage.size() - newSeq.size() << std::endl;
#endif
        int start2 = probReadPostDamage.size() - newSeq.size();

        // Ensure matrixIndex is within the bounds of scoringMatrix
        if (matrixIndex >= scoringMatrix.size()) {
            std::cerr << "Warning: matrixIndex " << matrixIndex << " is out of bounds for scoringMatrix with size " << scoringMatrix.size() << std::endl;
            continue;
        }

        // Assuming probReadPostDamage is aligned with newSeq and has the same length
        const std::array<double, 5>& probabilities = probReadPostDamage[start2 + i];

#ifdef DEBUGDAM
        std::cerr << "Updating index " << matrixIndex << " with new base " << newSeq[i] << " and probabilities ";
        for (double p : probabilities) {
            std::cerr << p << " ";
        }
        std::cerr << std::endl;

        // Debug: print current values before updating
        std::cerr << "Current values at index " << matrixIndex << ": A=" << std::exp(scoringMatrix[matrixIndex]['A']) 
                  << ", C=" << std::exp(scoringMatrix[matrixIndex]['C'])
                  << ", G=" << std::exp(scoringMatrix[matrixIndex]['G'])
                  << ", T=" << std::exp(scoringMatrix[matrixIndex]['T']) << std::endl;

        // Debug: print values to be added
        std::cerr << "Log probabilities to be added at index " << start2 + i << ": A=" << std::log(probabilities[0])
                  << ", C=" << std::log(probabilities[1])
                  << ", G=" << std::log(probabilities[2])
                  << ", T=" << std::log(probabilities[3]) << std::endl;
#endif
        // Update the scoring matrix for all bases by adding log probabilities
        scoringMatrix[matrixIndex]['A'] = std::log(probabilities[0]);
        scoringMatrix[matrixIndex]['C'] = std::log(probabilities[1]);
        scoringMatrix[matrixIndex]['G'] = std::log(probabilities[2]);
        scoringMatrix[matrixIndex]['T'] = std::log(probabilities[3]);
        scoringMatrix[matrixIndex]['-'] = log(probabilities[4]);

#ifdef DEBUGDAM
        // Debug: print updated values
        std::cerr << "Updated values at index " << matrixIndex << ": A=" << std::exp(scoringMatrix[matrixIndex]['A'])
                  << ", C=" << std::exp(scoringMatrix[matrixIndex]['C'])
                  << ", G=" << std::exp(scoringMatrix[matrixIndex]['G'])
                  << ", T=" << std::exp(scoringMatrix[matrixIndex]['T']) << std::endl;

        // Debug: check sum of probabilities
        double sumProbabilities = std::exp(scoringMatrix[matrixIndex]['A']) + 
                                  std::exp(scoringMatrix[matrixIndex]['C']) + 
                                  std::exp(scoringMatrix[matrixIndex]['G']) + 
                                  std::exp(scoringMatrix[matrixIndex]['T']);
        std::cerr << "Sum of probabilities at index " << matrixIndex << ": " << sumProbabilities << std::endl;

        if (std::abs(sumProbabilities - 1.0) > 1e-6) {
            std::cerr << "Warning: Sum of probabilities is not 1 at index in updating function " << matrixIndex << ". Sum: " << sumProbabilities << "\n";
        }
#endif
    }

    // Check for invalid values in the updated scoring matrix
    for (const auto& map : scoringMatrix) {
        for (const auto& [base, score] : map) {
            if (std::isinf(score) || std::isnan(score)) {
                std::cerr << "Error: Invalid value in scoring matrix during extension: " << score << "\n";
                return; // Exit the function to indicate an error
            }
        }
    }
}

void assembly::updateCountMatrix(std::vector<std::unordered_map<char, int>>& countMatrix,
                                 const std::string& contig,
                                 const std::string& newSeq,
                                 size_t startIdx,
                                 const std::vector<std::array<int, 5>>& countsRead,
                                 int best_overlap_length,
                                 bool contigfirst) {
    if (contigfirst) { // we add to the end of the contig
        size_t originalSize = countMatrix.size();

        if (countMatrix.size() < contig.size()) {
            countMatrix.resize(contig.size());
            for (size_t i = originalSize; i < contig.size(); ++i) {
                countMatrix[i]['A'] = 0;
                countMatrix[i]['C'] = 0;
                countMatrix[i]['G'] = 0;
                countMatrix[i]['T'] = 0;
                countMatrix[i]['-'] = 0;
            }
        }

        if (best_overlap_length + newSeq.size() != countsRead.size()) {
            std::cerr << "Incorrect sizing of the overlap region, count matrix could be wrong" << std::endl;
        }

        // Update the countMatrix completely from overlapping until new part of the contig
        for (size_t i = 0; i < best_overlap_length + newSeq.size(); ++i) {
            size_t matrixIndex = (startIdx - best_overlap_length) + i;

            // Ensure matrixIndex is within bounds
            if (matrixIndex >= countMatrix.size()) {
                std::cerr << "Warning: matrixIndex " << matrixIndex
                          << " is out of bounds for countMatrix with size " << countMatrix.size() << std::endl;
                continue;
            }

            const std::array<int, 5>& counts = countsRead[i];
            countMatrix[matrixIndex]['A'] += counts[0];
            countMatrix[matrixIndex]['C'] += counts[1];
            countMatrix[matrixIndex]['G'] += counts[2];
            countMatrix[matrixIndex]['T'] += counts[3];
            countMatrix[matrixIndex]['-'] += counts[4];
        }
    } else {
        // We add to the beginning of the contig
        int matrixStart = static_cast<int>(startIdx) - best_overlap_length;

        // Prepend space to the beginning if needed
        if (matrixStart < 0) {
            size_t prependCount = static_cast<size_t>(-matrixStart);
            std::vector<std::unordered_map<char, int>> newEntries(prependCount);
            for (auto& m : newEntries) {
                m['A'] = 0;
                m['C'] = 0;
                m['G'] = 0;
                m['T'] = 0;
                m['-'] = 0;
            }
            countMatrix.insert(countMatrix.begin(), newEntries.begin(), newEntries.end());
            matrixStart = 0;  // Reset after prepending
        }

        // Resize at the end if needed
        size_t requiredSize = contig.size();
        if (countMatrix.size() < requiredSize) {
			cerr << "This shouldn't be called" << endl;
            size_t originalSize = countMatrix.size();
            countMatrix.resize(requiredSize);
            for (size_t i = originalSize; i < requiredSize; ++i) {
                countMatrix[i]['A'] = 0;
                countMatrix[i]['C'] = 0;
                countMatrix[i]['G'] = 0;
                countMatrix[i]['T'] = 0;
                countMatrix[i]['-'] = 0;
            }
        }

        if (best_overlap_length + newSeq.size() != countsRead.size()) {
            std::cerr << "Incorrect sizing of the overlap region, count matrix could be wrong" << std::endl;
        }

        for (size_t i = 0; i < contig.size(); ++i) {
            size_t matrixIndex = matrixStart + i;
            if (matrixIndex >= countMatrix.size()) {
                std::cerr << "Warning: matrixIndex " << matrixIndex
                          << " is out of bounds for countMatrix with size " << countMatrix.size() << std::endl;
                continue;
            }

            const std::array<int, 5>& counts = countsRead[i];
            countMatrix[matrixIndex]['A'] += counts[0];
            countMatrix[matrixIndex]['C'] += counts[1];
            countMatrix[matrixIndex]['G'] += counts[2];
            countMatrix[matrixIndex]['T'] += counts[3];
            countMatrix[matrixIndex]['-'] += counts[4];
        }
    }

    // Check for invalid values in the updated scoring matrix
    for (const auto& map : countMatrix) {
        for (const auto& [base, score] : map) {
            if (std::isinf(score) || std::isnan(score)) {
                std::cerr << "Error: Invalid value in scoring matrix during extension: " << score << "\n";
                return; // Exit the function to indicate an error
            }
        }
    }
}





void assembly::adjustFinalCut(pair<bool, int>& finalcut, const pair<bool, int>& neighborCut) {
    if (neighborCut != finalcut) {
        if (finalcut.first == false) {
            if (neighborCut.second == 2 || neighborCut.second == 3) {
                finalcut = neighborCut;
            }
        } else if (neighborCut.first == false) {
            if (finalcut.second == 2) {
                finalcut = std::make_pair(false, 0);
            } else if (finalcut.second == 3) {
                finalcut.second = 1;
            }
        } else if (finalcut.first == true) {
            if (finalcut.second != neighborCut.second) {
                if ((finalcut.second == 1 && (neighborCut.second == 2 || neighborCut.second == 3)) ||
                    (finalcut.second == 3 && neighborCut.second == 1)) {
                    finalcut.second = 3;
                } else if (finalcut.second == 2 && neighborCut.second == 1) {
                    finalcut = std::make_pair(false, 0);
                } else if (finalcut.second == 3 && neighborCut.second == 1) {
                    finalcut.second = 1;
                }
            }
        }
    }
}

size_t assembly::findAlignmentPosition(const std::string& contig, const std::string& read, int allowedMismatches) {
    if (read.size() > contig.size()) return std::string::npos;

    int matchThreshold = static_cast<int>(SIMILARITY * read.size());
    for (size_t i = 0; i <= contig.size() - read.size(); ++i) {
        int matchCount = 0;
        int mismatchCount = 0;
        for (size_t j = 0; j < read.size(); ++j) {
            if (i + j < contig.size()) {
                if (isRYMatch(read[j], contig[i + j])) {
                    matchCount++;
                } else {
                    mismatchCount++;
                    if (mismatchCount > allowedMismatches) break;
                }
            }
        }
        if (matchCount >= matchThreshold) {
            return i; // Return the starting position of the match
        }
    }
    return std::string::npos;
}




bool assembly::isSubset(const std::string& a, const std::vector<int>& aNodeIds,
                        const std::string& b, const std::vector<int>& bNodeIds,
                        int allowedMismatches) {
    if (a.size() > b.size()) return false;

    // Check if all node IDs from a are in b
    bool allNodeIdsInB = std::all_of(aNodeIds.begin(), aNodeIds.end(), [&bNodeIds](int id) {
        return std::find(bNodeIds.begin(), bNodeIds.end(), id) != bNodeIds.end();
    });

    if (allNodeIdsInB) {
        return true;
    }

    // Check if the sequence a is a subset of sequence b
    int matchThreshold = static_cast<int>(SIMILARITY * a.size());
    for (size_t i = 0; i <= b.size() - a.size(); ++i) {
        int matchCount = 0;
        int mismatchCount = 0;
        for (size_t j = 0; j < a.size(); ++j) {
            if (i + j < b.size()) {
                if (isRYMatch(a[j], b[i + j])) {
                    matchCount++;
                } else {
                    mismatchCount++;
                    if (mismatchCount > allowedMismatches) break;
                }
            }
        }
        if (matchCount >= matchThreshold) {
            return true;
        }
    }
    return false;
}




std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> assembly::mergeAllPaths(
    const GraphAss& graph, const std::map<std::pair<int, int>, OverlapResult>& overlap_map, const std::vector<frags>& reads, int min_overlap_length, bool specifiedDeam, bool minimeta, GraphData new_graph_data, int lenMin) {

    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int> >> contigs;
    std::unordered_set<int> visited;

    // Check for invalid values in all initial scoring matrices
    for (const auto& read : reads) {
        const auto& probReadPostDamage = read.probReadPostDamage;
        for (const auto& probs : probReadPostDamage) {
            for (const auto& prob : probs) {
                if (std::isinf(prob) || std::isnan(prob)) {
                    std::cerr << "Error: Invalid value in read probabilities: " << prob << "\n";
                }
            }
        }
    }

    for (int startVertex = 0; startVertex < graph.adj_list.size(); ++startVertex) {
        if (visited.find(startVertex) != visited.end()) {
            continue; // Skip already visited vertices
        }

        // Initialize the contig with the sequence of the start vertex
        std::string contig = reads[startVertex].sequence;
        std::string RYcontig = reads[startVertex].rymer_sequence;
        auto scoringMatrix = createScoringMatrix(contig, reads[startVertex].probReadPostDamage);
        auto countMatrix = createCountMatrix(contig, reads[startVertex].countsRead);
        std::vector<int> nodeIDs = reads[startVertex].nodeIDs; // Start with the initial node IDs
        std::pair<bool, int> finalcut = reads[startVertex].cutbool;

#ifdef DEBUGMERGEALL
        std::cerr << "Start seq: " << contig << "\n";
#endif

        // Use a queue for BFS
        std::queue<int> vertexQueue;
        vertexQueue.push(startVertex);
        visited.insert(startVertex);

        while (!vertexQueue.empty())
        {
            int currentVertex = vertexQueue.front();
            vertexQueue.pop();

            // Sort the overlaps for the current vertex
            auto sortedNeighbors = sortOverlapsByScore(currentVertex, graph, overlap_map);

            for (const int neighbor : sortedNeighbors) {
                if (visited.find(neighbor) != visited.end()) {
                    continue; // Skip already visited vertices
                }
                std::vector<int> commonNodeIDs = {};  // Vector to store common node IDs
                std::vector<uint64_t> nodeIDDepths;
                std::vector<uint64_t> neighborDepths;
                std::vector<uint64_t> matchingDepths;

                // Assuming 'nodeIDs' and 'reads[neighbor].nodeIDs' are vectors of integers
                for (const auto& nodeID : nodeIDs) {
                    if (std::find(reads[neighbor].nodeIDs.begin(), reads[neighbor].nodeIDs.end(), nodeID) != reads[neighbor].nodeIDs.end()) {
                        commonNodeIDs.push_back(nodeID);
                    }
                }
                bool hasCommonNodeID = false;
                bool hasCommonDepth = false;
                // Now commonNodeIDs contains all the node IDs that are common to both vectors
                if (!commonNodeIDs.empty()) {
                    hasCommonNodeID = true;  // You have at least one common node ID
                    hasCommonDepth = true;
                }

                if(minimeta){
                    cerr << "We are entering mini meta mode ..." << endl;
                    
                    for (const auto& nodeID : nodeIDs) {
                        nodeIDDepths.push_back(new_graph_data.get_depth_by_node_id(nodeID));
                    }

                    // Step 2: Retrieve depths for all nodeIDs in reads[neighbor]
                   
                    for (const auto& neighborID : reads[neighbor].nodeIDs) {
                        neighborDepths.push_back(new_graph_data.get_depth_by_node_id(neighborID));
                    }

                    // Step 3: Check for matching depths and store the matching depths
                    
                    for (const auto& depth : nodeIDDepths) {
                        if (std::find(neighborDepths.begin(), neighborDepths.end(), depth) != neighborDepths.end()) {
                            matchingDepths.push_back(depth);  // Store matching depth
                        }
                    }

                    if (!matchingDepths.empty()){
                        hasCommonDepth = true;
                    }
                }
                int first_common_id;
                int pos_in_i;
                int pos_in_j;

                if (!hasCommonNodeID && !hasCommonDepth) {
                    continue; // Skip this neighbor as there are no common node IDs or depths
                }else if(hasCommonNodeID && hasCommonDepth){
                    first_common_id = commonNodeIDs.front();
                    pos_in_i = std::find(nodeIDs.begin(), nodeIDs.end(), first_common_id) - nodeIDs.begin();
                    pos_in_j = std::find(reads[neighbor].nodeIDs.begin(), reads[neighbor].nodeIDs.end(), first_common_id) - reads[neighbor].nodeIDs.begin();
                }else if(!hasCommonNodeID && hasCommonDepth){
                    cerr << "Overlaps are searched based on node depth." << endl;
                    first_common_id = matchingDepths.front();
                    pos_in_i = std::find(nodeIDDepths.begin(), nodeIDDepths.end(), first_common_id) - nodeIDDepths.begin();
                    pos_in_j = std::find(neighborDepths.begin(), neighborDepths.end(), first_common_id) - neighborDepths.begin();
                }
                // Use positions to decide merge direction
                bool is_i_to_j = (pos_in_i < pos_in_j);

                OverlapResult best_overlap = get_overlap_length_and_score(RYcontig, reads[neighbor].rymer_sequence, min_overlap_length, lenMin);
                OverlapResult best_overlap2 = get_overlap_length_and_score(reads[neighbor].rymer_sequence, RYcontig, min_overlap_length, lenMin);

                if(!specifiedDeam){
                     best_overlap = get_overlap_length_and_score(contig, reads[neighbor].sequence, min_overlap_length, lenMin);
                     best_overlap2 = get_overlap_length_and_score(reads[neighbor].sequence, contig, min_overlap_length, lenMin);
                }
                

#ifdef DEBUGMERGEALL
                std::cerr << "best_overlap length: " << best_overlap.length << " score: " << best_overlap.score << " seq size " << reads[neighbor].sequence.size() << "\n";
                std::cerr << "best_overlap2 length: " << best_overlap2.length << " score: " << best_overlap2.score << " contig size " << contig.size() << "\n";
                std::cerr << "bool is_i_to_j " << is_i_to_j << " pos in i " << pos_in_i << " pos in j " << pos_in_j <<  "\n";
#endif

                // Check for subset if the overlap length matches the sequence length
                bool isNeighborSubset = isSubset(reads[neighbor].sequence, reads[neighbor].nodeIDs, contig, nodeIDs, 0);
                bool isContigSubset = isSubset(contig, nodeIDs, reads[neighbor].sequence, reads[neighbor].nodeIDs, 0);
				
				// subset check
                if ((best_overlap.length == reads[neighbor].sequence.size() || best_overlap.length == contig.size()) &&
                    (isNeighborSubset || isContigSubset)) {

                    if (isContigSubset) {
						// back up before overwriting
						
						std::string oldContig = contig;
						
						auto oldCountMatrix = countMatrix;
						 // Replace contig with the read since contig is a subset of the read
                        contig = reads[neighbor].sequence;
                        RYcontig = reads[neighbor].rymer_sequence;
                        nodeIDs = reads[neighbor].nodeIDs;
						
                        scoringMatrix = createScoringMatrix(contig, reads[neighbor].probReadPostDamage); // Use the new scoring matrix
                        countMatrix = createCountMatrix(contig, reads[neighbor].countsRead); // create a new count matrix and add a +1 for every base in the read (new contig)

						size_t offset = findAlignmentPosition(reads[neighbor].sequence, oldContig, 1);

#ifdef DEBUGMERGEALL
						std::cerr << "Replacing contig with read " << reads[neighbor].sequence << " because contig is a subset of the read " << oldContig << endl;
						std::cerr << "alignmentPosition " << offset << "\n";
#endif

						if (offset != std::string::npos) {
					        // Add old contig support to countMatrix at the right offset
					        for (size_t i = 0; i < oldContig.size(); ++i) {
					            for (const auto& [base, count] : oldCountMatrix[i]) {
					                countMatrix[offset + i][base] += count;
#ifdef DEBUGMERGEALL
                                    std::cerr << "Updating countMatrix at contig position " << offset + i
                                              << " with base " << base
                                              << " from read " << oldContig
										      << " with count " << count << "\n";
#endif
					            }
					        }
								

                        } else {
#ifdef DEBUGMERGEALL
                            std::cerr << "Warning: Could not find alignment position for subset contig " << oldContig << ".\n";
#endif
                        }
                        // Insert only new node IDs
                        for (const auto& id : reads[neighbor].nodeIDs) {
                            if (std::find(nodeIDs.begin(), nodeIDs.end(), id) == nodeIDs.end()) {
                                nodeIDs.push_back(id);
                            }
                        }


                    } else {
                        // Find alignment position using rymer_sequence
                        size_t alignmentPosition = findAlignmentPosition(contig, reads[neighbor].sequence, 1);

#ifdef DEBUGMERGEALL
                        std::cerr << "alignmentPosition " << alignmentPosition << "\n";
#endif

                        if (alignmentPosition != std::string::npos) {
                            for (size_t i = 0; i < reads[neighbor].sequence.size(); ++i) {
                                size_t contigIndex = alignmentPosition + i;
                                if (contigIndex < contig.size()) {
                                    countMatrix[contigIndex][reads[neighbor].sequence[i]]++;
#ifdef DEBUGMERGEALL
                                    std::cerr << "Updating countMatrix at contig position " << contigIndex
                                              << " with base " << reads[neighbor].sequence[i]
                                              << " from read " << neighbor << "\n";
#endif
                                } else {
#ifdef DEBUGMERGEALL
                                    std::cerr << "Warning: Subset read " << neighbor
                                              << " has longer sequence than contig at index " << contigIndex << ".\n";
#endif
                                }
                            }
                        } else {
#ifdef DEBUGMERGEALL
                            std::cerr << "Warning: Could not find alignment position for subset read " << neighbor << ".\n";
#endif
                        }
                        for (const auto& id : reads[neighbor].nodeIDs) {
                            if (std::find(nodeIDs.begin(), nodeIDs.end(), id) == nodeIDs.end()) {
                                nodeIDs.push_back(id);
                            }
                        }
                    }

#ifdef DEBUGMERGEALL
                    // Print debug information
                    std::cerr << "Excluding read " << neighbor << " as it is a subset of the contig.\n";
                    std::cerr << "Contig: " << contig << "\n";
                    std::cerr << "Contig Node IDs: ";
                    for (const auto& id : nodeIDs) {
                        std::cerr << id << " ";
                    }
                    std::cerr << "\nRead: " << reads[neighbor].sequence << "\n";
                    std::cerr << "Read Node IDs: ";
                    for (const auto& id : reads[neighbor].nodeIDs) {
                        std::cerr << id << " ";
                    }
                    std::cerr << "\n";
#endif

                    visited.insert(neighbor);
                    continue; // Continue to the next neighbor after handling the subset
                } // subset check ends - 
				
				// Go on with overlaps (not subset)
                // Ignore reads with zero overlap lengths and scores, even if they have common node IDs
                if (best_overlap.length == 0 && best_overlap.score == 0 &&
                    best_overlap2.length == 0 && best_overlap2.score == 0 &&
                    hasCommonNodeID) {
#ifdef DEBUGMERGEALL
                    std::cerr << "Ignoring read " << neighbor << " due to zero overlap lengths and scores, despite common node IDs.\n";
#endif
                    continue;
                }

                // Ensure that the overlap length is not longer than the sequence lengths
                if (best_overlap.length > contig.size() || best_overlap.length > reads[neighbor].sequence.size()) {
                    continue;
                }
                if (best_overlap2.length > contig.size() || best_overlap2.length > reads[neighbor].sequence.size()) {
                    continue;
                }

                bool gowithoverlap2 = false;
                if(pos_in_i > pos_in_j){
                    gowithoverlap2 = false;
                }else if(pos_in_i < pos_in_j){
                    gowithoverlap2 = true;
                }else if(best_overlap2.score > best_overlap.score && best_overlap2.length > best_overlap.length && pos_in_i == pos_in_j){
                    gowithoverlap2 = true;
                }else if(best_overlap2.score > best_overlap.score || best_overlap2.length > best_overlap.length && pos_in_i == pos_in_j){
                    gowithoverlap2 = true;
                }
                // If best_overlap2 has a higher score, merge the read first, then the contig
                if (gowithoverlap2) {
                //if(is_i_to_j){
                    if (reads[neighbor].nodeIDs.front() != nodeIDs.front()) {
                        continue; // Skip this neighbor as the first node IDs do not match
						
                    }

                   // Check if all bases in the overlap region match with the specified conditions
                    bool ryMatch = true;
                    for (size_t i = 0; i < best_overlap2.length; ++i) {
                        if (i < 5 || i >= best_overlap2.length - 5) { // First and last 5 bases
                            if (!isRYMatch(contig[contig.size() - best_overlap2.length + i], reads[neighbor].sequence[i])) {
                                ryMatch = false;
                                break;
                            }
                        } else { // Middle bases must be a perfect match
                            if (contig[contig.size() - best_overlap2.length + i] != reads[neighbor].sequence[i]) {
                                ryMatch = false;
                                break;
                            }
                        }
                    }

                    if (!ryMatch) {
                        continue; // Skip this neighbor as the overlap in the RY sequence is not consistent
                    }

                    // Merge the non-overlapping part of the contig into the read
                    std::string newSeqPart = contig.substr(best_overlap2.length);
                    std::string newRYseq = RYcontig.substr(best_overlap2.length);
                    std::string updatedSequence = reads[neighbor].sequence + newSeqPart;
                    std::string updatedRYsequence = reads[neighbor].rymer_sequence + newRYseq;
#ifdef DEBUGMERGEALL
                    
                    std::cerr << "adding read first new seq part " << newSeqPart << " new rymer part " << newRYseq << "\n";
                    std::cerr << "read seq " << contig << " read ryseq " << RYcontig << "\n";
                    std::cerr << "new seq " << updatedSequence << "\n";
                    std::cerr << "new ry seq " << updatedRYsequence << "\n";
#endif
					

					// After that, you should handle updating the scoring matrix and any further calculations
					updateScoringMatrix(scoringMatrix, updatedSequence, newSeqPart, reads[neighbor].sequence.size(), reads[neighbor].probReadPostDamage);
					updateCountMatrix(countMatrix, updatedSequence, newSeqPart, reads[neighbor].sequence.size(), reads[neighbor].countsRead, best_overlap2.length, false);
					

                    if (scoringMatrix.size() != updatedSequence.size()) {
                        std::cerr << "Error: Scoring matrix and read size mismatch after merging.\n";
                        std::cerr << "Scoring matrix size: " << scoringMatrix.size() << " Read size: " << updatedSequence.size() << "\n";
                    }
                    if (countMatrix.size() != updatedSequence.size()) {
                        std::cerr << "Error: Count matrix and read size mismatch after merging.\n";
                        std::cerr << "Scoring matrix size: " << countMatrix.size() << " Read size: " << updatedSequence.size() << "\n";
                    }

                    // Update contig and RYcontig to reflect the merged sequence
                    contig = updatedSequence;
                    RYcontig = updatedRYsequence;

                    if(contig.size() != RYcontig.size()){
                        cerr << "error merging of reads wrong forward" << endl;
                    }

                    // Adjust finalcut based on the neighbor's cutbool
                    adjustFinalCut(finalcut, reads[neighbor].cutbool);
                    
                    // Insert only new node IDs
                    for (const auto& id : reads[neighbor].nodeIDs) {
                        if (std::find(nodeIDs.begin(), nodeIDs.end(), id) == nodeIDs.end()) {
                            nodeIDs.push_back(id);
                        }
                    }

                } else
				{ // Merge the read into the contig using best_overlap

				// Check if all bases in the overlap region match with the specified conditions
				    bool ryMatch = true;
				    for (size_t i = 0; i < best_overlap.length; ++i) {
				        if (i < 5 || i >= best_overlap.length - 5) { // First and last 5 bases
				            if (!isRYMatch(contig[contig.size() - best_overlap.length + i], reads[neighbor].sequence[i])) {
				                ryMatch = false;
				                break;
				            }
				        } else { // Middle bases must be a perfect match
				            if (contig[contig.size() - best_overlap.length + i] != reads[neighbor].sequence[i]) {
				                ryMatch = false;
				                break;
				            }
				        }
				    }

				    if (!ryMatch) {
				        continue; // Skip this neighbor as the overlap in the RY sequence is not consistent
				    }
					 
				    // Merge the non-overlapping part of the next read into the contig
				    std::string newSeqPart = reads[neighbor].sequence.substr(best_overlap.length);
				    std::string newRYseq = reads[neighbor].rymer_sequence.substr(best_overlap.length);
				    size_t startIdx = contig.size();
				    contig += newSeqPart;
				    RYcontig += newRYseq;

#ifdef DEBUGMERGEALL
				    if(contig.size() != RYcontig.size()){
				        cerr << "error merging of reads wrong forward" << endl;
				    }
				    std::cerr << "adding contig first new seq part " << newSeqPart << " new rymer part " << newRYseq << "\n";
				    std::cerr << "read seq " << reads[neighbor].sequence << " read ryseq " << reads[neighbor].rymer_sequence << "\n";
				    std::cerr << "new seq " << contig << "\n";
				    std::cerr << "new ry seq " << RYcontig << "\n";
#endif

				    // Update the scoring matrix with the new sequence part
				    updateScoringMatrix(scoringMatrix, contig, newSeqPart, startIdx, reads[neighbor].probReadPostDamage);
    
				    updateCountMatrix(countMatrix, contig, newSeqPart, startIdx, reads[neighbor].countsRead, best_overlap.length, true);

				    // Check consistency of scoring matrix and contig size
				    if (scoringMatrix.size() != contig.size()) {
				        std::cerr << "Error: Scoring matrix and contig size mismatch after merging.\n";
				        std::cerr << "Scoring matrix size: " << scoringMatrix.size() << " Contig size: " << contig.size() << "\n";
				    }

				    // Adjust finalcut based on the neighbor's cutbool
				    adjustFinalCut(finalcut, reads[neighbor].cutbool);

				    // Insert only new node IDs
				    for (const auto& id : reads[neighbor].nodeIDs) {
				        if (std::find(nodeIDs.begin(), nodeIDs.end(), id) == nodeIDs.end()) {
				            nodeIDs.push_back(id);
				        }
				    }
				}

				visited.insert(neighbor);
				vertexQueue.push(neighbor);

				}

				if(contig.size() != RYcontig.size()){
				    cerr << "error merging of reads wrong middle" << endl;
				    cerr << contig << " " << RYcontig << endl;
				}

			}
			if(contig.size() != RYcontig.size()){
			    cerr << "error merging of reads wrong end" << endl;
			    cerr << contig << " " << RYcontig << endl;
			}
			
			for (size_t i = 0; i < contig.size(); ++i) {
			    const auto& counts = countMatrix[i];
			    int maxCount = -1;
			    std::vector<char> maxBases;

			    // Find base(s) with the highest count
			    for (const auto& [base, count] : counts) {
			        if (count > maxCount) {
			            maxCount = count;
			            maxBases = {base};
			        } else if (count == maxCount) {
			            maxBases.push_back(base);
			        }
			    }

			    char originalBase = contig[i];

			    // If current base is not among the top ones, replace and print matrix
			    if (std::find(maxBases.begin(), maxBases.end(), originalBase) == maxBases.end()) {
			        contig[i] = maxBases[0];
#ifdef DEBUGMERGEALL
			        // Print full countMatrix for this position
			        std::cerr << "Corrected position " << i << " from '" << originalBase << "' to '" << contig[i] << "'\n";
			        std::cerr << "Count matrix at position " << i << ":\n";
			        for (const auto& [base, count] : counts) {
			        std::cerr << "  Base " << base << ": " << count << "\n";

			        }
#endif
			    }
			}
			

		contigs.push_back({contig, scoringMatrix, nodeIDs, RYcontig, countMatrix, finalcut});
                    
    }

    return contigs;
}



std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>>
assembly::removeSubsetContigs(const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> & contigs) {
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> filteredContigs;
    std::unordered_set<size_t> toRemove;

    for (size_t i = 0; i < contigs.size(); ++i) {
        if (toRemove.find(i) != toRemove.end()) {
            continue; // Skip already marked contigs
        }

        std::string seq_i = std::get<0>(contigs[i]);
        std::vector<int> nodeIds_i = std::get<2>(contigs[i]);

        for (size_t j = 0; j < contigs.size(); ++j) {
            if (i != j && toRemove.find(j) == toRemove.end()) {
                std::string seq_j = std::get<0>(contigs[j]);
                std::vector<int> nodeIds_j = std::get<2>(contigs[j]);

                bool isISubsetOfJ = isSubset(seq_i, nodeIds_i, seq_j, nodeIds_j, 2);
                bool isJSubsetOfI = isSubset(seq_j, nodeIds_j, seq_i, nodeIds_i, 2);

                if (isISubsetOfJ || isJSubsetOfI) {
                    if (isISubsetOfJ) {
                        toRemove.insert(i);
#ifdef DEBUGASS
                        std::cerr << "Contig " << i << " is a subset of Contig " << j << " (i in j).\n";
                        std::cerr << "seq_i: " << seq_i << " seq_j: " << seq_j << "\n";
#endif
                        break;
                    }
                    if (isJSubsetOfI) {
                        toRemove.insert(j);
#ifdef DEBUGASS
                        std::cerr << "Contig " << j << " is a subset of Contig " << i << " (j in i).\n";
                        std::cerr << "seq_i: " << seq_i << " seq_j: " << seq_j << "\n";
#endif
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < contigs.size(); ++i) {
        if (toRemove.find(i) == toRemove.end()) {
            filteredContigs.push_back(contigs[i]);
#ifdef DEBUGASS
            std::cerr << "Adding Contig " << i << "\n";
#endif
        } else {
#ifdef DEBUGASS
            std::cerr << "Excluding Contig " << i << " from the results.\n";
#endif
        }
    }

    return filteredContigs;
}


std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>> assembly::mergeContigs(
const std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>& contig1,
const std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>& contig2,
int overlapLength, bool isOverlapAtEnd, string mode)
{
    std::string seq1 = std::get<0>(contig1);
    std::string seq2 = std::get<0>(contig2);
    std::vector<std::unordered_map<char, double>> scoringMatrix1 = std::get<1>(contig1);
    std::vector<std::unordered_map<char, double>> scoringMatrix2 = std::get<1>(contig2);
    std::vector<int> nodeIDs1 = std::get<2>(contig1);
    std::vector<int> nodeIDs2 = std::get<2>(contig2);
    std::string ry1 = std::get<3>(contig1);
    std::string ry2 = std::get<3>(contig2);
    std::vector<std::unordered_map<char, int>> countMatrix1 = std::get<4>(contig1);
    std::vector<std::unordered_map<char, int>> countMatrix2 = std::get<4>(contig2);
    pair<bool, int> cutbool1 = get<5>(contig1);
    pair<bool, int> cutbool2 = get<5>(contig2);

    std::string mergedSeq;
    std::string mergedRySeq;
    std::vector<std::unordered_map<char, double>> mergedScoringMatrix;
    std::vector<int> mergedNodeIDs;
    std::vector<std::unordered_map<char, int>> mergedCountMatrix;
    pair<bool, int> finalcut = make_pair(false, 0); 

    if(seq1.size() != ry1.size() || seq2.size() != ry2.size()){
        cerr << "error different lengths between sequence and rymer sequence" << endl;

    } 
#ifdef DEBUGASS2	
	// Check consistency between sequence and count matrix for contig1
	if (seq1.size() != countMatrix1.size()) {
	    std::cerr << "Error: Sequence and count matrix size mismatch for contig1\n";
	} else {
	    for (size_t i = 0; i < seq1.size(); ++i) {
	        char expectedBase = seq1[i];
	        char mostFrequentBase = '\0';
	        int maxCount = -1;
	        int secondMaxCount = -1;

	        // Find the base with the highest count, and ensure it's distinct from others
	        for (const auto& [base, count] : countMatrix1[i]) {
	            if (count > maxCount) {
	                secondMaxCount = maxCount;
	                maxCount = count;
	                mostFrequentBase = base;
	            } else if (count > secondMaxCount) {
	                secondMaxCount = count;
	            }
	        }

	        // Only change the base if there's a significantly higher count for another base
	        if (mostFrequentBase != expectedBase && maxCount > secondMaxCount) {
	            // Swap the base in the sequence
	            seq1[i] = mostFrequentBase;
	            // std::cerr << "Updated base at position " << i
// 	                      << " in contig1: sequence had '" << expectedBase
// 	                      << "', but count matrix suggests '" << mostFrequentBase
// 	                      << "' with count = " << maxCount << "\n";
	        }
	    }
	}

	// Check consistency between sequence and count matrix for contig2
	if (seq2.size() != countMatrix2.size()) {
	    std::cerr << "Error: Sequence and count matrix size mismatch for contig2\n";
	} else {
	    for (size_t i = 0; i < seq2.size(); ++i) {
	        char expectedBase = seq2[i];
	        char mostFrequentBase = '\0';
	        int maxCount = -1;
	        int secondMaxCount = -1;

	        // Find the base with the highest count, and ensure it's distinct from others
	        for (const auto& [base, count] : countMatrix2[i]) {
	            if (count > maxCount) {
	                secondMaxCount = maxCount;
	                maxCount = count;
	                mostFrequentBase = base;
	            } else if (count > secondMaxCount) {
	                secondMaxCount = count;
	            }
	        }

	        // Only change the base if there's a significantly higher count for another base
	        if (mostFrequentBase != expectedBase && maxCount > secondMaxCount) {
	            // Swap the base in the sequence
	            seq2[i] = mostFrequentBase;
	            std::cerr << "Updated base at position " << i
	                      << " in contig2: sequence had '" << expectedBase
	                      << "', but count matrix suggests '" << mostFrequentBase
	                      << "' with count = " << maxCount << "\n";
	        }
	    }
	}
	
#endif
    if (isOverlapAtEnd)
    {
        if (overlapLength > seq2.size()) {
            std::cerr << "Error: Invalid overlap length at the end for sequences. Overlap length: " << overlapLength << ", seq1 size: " << seq1.size() << ", seq2 size: " << seq2.size() << "\n";
            return {seq1, scoringMatrix1, nodeIDs1, ry1, countMatrix1, cutbool1};
        }
        //cerr << "overlap length " << overlapLength << ", seq1 size: " << seq1.size() << ", seq2 size: " << seq2.size() << "\n";
        // Initialize merged matrices with the appropriate sizes
        mergedScoringMatrix.resize(seq1.size() + seq2.size() - overlapLength);
        mergedCountMatrix.resize(seq1.size() + seq2.size() - overlapLength);
        mergedSeq = merge_sequences(seq1, seq2, overlapLength);
        mergedRySeq = merge_sequences(ry1, ry2, overlapLength);
        mergedScoringMatrix = scoringMatrix1;
        mergedScoringMatrix.resize(mergedSeq.size());
        mergedCountMatrix = countMatrix1;
        mergedCountMatrix.resize(mergedSeq.size());

        int seq1Start = seq1.size() - overlapLength;
		// go through each base of the overlap to decide what to do. 
        for (int ov = 0; ov < overlapLength; ++ov) {
            int seq1Index = seq1Start + ov;
            int seq2Index = ov;
			// ANY MISMATCHING BASES IN THE OVERLAP
            if (seq1[seq1Index] != seq2[seq2Index])
			{
#ifdef DEBUGASS2
                std::cerr << "Mismatch between seq1[" << seq1Index << "]=" << seq1[seq1Index] << " and seq2[" << seq2Index << "]=" << seq2[seq2Index] << "\n";
                std::cerr << "Scoring matrix for seq1[" << seq1Index << "]: ";
                for (const auto& [base, score] : scoringMatrix1[seq1Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\nScoring matrix for seq2[" << seq2Index << "]: ";
                for (const auto& [base, score] : scoringMatrix2[seq2Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\n";

                std::cerr << "Count matrix for seq1[" << seq1Index << "]: ";
                for (const auto& [base, score] : countMatrix1[seq1Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                cerr << "\n";
                std::cerr << "Count matrix for seq2[" << seq2Index << "]: ";
                for (const auto& [base, score] : countMatrix2[seq2Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\n";
#endif
				// SPECIFIC CASE FOR GAPS - ALWAYS MAJORITY RULE
                if(seq1[seq1Index] == '-' || seq2[seq2Index] == '-'){
                   
                   if((mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]]) < (mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]])){
                        mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
                        mergedSeq[seq1Index] = seq2[seq2Index];
                        mergedRySeq[seq1Index] = ry2[seq2Index];
                   }
				   continue;
                }
				// TRANSITIONS
                if(seq1[seq1Index] == 'C' && seq2[seq2Index] == 'T' || seq1[seq1Index] == 'T' && seq2[seq2Index] == 'C' || seq1[seq1Index] == 'G' && seq2[seq2Index] == 'A' || seq1[seq1Index] == 'A' && seq2[seq2Index] == 'G')
		        {
			    if (mode == "reckless")
			    {
			        if ((mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]]) < (mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]]))
				{		
					

			            // reckless mode, everything is handle based on majority rule - the default is seq1 we only need to change to seq2 when the count is higher.
			            mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
			            mergedSeq[seq1Index] = seq2[seq2Index];
			            mergedRySeq[seq1Index] = ry2[seq2Index];
				}
				
			    } else if(mode == "strict"){
			        double countSeq1 = mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]];
			        double countSeq2 = mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]];

			        double total = countSeq1 + countSeq2;
			        if (countSeq1 / total >= 0.9) {
						continue;
  
             
		            } else if (countSeq2 / total >= 0.9) {
		            
		                mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
		                mergedSeq[seq1Index] = seq2[seq2Index];
		                mergedRySeq[seq1Index] = ry2[seq2Index];
	        
		      	    }else{
						
		        	    mergedSeq[seq1Index] = 'N';
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
		                mergedScoringMatrix[seq1Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
		                mergedRySeq[seq1Index] = 'N';
		      	    }

			    }else{
       
			        double countSeq1 = mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]];
			        double countSeq2 = mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]];

			        double total = countSeq1 + countSeq2;
			        if (countSeq1 / total >= 0.65) {
						continue;
  
             
		            } else if (countSeq2 / total >= 0.65) {
		            
		                mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
		                mergedSeq[seq1Index] = seq2[seq2Index];
		                mergedRySeq[seq1Index] = ry2[seq2Index];
	        
		      	    }else{
						
		        	    mergedSeq[seq1Index] = 'N';
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
						}
		                mergedScoringMatrix[seq1Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
		                mergedRySeq[seq1Index] = 'N';
		      	    }
			    }    
				
		     
		        } else{
			// TRANSVERSIONS
				    if (mode == "reckless")
				    {
						if ((mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]]) < (mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]]))
						{
							
			
					            // reckless mode, everything is handle based on majority rule - the default is seq1 we only need to change to seq2 when the count is higher.
					            mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
								for (const char base : {'A', 'C', 'G', 'T', '-'}) {
								    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
								}
					            mergedSeq[seq1Index] = seq2[seq2Index];
					            mergedRySeq[seq1Index] = ry2[seq2Index];
						}
			
				    } else if(mode == "strict"){
				    // In strict mode transversions are being masked as Ns. To avoid any missincopration based on ancient damage
				    // This is the normal mode: transistion have to have a 0.75 majority to be accepted
				        double countSeq1 = mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]];
				    	double countSeq2 = mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]];
		
						double total = countSeq1 + countSeq2;
                                
						if (countSeq2 / total >= 0.9) {

				                
			        	    mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
							}
			        	    mergedSeq[seq1Index] = seq2[seq2Index];
			        	    mergedRySeq[seq1Index] = ry2[seq2Index];
				    	} else if(countSeq1 / total >= 0.9){
							continue;
						}else{
				        	
			        	    mergedSeq[seq1Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
							}
			                mergedScoringMatrix[seq1Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
			                mergedRySeq[seq1Index] = 'N';
			    	
						}
				    }else{
				
				        double countSeq1 = mergedCountMatrix[seq1Index][seq1[seq1Index]] + countMatrix2[seq2Index][seq1[seq1Index]];
				    	double countSeq2 = mergedCountMatrix[seq1Index][seq2[seq2Index]] + countMatrix2[seq2Index][seq2[seq2Index]];
		
						double total = countSeq1 + countSeq2;
						if (countSeq2 / total >= 0.65) {

				            mergedScoringMatrix[seq1Index] = scoringMatrix2[seq2Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
							}
				            mergedSeq[seq1Index] = seq2[seq2Index];
				            mergedRySeq[seq1Index] = ry2[seq2Index];
			
				        } else if (countSeq1 / total >= 0.65){
				    	    continue;
				        }else{
				            mergedSeq[seq1Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq1Index][base] += countMatrix2[seq2Index][base];
							}
							
				            mergedScoringMatrix[seq1Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
				            mergedRySeq[seq1Index] = 'N';
				        }
				    }
				}	
				// MATCH we do nothing because it stays the same!	
		    }
		}

                
        int rest = seq2.size() - overlapLength;
        for (int r = 0; r < rest; ++r) {
            int index = r + overlapLength + seq1Start;
            if (index >= mergedScoringMatrix.size() || index >= mergedCountMatrix.size() || r + overlapLength >= scoringMatrix2.size() || r + overlapLength >= countMatrix2.size()) {
                std::cerr << "Error: Index out of bounds during matrix update\n";
                continue;
            }
            mergedScoringMatrix[index] = scoringMatrix2[r + overlapLength];
            mergedCountMatrix[index] = countMatrix2[r + overlapLength];
        }

        if (mergedScoringMatrix.size() != mergedSeq.size() || mergedCountMatrix.size() != mergedSeq.size()) {
            std::cerr << "Error: scoring matrix or count matrix and sequence are not the same size.\n";
            std::cerr << "Sequence size: " << mergedSeq.size() << ", Scoring matrix size: " << mergedScoringMatrix.size() << ", Count matrix size: " << mergedCountMatrix.size() <<"\n";
        }

        std::unordered_set<int> seenNodeIDs(nodeIDs1.begin(), nodeIDs1.end());
        mergedNodeIDs = nodeIDs1;
        for (const auto& id : nodeIDs2) {
            if (seenNodeIDs.find(id) == seenNodeIDs.end()) {
                mergedNodeIDs.push_back(id);
                seenNodeIDs.insert(id);
            }
        }

        if (cutbool2 != cutbool1){
            if(cutbool1.first == false){
                if(cutbool2.second == 2 || cutbool2.second == 3){
                    cutbool1 = cutbool2;
                }
                
            }
            else if(cutbool2.first == false){
                if(cutbool1.second == 2){
                    cutbool1 = make_pair(false, 0);
                }else if(cutbool1.second == 3){
                    cutbool1.second = 1;
                }

            }
            else if(cutbool1.first == true){
                if(cutbool1.second != cutbool2.second){
                    if(cutbool1.second == 1 && cutbool2.second == 2 || cutbool1.second == 1 && cutbool2.second == 3){
                        cutbool1.second = 3;
                    }else if (cutbool1.second == 2 && cutbool2.second == 1){
                        cutbool1 = make_pair(false, 0);
                    }else if (cutbool1.second == 3 && cutbool2.second == 1){
                        cutbool1.second = 1;
                    }
                }
            }
        }
        finalcut = cutbool1;
	// OVERLAP IS THE OTHER WAY AROUND - SENARIO 2 
    } else {
        if (overlapLength > seq1.size()) {
            std::cerr << "Error: Invalid overlap length at the beginning for sequences. Overlap length: " << overlapLength << ", seq1 size: " << seq1.size() << ", seq2 size: " << seq2.size() << "\n";
            return {seq2, scoringMatrix2, nodeIDs2, ry2, countMatrix2, cutbool2};
        }
        // Initialize merged matrices with the appropriate sizes
        mergedScoringMatrix.resize(seq2.size() + seq1.size() - overlapLength);
        mergedCountMatrix.resize(seq2.size() + seq1.size() - overlapLength);

        mergedSeq = merge_sequences(seq2, seq1, overlapLength);
        mergedRySeq = merge_sequences(ry2, ry1, overlapLength);
        mergedScoringMatrix = scoringMatrix2;
        mergedScoringMatrix.resize(mergedSeq.size());
        mergedCountMatrix = countMatrix2;
        mergedCountMatrix.resize(mergedSeq.size());

        int seq2Start = seq2.size() - overlapLength;
        for (int ov = 0; ov < overlapLength; ++ov) {
            int seq2Index = seq2Start + ov;
            int seq1Index = ov;


            if (seq2[seq2Index] != seq1[seq1Index]) {
#ifdef DEBUGASS2
                std::cerr << "Mismatch between seq2[" << seq2Index << "]=" << seq2[seq2Index] << " and seq1[" << seq1Index << "]=" << seq1[seq1Index] << "\n";
                std::cerr << "Scoring matrix for seq2[" << seq2Index << "]: ";
                for (const auto& [base, score] : scoringMatrix2[seq2Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\nScoring matrix for seq1[" << seq1Index << "]: ";
                for (const auto& [base, score] : scoringMatrix1[seq1Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\n";

                std::cerr << "Count matrix for seq1[" << seq1Index << "]: ";
                for (const auto& [base, score] : countMatrix1[seq1Index]) {
                    std::cerr << base << ": " << score << " ";
                }

                std::cerr << "Count matrix for seq2[" << seq2Index << "]: ";
                for (const auto& [base, score] : countMatrix2[seq2Index]) {
                    std::cerr << base << ": " << score << " ";
                }
                std::cerr << "\n";
                if(seq1[seq1Index] == '-' || seq2[seq2Index] == '-'){
                    //cerr << "Error we have a dash" << endl;
                }
                std::cerr << "\n";
#endif
                if(seq1[seq1Index] == '-' || seq2[seq2Index] == '-'){
                   
                   if((mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]]) < (countMatrix1[seq1Index][seq1[seq1Index]] + mergedCountMatrix[seq2Index][seq1[seq1Index]])){
                        mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
						for (const char base : {'A', 'C', 'G', 'T', '-'}) {
						    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
						}
						
                        mergedSeq[seq2Index] = seq2[seq1Index];
                        mergedRySeq[seq2Index] = ry2[seq1Index];
                   }
				   continue;
                }

				// TRANSITIONS
                if(seq1[seq1Index] == 'C' && seq2[seq2Index] == 'T' || seq1[seq1Index] == 'T' && seq2[seq2Index] == 'C' || seq1[seq1Index] == 'G' && seq2[seq2Index] == 'A' || seq1[seq1Index] == 'A' && seq2[seq2Index] == 'G')
					
				{				
					
					if (mode == "reckless")
					{
						if ((mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]]) < (mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]]))
							
						{
				
						
		        			// reckless mode, everything is handle based on majority rule - the default is seq1 we only need to change to seq2 when the count is higher.
		        			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
		        			mergedSeq[seq2Index] = seq1[seq1Index];
		        			mergedRySeq[seq2Index] = ry1[seq1Index];
						}
						
				    } else if(mode == "strict"){
						// This is the normal mode: transistion have to have a 0.75 majority to be accepted
		    			double countSeq1 = mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]];
		    			double countSeq2 = mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]];

						double total = countSeq1 + countSeq2;
						if (countSeq2 / total >= 0.9) {
              
	            			cerr << "switching " << countSeq2 / total << endl;
	            			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	            			mergedSeq[seq2Index] = seq1[seq1Index];
	            			mergedRySeq[seq2Index] = ry1[seq1Index];
	        			} else if(countSeq1 / total >= 0.9) {
	    			        continue;
	        			}else{
	    			  
	            			mergedSeq[seq2Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	                        mergedScoringMatrix[seq2Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
	                        mergedRySeq[seq2Index] = 'N';
	        			}
					}else{
						// This is the normal mode: transistion have to have a 0.75 majority to be accepted
		    			double countSeq1 = mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]];
		    			double countSeq2 = mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]];

						double total = countSeq1 + countSeq2;
						if (countSeq2 / total >= 0.65) {
              
	            			cerr << "switching " << countSeq2 / total << endl;
	            			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	            			mergedSeq[seq2Index] = seq1[seq1Index];
	            			mergedRySeq[seq2Index] = ry1[seq1Index];
	        			} else if(countSeq1 / total >= 0.65) {
	    			        continue;
	        			}else{
	    			  
	            			mergedSeq[seq2Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	                        mergedScoringMatrix[seq2Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
	                        mergedRySeq[seq2Index] = 'N';
	        			}
					}
						
				     
				} else{
					// TRANSVERSIONS
					if (mode == "reckless")
					{
						if ((mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]]) < (mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]]))
							
						{
		        			// reckless mode, everything is handle based on majority rule - the default is seq2 we only need to change to seq1 when the count is higher.
		        			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
		        			mergedSeq[seq2Index] = seq1[seq1Index];
		        			mergedRySeq[seq2Index] = ry1[seq1Index];
						}
						
				    } else if(mode == "strict"){
						// In strict mode transversions are being masked as Ns. To avoid any missincopration based on ancient damage
						// This is the normal mode: transistion have to have a 0.75 majority to be accepted
		    			double countSeq1 = mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]];
		    			double countSeq2 = mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]];
						
						double total = countSeq1 + countSeq2;
						if (countSeq2 / total >= 0.9) {
	            			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	            			mergedSeq[seq2Index] = seq1[seq1Index];
	            			mergedRySeq[seq2Index] = ry1[seq1Index];
	        			} else if (countSeq1 / total >= 0.9) {
							continue;
	        			}else{
	            			mergedSeq[seq2Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	                        mergedScoringMatrix[seq2Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
	                        mergedRySeq[seq2Index] = 'N';
	        			}
					}else{
		    			double countSeq1 = mergedCountMatrix[seq2Index][seq2[seq2Index]] + countMatrix1[seq1Index][seq2[seq2Index]];
		    			double countSeq2 = mergedCountMatrix[seq2Index][seq1[seq1Index]] + countMatrix1[seq1Index][seq1[seq1Index]];
						
						double total = countSeq1 + countSeq2;
						if (countSeq2 / total >= 0.65) {
						
	            			mergedScoringMatrix[seq2Index] = scoringMatrix1[seq1Index];
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	            			mergedSeq[seq2Index] = seq1[seq1Index];
	            			mergedRySeq[seq2Index] = ry1[seq1Index];
	        			} else if (countSeq1 / total >= 0.65) {
	    			        
							continue;
	        			}else{
							//cerr << "Masking " << countSeq2 / total << endl;
	            			mergedSeq[seq2Index] = 'N';
							for (const char base : {'A', 'C', 'G', 'T', '-'}) {
							    mergedCountMatrix[seq2Index][base] += countMatrix1[seq1Index][base];
							}
							
	                        mergedScoringMatrix[seq2Index] = {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}};
	                        mergedRySeq[seq2Index] = 'N';
	        			}
					}
				}	
			// MATCH we do nothing because it stays the same!	
			}
		}


        int rest = seq1.size() - overlapLength;
        for (int r = 0; r < rest; ++r) {
            int index = r + overlapLength + seq2Start;
            if (index >= mergedScoringMatrix.size() || index >= mergedCountMatrix.size() || r + overlapLength >= scoringMatrix1.size() || r + overlapLength >= countMatrix1.size()) {
                std::cerr << "Error: Index out of bounds during matrix update\n";
                continue;
            }
            mergedScoringMatrix[index] = scoringMatrix1[r + overlapLength];
            mergedCountMatrix[index] = countMatrix1[r + overlapLength];
        }

        if (mergedScoringMatrix.size() != mergedSeq.size() || mergedCountMatrix.size() != mergedSeq.size()) {
            std::cerr << "Error: scoring matrix or count matrix and sequence are not the same size.\n";
            std::cerr << "Sequence size: " << mergedSeq.size() << ", Scoring matrix size: " << mergedScoringMatrix.size() << ", Count matrix size: " << mergedCountMatrix.size() <<"\n";
        }

        std::unordered_set<int> seenNodeIDs(nodeIDs2.begin(), nodeIDs2.end());
        mergedNodeIDs = nodeIDs2;
        for (const auto& id : nodeIDs1) {
            if (seenNodeIDs.find(id) == seenNodeIDs.end()) {
                mergedNodeIDs.push_back(id);
                seenNodeIDs.insert(id);
            }
        }

        if (cutbool2!= cutbool1){
            if(cutbool2.first == false){
                if(cutbool1.second == 2 || cutbool1.second == 3){
                    cutbool2 = cutbool1;
                }
                
            }
            else if(cutbool1.first == false){
                if(cutbool2.second == 2){
                    cutbool2 = make_pair(false, 0);
                }else if(cutbool2.second == 3){
                    cutbool2.second = 1;
                }

            }
            else if(cutbool2.first == true){
                if(cutbool2.second != cutbool1.second){
                    if(cutbool2.second == 1 && cutbool1.second == 2 || cutbool2.second == 1 && cutbool1.second == 3){
                        cutbool2.second = 3;
                    }else if (cutbool2.second == 2 && cutbool1.second == 1){
                        cutbool2 = make_pair(false, 0);
                    }else if (cutbool2.second == 3 && cutbool1.second == 1){
                        cutbool2.second = 1;
                    }
                }
            }
        }
        finalcut = cutbool2;
	    
	}


	for (size_t i = 0; i < mergedSeq.size(); ++i) {
	    const auto& counts = mergedCountMatrix[i];
	    int maxCount = -1;
	    int totalCount = 0;
	    std::vector<char> maxBases;
	    int nonZeroBases = 0;

	    // Determine totalCount, maxCount, and number of bases with non-zero counts
	    for (const auto& [base, count] : counts) {
	        if (count > 0) {
	            ++nonZeroBases;
	        }
	        totalCount += count;
	        if (count > maxCount) {
	            maxCount = count;
	            maxBases = {base};
	        } else if (count == maxCount) {
	            maxBases.push_back(base);
	        }
	    }

	    char originalBase = mergedSeq[i];

	    if (mode == "reckless") {
	        // Handle case where no base has any support
	        if (maxCount <= 0 || maxBases.empty()) {
	            mergedSeq[i] = 'N';
	            std::cerr << "Reckless mode: All counts are zero at position " << i
	                      << ", masking as 'N'\n";
	            continue;
	        }

	        // Only correct if originalBase is not among the most frequent bases
	        if (std::find(maxBases.begin(), maxBases.end(), originalBase) == maxBases.end()) {
	            mergedSeq[i] = maxBases[0];
	            std::cerr << "Reckless mode: Corrected position " << i
	                      << " from '" << originalBase << "' to '" << mergedSeq[i] << "'\n";
	            std::cerr << "Count matrix at position " << i << ":\n";
	            for (const auto& [base, count] : counts) {
	                std::cerr << "  Base " << base << ": " << count << "\n";
	            }
	        }

	    } else if (mode == "normal") {
	        if (nonZeroBases > 1) {
	            double freq = static_cast<double>(maxCount) / totalCount;
	            if (freq >= 0.65) {
	                mergedSeq[i] = maxBases[0];
	            } else {
	                mergedSeq[i] = 'N';
	                std::cerr << "Normal mode: Masked position " << i << " due to low frequency\n";
	            }
	        }

	    } else if (mode == "strict") {
	        if (nonZeroBases > 1) {
	            double freq = static_cast<double>(maxCount) / totalCount;

	            if (maxBases.size() == 1) {
	                char topBase = maxBases[0];

	                bool isTransition =
	                    (originalBase == 'A' && topBase == 'G') ||
	                    (originalBase == 'G' && topBase == 'A') ||
	                    (originalBase == 'C' && topBase == 'T') ||
	                    (originalBase == 'T' && topBase == 'C');

	                if ((isTransition && freq >= 0.9) || (!isTransition && freq >= 0.65)) {
	                    mergedSeq[i] = topBase;
	                } else {
	                    mergedSeq[i] = 'N';
	                    std::cerr << "Strict mode: Masked position " << i << " due to insufficient confidence\n";
	                }
	            } else {
	                mergedSeq[i] = 'N';  // Tie between top bases
	                std::cerr << "Strict mode: Masked position " << i << " due to tie between top bases\n";
	            }
	        }
	    }
	}
	
	
	
	
    return {mergedSeq, mergedScoringMatrix, mergedNodeIDs, mergedRySeq, mergedCountMatrix, finalcut};
}



bool assembly::basesMatchWithDamage(char a, char b) {
    // Consider ancient DNA damage: C -> T and G -> A
    if (a == b) return true;
    if ((a == 'C' && b == 'T') || (a == 'T' && b == 'C')) return true;
    if ((a == 'G' && b == 'A') || (a == 'A' && b == 'G')) return true;
    return false;
}

MergeResult assembly::tryMergeContigs(
    const std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>>& contigs,
    int i,
    int j,
    int min_overlap_length,
    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap,
    std::unordered_set<int>& mergedIndices, bool specifiedDeam, int lenMin
) {
    MergeResult result = {false, OverlapResult(), false, -1};

    // Extract node IDs
    std::vector<int> nodeIDs_i = std::get<2>(contigs[i]);
    std::vector<int> nodeIDs_j = std::get<2>(contigs[j]);
    // Extract sequences


    // Use an unordered_set for efficient lookup of node IDs in contig j
    std::unordered_set<int> nodeIDs_set_j(nodeIDs_j.begin(), nodeIDs_j.end());

    // Find common node IDs without sorting
    std::vector<int> common_nodeIDs;
    for (int id : nodeIDs_i) {
        if (nodeIDs_set_j.find(id) != nodeIDs_set_j.end()) {
            common_nodeIDs.push_back(id);
        }
    }

    if (common_nodeIDs.empty())
    {
        //std::cerr << "No common node IDs found between contigs " << i << " and " << j << "\n";
        return result; // Skip merge if no common node ID
        
    } else
    {
#ifdef DEBUGASS2
        std::cerr << "Node IDs in contig " << i << ": ";
        for (int id : nodeIDs_i) {
            std::cerr << id << " ";
        }
        std::cerr << std::endl;

        std::cerr << "Node IDs in contig " << j << ": ";
        for (int id : nodeIDs_j) {
            std::cerr << id << " ";
        }
        std::cerr << std::endl;

        std::cerr << "Common Node IDs are: ";
        for (int id : common_nodeIDs) {
            std::cerr << id << " ";
        }
        std::cerr << std::endl;
#endif

    }
    // Later usage of common_nodeIDs
    for (int common_id : common_nodeIDs) {
        // Confirm the common_id is indeed in both contig node lists
        if (std::find(nodeIDs_i.begin(), nodeIDs_i.end(), common_id) == nodeIDs_i.end() ||
            std::find(nodeIDs_j.begin(), nodeIDs_j.end(), common_id) == nodeIDs_j.end()) {
            std::cerr << "Error: Common node ID " << common_id << " not found in both lists!" << "\n";
        } else {
            //std::cerr << "Confirmed: Common node ID " << common_id << " is present in both lists." << "\n";
        }
    }

    std::string seq_i = std::get<0>(contigs[i]);
    std::string seq_j = std::get<0>(contigs[j]);

    OverlapResult best_overlap;
    bool isOverlapAtEnd;

    std::string ry_i = std::get<3>(contigs[i]);
    std::string ry_j = std::get<3>(contigs[j]);

    pair<bool, int> cutbool_i = get<5>(contigs[i]);
    pair<bool, int> cutbool_j = get<5>(contigs[j]);
    
    OverlapResult overlapResult_end;
    OverlapResult overlapResult_begin;
    // Check for suitable overlap
    overlapResult_end = get_overlap_length_and_score(ry_i, ry_j, min_overlap_length, lenMin);
    overlapResult_begin = get_overlap_length_and_score(ry_j, ry_i, min_overlap_length, lenMin);

    OverlapResult seqTE = get_overlap_length_and_score(seq_i, seq_j, min_overlap_length, lenMin);
    OverlapResult seqTS = get_overlap_length_and_score(seq_j, seq_i, min_overlap_length, lenMin);

    if((seqTE.length > overlapResult_end.length && seqTE.score > overlapResult_end.score) || (seqTS.length > overlapResult_begin.length && seqTS.score > overlapResult_begin.score)){
        overlapResult_end = seqTE;
        overlapResult_begin = seqTS;
    }

    // if(!specifiedDeam){
//          overlapResult_end = get_overlap_length_and_score(seq_i, seq_j, min_overlap_length, lenMin);
//          overlapResult_begin = get_overlap_length_and_score(seq_j, seq_i, min_overlap_length, lenMin);
//     }
#ifdef DEBUGASS2
    // // Debug output for overlap results
    std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;
#endif
    // Determine positions of the first common node ID in both contig node lists
    int first_common_id = common_nodeIDs.front();
    int pos_in_i = std::find(nodeIDs_i.begin(), nodeIDs_i.end(), first_common_id) - nodeIDs_i.begin();
    int pos_in_j = std::find(nodeIDs_j.begin(), nodeIDs_j.end(), first_common_id) - nodeIDs_j.begin();

    // Use positions to decide merge direction
    bool is_i_to_j = (pos_in_i < pos_in_j);
    if(pos_in_i != 0 && pos_in_j != 0 && pos_in_i != pos_in_j){
        if(overlapResult_end.length > overlapResult_begin.length && overlapResult_end.score > overlapResult_begin.score){
            is_i_to_j = false;
        }else{
            is_i_to_j = true;
        }
    }
#ifdef DEBUGASS2
    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
#endif
    

    // Special case both contigs begin on the same node id
    if (((overlapResult_end.length > common_nodeIDs.size() && overlapResult_end.score > 0) || (overlapResult_begin.score > 0 && overlapResult_begin.length > common_nodeIDs.size() ))) {
        if (pos_in_i == pos_in_j) {
            //std::cerr << "Special case: both contigs begin on the same node id\n";
            if (overlapResult_end.score > overlapResult_begin.score && overlapResult_end.length > overlapResult_begin.length) {
                best_overlap = overlapResult_end;
                isOverlapAtEnd = true;
            } else if (overlapResult_end.score < overlapResult_begin.score && overlapResult_end.length < overlapResult_begin.length) {
                best_overlap = overlapResult_begin;
                isOverlapAtEnd = false;
            } else {
#ifdef DEBUGASS4
        std::cerr << "last else statement for pos equal \n";
         std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

                return result; // no way to turst the results otheriwse
                // if (overlapResult_end.length > overlapResult_begin.length) {
                //     best_overlap = overlapResult_end;
                //     isOverlapAtEnd = true;
                // } else if (overlapResult_begin.length > overlapResult_end.length) {
                //     best_overlap = overlapResult_begin;
                //     isOverlapAtEnd = false;
                // }
            }
        } else {
            // Compare scores and decide tie-breaking rule
            if ((overlapResult_end.score > overlapResult_begin.score && overlapResult_end.length > overlapResult_begin.length) && !is_i_to_j) {
                best_overlap = overlapResult_end;
                isOverlapAtEnd = true; 
            } else if ((overlapResult_end.score < overlapResult_begin.score && overlapResult_end.length < overlapResult_begin.length) && is_i_to_j) {
                best_overlap = overlapResult_begin;
                isOverlapAtEnd = false;
            } else {

                if(cutbool_i.first == true || cutbool_j.first == true)
                {
#ifdef DEBUGASS2
                    // Debug output for overlap results
                    cerr << "cutbool_i " << cutbool_i << " cutbool_j " << cutbool_j << endl;
                    std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif
                    if (overlapResult_end.score > overlapResult_begin.score && overlapResult_end.length > overlapResult_begin.length) {
                        best_overlap = overlapResult_end;
                        isOverlapAtEnd = true;
                    } else if (overlapResult_end.score < overlapResult_begin.score && overlapResult_end.length < overlapResult_begin.length) {
                        best_overlap = overlapResult_begin;
                        isOverlapAtEnd = false;
                    }else{
#ifdef DEBUGASS4
                    std::cerr << "Else statement for cutbool true \n";
                    std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

                        return result;
                        
                    }
                }else{
#ifdef DEBUGASS2
                    std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif
                    if(!is_i_to_j){
                        best_overlap = overlapResult_end;
                        isOverlapAtEnd = true;
                    }else if(is_i_to_j){
                        best_overlap = overlapResult_begin;
                        isOverlapAtEnd = false;   
                    }else{
#ifdef DEBUGASS4
                    std::cerr << "Most inner if statement after onlty node id direction is asked \n";
                    std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

                        return result;
                    }
                }
            }

                    
            
        }
    }


    // Attempt merging based on node IDs and sequences
    else if (common_nodeIDs.size() == 1) {
        //std::cerr << "one common node id and no overlap length" << std::endl;
        int common_node_id = common_nodeIDs.front();
        std::string common_node_seq = std::get<0>(nodeSequenceMap.at(common_node_id));

        auto findMatch = [this](const std::string& contig, const std::string& node_seq, bool from_start) -> size_t {
            size_t max_match_length = std::min(contig.size(), node_seq.size());
            size_t len = 0;

            if (from_start) {
                // Start from the beginning of the contig and the end of the node_seq
                for (size_t i = 0; i < max_match_length; ++i) {
                    if (!basesMatchWithDamage(contig[i], node_seq[node_seq.size() - max_match_length + i])) {
                        break; // Stop at the first mismatch considering DNA damage
                    }
                    len++;
                }
            } else {
                // Start from the end of the contig and the beginning of the node_seq
                for (size_t i = 0; i < max_match_length; ++i) {
                    if (!basesMatchWithDamage(contig[contig.size() - max_match_length + i], node_seq[i])) {
                        break; // Stop at the first mismatch considering DNA damage
                    }
                    len++;
                }
            }

            return len;
        };


        if (nodeIDs_i.back() == common_node_id && nodeIDs_j.front() == common_node_id && cutbool_i.second < 2 && cutbool_j.second != 1 && cutbool_j.second != 3) {
            size_t match_length_i = findMatch(seq_i, common_node_seq, false);
            size_t match_length_j = findMatch(seq_j, common_node_seq, true);
            //cerr << "match_length_i " << match_length_i << " match_length_j " << match_length_j << endl;
            if (match_length_i > 0 && match_length_j > 0) {
                int numN = common_node_seq.size() - (match_length_i + match_length_j);
                result.canMerge = true;
                result.isOverlapAtEnd = true;
                result.numN = numN;
            } else {
                result.canMerge = true;
                result.isOverlapAtEnd = true;
                result.numN = 0;
            }
        }

        if (nodeIDs_j.back() == common_node_id && nodeIDs_i.front() == common_node_id&& cutbool_j.second < 2 && cutbool_i.second != 1 && cutbool_i.second != 3) {
            size_t match_length_j = findMatch(seq_j, common_node_seq, false);
            size_t match_length_i = findMatch(seq_i, common_node_seq, true);
            //cerr << "2match_length_i " << match_length_i << " match_length_j " << match_length_j << endl;
            if (match_length_j > 0 && match_length_i > 0) {
                int numN = common_node_seq.size() - (match_length_i + match_length_j);
                result.canMerge = true;
                result.isOverlapAtEnd = false;
                result.numN = numN;
            } else {
                result.canMerge = true;
                result.isOverlapAtEnd = false;
                result.numN = 0;
            }
        }
        //cerr << "return after node id length == 1" << endl;
        return result;

    }

    // Handle case with multiple common nodes but no overlap
    else if (overlapResult_end.length == 0 && overlapResult_begin.length == 0 && common_nodeIDs.size() > 1)
    {
#ifdef DEBUGASS2
        std::cerr << "No overlap but common node IDs exist\n";
         std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

        
        OverlapResult overlapResult_end;
        OverlapResult overlapResult_begin;
        // Check for suitable overlap
        overlapResult_end = get_overlap_length_and_score(ry_i, ry_j, min_overlap_length, lenMin);
        overlapResult_begin = get_overlap_length_and_score(ry_j, ry_i, min_overlap_length, lenMin);

        OverlapResult seqTE = get_overlap_length_and_score(seq_i, seq_j, min_overlap_length, lenMin);
        OverlapResult seqTS = get_overlap_length_and_score(seq_j, seq_i, min_overlap_length, lenMin);

        if((seqTE.length > overlapResult_end.length && seqTE.score > overlapResult_end.score) || (seqTS.length > overlapResult_begin.length && seqTS.score > overlapResult_begin.score)){
            overlapResult_end = seqTE;
            overlapResult_begin = seqTS;
        }

        // if(!specifiedDeam){
//              overlapResult_end = get_overlap_length_and_score(seq_i, seq_j, min_overlap_length, lenMin);
//              overlapResult_begin = get_overlap_length_and_score(seq_j, seq_i, min_overlap_length, lenMin);
//         }
        if ((overlapResult_end.length >= common_nodeIDs.size()  && overlapResult_end.score > 0) || (overlapResult_begin.score > 0 && overlapResult_begin.length >= common_nodeIDs.size() )){

            if (pos_in_i == pos_in_j) {
                if (overlapResult_end.score > overlapResult_begin.score && overlapResult_end.length > overlapResult_begin.length) {
                    best_overlap = overlapResult_end;
                    isOverlapAtEnd = true;
                } else if (overlapResult_end.score < overlapResult_begin.score && overlapResult_end.length < overlapResult_begin.length) {
                    best_overlap = overlapResult_begin;
                    isOverlapAtEnd = false;
                } else {
#ifdef DEBUGASS4
                std::cerr << "No overlap but common node IDs when pos equal \n";
                 std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                            std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

                    // if (overlapResult_end.length > overlapResult_begin.length) {
                    //     best_overlap = overlapResult_end;
                    //     isOverlapAtEnd = true;
                    // } else if (overlapResult_begin.length > overlapResult_end.length) {
                    //     best_overlap = overlapResult_begin;
                    //     isOverlapAtEnd = false;
                    // }
                    return result;
                }
            }

            bool is_i_to_j = (pos_in_i < pos_in_j);
            // cerr << "is_i_to_j " << is_i_to_j << endl;
            // cerr << "cutbool_i " << cutbool_i << " cutbool_j " << cutbool_j << endl;
            if (cutbool_j.first == true || cutbool_i.first == true){
                if (overlapResult_end.score > overlapResult_begin.score && overlapResult_end.length > overlapResult_begin.length) {
                    best_overlap = overlapResult_end;
                    isOverlapAtEnd = true;
                } else if (overlapResult_end.score < overlapResult_begin.score && overlapResult_end.length < overlapResult_begin.length) {
                    best_overlap = overlapResult_begin;
                    isOverlapAtEnd = false;
                } else {
#ifdef DEBUGASS4
        std::cerr << "No overlap but common node IDs exist when cutbool true \n";
         std::cerr << "OverlapResult_end: length = " << overlapResult_end.length << ", score = " << overlapResult_end.score << std::endl;
                    std::cerr << "OverlapResult_begin: length = " << overlapResult_begin.length << ", score = " << overlapResult_begin.score << std::endl;

                    std::cerr << "Position in contig " << i << ": " << pos_in_i << ", Position in contig " << j << ": " << pos_in_j << std::endl;
                    std::cerr << "Merge direction is i to j: " << is_i_to_j << std::endl;
                    std::cerr << "Sequence for contig " << i << ": " << seq_i << std::endl;
                    std::cerr << ry_i << std::endl;
                    std::cerr << "Sequence for contig " << j << ": " << seq_j << std::endl;
                    std::cerr << ry_j << std::endl;
                        std::cerr << "Node IDs in contig " << i << ": ";
                        for (int id : nodeIDs_i) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;

                        std::cerr << "Node IDs in contig " << j << ": ";
                        for (int id : nodeIDs_j) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                            }
                        std::cerr << std::endl;

                        std::cerr << "Common Node IDs are: ";
                        for (int id : common_nodeIDs) {
                            const auto& node_tuple = nodeSequenceMap.at(id);
                            std::cerr << id << " (" 
                                      << std::get<0>(node_tuple) << ", " 
                                      << std::get<1>(node_tuple) << ", " 
                                      << std::get<2>(node_tuple) << ") ";
                        }
                        std::cerr << std::endl;
#endif

                    // if (overlapResult_end.length > overlapResult_begin.length) {
                    //     best_overlap = overlapResult_end;
                    //     isOverlapAtEnd = true;
                    // } else if (overlapResult_begin.length > overlapResult_end.length) {
                    //     best_overlap = overlapResult_begin;
                    //     isOverlapAtEnd = false;
                    // }
                    return result;
                }
            }else{
                if (!is_i_to_j) {
                    best_overlap = overlapResult_end;
                    isOverlapAtEnd = true;
                } else if (is_i_to_j) {
                    best_overlap = overlapResult_begin;
                    isOverlapAtEnd = false;
                } else {
                    //cerr << "No overlap but node ids when only node id direction is  there" << endl;
                    return result;
                    
                }
            }
            
            //std::cerr << "For no overlap, scores are: overlapResult_begin.score = " << overlapResult_begin.score << ", overlapResult_end.score = " << overlapResult_end.score << " overlapResult_begin.length = " << overlapResult_begin.length << ", overlapResult_end.length = " << overlapResult_end.length << std::endl;
        } else {
            bool subsetTest = isSubset(seq_i, nodeIDs_i, seq_j, nodeIDs_j, 0);
                    if (subsetTest == true){
                        //cerr << "Error: This actually happens" << endl;
                        OverlapResult ovSign;
                        ovSign.length = -1;
                        ovSign.score = -1;
                        ovSign.sequence = "subset";
                        best_overlap = ovSign;
                        isOverlapAtEnd = true;
                    }
                    
            return result; // Skip merge if no suitable overlap found
        }
    } else {
        //std::cerr << "there is no overlap here " << std::endl;
        return result;
    }

    mergedIndices.insert(i);
    mergedIndices.insert(j);
    result.canMerge = true;
    result.best_overlap = best_overlap;
    result.isOverlapAtEnd = isOverlapAtEnd;

#ifdef DEBUGASS2
    std::cerr << "Final merge decision for contigs " << i << " and " << j << ": canMerge = " << result.canMerge << ", isOverlapAtEnd = " << result.isOverlapAtEnd << ", numN = " << result.numN << ", best_overlap.length = " << best_overlap.length << ", best_overlap.score = " << best_overlap.score << "\n";
    cerr << "Sequence i " << seq_i << " sequence j " << seq_j << endl;
#endif

    return result;
}



std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> assembly::mergeOverlappingContigs(
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int> >> contigs,
    int min_overlap_length,
    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap, bool specifiedDeam, string mode, int lenMin
) {
    // Step 1: Remove subset contigs
    contigs = removeSubsetContigs(contigs);

    // Check for invalid values in all contigs' scoring matrices
    for (const auto& contig : contigs) {
        const auto& scoringMatrix = std::get<1>(contig);
        for (const auto& map : scoringMatrix) {
            for (const auto& [base, score] : map) {
                if (std::isinf(score) || std::isnan(score)) {
                    std::cerr << "Error: Invalid value in scoring matrix: " << score << "\n";
                }
            }
        }
    }
    

    struct OverlapInfo {
        size_t i;
        size_t j;
        MergeResult mergeResult;
    };

    std::vector<OverlapInfo> potentialOverlaps;

    // Identify potential overlaps
    for (size_t i = 0; i < contigs.size(); ++i) {
        for (size_t j = i + 1; j < contigs.size(); ++j) {
            std::unordered_set<int> tempMergedIndices; // Create a temporary unordered_set
            MergeResult mergeResult = tryMergeContigs(contigs, i, j, min_overlap_length, nodeSequenceMap, tempMergedIndices, specifiedDeam, lenMin);
            if (mergeResult.canMerge) {
                if(mergeResult.best_overlap.length == -1 && mergeResult.best_overlap.sequence == "subset"){
                    cerr << "Error: undetected subset. " << endl;
                    continue;
                }
                if((mergeResult.isOverlapAtEnd && mergeResult.best_overlap.length > get<0>(contigs[j]).size()) || !mergeResult.isOverlapAtEnd && mergeResult.best_overlap.length > get<0>(contigs[i]).size()){
                        std::cerr << "Overlap length exceeds or matches sequence length; cannot merge.\n";
                    continue;
                }
#ifdef DEBUGASS2

                std::cerr << "Potential overlap found between contig " << i << " and contig " << j << "\n";
                std::cerr << "Attempting to merge " << i << " with " << j << " using overlap length " << mergeResult.best_overlap.length << " on sequences of sizes " << get<0>(contigs[i]).size() << " and " << get<0>(contigs[j]).size() << std::endl;

#endif
                potentialOverlaps.push_back({i, j, mergeResult});
            }
        }
    }

    std::unordered_map<size_t, std::vector<OverlapInfo>> overlapCount;
    for (const auto& overlap : potentialOverlaps) {
        overlapCount[overlap.i].push_back(overlap);
    }

    std::vector<OverlapInfo> filteredOverlaps;

    for (const auto& entry : overlapCount) {
        const auto& overlaps = entry.second;
        std::unordered_map<bool, OverlapInfo> bestOverlap; // Key is true for end, false for beginning

        for (const auto& overlap : overlaps) {
            bool isEnd = overlap.mergeResult.isOverlapAtEnd;

            if (bestOverlap.find(isEnd) == bestOverlap.end() ||
                overlap.mergeResult.best_overlap.score > bestOverlap[isEnd].mergeResult.best_overlap.score ||
                (overlap.mergeResult.best_overlap.score == bestOverlap[isEnd].mergeResult.best_overlap.score &&
                 overlap.mergeResult.best_overlap.length > bestOverlap[isEnd].mergeResult.best_overlap.length)) {
                bestOverlap[isEnd] = overlap;
            }
        }

        for (const auto& best : bestOverlap) {
            filteredOverlaps.push_back(best.second);
        }
    }

    std::cerr << "Total filtered overlaps: " << filteredOverlaps.size() << "\n";

    std::unordered_set<int> mergedIndices;
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, pair<bool, int>>> mergedContigs;

    for (const auto& overlap : filteredOverlaps)
    {
        size_t i = overlap.i;
        size_t j = overlap.j;

        if (mergedIndices.find(i) != mergedIndices.end() || mergedIndices.find(j) != mergedIndices.end()) {
            //std::cerr << "Skipping already merged contigs " << i << " and/or " << j << "\n";
            continue; // Skip already merged contigs
        }

        const MergeResult& mergeResult = overlap.mergeResult;
        if (mergeResult.numN > -1) // We have Ns inbetween the merging contigs 
        {
            //cerr << "We have Ns between the contigs, special case merging starting" << endl;
            std::string Ns(mergeResult.numN, 'N');

            // SPECIAL CASE: Merge the sequences with Ns in between //
            std::string mergedSeq;
            std::string mergedRySeq;
            std::vector<std::unordered_map<char, double>> mergedScoringMatrix;
            vector<unordered_map<char, int>> mergedCountMatrix;
            pair<bool, int> mergedcut = make_pair(false, 0);

            std::vector<int> mergedNodeIDs;

            if (mergeResult.isOverlapAtEnd) {

#ifdef DEBUGASS2
                std::cerr << "Merging contig " << i << " (end) with contig " << j << " with overlap length: " << mergeResult.best_overlap.length << "\n";
                cerr << "contig i seq " << std::get<0>(contigs[i]) << " contitg j seq " << get<0>(contigs[j]) << endl;
                cerr << "contig i seq " << std::get<3>(contigs[i]) << " contitg j seq " << get<3>(contigs[j]) << endl;

            
#endif
                mergedSeq = std::get<0>(contigs[i]) + Ns + std::get<0>(contigs[j]);
                mergedRySeq = std::get<3>(contigs[i]) + Ns + std::get<3>(contigs[j]);
                mergedScoringMatrix = std::get<1>(contigs[i]);
                mergedScoringMatrix.resize(mergedSeq.size());
                mergedCountMatrix = get<4>(contigs[i]);
                mergedCountMatrix.resize(mergedSeq.size());

                int rest = std::get<0>(contigs[i]).size();
                for (int sp = 0; sp < Ns.size(); ++sp) {
                    mergedScoringMatrix[rest + sp] = std::unordered_map<char, double>(); // Empty map for 'N's
                    mergedCountMatrix[rest + sp] = unordered_map<char, int>();
                }

                for (int sp = 0; sp < std::get<0>(contigs[j]).size(); ++sp) {
                    mergedScoringMatrix[rest + Ns.size() + sp] = std::get<1>(contigs[j])[sp];
                    mergedCountMatrix[rest + Ns.size() + sp] = std::get<4>(contigs[j])[sp];
                }


                mergedNodeIDs = std::get<2>(contigs[i]);
                mergedNodeIDs.insert(mergedNodeIDs.end(), std::get<2>(contigs[j]).begin(), std::get<2>(contigs[j]).end());
                pair<bool, int> cutbool1 = get<5>(contigs[j]);
                pair<bool, int> cutbool2 = get<5>(contigs[i]);
                if (cutbool2!= cutbool1){
                    if(cutbool2.first == false){
                        if(cutbool1.second == 2 || cutbool1.second == 3){
                            cutbool2 = cutbool1;
                        }
                        
                    }
                    else if(cutbool1.first == false){
                        if(cutbool2.second == 2){
                            cutbool2 = make_pair(false, 0);
                        }else if(cutbool2.second == 3){
                            cutbool2.second = 1;
                        }

                    }
                    else if(cutbool2.first == true){
                        if(cutbool2.second != cutbool1.second){
                            if(cutbool2.second == 1 && cutbool1.second == 2 || cutbool2.second == 1 && cutbool1.second == 3){
                                cutbool2.second = 3;
                            }else if (cutbool2.second == 2 && cutbool1.second == 1){
                                cutbool2 = make_pair(false, 0);
                            }else if (cutbool2.second == 3 && cutbool1.second == 1){
                                cutbool2.second = 1;
                            }
                        }
                    }
                }
                mergedcut = cutbool2;

            } else {

#ifdef DEBUGASS2
                std::cerr << "Merging contig N cases " << j << " (beginning) with contig " << i << " with overlap length: " << mergeResult.best_overlap.length << "\n";
                cerr << "contig i seq " << std::get<0>(contigs[i]) << " contitg j seq " << get<0>(contigs[j]) << endl;
#endif
                mergedSeq = std::get<0>(contigs[j]) + Ns + std::get<0>(contigs[i]);
                mergedRySeq = std::get<3>(contigs[j]) + Ns + std::get<3>(contigs[i]);
                mergedScoringMatrix = std::get<1>(contigs[j]);
                mergedScoringMatrix.resize(mergedSeq.size());
                mergedCountMatrix = get<4>(contigs[j]);
                mergedCountMatrix.resize(mergedSeq.size());

                int rest = std::get<0>(contigs[j]).size();
                for (int sp = 0; sp < Ns.size(); ++sp) {
                    mergedScoringMatrix[rest + sp] = std::unordered_map<char, double>(); // Empty map for 'N's
                    mergedCountMatrix[rest + sp] = unordered_map<char, int>();
                }

                for (int sp = 0; sp < std::get<0>(contigs[i]).size(); ++sp) {
                    mergedScoringMatrix[rest + Ns.size() + sp] = std::get<1>(contigs[i])[sp];
                    mergedCountMatrix[rest + Ns.size() + sp] = std::get<4>(contigs[i])[sp];
                }

                mergedNodeIDs = std::get<2>(contigs[j]);
                mergedNodeIDs.insert(mergedNodeIDs.end(), std::get<2>(contigs[i]).begin(), std::get<2>(contigs[i]).end());

                pair<bool, int> cutbool2 = get<5>(contigs[j]);
                pair<bool, int> cutbool1 = get<5>(contigs[i]);
                if (cutbool2!= cutbool1){
                    if(cutbool2.first == false){
                        if(cutbool1.second == 2 || cutbool1.second == 3){
                            cutbool2 = cutbool1;
                        }
                        
                    }
                    else if(cutbool1.first == false){
                        if(cutbool2.second == 2){
                            cutbool2 = make_pair(false, 0);
                        }else if(cutbool2.second == 3){
                            cutbool2.second = 1;
                        }

                    }
                    else if(cutbool2.first == true){
                        if(cutbool2.second != cutbool1.second){
                            if(cutbool2.second == 1 && cutbool1.second == 2 || cutbool2.second == 1 && cutbool1.second == 3){
                                cutbool2.second = 3;
                            }else if (cutbool2.second == 2 && cutbool1.second == 1){
                                cutbool2 = make_pair(false, 0);
                            }else if (cutbool2.second == 3 && cutbool1.second == 1){
                                cutbool2.second = 1;
                            }else{
                                continue;
                            }
                        }
                    }else{
                        continue;
                    }
                }
                mergedcut = cutbool2;
            }

            // Debug: Ensure scoring matrix is not empty
            for (size_t idx = 0; idx < mergedScoringMatrix.size(); ++idx) {
                if (mergedScoringMatrix[idx].empty()) {
                    //std::cerr << "Warning: Ns thing Empty scoring matrix at index " << idx << " in merged contig.\n";
                }
            }
            // bool mergedcut = false;
            // if(get<5>(contigs[i]) == true || get<5>(contigs[j]) == true){
            //    mergedcut = true;
            //}
            

            auto con = std::make_tuple(mergedSeq, mergedScoringMatrix, mergedNodeIDs, mergedRySeq, mergedCountMatrix, mergedcut);
            mergedContigs.push_back(con);
        } else { /// we have no Ns between the two contigs - we can merge directly

#ifdef DEBUGASS2
            std::cerr << "Directly merging contig " << i << " with contig " << j << " with overlap length: " << mergeResult.best_overlap.length << "\n";
            cerr << "at end " << mergeResult.isOverlapAtEnd << endl;
            cerr << "contig i seq " << std::get<0>(contigs[i]) << " contitg j seq " << get<0>(contigs[j]) << endl;
            cerr << "contig i seq " << std::get<3>(contigs[i]) << " contitg j seq " << get<3>(contigs[j]) << endl;
            cerr << "contig i size " << std::get<0>(contigs[i]).size() <<" contitg j size " << get<0>(contigs[j]).size() << endl;
            cerr << "contig i size " << std::get<3>(contigs[i]).size() <<" contitg j size " << get<3>(contigs[j]).size() << endl;
#endif

            auto con = mergeContigs(contigs[i], contigs[j], mergeResult.best_overlap.length, mergeResult.isOverlapAtEnd, mode);
            mergedContigs.push_back(con);
        }

        mergedIndices.insert(i);
        mergedIndices.insert(j);
    }

    // Add any contigs that were not merged
    for (size_t i = 0; i < contigs.size(); ++i) {
        if (mergedIndices.find(i) == mergedIndices.end()) {
            mergedContigs.push_back(contigs[i]);
        }
    }

    return mergedContigs;
}



// Function to print the scoring matrix in probability space
void assembly::printScoringMatrix(const std::vector<std::unordered_map<char, double>>& scoringMatrix, std::string &outfile, size_t startAt) {
    std::ofstream scoreMat;
    scoreMat.open(outfile);
    size_t index = startAt;
    for (const auto& scoreMap : scoringMatrix) {
        scoreMat << "Position " << index + 1 << ":\t";
        if (scoreMap.empty()) {
            //scoreMat << "A=0\tC=0\tG=0\tT=0\t"; // Output zeros for Ns
        } else {
            for (const auto& score : scoreMap) {
                scoreMat << score.first << "=" << std::exp(score.second) << "\t";
            }
        }
        scoreMat << "\n";
        ++index;
    }

    scoreMat.close();
}

// Function to print the scoring matrix in probability space
void assembly::printCountMatrix(const std::vector<std::unordered_map<char, int>>& countMatrix, std::string &outfile, size_t startAt) {
    std::ofstream scoreMat;
    scoreMat.open(outfile);
    size_t index = startAt;
    for (const auto& scoreMap : countMatrix) {
        scoreMat << "Position " << index + 1 << ":\t";
        if (scoreMap.empty()) {
            //scoreMat << "A=0\tC=0\tG=0\tT=0\t"; // Output zeros for Ns
        } else {
            for (const auto& score : scoreMap) {
                scoreMat << score.first << "=" << score.second << "\t";
            }
        }
        scoreMat << "\n";
        ++index;
    }

    scoreMat.close();
}

std::optional<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>>
assembly::checkAndMergeContigs(
    std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>& contig_i,
    std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>& contig_j,
    const std::map<int, std::tuple<std::string, size_t, int>>& nodeSequenceMap, int lenMin, int scoreMin, bool specifiedDeam, string mode) 
{
    std::string seq_i = std::get<3>(contig_i);
    std::string seq_j = std::get<3>(contig_j);
    std::vector<int> nodeIds_i = std::get<2>(contig_i);
    std::vector<int> nodeIds_j = std::get<2>(contig_j);

    std::string s_i = std::get<0>(contig_i);
    std::string s_j = std::get<0>(contig_j);
    auto cutbool_i = std::get<5>(contig_i);
    auto cutbool_j = std::get<5>(contig_j);

#ifdef DEBUGASS2
    std::cerr << "Checking overlap between contig_i and contig_j\n";
    std::cerr << "Contig_i sequence: " << s_i << "\n";
    std::cerr << "Contig_j sequence: " << s_j << "\n";
#endif

    // Check for common node IDs between contig_i and contig_j using a loop
    std::vector<int> commonNodeIds;
    for (const auto& id_i : nodeIds_i) {
        if (std::find(nodeIds_j.begin(), nodeIds_j.end(), id_i) != nodeIds_j.end()) {
            commonNodeIds.push_back(id_i);
        }
    }

    if (!commonNodeIds.empty()) {
#ifdef DEBUGASS2
        std::cerr << "Common node IDs found between contig_i and contig_j: ";
        for (const auto& id : commonNodeIds) {
            std::cerr << id << " ";
        }
        std::cerr << "\n";
#endif
    }

    // Check for suitable overlap in both directions
    OverlapResult overlapResult = get_overlap_length_and_score(s_i, s_j, lenMin, scoreMin);
    OverlapResult overlapResult2 = get_overlap_length_and_score(s_j, s_i, lenMin, scoreMin);
    if (!specifiedDeam) {
        overlapResult = get_overlap_length_and_score(seq_i, seq_j, lenMin, scoreMin);
        overlapResult2 = get_overlap_length_and_score(seq_j, seq_i, lenMin, scoreMin);
    }

#ifdef DEBUGASS2
    std::cerr << "OverlapResult: length = " << overlapResult.length << ", score = " << overlapResult.score << " length of common node ids " << commonNodeIds.size() << "\n";
    std::cerr << "OverlapResult2: length = " << overlapResult2.length << ", score = " << overlapResult2.score << "\n";
    cerr << "seq size " << s_i.size() << " " << s_j.size() << endl;
#endif


    bool canMergeWithOverlapResult = (cutbool_i.first && (cutbool_i.second == 2 || cutbool_i.second == 3)) || (cutbool_j.first && (cutbool_j.second == 1 || cutbool_j.second == 3));
    bool canMergeWithOverlapResult2 = (cutbool_j.first && (cutbool_j.second == 2 || cutbool_j.second == 3)) || (cutbool_i.first && (cutbool_i.second == 1 || cutbool_i.second == 3));

    // Implement check to ignore if overlap length is larger than the sequence length
    if (canMergeWithOverlapResult && overlapResult.length > lenMin && overlapResult.score > scoreMin &&
        overlapResult.length < s_i.size() && overlapResult.length < s_j.size() &&
        (overlapResult.score > overlapResult2.score || overlapResult.length > overlapResult2.length)) {
        
        bool isOverlapAtEnd = true;
        auto mergedContig = mergeContigs(contig_i, contig_j, overlapResult.length, isOverlapAtEnd, mode);

#ifdef DEBUGASS2
        std::cerr << "Merged contigs: new sequence overlap  = " << std::get<0>(mergedContig) << "\n";
#endif

        return mergedContig;
    }

    // Check cutbool conditions for overlapResult2
    if (canMergeWithOverlapResult2 && overlapResult2.length > lenMin && overlapResult2.score > scoreMin &&
        overlapResult2.length < s_i.size() && overlapResult2.length < s_j.size() &&
        (overlapResult.score < overlapResult2.score || overlapResult.length < overlapResult2.length)) {
        
        bool isOverlapAtEnd = false;
        auto mergedContig = mergeContigs(contig_j, contig_i, overlapResult2.length, isOverlapAtEnd, mode);

#ifdef DEBUGASS2
        std::cerr << "Merged contigs: new sequence overlap2 = " << std::get<0>(mergedContig) << "\n";
#endif

        return mergedContig;
    }

    // If no suitable overlap is found but there are common node IDs, calculate forced overlap
    if (!commonNodeIds.empty()) {
#ifdef DEBUGASS2
        std::cerr << "Forcing overlap due to common node IDs\n";
#endif
        auto findMatch = [this](const std::string& contig, const std::string& node_seq, bool from_start) -> size_t {
            size_t max_match_length = std::min(contig.size(), node_seq.size());
            size_t len = 0;

            if (from_start) {
                // Start from the beginning of the contig and the end of the node_seq
                for (size_t i = 0; i < max_match_length; ++i) {
                    if (!basesMatchWithDamage(contig[i], node_seq[node_seq.size() - max_match_length + i])) {
                        break; // Stop at the first mismatch considering DNA damage
                    }
                    len++;
                }
            } else {
                // Start from the end of the contig and the beginning of the node_seq
                for (size_t i = 0; i < max_match_length; ++i) {
                    if (!basesMatchWithDamage(contig[contig.size() - max_match_length + i], node_seq[i])) {
                        break; // Stop at the first mismatch considering DNA damage
                    }
                    len++;
                }
            }

            return len;
        };

        std::set<int> uniqueNodeIds(commonNodeIds.begin(), commonNodeIds.end());

        size_t forcedOverlapLength = 0;

        // Calculate forced overlap length based on unique node IDs and nodeSequenceMap
        for (int nodeID : uniqueNodeIds) {
            auto it = nodeSequenceMap.find(nodeID);
            if (it != nodeSequenceMap.end()) {
                const auto& node_tuple = it->second;
                const std::string& node_seq = std::get<0>(node_tuple);
                size_t node_length = std::get<1>(node_tuple);

                if (nodeID == *uniqueNodeIds.begin()) {
                    // For the first unique node ID, use findMatch
                    size_t matchLengthStart = findMatch(s_i, node_seq, true);
                    forcedOverlapLength += matchLengthStart;
                } else if (nodeID == *uniqueNodeIds.rbegin()) {
                    // For the last unique node ID, use findMatch
                    size_t matchLengthEnd = findMatch(s_j, node_seq, false);
                    forcedOverlapLength += matchLengthEnd;
                } else {
                    // For other node IDs, add the node length from the map
                    forcedOverlapLength += node_length;
                }
            }
        }

        if (forcedOverlapLength >= s_i.size()) {
            return contig_i;
        }

        bool isOverlapAtEnd = true;
        auto mergedContig = mergeContigs(contig_i, contig_j, forcedOverlapLength, isOverlapAtEnd, mode);

#ifdef DEBUGASS2
        std::cerr << "Forced overlap length " << forcedOverlapLength << "\n";
        std::cerr << "Forced merged contigs: new sequence = " << std::get<0>(mergedContig) << "\n";
#endif

        return mergedContig;
    }

    return std::nullopt;
}







const int assembly::run(int argc, char *argv[] , const string & cwdProg)
{

  assemblySetup();

  int lastOpt=1;
  int n_threads = 1;
  bool interleaved=false;
  string fastq1filename, fastq2filename, gamfilename, samplename;
  string tmpdir = "/tmp/";
  bool sbdirspecified = false;
  bool specifiedDeam = false;
  string sbdir = getFullPath(cwdProg+"../share/vgan/keelime_dir/");
  string dbprefixS = "keelime_db";
  string outputfilename = "keeOut";
  bool dbprefixFound = false; 
  int lengthToProf = 5;
  string pathName;
  bool pathNameFound = false;
  int lenMin = 10;
  int scoreMin = 15;
  bool unknownRef = false;
  bool unknown = false;
  string altMin = dbprefixS;
  bool minM = false;
  bool minimeta = false;
  bool useRemaining = false;
  string mode = "normal";
  int covBase = 1;

    string deam5pfreqE;
    string deam3pfreqE;

    ifstream file(cwdProg + "../share/vgan/damageProfiles/none.prof");

    if (file) {
        // Load deamination profiles
        deam5pfreqE = getFullPath(cwdProg + "../share/vgan/damageProfiles/none.prof");
        deam3pfreqE = getFullPath(cwdProg + "../share/vgan/damageProfiles/none.prof");
    } else {
        // Handle error when getFullPath fails
        //cerr << "Warning: deamination profile not found. Proceeding without the assumption of ancient damage." << endl;
        deam5pfreqE = ""; // Set to empty to handle it later
        deam3pfreqE = "";
    }

  for(int i=1;i<(argc);i++)
  {

    if(string(argv[i]) == "--keelime_dir" || string(argv[i]) == "--keelime-dir"){
        sbdir = argv[i+1];
        sbdirspecified=true; 
        if(sbdir.back() != '/'){sbdir += '/';}
        continue;
    }

    if(string(argv[i]) == "--dbprefix"){
        dbprefixS = argv[i+1];
        dbprefixFound=true;
        continue;
    }

    if(string(argv[i]) == "-fq1"){
        fastq1filename = argv[i+1];
        const int idx = fastq1filename.find_last_of("/");
        samplename=fastq1filename.substr(idx + 1);
        samplename=fastq1filename;
        if (fastq1filename.ends_with(".fa") || fastq1filename.ends_with(".fasta") || fastq1filename.ends_with(".fa.gz") || fastq1filename.ends_with(".fasta.gz"))
            {throw runtime_error("[keelime] Input file must be FASTQ, not FASTA");}
        continue;
                            }

    if(string(argv[i]) == "-fq2"){
        fastq2filename = argv[i+1];
        if (fastq2filename.ends_with(".fa") || fastq2filename.ends_with(".fasta") || fastq2filename.ends_with(".fa.gz") || fastq2filename.ends_with(".fasta.gz"))
            {throw runtime_error("[keelime] Input file must be FASTQ, not FASTA");}
        continue;
    }

    if(string(argv[i]) == "-g"){
        gamfilename = argv[i+1];
        samplename = gamfilename;
        continue;
                            }
    if(string(argv[i]) == "--deam5p"  ){
            deam5pfreqE=string(argv[i+1]);
        specifiedDeam=true;
            continue;
        }

    if(string(argv[i]) == "--deam3p"  ){
        deam3pfreqE=string(argv[i+1]);
    specifiedDeam=true;
        continue;
    }

    if(string(argv[i]) == "-t"){
      if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw runtime_error("[keelime] Error, invalid number of threads");}
      if (stoi(argv[i+1]) == -1) {n_threads = thread::hardware_concurrency();}
      else if (stoi(argv[i+1]) <= thread::hardware_concurrency()) {
              n_threads = stoi(argv[i+1]);
                                                                       }
    else {
           cerr << "[keelime] Warning, specified number of threads is greater than the number available. Using " << n_threads << " threads\n";
           n_threads = thread::hardware_concurrency();
         }
        continue;
                                   }

      if(string(argv[i]) == "-z"){
          tmpdir = argv[i+1];
          if (tmpdir.back() != '/') {tmpdir += '/';}
          continue;
                             }
        if(string(argv[i]) == "-M"){
            altMin = argv[i+1];
            continue;
                                }

      if(string(argv[i]) == "-i"){
            interleaved = true;
            if (fastq2filename != ""){throw runtime_error("[keelime] If interleaved option chosen, keelime expects only one FASTQ file");}
            continue;
                               }
        if(string(argv[i]) == "-o"){
            outputfilename = argv[i+1];
            continue;
        }
        if(string(argv[i]) == "-p" || string(argv[i]) == "--path"){
            pathName = argv[i+1];
            pathNameFound=true;
            continue;
        }
        if(string(argv[i]) == "-mL" || string(argv[i]) == "--minLength"){
            lenMin = stoi(argv[i+1]);
            continue;
        }
        if(string(argv[i]) == "-mS" || string(argv[i]) == "--lenMin"){
            scoreMin = stoi(argv[i+1]);
            continue;
        }
        if(string(argv[i]) == "-uR" || string(argv[i]) == "--unknownRef"){
            unknown = argv[i+1];
            unknownRef = true;
            continue;
        }
        if(string(argv[i]) == "-mM" || string(argv[i]) == "--miniMeta"){
            minM = argv[i+1];
            minimeta = true;
            continue;
        }
        if(string(argv[i]) == "-m" || string(argv[i]) == "--mode"){
            mode = argv[i+1];
            continue;
        }
        if(string(argv[i]) == "-uC" || string(argv[i]) == "--useContigs"){
            unknown = argv[i+1];
            useRemaining = true;
            continue;
        }
        if(string(argv[i]) == "-c" || string(argv[i]) == "--coverageBase"){
            covBase = stoi(argv[i+1]);
            continue;
        }

        
    }


    Euka ek;

    string altMinFull;
    string dbprefix              = sbdir + dbprefixS;
    if (altMin == "keelime_db"){
        altMinFull            = sbdir + dbprefixS;
    }else{
        altMinFull            = sbdir + altMin;
    }
    

    string ogfilename     = dbprefix+".og";
    string vgfilename     = dbprefix+".vg";
    string cladefilename  = sbdir + "keelime_db.clade";
    //string binsfilename   = sbdir + "keelime_db.bins";
    string gbwtfilename = dbprefix+".gbwt";

    if (specifiedDeam){
        cerr << "Computing probability of ancient damage based on provided damage profiles." << endl;
        
    }else{
        cerr << "[keelime] Warning: No damage profiles provided. Assuming no ancient damage has occured." << endl;
    }
    cerr << "[keelime] running in " << mode << " assembly mode." << endl;
    Damage dmg;
    dmg.initDeamProbabilities(deam5pfreqE,deam3pfreqE);
    
    if (!dbprefixFound){
        throw runtime_error("No database specified. Please choose a taxon of interest and specify it with the --dbprefix option. You can create a graph for you taxon of interested by using the make_graph_files.sh script.");
    }
    if(pathName.empty()){
        throw std::runtime_error("No reference path was specified. Please provide a refernce name from the list (use 'vg paths -L -x [your vg graph] with the -p option. Aboorting.");
    }

    vector<Clade *> * clade_vec = ek.load_clade_info(cladefilename, lengthToProf);
    if (clade_vec->empty()) {
    throw std::runtime_error("Error: The clade vector is empty. Unable to proceed. Check if the soibean.clade file is not empty.");
    }

    const vector<double> qscore_vec = ek.get_qscore_vec();
    cerr << "Reading in variation graph ..." << endl;

  
    auto [nodevector, minid, graph, node_path_matrix, path_names] = ek.readPathHandleGraph(ogfilename, 1, gbwtfilename, dbprefixS, clade_vec);
    if (nodevector.empty()) {
    throw std::runtime_error("Error: The nodevector is empty. Unable to proceed.");
    }
    if (path_names.empty()) {
    throw std::runtime_error("Error: The path_names vector is empty. Unable to proceed.");
    }

    //////// GIRAFFE MAPPING ///////
    string first_fifo = tmpdir + random_string(9);
    const char * fifo_A = first_fifo.c_str();
    mkfifo(fifo_A, 0666);
    pid_t wpid;
    const vg::subcommand::Subcommand* sc = NULL;

    pid_t pid1 = fork();
    
    if (pid1 == -1) {
    throw runtime_error("Error in fork");
    }
    
    if(pid1 == 0) {    // Child process
    // Code for writing to the FIFO
        if (gamfilename == "" ) {
            //cerr << "Mapping reads..." << endl;
            ek.map_giraffe(fastq1filename, fastq2filename, n_threads, interleaved, fifo_A, sc, tmpdir, sbdir, dbprefix, altMinFull);
            //cerr << "and done" << endl;
            exit(0);
        
        } else {
            // Redirect buffer in case of GAM input
            ifstream src(gamfilename);
            ofstream dst(fifo_A);
            dst << src.rdbuf();
            exit(0);
        }
    }
    
    int status;
    
    cerr << "...done!" << endl;
    
    ////// EXTRACTING GRAPH INFORMATION //////
    int startNode = 0;
    int endNode = 0;
    for (int i = 1; i<clade_vec->size(); i+=6){
        
            if(clade_vec->at(i)->name == dbprefixS){
                
                startNode = clade_vec->at(i)->snode;
                endNode = clade_vec->at(i)->enode;
            }else{
                startNode = minid;
                endNode = nodevector.size();
            }
        }


    GraphData new_graph_data;
    reindex_odgi_graph(graph, new_graph_data, startNode, endNode);
    


    
    ////////// ANALYSING GAM FILE ////////
    std::ifstream gam_file(fifo_A, std::ios::binary);
    vg::io::ProtobufIterator<vg::Alignment> iter(gam_file);

    cerr << "Analysing data" << endl;

    std::vector<frags> sequences;

    for (; iter.has_current(); iter.advance()) {
        const vg::Alignment& a = *iter;
        if (a.identity() != 0 && a.sequence().size() > 25) {
            std::vector<int> nodeIDs;
            std::vector<int> offsets;
            std::vector<int> baseQ;
            std::vector<int> coverage;

            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());
            // cerr << "read " << read_seq << " ";
            // cerr << "graph " << graph_seq << endl;
            for (int n = 0; n < a.path().mapping().size(); ++n) {


                const auto& mapping = a.path().mapping(n);
                nodeIDs.emplace_back(mapping.position().node_id());
                //cerr << mapping.position().node_id() << " ";
                offsets.emplace_back(mapping.position().offset());
                

                int node_coverage = 0;
                for (const auto& edit : mapping.edit()) {
                    node_coverage += edit.from_length();
                }
                coverage.emplace_back(node_coverage);
            }
            //cerr << endl;

            std::string sequence;
            string graphSeq;
            if (a.path().mapping(0).position().is_reverse()) {
                std::reverse(nodeIDs.begin(), nodeIDs.end());
                std::reverse(offsets.begin(), offsets.end());
                std::reverse(coverage.begin(), coverage.end());
                sequence = reverse_complement(read_seq);
                graphSeq = reverse_complement(graph_seq);
            } else {
                sequence = read_seq;
                graphSeq = graph_seq;
            }

            std::vector<std::array<double, 5>> probReadPostDamage;
			vector<array<int,5>>countsRead;

            for (int s = 0; s < sequence.size(); ++s) {
                baseQ.emplace_back(a.quality()[s]);
            }
            if (a.path().mapping(0).position().is_reverse()) {
                reverse(baseQ.begin(), baseQ.end());
            }
            int counter = 0;
            for(int s = 0; s <sequence.size(); ++s){


                if(graphSeq[s] == '-' || graphSeq[s] == 'N' || graphSeq[s] == 'S'){
                     counter++;
                }



                std::array<double, 5> probBasePostDamage = {0.0, 0.0, 0.0, 0.0, 0.0};
				array<int,5> countBase = {0,0,0,0,0};

                if(sequence[s] == 'N' || sequence[s] == 'S'){
                   
                    probBasePostDamage = {0.20, 0.20, 0.20, 0.20, 0.20};
					//countBase = {0,0,0,0,0};
					
                    probReadPostDamage.push_back(probBasePostDamage);
					countsRead.push_back(countBase);
					
				}else if(sequence[s] == '-'){
				    probBasePostDamage = {INDELERRORPROB/4, INDELERRORPROB/4, INDELERRORPROB/4, INDELERRORPROB/4, 1-INDELERRORPROB};
				    probReadPostDamage.push_back(probBasePostDamage);
			
					countBase = {0,0,0,0,1};
					countsRead.push_back(countBase);
                
                }else{

                    double probBasePreDamage[5] = {0.0};
					std::array<int, 5> countBase = {0, 0, 0, 0, 0};
					

                    for (int bpo = 0; bpo < 5; ++bpo) {
                        if ("ACGT-"[bpo] == sequence[s]) {
                            probBasePreDamage[bpo] = 1 - qscore_vec[static_cast<unsigned char>(baseQ[s])];
							countBase[bpo]++;
                        } else {
                            probBasePreDamage[bpo] = qscore_vec[static_cast<unsigned char>(baseQ[s])] / 4;
                        }
                    }
#ifdef DEBUGDAM
                    std::cerr << "Pre damage matrix" << std::endl;
                    for (int bpo = 0; bpo < 5; ++bpo) {
                        std::cerr << std::setprecision(14) << bpo << "\t" << probBasePreDamage[bpo] << std::endl;
                    }
#endif
                    
                   
                    for (int bpd = 0; bpd < 4; ++bpd) {
                        for (int bpo = 0; bpo < 4; ++bpo) {
                            
                            probBasePostDamage[bpd] += probBasePreDamage[bpo]*dmg.subDeamDiNuc[sequence.size()][s].p[bpo][bpd];
                            
                            
                        }
                    }
                    
                    // Handle indel error separately
                    for (int bpo = 0; bpo < 5; ++bpo) {
                        probBasePostDamage[4] += probBasePreDamage[bpo] * INDELERRORPROB; // This is separate to ensure proper handling
                        
                    }
                    for (int bpo = 0; bpo < 4; ++bpo) {
                        if ("ACGT"[bpo] == sequence[s]) {
                                probBasePostDamage[bpo] -= INDELERRORPROB;
                        } else {
                           continue;
                        }
                    }
                    


#ifdef DEBUGDAM
                    std::cerr << "After damage matrix " << std::endl;
                    for (int bpd = 0; bpd < 5; ++bpd) {
                        std::cerr << std::setprecision(14) << bpd << "\t" << probBasePostDamage[bpd] << std::endl;
                    }
#endif


                    countsRead.push_back(countBase);
                    probReadPostDamage.push_back(probBasePostDamage);
                }
            }
            pair<bool,int> cutbool = make_pair(false, 0);
            int dec = 0;
            if (counter >= 3){
                if (graph_seq[0] == '-'){
                    dec = 1;
                }else {
                    dec = 2;
                }
                cutbool = make_pair(true, dec);
            }

            std::string name = a.name();

            sequences.push_back({name, sequence, graphSeq, nodeIDs, offsets, baseQ, coverage, probReadPostDamage, countsRead, cutbool});
        }
    }

        
    /////////////////////////////////////////////


    // Sort the reads by the first node ID in each struct
    // Sorting using the compareByFirstNodeID with GraphData
    std::stable_sort(sequences.begin(), sequences.end(), [&](const frags& a, const frags& b) {
        return assembly::compareByFirstNodeID(a, b, new_graph_data);
    });
    convertToRYmerSpace(sequences);

    std::map<int, int> offsetsMap;
    map<int, int> covMap;
    

    // Fill offsetsMap and nodeCoverageMap
    for (const auto& seq : sequences) {
        for (size_t i = 0; i < seq.nodeIDs.size(); ++i) {
            int nodeID = seq.nodeIDs[i];
            
            offsetsMap[nodeID] = seq.offsets[i];
            covMap[nodeID] = seq.coverage[i];
        }
    }

    if (sequences.empty()) {
    throw std::runtime_error("Error: No reads are mapped. Unable to proceed.");
    }
    cerr << "Building overlap graph ..." << endl; 
    int n = sequences.size(); 
    int min_overlap_length = lenMin;
    
    GraphAss g(n);
    std::map<std::pair<int, int>, OverlapResult> overlap_map; // To store overlap lengths

    // Find overlaps and build the graph
    pair<unordered_map<int,vector<int>>, int> density_vec = initial_overlap(sequences, new_graph_data.node_depths);

    // FIND INITIAL OVERLAPS 
    find_overlaps(sequences, lenMin, g, overlap_map, density_vec.first, scoreMin, density_vec.second, lenMin, specifiedDeam);
    g.print_adjacency_list();
    cerr << "... done!" << endl;
    cerr << "Building contigs ..." << endl;
    int startVertex = 0;

    // CONSTRUCT FIRST SET OF CONTIGS BASED ON THE BREADTH FIRST SEARCH
    auto contigs = mergeAllPaths(g, overlap_map, sequences, 1, specifiedDeam, minimeta, new_graph_data, lenMin);

    auto filteredContigs = removeSubsetContigs(contigs);
    cerr << "... done!" << endl;
    

    std::vector<std::pair<int, size_t>> nodeLengths;
    std::map<int, std::tuple<std::string, size_t, int>> nodeSequenceMap; 

    auto path_handle = graph.get_path_handle(pathName);
    int order = 0; // Start numbering from 0
    nid_t last_node_id = 0; // Initialize to zero or an appropriate invalid value

    graph.for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
        handle_t handle = graph.get_handle_of_step(step);
        nid_t node_id = graph.get_id(handle);

        std::string node_sequence = graph.get_sequence(handle);
        size_t node_length = node_sequence.length();
        nodeSequenceMap[node_id] = std::make_tuple(node_sequence, node_length, order); // Store sequence, length, and order
        nodeLengths.push_back({node_id, node_length});

        order++; // Increment order for the next node
    });
    

    std::vector<std::pair<int, size_t>> nodeLengths2;
    std::map<int, std::tuple<std::string, size_t, int>> nodeSequenceMap2; 

    // Iterate over all paths in the graph
    graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        int order = 0; // Start numbering from 0 for each path
        nid_t last_node_id = 0; // Initialize to zero or an appropriate invalid value
        
        graph.for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
            handle_t handle = graph.get_handle_of_step(step);
            nid_t node_id = graph.get_id(handle);

            std::string node_sequence = graph.get_sequence(handle);
            size_t node_length = node_sequence.length();
            nodeSequenceMap2[node_id] = std::make_tuple(node_sequence, node_length, order); // Store sequence, length, and order
            nodeLengths2.push_back({node_id, node_length});

            order++; // Increment order for the next node
        });
    });

    
    /////////// RUNNING MERGING OF CONTIGS ///////////////////////
    cerr << "Merging overlapping contigs ..." << endl;
    auto mergedContigs = mergeOverlappingContigs(filteredContigs, 1, nodeSequenceMap2, specifiedDeam, mode, lenMin);

    bool changesMade;
    int iteration = 0; 

    do {
        auto previousSize = mergedContigs.size();
        mergedContigs = mergeOverlappingContigs(mergedContigs, 1, nodeSequenceMap2, specifiedDeam, mode, lenMin);
        changesMade = (mergedContigs.size() < previousSize);


    } while (changesMade);


    if (!mergedContigs.empty()) {
        const auto& finalContig = mergedContigs[0];
        const auto& finalScoringMatrix = std::get<1>(finalContig);
		const auto& finalCountMatrix = get<4>(finalContig);
    
    }
    cerr << "... done!" << endl;
   /////////////////// ORDERING CONTIGS ////////////////////////
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> orderedContigs;

    std::cerr << "\nOrdering contigs...\n";

    // Primary ordering based on the front node ID
    for (size_t n = 0; n < nodeLengths.size(); ++n) {
        int key = nodeLengths[n].first;
        //std::cerr << "Looking for contig with front node ID: " << key << "\n";

        auto it_contig = std::find_if(mergedContigs.begin(), mergedContigs.end(), [key](auto& contig) {
            //std::cerr << "Checking contig with front node ID: " << std::get<2>(contig).front() << "\n";
            return std::get<2>(contig).front() == key;
        });

        if (it_contig != mergedContigs.end()) {
            orderedContigs.push_back(*it_contig);
            //std::cerr << "Found and added contig with front node ID: " << key << "\n";
        } else {
            //std::cerr << "Did not find contig with front node ID: " << key << "\n";
        }
    }

    std::unordered_set<int> addedKeys;
    for (const auto& contig : orderedContigs) {
        addedKeys.insert(std::get<2>(contig).front());
    }

    // Collect unmatched contigs
    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> unmatchedContigs;
    for (const auto& contig : mergedContigs) {
        if (addedKeys.find(std::get<2>(contig).front()) == addedKeys.end()) {
            unmatchedContigs.push_back(contig);
            //std::cerr << "Contig with front node ID " << std::get<2>(contig).front() << " was not added to orderedContigs.\n";
        }
    }

    // Insert unmatched contigs based on their range
    for (const auto& unmatchedContig : unmatchedContigs) {
        int rangeUnmatched = std::get<2>(nodeSequenceMap2.at(std::get<2>(unmatchedContig).front()));
        //std::cerr << "Inserting unmatched contig with range: " << rangeUnmatched << "\n";

        auto it = std::find_if(orderedContigs.begin(), orderedContigs.end(), [&nodeSequenceMap2, rangeUnmatched](const auto& contig) {
            int range = std::get<2>(nodeSequenceMap2.at(std::get<2>(contig).front()));
            return range > rangeUnmatched;
        });

        orderedContigs.insert(it, unmatchedContig);
    }

    std::cerr << "Contigs ordered: " << orderedContigs.size() << "\n";
    

   ////////////////// FINAL MERGING CHECK ////////////////////////
    std::cerr << "\nPerforming final merging check...\n";

    std::vector<std::tuple<std::string, std::vector<std::unordered_map<char, double>>, std::vector<int>, std::string, std::vector<std::unordered_map<char, int>>, std::pair<bool, int>>> finalMergedContigs;
    std::unordered_set<size_t> mergedIndexSet;

    bool mergeOccurred;
    do {
        mergeOccurred = false;
        finalMergedContigs.clear();
        mergedIndexSet.clear();

        for (size_t i = 0; i < orderedContigs.size(); ++i) {
            // Access the vector of ints (the 3rd element of the tuple)
            std::vector<int>& intVector = std::get<2>(orderedContigs[i]);

            if (mergedIndexSet.find(i) != mergedIndexSet.end()) continue; // Skip already merged contigs

            auto& contig = orderedContigs[i];

            // Check for the next contig in sequence
            size_t next_i = i + 1;
            if (next_i < orderedContigs.size()) {
                auto& next_contig = orderedContigs[next_i];

                if (mergedIndexSet.find(next_i) == mergedIndexSet.end()) {
                    auto mergedContigOpt = checkAndMergeContigs(contig, next_contig, nodeSequenceMap, lenMin, scoreMin, specifiedDeam, mode);
                    if (mergedContigOpt) {
                        contig = *mergedContigOpt;
                        mergedIndexSet.insert(next_i);
                        mergeOccurred = true;
                    }
                }
            }

            finalMergedContigs.push_back(contig);
            mergedIndexSet.insert(i);
        }

        orderedContigs = finalMergedContigs;

    } while (mergeOccurred);

    std::cerr << "Final merging check completed. Total contigs after merge: " << finalMergedContigs.size() << "\n";

    cerr << "\nWriting contigs to fasta: " << outputfilename + "Contig.fa.gz" << endl;
    string contig = outputfilename + "Contig.fa";
    writeContigsToFasta(finalMergedContigs, contig);

    /////////////////// CONSENSUS CALLING ////////////////////////
    std::cerr << "\nCalling consensus sequence...\n";
    std::string fasta;
    int counter = 0, con = 0;
    std::unordered_set<size_t> usedContigs;

    bool isInContig = false;
    int lastNodeId = -1;
    std::string contig_seq;
    
    // extend the existing contig by checking where the best overlap is. This is assuming that the normal overlap function didnt have results ( this function is stored here for convinience)
    auto findMatch = [this](const std::string& contig, const std::string& node_seq, bool from_start) -> int {
        size_t max_match_length = std::min(contig.size(), node_seq.size());
        size_t len = 0;
	// checking witch way to merge 
        if (from_start) {
            // Start from the beginning of the contig and the end of the node_seq
            for (size_t i = 0; i < max_match_length; ++i) {
                if (!basesMatchWithDamage(contig[i], node_seq[node_seq.size() - max_match_length + i])) {
                    break; // Stop at the first mismatch considering DNA damage
                }
                len++;
            }
        } else {
            // Start from the end of the contig and the beginning of the node_seq
            for (size_t i = 0; i < max_match_length; ++i) {
                if (!basesMatchWithDamage(contig[contig.size() - max_match_length + i], node_seq[i])) {
                    break; // Stop at the first mismatch considering DNA damage
                }
                len++;
            }
        }

        return (len == 0) ? -1 : static_cast<int>(len);
    };


    size_t currentPosition = 0;  // Tracks the continuous position in the final FASTA sequence - This should have a maximum seq length 
    std::vector<std::unordered_map<char, double>> finalScoringMatrix;  // Full scoring matrix for the entire sequence
	std::vector<std::unordered_map<char, int>> finalCountMatrix; 
    
    // this is going through each node of the nodeLengths vector (all node IDs for the chosen path) 
    for (size_t n = 0; n < nodeLengths.size(); ++n) {
        int key = nodeLengths[n].first; // node ID
        size_t value = nodeLengths[n].second; // node length 
        // looking for a contig based on the node ID key (again not optimal function storage)
        auto findContig = [&](int key) {
            return std::find_if(finalMergedContigs.begin(), finalMergedContigs.end(), [key](auto& contig) {
                return std::find(std::get<2>(contig).begin(), std::get<2>(contig).end(), key) != std::get<2>(contig).end();
            });
        }; 

        auto it_contig = findContig(key);
        // contig is valid and not used before!
        if (it_contig != finalMergedContigs.end() && usedContigs.find(it_contig - finalMergedContigs.begin()) == usedContigs.end()) {
            auto& contig = *it_contig;
            contig_seq = std::get<0>(contig); // Sequence of the contig
            const auto& scoringMatrix = std::get<1>(contig); //scoring matrix of the contig
			const auto& countMatrix = get<4>(contig); // count matrix of the contig
            std::pair<bool, int> cutbool = std::get<5>(contig); // do we have trustworthy endings on node IDs?

            int nodeID_to_use = key;
            for (int nodeID : std::get<2>(contig)) {
                if (std::find_if(nodeLengths.begin(), nodeLengths.end(), [nodeID](const std::pair<int, size_t>& p) { return p.first == nodeID; }) != nodeLengths.end()) {
                    nodeID_to_use = nodeID;
                    break;
                }
            }

            std::string node_seq = std::get<0>(nodeSequenceMap.at(nodeID_to_use));
            int overlap_start = findMatch(contig_seq, node_seq, true); // uses the findMatch function to check where on the node sequence we place the existing contig and the new contig 
            
            int NsToAdd_start = 0; // how many Ns inbetween the two contigs
            if (overlap_start == -1) {
                NsToAdd_start = 0; // if all fails and damage or sth else results in no matches just don't add Ns
            }else{
                NsToAdd_start = node_seq.size() - overlap_start; 
                if(NsToAdd_start < 0 ){
                    NsToAdd_start = 0;
                }
            }

            // Add Ns to fasta and zeros to scoring matrix if not restricted by cutbool
            if (cutbool.first) {
                NsToAdd_start = 0; // if we don't trust the node ID in the first place just put no N inbetween. 
            }
            //std::cerr << "Node ID: " << key << ", Adding Ns (start): " << NsToAdd_start << std::endl;

            fasta += std::string(NsToAdd_start, 'N');
            finalScoringMatrix.resize(currentPosition + NsToAdd_start, {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}});
			finalCountMatrix.resize(currentPosition + NsToAdd_start, {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'-', 0}});
            currentPosition += NsToAdd_start;

            // Add contig sequence to fasta and extend the scoring matrix with contig's matrix
            fasta += contig_seq;
            finalScoringMatrix.insert(finalScoringMatrix.end(), scoringMatrix.begin(), scoringMatrix.end());
			finalCountMatrix.insert(finalCountMatrix.end(), countMatrix.begin(), countMatrix.end());
            currentPosition += contig_seq.size();

            // Mark the contig as used
            usedContigs.insert(it_contig - finalMergedContigs.begin());

            lastNodeId = std::get<2>(contig).back();
            isInContig = true;
        } else { // if the node ID is not part of the contig
            if (isInContig) { // checking for overlaps again (did we miss an overlap or a subset)
                if (key == lastNodeId) { // find the ending of the contig and merge on the node ID (safty measure)
                    std::string node_seq = std::get<0>(nodeSequenceMap.at(key));
                    int overlap = findMatch(contig_seq, node_seq, false);
                    
                    int NsToAdd = 0;
                    if (overlap == -1) {
                        NsToAdd = 1;
                    }else{
                        NsToAdd = node_seq.size() - overlap;
                        if(NsToAdd < 0 ){
                            NsToAdd = 0;
                        }
                    }
                    if (std::pair<bool, int> cutbool = std::get<5>(*it_contig); cutbool.first) {
                        NsToAdd = 1;
                    }
                    //std::cerr << "Node ID: " << key << ", Adding Ns (overlap): " << NsToAdd << std::endl;

                    fasta += std::string(NsToAdd, 'N');
                    finalScoringMatrix.resize(currentPosition + NsToAdd, {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}});
					finalCountMatrix.resize(currentPosition + NsToAdd, {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'-', 0}});
                    currentPosition += NsToAdd;
                    isInContig = false;
                }
            } else {
                //std::cerr << "Node ID: " << key << ", Adding Ns (not in contig): " << value << std::endl;
                if(unknownRef){
		    cerr << "Unknown reference specifed: The progam will not try to bridge unknown Node IDs between contigs with 'N'. " << endl;
                    //fasta += std::string(1, 'N');
                }else{
                    finalScoringMatrix.resize(currentPosition + value, {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}});
					finalCountMatrix.resize(currentPosition + value, {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'-', 0}});
                    currentPosition += value;
                    fasta += string(value, 'N');
                }
            }
        }
    }

    // Add any remaining unused contigs - this is a senario that would only happen if we have a large divergence and node IDs couldnt bridge all the contigs. 

    if(useRemaining){
    for (size_t i = 0; i < finalMergedContigs.size(); ++i) {
        if (usedContigs.find(i) == usedContigs.end()) {
            auto& contig = finalMergedContigs[i];
            contig_seq = std::get<0>(contig);
            const auto& scoringMatrix = std::get<1>(contig);
			const auto& countMatrix = get<4>(contig);
            std::pair<bool, int> cutbool = std::get<5>(contig);

            // Add Ns to the start if necessary and not restricted by cutbool
            int NsToAdd_start = 0;
            if (isInContig) {
                auto iter = nodeSequenceMap.find(std::get<2>(contig).front());
                if (iter != nodeSequenceMap.end()) {
                    std::string node_seq = std::get<0>(iter->second);  // Access the found element
                    NsToAdd_start = findMatch(contig_seq, node_seq, true);
                    if(NsToAdd_start < 0){
                        NsToAdd_start = 0;
                    }
                } else {
                    NsToAdd_start = 0;  // Default to adding one 'N' if not found
                }
            }
            if (cutbool.first) {
                NsToAdd_start = 1;
            }
            std::cerr << "Adding remaining contig with Ns (start): " << NsToAdd_start << std::endl;

            fasta += std::string(NsToAdd_start, 'N');
            finalScoringMatrix.resize(currentPosition + NsToAdd_start, {{'A', log(0.20)}, {'C', log(0.20)}, {'G', log(0.20)}, {'T', log(0.20)}, {'-', log(0.20)}});
			finalCountMatrix.resize(currentPosition + NsToAdd_start, {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}, {'-', 0}});
            currentPosition += NsToAdd_start;

            // Add contig sequence to fasta and extend the scoring matrix with contig's matrix
            fasta += contig_seq;
            finalScoringMatrix.insert(finalScoringMatrix.end(), scoringMatrix.begin(), scoringMatrix.end());
			finalCountMatrix.insert(finalCountMatrix.end(), countMatrix.begin(), countMatrix.end());
            currentPosition += contig_seq.size();
        }
    }
    }
	
	for (size_t i = 0; i < fasta.size(); ++i) {
	    char base = fasta[i];
	    int maxCount = -1;
	    char mostFrequentBase = 'N';
	    int totalCount = 0;

	    // Sum counts and determine most frequent base
	    for (const auto& [b, count] : finalCountMatrix[i]) {
	        totalCount += count;
	        if (count > maxCount) {
	            maxCount = count;
	            mostFrequentBase = b;
	        }
	    }

	    // Apply conBase threshold: mask with 'N' if total count is too low
	    if (totalCount < covBase) {
	        std::cerr << "Masking position " << i << " as 'N' due to low coverage ("
	                  << totalCount << " < " << covBase << ")\n";
	        fasta[i] = 'N';  // <-- Update the fasta sequence directly
	        continue;        // Skip mismatch check if masked
	    }

	    // Check for mismatch
	    if (base != mostFrequentBase) {
	        std::cerr << "Mismatch at position " << i << ": "
	                  << "sequence has '" << base
	                  << "', but most frequent base is '" << mostFrequentBase << "'\n";
	        // Optionally update to most frequent base:
	        fasta[i] = mostFrequentBase;
	    }

	    std::cerr << "Total count at position " << i << " = " << totalCount << "\n";
	}
	
	

    std::cerr << "... consensus sequence saved to " << outputfilename + "Consensus.fa.gz" << std::endl;
    std::string cons = outputfilename + "Consensus.fa";
    saveToFastaGz(fasta, cons);
    std::string scores = outputfilename + "ProbabilityMatrix.tsv";
	string counts = outputfilename + "CountMatrix.tsv";
    printScoringMatrix(finalScoringMatrix, scores);
	printCountMatrix(finalCountMatrix, counts);
    std::cerr << "Probability matrix saved to " << scores << std::endl;
	cerr << "Count matrix saved to " << counts << endl;




return 0;
}
