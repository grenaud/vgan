#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <gzstream.h>
#include "libgab.h"
#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "readVG.h"
#include "readGAM.h"
#include "subcommand/subcommand.hpp"
#include "bdsg/odgi.hpp"

using namespace std;




class Haplocart{
private:

public:

const double getBaseFrequency(const char base);

/** @brief Read in the graph.
 *
 *  Deserialize the graph, which is in ODGI format.
 */

void readPathHandleGraph (shared_ptr<Trailmix_struct> &dta);

/** @brief Update likelihood_vec for a given read.
 *
 *  Given a read, update our vector of log-likelihoods appropriately.
 */

inline const vector<double> update_likelihood(shared_ptr<Trailmix_struct> &dta, const AlignmentInfo* read_info,
                                 vector<double> log_likelihood_vec) noexcept;

/** @brief Update log likelihood vector for each thread.
 *
 *  Same as above (used for multithreading).
 */

const vector<double> update(shared_ptr<Trailmix_struct> &dta, const int i, vector<double> log_likelihood_vec) noexcept;

// Process mapping

/**
 * @brief Penalize unsupporting haplotypes
 *
 * If a read maps to a node unsupported by a putative haplotype, it incurs a penalty commensurate with 3/4 mismatch and 1/4 match
 *
 */

inline const double get_log_lik_if_unsupported(const vg::Mapping &mppg, const vector<int> &quality_scores);


/**
 * @brief Update log likelihood vector upon observing a mapping to a single node
 *
 * Given a mapping to a node, we do a Bayesian update on our haplotypes
 */


inline const vector<double> process_mapping(shared_ptr<Trailmix_struct> &dta, const AlignmentInfo* read_info, vector<double> log_likelihood_vec,
const vg::Mapping &mppg, string &mapping_seq,
const vector<int> &quality_scores, const string &graph_seq, unsigned int path_idx) noexcept;

// Load

/**
 * @brief Load pangenome map.
 *
 * Load precomputed map of node -> pangenome base for our mutational model.
 *
 */

void load_pangenome_map(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Load mappabilities
 *
 * Load GenMap mappability scores
 *
 */

void load_mappabilities(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Load a vector of path names, in order
 *
 */

void load_path_names(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Load path supports into memory.
 *
 * Load a boolean matrix of dimension n_paths, n_nodes
 * where each indicator variable tells us where a given path
 * goes traverses a given node.
 *
 */

// const vector<vector<bool>> load_path_supports(shared_ptr<Trailmix_struct> &dta);
//void load_path_supports(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Compute per-base probability of no sequencing error.
 *
 * For each base in the given read, compute the probability that no sequencing error has occurred.
 * This is inferred from the PHRED-encoded quality scores.
 */

const vector<double> get_p_no_seq_error_mapping(shared_ptr<Trailmix_struct> &dta, string &mapping_seq, const vector<int> &quality_scores,
            const string &graph_seq);

/**
 * @brief Check if mutation is a transversion
 */

inline constexpr bool transversion(const char &base1, const char &base2);

/**
 * @brief Check if value is in a given range
 */

const bool inRange(const unsigned &low, const unsigned &high, const unsigned &x);

/**
 * @brief Compute per-base probability of mutation.
 *
 * For each base in the given read, compute the probability that a mutation has occurred.
 * This is based solely upon the estimated mutation rate of the considered locus of the mitogenome.
 *
 */

const double get_p_obs_base(const size_t pangenome_base, const char mapping_base, const char graph_base, \
                                 const double epsilon, const size_t generations, const size_t base_index_on_read, const size_t Lseq, shared_ptr<Trailmix_struct> &dta);

const pair<vector<string>, vector<string>> read_fasta(const string &fastafilename);

/**
 * @brief If we don't have GAM input, map whatever we have
 *
 * For FASTA, FASTQ, BAM, or CRAM input, map with vg giraffe
 * to generate a GAM file.
 *
 */


void map_giraffe(string &fastaseq, string &fastq1filename, string &fastq2filename, const int n_threads, bool interleaved,
                        double background_error_prob, const string & samplename, const char * fifo_A, const vg::subcommand::Subcommand* sc,
                        const string & tmpdir, const string &graph_dir_path, const bool quiet, const string &deam3pfreqE, const string &deam5pfreqE, const string &posterior_threshold="1.01");

/**
 * @brief Get dummy quality score from background error probability
 *
 * For consensus FASTA input we need to create temporary FASTQ files with dummy quality scores.
 * This function returns the quality score corresponding to a given background error probability
 *
 */

const char get_dummy_qual_score(const double &background_error_prob);


/**
 * @brief Remove PCR duplicates from a GAM file.
 *
 * Remove PCR duplicates from a GAM file.
 *
 */

shared_ptr<vector<AlignmentInfo*>> remove_duplicates(shared_ptr<vector<AlignmentInfo*>> sorted_algnvector, const int n_threads, const bool &quiet);

/**
 * @brief Helper function to update vector of booleans when encountering a paired forward-end read.
 *
 * Helper function to update vector of booleans when encountering a paired forward-end read.
 *
 */

const vector<bool> update_for_paired_forward(vector<bool> &is_dup, shared_ptr<vector<shared_ptr<AlignmentInfo>>> sorted_algnvector,
const long unsigned int i, const int &n_threads);

/**
 * @brief  Helper function to update vector of booleans when encountering a single-end read.
 *
 * Helper function to update vector of booleans when encountering a single-end read.
 *
 */

vector<bool> update_for_single(vector<bool> &is_dup, shared_ptr<vector<shared_ptr<AlignmentInfo>>> sorted_algnvector, long unsigned int i, const int n_threads);

/**
 * @brief Sort GAM file
 *
 * Call vg gamsort to sort the GAM file before performing duplicate removal.
 *
 */

void gamsort(const int n_threads, const bool interleaved, char const * fifo_B, char const * fifo_C, const string & tmpdir);

/**
 * @brief Filter sorted GAM file
 *
 * Call vg filter to filter reads
 *
 */


void filter(const int n_threads, const bool interleaved, char const * fifo_A, char const * fifo_B);


/**
 * @brief Get the background frequency of a base.
 *
 * Return the frequency of a given nucleobase within the graph, i.e. the human mitogenome.
 *
 */

const double get_background_freq(const char &base);

/**
 * @brief Autodetect input file format
 *
 *
 */

const string autodetect(const string & input_file);

/**
 * @brief Get the probability that an alignment was incorrectly mapped.
 *
 * Given a mapping quality score (in PHREDs), return the corresponding probability of incorrect mapping.
 *
 */

constexpr double get_p_incorrectly_mapped(const int &Q);

/**
 * @brief Precompute incorrect mapping probabilities
 *
 * Precompute incorrect mapping probabilities for each possible mapping quality score.
 *
 */

void precompute_incorrect_mapping_probs(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Compute confidence in our predicted haplogroup assignment
 *
 * For the predicted haplotype as well as all assignments ancestral to the predicted haplotype, compute the posterior probability
 * of that haplotype and all other haplotypes subsumed within its clade.
 *
 */

void get_posterior(const vector<double> &final_vec, const vector<string> path_names,
                               map<string, vector<string>> parents, map<string, vector<string>> children,
                               const string &samplename, const string &predicted_haplotype, const string &posteriorfilename, const bool webapp);

/**
 * @brief Write posterior output file
 *
 */

void write_posterior_log(const string &samplename, const string &posteriorfilename, vector<string> clade_vec, vector<double> confidence_vec,
                         const bool &webapp);


/**
 * @brief Load file of parents
 *
 * Load file of parents per haplotype, going back to the mt-MRCA. The file was creating via modifying a script in the mixemt github repository
 *
 */

void load_parents(shared_ptr<Trailmix_struct> &dta);

/**
 * @brief Load file of children
 *
 * Load file of children for each haplotype. The file was creating via modifying a script in the mixemt github repository
 *
 */

void load_children(shared_ptr<Trailmix_struct> &dta);


/**
 * @brief Find total posterior probability that the true haplotype lies within a given clade.
 *
 * Find total posterior probability that the true haplotype lies within a given clade.
 */

const vector<double> get_posterior_of_clade(vector<double> all_tops, vector<double> final_vec, set<string> haplotypes,
                                                 map<string, vector<string>> children, const double total_ll,
                                                 const vector<string> path_names, bool initial);

const set<string> get_children(set<string> haplotypes, map<string, vector<string>> children);

/**
 * @brief Add log probabilities in a numerically safe way.
 *
 * Compute log(exp(x) + exp(y) + ...) while avoiding underflows or loss of precision.
 */

const double sum_log_likelihoods(vector<double> all_top);

/**
 * @brief Convert FASTA to FASTQ
 *
 * Since Giraffe only takes FASTQ input we synthetically create such input with dummy quality scores.
 */

const string fa2fq(const string & fastaseq, const char & dummyqualscore, const string & tmpdir);

/**
 * @brief Write a single FASTQ read to the FIFO
 */


void write_fq_read(auto & dummyFASTQFile, int offset, const int window_size, const string &fastaseq, char dummyqualscore);


bool contains_no_inf(const std::vector<double>& v);

void load_path_supports(shared_ptr<Trailmix_struct> &dta);

void modifyPathNameInPlace(shared_ptr<Trailmix_struct> &dta, string &path_name, bool first=false){
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

    if(first){
    dta->originalPathNames[path_name] = original_path_name;
             }

}

void readPHG(shared_ptr<Trailmix_struct> &dta) {

    cerr << "[TrailMix] Deserializing graph: " << dta->graphfilename << endl;
    cerr << "GRAPH DIR: " << dta->graph_dir << endl;
    dta->graph.deserialize(dta->graph_dir+dta->graphfilename);
    cerr << "[TrailMix] Done deserializing" << endl;
    dta->minid = dta->graph.min_node_id();
    dta->maxid = dta->graph.max_node_id();

    // create path support file
    size_t N_nodes = dta->graph.get_node_count();

    unsigned int p_index = 0;
    dta->graph.for_each_path_handle([&](const handlegraph::path_handle_t &path_handle) {
    string path_name = dta->graph.get_path_name(path_handle);
    modifyPathNameInPlace(dta, path_name, true);

size_t pos = 0;
        dta->path_names.emplace_back(path_name);
        p_index++;
    });

    const size_t N_paths = dta->path_names.size();

    vector<vector<bool>> node_path_matrix(N_nodes, vector<bool>(N_paths, false));

    for (size_t path_id = 0; path_id < N_paths; ++path_id) {
        // Get the nodes in the path
        gbwt::vector_type path_nodes = dta->gbwt->extract(path_id);
        for (const auto& node_id : path_nodes) {
            //cout << node_id << " ";
            // Convert the GBWT node ID to the node handle in the graph
            bdsg::handle_t handle = dta->graph.get_handle(node_id);
            // Determine the index in the node_path_matrix
            int64_t index = dta->graph.get_id(handle) - 1;
            //cout << index << endl;
            if (index >= 0 && index < N_nodes) {
                node_path_matrix[index][path_id] = true;
                                //cout << "I get positive too " << endl;
            }

        }
        //cout << endl;

    }

    int nbpaths = dta->path_names.size();

    for(int64 i=dta->minid;i<=dta->maxid;++i){
        //cerr << i << endl;
        NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,0);
        const auto nodehandle = dta->graph.get_handle(i);
        string seqtoadd = dta->graph.get_sequence(nodehandle);
        nodetoadd->seq = seqtoadd;
        long unsigned int j;
        for (j = 0; j < nbpaths; ++j) {
            //nodetoadd->pathsgo[j] = node_path_matrix[i - 1][j];
        }
        dta->nodevector.emplace_back(move(nodetoadd));

        //std::this_thread::sleep_for (std::chrono::milliseconds(10));
    }

    assert(dta->nodevector.size() != 0);
    assert(dta->minid != 0);
    // Check if at least one true value exists
    bool hasTrueValue = false;
    for (const auto& row : node_path_matrix) {
        for (bool value : row) {
            if (value) {
                hasTrueValue = true;
                                break;
            }
        }
        if (hasTrueValue) {
            break;
        }
    }

    assert(dta->path_names.size() != 0);

   return;
}



Haplocart();
Haplocart(const Haplocart & other) = delete;
~Haplocart();
Haplocart & operator= (const Haplocart & other) = delete;

const string usage() const;
void run(int argc, char *argv[], shared_ptr<Trailmix_struct> &dta);

};


