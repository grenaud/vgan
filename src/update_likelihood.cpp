#include "HaploCart.h"
#include "vgan_utils.h"
#include "process_mapping.cpp"

const vector<long double> Haplocart::update(const int i, const map<const string, int> &pangenome_map, const vector<NodeInfo *> &nodevector,
                                       const AlignmentInfo* read_info, vector<long double> log_likelihood_vec,
                                       const vector<double> &qscore_vec, const vector<double> &mappabilities, const int nbpaths, const bool quiet,
                                       bool use_background_error_prob, double background_error_prob, const vector<double> &incorrect_mapping_vec,
                                       int n_reads, const int minid, const bool is_consensus_fasta, const int n_threads,
                                       const bdsg::ODGI &graph) noexcept {

const vector<long double> ret = Haplocart::update_likelihood(pangenome_map, read_info, nodevector, log_likelihood_vec, qscore_vec, mappabilities,
                                                       nbpaths, quiet, use_background_error_prob, background_error_prob,
                                                       incorrect_mapping_vec, minid, is_consensus_fasta, n_threads, graph);

return ret;
}

inline const vector<long double> Haplocart::update_likelihood(const map<const string, int> &pangenome_map, const AlignmentInfo* read_info,
                                 const vector<NodeInfo *> &nodevector,
                                 vector<long double> log_likelihood_vec, const vector<double> &qscore_vec,
                                 const vector<double> &mappabilities, const int nbpaths, const bool verbose, bool use_background_error_prob,
                                 const double &background_error_prob, const vector<double> &incorrect_mapping_vec,
                                 const int minid,
                                 const bool is_consensus_fasta, const int n_threads, const bdsg::ODGI &graph) noexcept
{

const auto path = read_info->path;
const tuple<string, string, vector<int>> seq_tuple = reconstruct_graph_sequence(graph, path, read_info->seq);
string algnseq = get<1>(seq_tuple);
int position_in_read = 0;

for(int i=0;i<path.mapping().size();++i){
     vector<int> quality_scores;
     const vector<int> mppg_sizes = get<2>(seq_tuple);
     string graph_seq = get<0>(seq_tuple).substr(position_in_read, mppg_sizes[i]);
     string read_seq = algnseq.substr(position_in_read, mppg_sizes[i]);
     long unsigned int j;
     //#pragma omp parallel for num_threads(n_threads) private(j, log_likelihood_vec) schedule(static)
     for(long unsigned int j=position_in_read;j<position_in_read + algnseq.size();++j) {
         const int qscore = int(read_info->quality_scores[j]);
         if (qscore >= 90) {use_background_error_prob = true;}
         quality_scores.emplace_back(qscore);
                                                                                       }
     position_in_read += read_seq.size();
     log_likelihood_vec = process_mapping(pangenome_map, read_info, log_likelihood_vec, path.mapping()[i], algnseq, nodevector, quality_scores,
                                          graph_seq, qscore_vec, mappabilities, nbpaths, use_background_error_prob, background_error_prob,
                                          incorrect_mapping_vec, minid, is_consensus_fasta, n_threads);
                               }


     return log_likelihood_vec;
 }
