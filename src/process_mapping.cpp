#include "HaploCart.h"
#include "miscfunc.h"

inline const double Haplocart::get_log_lik_if_unsupported(const vg::Mapping &mppg, const vector<int> &quality_scores)
{
    // If a node is unsupported by a path, penalize with approximately 3/4 mismatches and 1/4 matches.
    bool mismatch = true;
    int counter = 0;
    double ret = 0;
    for (const int &Q: quality_scores) {
        if (mismatch) {
              ret += log(get_p_seq_error(Q));
                      }
        else {
              ret += log(1 - get_p_seq_error(Q));
             }

        counter += 1;
        if (counter % 4 == 4) {mismatch = false;}
        else {mismatch = true;}
                                       }

    return ret;
}

inline const vector<long double> Haplocart::process_mapping(const map<const string, int> &pangenome_map, const AlignmentInfo * read_info,
vector<long double> log_likelihood_vec, const vg::Mapping &mppg, string &mapping_seq, const vector<NodeInfo *> &nodevector,
const vector<int> &quality_scores, string &graph_seq, const vector<double> &qscore_vec, const vector<double> &mappabilities,
const int &nbpaths, bool &use_background_error_prob, const double &background_error_prob, const vector<double> &incorrect_mapping_vec,
const int &minid, const bool &is_consensus_fasta, const int & n_threads) noexcept
{

const int pangenome_base = pangenome_map.at(to_string(mppg.position().node_id()));

const double mappability = mappabilities[pangenome_base];
vector<double> background_freqs;

long double p_correctly_mapped;

if (is_consensus_fasta == false) {
   p_correctly_mapped = (1-incorrect_mapping_vec[read_info->mapping_quality]) * mappability;
                                 }

unsigned int i;
//#pragma omp parallel for private(i) num_threads(n_threads)
for (i=0; i != mapping_seq.size(); ++i) {
   background_freqs.emplace_back(Haplocart::get_background_freq(mapping_seq[i]));
}

const vector<long double> p_no_seq_error = Haplocart::get_p_no_seq_error_mapping(mapping_seq, quality_scores, graph_seq,
                                                                                 qscore_vec, use_background_error_prob,
                                                                                 background_error_prob, n_threads);

	for(int i=0;i!=nbpaths;++i){
          const bool is_supported = nodevector.at(mppg.position().node_id()-minid)->pathsgo[i];

	  if (is_supported) {
                  long double log_lik_if_mapped = 0;
		  for (long unsigned int j = 0; j != graph_seq.size(); ++j)
		      {
                                // Discard ambiguous bases
                                if (graph_seq[j] == 'N' || mapping_seq[j] == 'N') {continue;}
                                if ( !isValidDNA(graph_seq[j]) || !isValidDNA(mapping_seq[j])) {continue;}

                                // Get the probability of observing the base given that it was generated from the putative haplotype
                                const long double p_obs_base = get_p_obs_base(pangenome_base, mapping_seq[j], graph_seq[j],
                                                                                          p_no_seq_error[j], 8);
                                long double log_prod;
                                if (is_consensus_fasta == false) {
                                       log_prod = log(((1-p_correctly_mapped) * background_freqs[j]) + (p_correctly_mapped * p_obs_base));
                                                                 }
                                else {
                                       log_prod = log((1-background_error_prob) * p_obs_base);
                                     }

				log_lik_if_mapped += log_prod;
		      }

		  log_likelihood_vec[i] += log_lik_if_mapped;

		  }

	  else if (!is_supported)           {
                  const double log_lik_if_unsupported = get_log_lik_if_unsupported(mppg, quality_scores);
                  //#pragma omp critical
                  log_likelihood_vec[i] += log_lik_if_unsupported;
		                            }
	}

return log_likelihood_vec;
}



