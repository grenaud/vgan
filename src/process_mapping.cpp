#pragma once
#include "HaploCart.h"
#include "miscfunc.h"

inline const double Haplocart::getBaseFrequency(const char base) {
    switch(base) {
        case 'A': return 0.27532;
        case 'C': return 0.30044;
        case 'G': return 0.25780;
        case 'T': return 0.16644;
        default:
            // Handle error - unknown base
            throw std::invalid_argument("Invalid base");
    }
}

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

inline const vector<double> Haplocart::process_mapping(shared_ptr<Trailmix_struct> &dta, const AlignmentInfo* read_info,
vector<double> log_likelihood_vec, const vg::Mapping &mppg, string &mapping_seq,
const vector<int> &quality_scores, const string &graph_seq, unsigned int path_idx) noexcept
{

int base_index_on_read;
const bool isrev = read_info->path.mapping()[0].position().is_reverse();
if (isrev){base_index_on_read=0;}
else{base_index_on_read=read_info->seq.size()-1;}

const int pangenome_base = dta->pangenome_map.at(to_string(mppg.position().node_id()));

const double mappability = dta->mappabilities[pangenome_base];
vector<double> background_freqs;

double p_correctly_mapped;

if (dta->is_consensus_fasta == false) {
   p_correctly_mapped = (1-dta->incorrect_mapping_vec[read_info->mapping_quality]) * mappability;
                                      }

for (unsigned int i=0; i != mapping_seq.size(); ++i) {
   background_freqs.emplace_back(Haplocart::get_background_freq(mapping_seq[i]));
}


const vector<double> p_no_seq_error = Haplocart::get_p_no_seq_error_mapping(dta, mapping_seq, quality_scores, graph_seq);

          const bool is_supported = dta->nodevector.at(mppg.position().node_id()-dta->minid)->pathsgo[path_idx];
          //const bool is_supported = dta->path_supports[mppg.position().node_id()-dta->minid][i];

	  if (is_supported) {


                  double log_lik_if_supported = 0;
		  for (long unsigned int j = 0; j != graph_seq.size(); ++j)
		      {
                                // Discard ambiguous bases
                                if (graph_seq[j] == 'N' || mapping_seq[j] == 'N') {continue;}
                                if ( !isValidDNA(graph_seq[j]) || !isValidDNA(mapping_seq[j])) {continue;}

                                if (!isrev){
                                    base_index_on_read++;
                                           }
                                else{
                                    base_index_on_read--;
                                    }

                                const double p_obs_base = get_p_obs_base(pangenome_base, mapping_seq[j], graph_seq[j], \
                                               p_no_seq_error[j], 8, base_index_on_read, read_info->seq.size(), dta);

                                double log_prod;
                                if (dta->is_consensus_fasta == false) {
                                           double base_frequency = Haplocart::getBaseFrequency(mapping_seq[j]);
                                           log_prod = log(((1-p_correctly_mapped) * base_frequency) + (p_correctly_mapped * p_obs_base));
                                                                      }
                                else {
                                       log_prod = log((1-dta->background_error_prob) * p_obs_base);
                                     }

				log_lik_if_supported += log_prod;
		      }
                          log_likelihood_vec[path_idx] += 0; //log_lik_if_supported;
		  }

	  else if (!is_supported)           {

                  for (unsigned int b = 0; b!= graph_seq.size(); ++b){
                      if (graph_seq[b] == 'N' || mapping_seq[b] == 'N') {continue;}
                      if ( !isValidDNA(graph_seq[b]) || !isValidDNA(mapping_seq[b])) {continue;}
                         const double llu = get_log_lik_if_unsupported(mppg, quality_scores);
                                                                     }

                  log_likelihood_vec[path_idx] += -100.0; //get_log_lik_if_unsupported(mppg, quality_scores);
		                            }

return log_likelihood_vec;
}



