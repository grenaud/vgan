#include "HaploCart.h"

const vector<long double> Haplocart::get_p_no_seq_error_mapping(const string &mapping_seq, const vector<int> &quality_scores,
               const string &graph_seq, const vector<double> &qscore_vec, const bool use_background_error_prob, const double &background_error_prob,
               const int n_threads)
{

vector<long double> ret;
long unsigned int i;
for (i=0; i<graph_seq.size(); ++i)
    {
        char graph_base = graph_seq[i];
        char read_base = mapping_seq[i];
        if (use_background_error_prob) {
            if (graph_base == read_base) {ret.emplace_back(background_error_prob);}
            else {ret.emplace_back(1-background_error_prob);}
                                       }
        else {
            if (graph_base == read_base) {ret.emplace_back(qscore_vec[quality_scores[i]]);}
            else {
                   ret.emplace_back(1-qscore_vec[quality_scores[i]]);
                 }
             }
    }

return ret;
}


inline constexpr bool Haplocart::transversion(const char base1, const char base2) {
    if ((base1=='A' || base1=='G') && (base2=='C' || base2=='T')) {return true;}
    else if ((base1=='C' || base1=='T') && (base2=='A' || base2=='G')) {return true;}
    else {return false;}
}

const bool Haplocart::inRange(const unsigned &low, const unsigned &high, const unsigned &x) {return (low <= x && x <= high);}

const long double Haplocart::get_p_obs_base(const int pangenome_base, const char mapping_base, const char graph_base,
                                       const double epsilon, const int generations)
{
  // TODO: Change this to a switch statement

  double mu;
  if(inRange(57, 372, pangenome_base)) {mu = 1.64273e-7;} // HVS I
  else if(inRange(1, 56, pangenome_base) || inRange(373, 576, pangenome_base)) {mu = 2.29640e-8;}   // HVS II
  else if(inRange(16384, 16569, pangenome_base)) { mu = 1.54555e-8;} // Control region, remaining
  else if (inRange(3307, 4262, pangenome_base) || inRange(4470, 5511, pangenome_base) || inRange(5904, 7445, pangenome_base)
        || inRange(7586, 8269, pangenome_base) || inRange(8366, 9990, pangenome_base) || inRange(10059, 10403, pangenome_base)
        || inRange(10470, 12137, pangenome_base) || inRange(12337, 14673, pangenome_base) || inRange(14747, 15886, pangenome_base)) // Protein coding
  {mu = 8.87640e-9 * (2/3) * 1.92596e-8 * (1/3);}
  else if (inRange(577, 647, pangenome_base) || inRange(1602, 1670, pangenome_base) || inRange(3230, 3304, pangenome_base)
        || inRange(4263, 4400, pangenome_base) || inRange(4402, 4469, pangenome_base) || inRange(5512, 5579, pangenome_base)
        || inRange(5587, 5654, pangenome_base) || inRange(5657, 5728, pangenome_base) || inRange(5761, 5891, pangenome_base)
        || inRange(7446, 7514, pangenome_base) || inRange(7518, 7585, pangenome_base) || inRange(8295, 8364, pangenome_base)
        || inRange(15888, 15953, pangenome_base) || inRange(15956, 16023, pangenome_base)) // tRNA
  {mu = 6.91285e-9;}
  else if (inRange(648, 1601, pangenome_base) || inRange(1671, 3229, pangenome_base)) //rRNA
  {mu = 6.91285e-9;}
  else{mu = 2.48537e-8;}

  // In the paper this is given in years, so we multiply by 30 to get mutations per generation
  mu *= 30;

  double P_no_mutation = pow(1 - mu, generations);
    double P_mutation = 1 - P_no_mutation;

    double R = 22.0 / 23.0; // Transition/transversion ratio
    double P_transition = 22.0 / 23.0;
    double P_transversion = (1.0-23.0) * 2;

    double P_no_error = 1 - epsilon;
    double P_error = epsilon / 3.0; // Equal probability for other bases

auto isTransition = [](char base1, char base2) {
    return (base1 == 'A' && base2 == 'G') || 
           (base1 == 'G' && base2 == 'A') ||
           (base1 == 'C' && base2 == 'T') || 
           (base1 == 'T' && base2 == 'C');
};


    double total_prob = 0.0;
    for (char true_base : {'A', 'C', 'G', 'T'}) {
        // Probability of true base given reference base
        double P_true_given_ref;
        if (true_base == graph_base) {
            P_true_given_ref = P_no_mutation;
        } else if (isTransition(graph_base, true_base)) {
            P_true_given_ref = P_mutation * P_transition;
        } else {
            P_true_given_ref = P_mutation * P_transversion;
        }

        // Probability of read base given true base
        double P_read_given_true = (mapping_base == true_base) ? P_no_error : P_error;

        total_prob += P_true_given_ref * P_read_given_true;
    }

  return total_prob;
}





