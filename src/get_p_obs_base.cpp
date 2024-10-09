
#include "HaploCart.h"

// Define the isTransition function
bool isTransition(char base1, char base2) {
    return (base1 == 'A' && base2 == 'G') || (base1 == 'G' && base2 == 'A') ||
           (base1 == 'C' && base2 == 'T') || (base1 == 'T' && base2 == 'C');
}

const vector<double> Haplocart::get_p_no_seq_error_mapping(shared_ptr<Trailmix_struct> &dta, string &mapping_seq, const vector<int> &quality_scores,
               const string &graph_seq)
{

vector<double> ret;

long unsigned int i;
//#pragma omp parallel for private(i) num_threads(1)
for (size_t i=0; i<graph_seq.size(); ++i)
    {
        char graph_base = graph_seq[i];
        char read_base = mapping_seq[i];
        if (dta->use_background_error_prob) {
            if (graph_base == read_base) {ret.emplace_back(dta->background_error_prob);}
            else {ret.emplace_back(1-dta->background_error_prob);}
                                       }
        else {
            if (graph_base == read_base) {ret.emplace_back(dta->qscore_vec[quality_scores[i]]);}
            else {
                   ret.emplace_back(1-dta->qscore_vec[quality_scores[i]]);
                 }
             }
    }

return ret;
}

const bool Haplocart::inRange(const unsigned &low, const unsigned &high, const unsigned &x) {return (low <= x && x <= high);}

const double Haplocart::get_p_obs_base(const size_t pangenome_base, const char mapping_base, const char graph_base,
                                       const double epsilon, const size_t generations, \
                                       const size_t base_index_on_read, const size_t Lseq, shared_ptr<Trailmix_struct> &dta)
{

double mu;
switch(pangenome_base){
    case 57 ... 372:
        mu = 1.64273e-7; // HVS I
        break;
    case 1 ... 56:
    case 373 ... 576:
        mu = 2.29640e-8; // HVS II
        break;
    case 16384 ... 16569:
        mu = 1.54555e-8; // Control region, remaining
        break;
    case 3307 ... 4262:
    case 4470 ... 5511:
    case 5904 ... 7445:
    case 7586 ... 8269:
    case 8366 ... 9990:
    case 10059 ... 10403:
    case 10470 ... 12137:
    case 12337 ... 14673:
    case 14747 ... 15886:
        mu = 8.87640e-9 * (2/3) * 1.92596e-8 * (1/3); // Protein coding
        break;
    case 577 ... 647:
    case 1602 ... 1670:
    case 3230 ... 3304:
    case 4263 ... 4400:
    case 4402 ... 4469:
    case 5512 ... 5579:
    case 5587 ... 5654:
    case 5657 ... 5728:
    case 5761 ... 5891:
    case 7446 ... 7514:
    case 7518 ... 7585:
    case 8295 ... 8364:
    case 15888 ... 15953:
    case 15956 ... 16023:
        mu = 6.91285e-9; // t-RNA
        break;
    case 648 ... 1601:
    case 1671 ... 3229:
        mu = 6.91285e-9; // r-RNA
        break;
    default:
        mu = 2.48537e-8;
                   }


  // In the paper this is given in years, so we multiply by 30 to get mutations per generation
  mu *= 30;

double P_no_mutation = pow(1 - mu, generations);
    double P_mutation = 1 - P_no_mutation;

    double R = 22.0 / 23.0; // Transition/transversion ratio
    double P_transition = 22.0 / 23.0;
    double P_transversion = (1.0-23.0) * 2;

    double P_no_error = 1 - epsilon;
    double P_error = epsilon / 3.0; // Equal probability for other bases

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





