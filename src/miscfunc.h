#pragma once
#include "libgab.h"
#include <random>
#include <limits>
#include <unordered_set>

#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << setprecision(10) << v[i] << '\t';}cerr << endl << endl;
#define MIN2(a,b) (((a)<(b))?(a):(b))
#define MAX2(a,b) (((a)>(b))?(a):(b))
#define probably_true(x) __builtin_expect(!!(x), 1)

// Compute the mean of a sequence
inline double mean(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Compute the variance of a sequence
inline double variance(const std::vector<double>& v, double mean) {
    double sum = 0.0;
    for (const auto& i : v) {
        double diff = i - mean;
         sum += diff * diff;
    }
    return sum / (v.size() - 1);
}

inline double autocorrelation(const std::vector<double>& v, int k, double m, double var) {
    double numer = 0.0;
    const size_t n = v.size();
    for (size_t i = 0; i < n - k; ++i) {
        numer += ((v[i] - m) * (v[i + k] - m));
    }
    return numer / ((n - k) * var);
}

inline double autocorrelation(const std::vector<double>& v, int k) {
    const size_t n = v.size();
    const double m = mean(v);
    const double var = variance(v, m);

    double numer = 0.0;
    for (size_t i = 0; i < n - k; ++i) {
        numer += ((v[i] - m) * (v[i + k] - m));
    }
    return numer / ((n - k) * var);
}

inline double effectiveSampleSize(const std::vector<double>& v) {
    const size_t n = v.size();
    const int max_lag = std::min(static_cast<size_t>(100), n / 2);  // Use fixed max lag
    const double m = mean(v);
    const double var = variance(v, m);
    
    double rho_prev = 1.0;
    double rho_current = autocorrelation(v, 1, m, var);
    double rho_hat_tot = 2.0 * rho_current;  // Initialize with 2 * rho[1]
    
    for (int t = 2; t <= max_lag; ++t) {
        double rho_next = autocorrelation(v, t, m, var);
        
        if (rho_next + rho_current <= 0) break;  // Early termination
        
        rho_hat_tot += 2.0 * rho_next;
        
        // Early termination if autocorrelation becomes very small
        if (std::abs(rho_next) < 0.01) break;
        
        rho_prev = rho_current;
        rho_current = rho_next;
    }
    
    return n / (1 + rho_hat_tot);
}


typedef struct {
    double s[12];
 } substitutionRates;


//  model->obs
//  0  A->A
//  1  A->C
//  2  A->G
//  3  A->T
//  4  C->A
//  5  C->C
//  6  C->G
//  7  C->T
//  8  G->A
//  9  G->C
//  10 G->G
//  11 G->T
//  12 T->A
//  13 T->C
//  14 T->G
//  15 T->T

typedef struct {
    double s[16];
} probSubstition;


typedef struct {
    double p[4][4];
} diNucleotideProb;


static void readNucSubstitionRatesFreq(const string filename,vector<substitutionRates> & subVec){
    igzstream subFP;

    subFP.open(filename.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (subFP.good()){
	vector<string> fields;
	string line;

	//header
	if ( !getline (subFP,line)){
	    throw std::runtime_error("Unable to open file "+filename);
	}
	fields = allTokens(line,'\t');
        if (fields.size() == 13){fields.pop_back();}

	if(fields.size() != 12){
	    throw std::runtime_error("line from error profile has " + to_string(fields.size()) + " fields rather than 12");
	}


	//probs
	while ( getline (subFP,line)){

	    fields = allTokens(line,'\t');
            if (fields.size() == 13){fields.pop_back();}

	    if(fields.size() != 12){
		throw std::runtime_error("line from error profile has " + to_string(fields.size()) + " fields rather than 12");
	    }

	    substitutionRates tempFreq;


	    for(unsigned int k=0;k<12;k++){
		//for(unsigned int t=0;t<=2;t++){
		tempFreq.s[k]=destringify<double>(fields[k]);
		//}
	    }



	    subVec.emplace_back( tempFreq );
	}
	subFP.close();
    }else{
	throw std::runtime_error("Unable to open file "+filename);
    }



}


void write_path_supports(const string & pathSupportFile, auto nodevector){

cerr << "WRITING PATH SUPPORTS" << endl;
cerr << "Number of nodes: " << nodevector->size() << endl;
ofstream fout(pathSupportFile,ios::out | ios::binary);

for (int i=0; i<nodevector->size(); ++i)
{
  if (i % 100 == 0) {cerr << i << endl;}
    for (int j = 0; j < nodevector->at(i)->nbpaths; ++j) {
        fout << nodevector->at(i)->pathsgo[j]?"1":"0";
                                                         }
fout << "\n";
}
                                                                         };


static const inline string random_string(std::string::size_type length)
//Generate random string of a specified length
// https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
{
    static auto& chrs = "0123456789"
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    thread_local static std::mt19937 rg{std::random_device{}()};
    thread_local static std::uniform_int_distribution<std::string::size_type> pick(0, sizeof(chrs) - 2);

    std::string s;

    s.reserve(length);

    while(length--)
        s += chrs[pick(rg)];

    return s;
}

/**
 * @brief Given a quality score, return the probability of sequencing error
 */
static inline const double get_p_seq_error(const int &Q)
{
if (Q > 2)
{
    return pow(10, ((-1 * Q)*0.1));
}
else{return 0.25;}
}


/**
 * get_qscore_vec gives a vector of
 *
 * This function creates a vector of probabilities for each phred score. It is used to look up the probability of a sequenceing error for each base
 *
 * @return vector of probabilities for phred scores 0-100
 *
 */

static inline const vector<double> get_qscore_vec() {

    vector<double> qscore_vec;
    for (int Q=0; Q<100;++Q) {
        if (Q >= 2) {
           qscore_vec.emplace_back(pow(10, ((-1 * Q)*0.1)));
                              }
        else{
           qscore_vec.emplace_back(0.25);
            }

    }
    return qscore_vec;
}


inline constexpr double get_p_incorrectly_mapped(const int Q)
{return pow(10, ((-1 * Q)*0.1));}

int main_filter(int argc, char** argv);
int main_gamsort(int argc, char** argv);
int main_giraffe(int argc, char** argv);

// void update_clade_dis(vector <double> clade_like, vector<double> clade_not_like, double frac, int vec_len){

//     std::transform(clade_like.begin(), clade_like.end(), clade_like.begin(), [&frac](auto & c){return c*frac;});

//     std::transform(clade_not_like.begin(), clade_not_like.end(), clade_not_like.begin(),[&vec_len](auto & c){return c/vec_len;});

// }

constexpr double logit(const double like_rato){

    const double p = 1 / (1 + (exp(-(like_rato -1))/1));

    return p;

}

//Taken from: https://stackoverflow.com/questions/11964552/finding-quartiles
template <typename T1, typename T2> typename T1::value_type quant(const T1 &x, T2 q)
{
    assert(q >= 0.0 && q <= 1.0);

    const auto n  = x.size();
    const auto id = (n - 1) * q;
    const auto lo = floor(id);
    const auto hi = ceil(id);
    const auto qs = x[lo];
    const auto h  = (id - lo);

    return (1.0 - h) * qs + h * x[hi];
}

