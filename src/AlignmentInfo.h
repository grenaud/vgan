/*
 * AlignmentInfo
 * Date: Feb-07-2022
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef AlignmentInfo_h
#define AlignmentInfo_h
#include <string>
#include "alignment.hpp"
#include "AlignmentInfo.h"

using namespace std;

class AlignmentInfo{
private:

public:

    struct BaseInfo {
      char readBase;
      char referenceBase;
      bool pathSupport = true;
      double logLikelihood = log(0.00000000001);
      double logLikelihoodNoDamage = log(0.00000000001);
      // P(readBase | b_s), for each hypothetical post-mutation base b_s
      // (A,C,G,T), combining deamination then sequencing error (Figure 2
      // of the TrailMix manuscript, stages 2+3 of the Markov chain). This
      // is independent of branch time t, so it is precomputed once per
      // alignment; the HKY mutation step (stage 1, t-dependent) is applied
      // at runtime by weighting this vector and marginalizing over b_s.
      double damageSeqErrProb[4] = {0.0, 0.0, 0.0, 0.0};
      double damageSeqErrProbNoDamage[4] = {0.0, 0.0, 0.0, 0.0};
      BaseInfo() : pathSupport(false) {}  // Default constructor
    };


    string seq;
    string name;
    vg::Path path;
    int32_t mapping_quality;
    string quality_scores;
    bool is_paired;
    double identity;
    int n_reads;
    vector<string> mostProbPath;
    unordered_map <string,double> pathMap;
    unordered_map <string, bool> supportMap;
    unordered_map <string, vector<vector<BaseInfo>>> detailMap;

    AlignmentInfo();
    AlignmentInfo(const AlignmentInfo & other);
    ~AlignmentInfo();
    AlignmentInfo & operator= (const AlignmentInfo & other);

};
#endif
