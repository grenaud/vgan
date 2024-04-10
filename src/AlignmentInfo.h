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
      bool pathSupport;
      double logLikelihood;
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
