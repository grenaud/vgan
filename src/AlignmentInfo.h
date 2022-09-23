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
    string seq;
    string name;
    vg::Path path;
    int32_t mapping_quality;
    string quality_scores;
    bool is_paired;
    double identity;

    AlignmentInfo();
    AlignmentInfo(const AlignmentInfo & other);
    ~AlignmentInfo();
    AlignmentInfo & operator= (const AlignmentInfo & other);

};
#endif
