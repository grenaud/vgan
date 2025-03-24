#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gzstream.h>
#include "libgab.h"
#include "miscfunc.h"

using namespace std;


class Damage{
private:

public:

    Damage();
    Damage(const Damage & other);
    ~Damage();
    Damage & operator= (const Damage & other);

    //deamination functions
    void initDeamProbabilities(const string & deam5pfreqE,const string & deam3pfreqE);
    void combineDeamRates(long double f1[4],long double f2[4],long double f[4],int b);

    //Substitution rates due to deamination
    vector<probSubstition> sub5p;
    vector<probSubstition> sub3p;
    vector<diNucleotideProb> sub5pDiNuc;
    vector<diNucleotideProb> sub3pDiNuc;

    //first dimension is fragment length, second is position
    vector< vector<probSubstition> > subDeam;
    //first dimension is fragment length, second is position
    vector <vector<diNucleotideProb> > subDeamDiNuc;

    probSubstition   defaultSubMatch;
    diNucleotideProb defaultSubMatchMatrix;

    unsigned int MINLENGTHFRAGMENT  =    15;      // minimal length for fragment
    unsigned int MAXLENGTHFRAGMENT  =    1000;    //  maximal length for fragment
    double base_freq [int ('T') +1];
    double t_T_ratio [int ('T') +1][int ('T') +1];
    bool rare_bases [int ('Y') ];
    float probBasePreDamage[4];
    double probBasePostDamage [4];

};
