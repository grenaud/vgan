#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gzstream.h>
//#include <variant>

#include "libgab.h"

#include "NodeInfo.h"
#include "AlignmentInfo.h"
#include "Clade.h"

#include "readGAM_Euka.h"
#include "readVG.h"

//#include "miscfunc.h"
#include "subcommand/subcommand.hpp"
#include "bdsg/odgi.hpp"
#include "MCMC.h"

using namespace std;

class Euka{
private:


    // //Substitution rates due to deamination
    // vector<probSubstition> sub5p;
    // vector<probSubstition> sub3p;
    // vector<diNucleotideProb> sub5pDiNuc;
    // vector<diNucleotideProb> sub3pDiNuc;

    // //first dimension is fragment length, second is position
    // vector< vector<probSubstition> > subDeam;
    // //first dimension is fragment length, second is position
    // vector <vector<diNucleotideProb> > subDeamDiNuc;

    // probSubstition   defaultSubMatch;
    // diNucleotideProb defaultSubMatchMatrix;

    // unsigned int MINLENGTHFRAGMENT  =    25;      // minimal length for fragment
    // unsigned int MAXLENGTHFRAGMENT  =    1000;    //  maximal length for fragment
    // //const double ENTROPY_SCORE_THRESHOLD = 1.17;
    //const unsigned int MINNUMOFBINS = 6; //minimum number of bins that have to be filled. Prevents clades with multiple bins below the entropy threshold to be detected by accident. 
    //const unsigned int MINNUMOFREADS = 10; // This does not represent the number of reads, but rather the number of hits each bin needs to have to be considered.
    //const unsigned int MINIMUMMQ = 29;
    double base_freq [int ('T') +1];
    double t_T_ratio [int ('T') +1][int ('T') +1];
    bool rare_bases [int ('Y') ];

    constexpr double get_p_seq_error(const int &Q);
    
    constexpr double base_frequency(char &base);
    vector <double> get_avg(vector<double> dam, int l, int no_clades);

    // //deamination functions
    // void initDeamProbabilities(const string & deam5pfreqE,const string & deam3pfreqE);
    // void combineDeamRates(long double f1[4],long double f2[4],long double f[4],int b);
    


public:
    vector<double> get_qscore_vec();
    vector<vector<tuple<int, int, double, double > > > load_clade_chunks(string clade_chunk_path);
    const vector<vector<bool>> load_path_supports_Euka(const string & pathsupportfile);
    vector<Clade *> *  load_clade_info(const string clade_info_path, int lengthToProf); //Mikkel last argument
    const vector<long double> compute_init_vec(vector<Clade *> * clade_vec, vector<int> &clade_id_list);
    // mapping function
    void map_giraffe(string fastq1filename, string fastq2filename, const int n_threads, bool interleaved,
                            const char * fifo_A, const vg::subcommand::Subcommand* sc,
                            const string & tmpdir, const string & cwdProg, const string &prefix, const string &minprefix);
    
    tuple<vector<NodeInfo *>, int, bdsg::ODGI, vector<vector<bool>>, vector<string>> readPathHandleGraph(string &ogfilename, int n_threads, string &gbtwfilename, string &db_prefix,  vector<Clade *> *&  clade_vec);
    // get the average damage profile of all clades for 2.round euka input

    Euka();
    Euka(const Euka & other);
    ~Euka();
    Euka & operator= (const Euka & other);

    const string usage() const;
    const int run(int argc, char *argv[] , const string cwdProg);

};
