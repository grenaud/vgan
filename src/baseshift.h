#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include "readGAM.h"
#include "Trailmix_struct.h"

using namespace std;

class Baseshift{

private:
    unsigned int** baseshift_data_array; // store the entire 2d data array
    int lengthToProf; // Store how many bases we want to look at from each ends.
    int dna2int ['T'+1]; // dna to int translater

    int toIndex[4][4]={
    {0,1,2,3},
    {4,5,6,7},
    {8,9,10,11},
    {12,13,14,15}
    };
    // Graph base first, read base last
    // {AA,AC,AG,AT},
    // {CA,CC,CG,CT,
    // {GA,GC,GG,GT},
    // {TA,TC,TG,TT}

public:
    Baseshift(int lengthMax, unsigned int** baseshift_data_array);
    Baseshift(const int lengthMax, const shared_ptr<Trailmix_struct> &dta);
    //~baseshift();

    void baseshift_calc(const string &graph_seq, const string &read_seq);
    void baseshift_calc(const string &graph_seq, const string &read_seq, shared_ptr<Trailmix_struct> &dta);

    void display_counts(const string &prof_out_file) const;

    vector<vector<double>> display_prof(string prof_out_file, string ends = "both") const;
    void display_prof(const string &prof_out_file, string ends, shared_ptr<Trailmix_struct> &dta) const;



    //Baseshift(const Baseshift & other);
    //~Baseshift();
    //Baseshift & operator= (const Baseshift & other);
};