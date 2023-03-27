#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

using namespace std;

class Baseshift{

private:
    int** baseshift_data_array; // store the entire 2d data array
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
    Baseshift(int lengthMax, int** baseshift_data_array);
    //~Baseshift();

    void baseshift_calc(string graph_seq, string read_seq);

    void print_counts(string prof_out_file) const;

    vector<vector<double>> print_prof(string prof_out_file, string ends = "both") const;

    //baseshift(const baseshift & other);
	//~baseshift();
	//baseshift & operator= (const baseshift & other);
};
