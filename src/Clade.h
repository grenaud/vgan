/*
 * Clade
 * Date: Feb-04-2022 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Clade_h
#define Clade_h

#include <string>
#include <iostream>
#include <fstream>

#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
#include "utility.hpp"
#include "alignment.hpp"


using namespace std;
using namespace google::protobuf;

class Clade{
private:
    
public:

    int64   id;
    string name;
    double dist;
    int count;
    vector<double> clade_like;
    vector<double> clade_not_like; 
    vector<int> inSize; 
    int** baseshift_clade_array; //Mikkel code
    vector<string>nameStorage; 
    Clade(const int64 id,const string name,const double dist, const int count, const vector<double> clade_like, const vector<double> clade_not_like, const vector<int> inSize, int** baseshift_clade_array, const vector<string> nameStorage); //Mikkel code last argument
    Clade(const Clade & other);
    ~Clade();
    Clade & operator= (const Clade & other);

};
#endif
