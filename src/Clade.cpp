/*
 * Clade
 * Date: Feb-04-2022
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "Clade.h"

//using namespace google::protobuf;


Clade::Clade(int64 id, string name, double dist, int count, std::vector<double> clade_like, vector<double> clade_not_like,const vector<int> inSize, int** baseshift_clade_array, const vector<string> nameStorage){ // Mikkel code last argument
    this->id = id;
    this->name = name;
    this->dist = dist;
    this->count = count;
    this->clade_like = clade_like;
    this->clade_not_like = clade_not_like;
    this->inSize = inSize;
    this->baseshift_clade_array = baseshift_clade_array; // Mikkel code
    this->nameStorage = nameStorage;


}

Clade::~Clade(){

}


