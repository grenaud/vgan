/*
 * NodeInfo
 * Date: Feb-04-2022 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef NodeInfo_h
#define NodeInfo_h


//include
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"


//VG src/
#include "utility.hpp"
#include "alignment.hpp"

using namespace std;
using namespace google::protobuf;

class NodeInfo{
private:

public:
    string  seq;
    int64   id; //TODO this is completely wasteful if it is same idx as the vector, should be removed at some point
    bool *  pathsgo;
    int cladeid; //TODO this is also wasteful, we should store it per region of the vector or have a vector of vectors
    int nbpaths; //TODO same for cladeid

    //constructors
    NodeInfo(const int64 id,const int nbpaths,const int cladeid);
    NodeInfo(const NodeInfo & other);
    //NodeInfo(); //NodeInfo ni; 

    //destructor
    ~NodeInfo();
    
    NodeInfo & operator= (const NodeInfo & other);

};
#endif
