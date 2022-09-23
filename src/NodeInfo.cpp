#include "NodeInfo.h"


//Default constructor
// NodeInfo::NodeInfo(){
//     //cerr<<"Default constructor"<<endl;
//     this->nbpaths   = -1;
// }


//constructor we are normally using
NodeInfo::NodeInfo(const int64 id,const int nbpaths,const int cladeid){
    this->seq       = "";
    this->nbpaths   = nbpaths;
    this->pathsgo   = new bool [nbpaths]; //1 bool for each path
    for(int j=0;j<nbpaths;j++){
    	this->pathsgo[j]=false;
    }

    this->cladeid=cladeid;
    //cerr<<"END constructor ID:"<<this->id<<" "<<this->nbpaths<<" "<<this->pathsgo<<endl;
}

NodeInfo::~NodeInfo(){
    //    cerr<<"BEGIN destructor ID:"<<this->id<<" "<<this->nbpaths<<" "<<this->pathsgo<<endl;
    if(nbpaths!=-1)
	//delete(this->pathsgo);  // delete command commented out for running purpose -- Gabriel pls come back to this!
	delete [] this->pathsgo;  // delete command commented out for running purpose -- Gabriel pls come back to this!
    //cerr<<"END destructor ID:"<<this->id<<endl;
}


