#ifndef readVG_h
#define readVG_h
#include <chrono>

//include
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
#include "bdsg/odgi.hpp"
#include "handlegraph/path_handle_graph.hpp"

//VG src/
#include "utility.hpp"
#include "alignment.hpp"


#include "NodeInfo.h"

  //#define VERBOSE
  #define ADDDATA
  //#define DEBUGREADGRAPH

using namespace std;
using namespace google::protobuf;


//int main (int argc, char *argv[]) {
static auto readVG(const string & vgfilename){
    //    string vgfilename = string(argv[1]);
#ifdef VERBOSE
    cerr<<"reading VG file "<<vgfilename<<endl;
#endif

    vector<NodeInfo *> * nodevector=NULL;
    bool foundNot1k=false;
    int cladeid=0;
    int nbpaths=0;
    function<void(vg::Graph&)> lambda = [&nodevector,&cladeid,&foundNot1k,&nbpaths](vg::Graph& g) {
#ifdef DEBUGREADGRAPH
					    cout<<"function "<<endl;
#endif
					    //checking the # of paths
					    //if(nbpaths == -1 ){
                                            auto g_path = g.path();
					    nbpaths=g_path.size();
                                            //cout << "nbpaths: " << nbpaths << endl;
						// }else{
					    // 	if(nbpaths != g.path().size()){
					    // 	     cerr<<"The number of path is not identical between chunks of the graph "<<nbpaths<<" vs "<<g.path().size()<<endl;
					    // 	     exit(1);
					    // 	}
					    // }

					    RepeatedPtrField<vg::Node> nd= g.node();

					    int64 maxid=0;
					    for(int i=0;i<nd.size();++i){
						if(nd[i].id() > maxid) maxid = nd[i].id();
					    }
//#ifdef DEBUGREADGRAPH
					    cout<<"nbpaths "<< g.path().size()<<" node size "<<nd.size()<<" max ID "<<maxid<<" "<<foundNot1k<<endl;
//#endif


					    //we need to initialize the vector and add NodeInfo
					    if(nodevector == NULL){
						nodevector =  new vector<NodeInfo *> ();
#ifdef DEBUGREADGRAPH
						cout<<"ADDINGnodes "<<maxid<<endl;
#endif
						for(int64 i=0;i<=maxid;i++){
						    //cout<<"i "<<i<<endl;
						    NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
						    nodevector->emplace_back(nodetoadd);
						}

					    }else{
						//new clade
						if(foundNot1k){ //need to add some nodes
                                                    cladeid++;
						    unsigned int orgSize = nodevector->size();
						    for(int64 i=orgSize;i<=maxid;++i){
							//cout<<"i "<<i<<endl;
#ifdef DEBUGREADGRAPH
							cout<<"adding nodes1 i "<<i<<endl;
#endif
							NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
							nodevector->emplace_back(nodetoadd);
						    }
						    foundNot1k=false;
						}else{
						    //we are still in the same connected component but we might have found a new max ID
						    //add the new nodes if we have found a new maxID
						    unsigned int orgSize = nodevector->size();
						    for(int64 i=orgSize;i<=maxid;++i){
#ifdef DEBUGREADGRAPH
							cout<<"adding nodes2 i "<<i<<endl;
#endif
							NodeInfo * nodetoadd = new NodeInfo(i,nbpaths,cladeid);
							nodevector->emplace_back(nodetoadd);
						    }
						}

					    }

					    if(nd.size() != 1000){//last chunk of current graph
						//foundNot1k=true;
					    }

					    //populating the sequence field
					    for(int i=0;i<nd.size();i++){
#ifdef DEBUGREADGRAPH
						cout<<"nd["<<i<<"] "<<nd[i].name()<<" seq="<<nd[i].sequence()<<" ID="<<nd[i].id()<<endl;
#endif
						nodevector->at(nd[i].id())->seq  = nd[i].sequence();
					    }

					    //cout << pb2json(g) << endl;
#ifdef DEBUGREADGRAPH
					    cout<<"------paths--------"<<endl;
#endif
					    RepeatedPtrField<vg::Path> pth= g.path();
#ifdef DEBUGREADGRAPH
					    cout<<"pth "<<pth.size()<<endl;
#endif

					    for(int i=0;i<pth.size();++i){ //for every path
#ifdef DEBUGREADGRAPH
						cout<<"pth["<<i<<"] "<<pth[i].name()<<endl;
#endif
						RepeatedPtrField<vg::Mapping> mpg= pth[i].mapping();
#ifdef DEBUGREADGRAPH
						cout<<"mpg "<<mpg.size()<<endl;
						cout<<"mpg[j] rank node_id offset is_reverse name"<<endl;
#endif
						for(int j=0;j<mpg.size();j++){
						    vg::Position ps = mpg[j].position();
#ifdef DEBUGREADGRAPH
						    cout<<"path["<<i<<"] mpg["<<j<<"] "<<mpg[j].rank()<<" "<<ps.node_id()<<" "<<ps.offset()<<" "<<ps.is_reverse()<<" "<<ps.name()<<endl;
#endif
						    nodevector->at( ps.node_id() )->pathsgo[ i ]  = true;
						}
					    }

#ifdef DEBUGREADGRAPH
					    cout<<"------nodes--------"<<endl;
#endif

					}; //end lambda


    //cerr<<"reading2 "<<vgfilename<<endl;

    vg::get_input_file(vgfilename, [&](istream& in) {
				       vg::io::for_each(in, lambda);
				   });
#ifdef DEBUGREADGRAPH
    cout<<"data struct"<<endl;
#endif

    //TODO maybe have a vector of vectors
    // for(unsigned int i=0;i<nodevector->size();i++){
    // 	//building a small string
    // 	string s="";
    // 	for(int j=0;j<nodevector->at(i)->nbpaths;j++){
    // 	    s=s+(nodevector->at(i)->pathsgo[ j ]?"1":"0");
    // 	}
    // 	//Idx NodeID  SEQ paths?
    // 	cout<<"idx "<<i<<" "<<nodevector->at(i)->cladeid<<" "<<nodevector->at(i)->id<<" "<<nodevector->at(i)->seq<<" "<<s<<endl;
    // }

#ifdef VERBOSE
    cerr<<"done reading VG file "<<vgfilename<<endl;
#endif

cerr << "NUMBER OF PATHS: " << nbpaths << endl;
cerr << "WRITING PATH SUPPORTS" << endl;
ofstream fout("path_supports",ios::out | ios::binary);

for (int i=0; i<nodevector->size(); ++i)
{
  if (i % 100 == 0) {cout << i << endl;}
    for (int j = 0; j < nbpaths; ++j) {
        fout << nodevector->at(i)->pathsgo[j]?"1":"0";
                                      }
fout << "\n";
}

cout << "DONE" << endl;
   return nodevector;
}

#endif
