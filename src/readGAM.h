#pragma once
#include <iostream>
#include <fstream>
#include <filesystem>
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
#include "utility.hpp"
#include "alignment.hpp"
#include "AlignmentInfo.h"
#define VERBOSE

using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

static vector<AlignmentInfo *> * readGAM(const char * fifo_C, const bool quiet, const bool dump_json, const string &jsonfilename, const string &fastafilename){

#ifdef VERBOSE
    const string strn(fifo_C);
    const int idx = strn.find_last_of("/");
    const string to_print= strn.substr(idx + 1);
    if (!quiet) {cerr<<"Reading GAM file "<<to_print<<endl;}
#endif

    auto normal_cin = std::cin.rdbuf();
    ifstream cin_(fifo_C);
    std::cin.rdbuf(cin_.rdbuf());
    int n_reads=0;
    vector<AlignmentInfo *> * read_vec=NULL;
    read_vec = new vector<AlignmentInfo *>();
    ofstream json_out(jsonfilename);

    function<void(vg::Alignment&)> lambda = [&n_reads,&read_vec, &dump_json, &json_out](vg::Alignment& a) {
                                                if (dump_json){json_out << pb2json(a) << "\n";}
						++n_reads;
						    AlignmentInfo * ai = new AlignmentInfo();
						    ai->seq = a.sequence();
							ai->path = a.path();
							ai->mapping_quality = a.mapping_quality();
							ai->quality_scores = a.quality();
                                                        ai->is_paired = a.read_paired();
                                                        ai->identity = a.identity();
                                                    if (ai->identity != 0) // Discard unmapped reads
						        {read_vec->emplace_back(ai);}
					    };

    vg::get_input_file("-", [&](istream& in) {vg::io::for_each(in, lambda);});
    std::cin.rdbuf(normal_cin);

#ifdef VERBOSE
    if (!quiet) {
    cerr<<"Done reading GAM file "<< fifo_C <<'\n';
    if (fastafilename == "") {
        cerr<<"Found " << read_vec->size() << " reads." <<'\n';
                             }
    else {
        cerr<<"Found 1 read." <<'\n';
         }
                }
#endif

    if (std::filesystem::is_fifo(fifo_C) && !dump_json) {remove(fifo_C);}
    return read_vec;
}


// For vg dup_rm

static vector<vg::Alignment> readGAM(const char* gamfile){
    vector<vg::Alignment> algnvec;
    function<void(vg::Alignment&)> duprm_lambda = [&algnvec](vg::Alignment& a) {
        algnvec.emplace_back(a);
                                                                               };
    vg::get_input_file(gamfile, [&](istream& in) {vg::io::for_each(in, duprm_lambda);});
    return algnvec;
                                                          }

