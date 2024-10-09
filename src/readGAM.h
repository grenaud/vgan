#pragma once
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
#include "Trailmix_struct.h"
#include "miscfunc.h"
//#include "TrailMix.h"
#define VERBOSE

using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

static shared_ptr<vector<AlignmentInfo*>> readGAM(shared_ptr<Trailmix_struct> &dta){

const bool dump_json = dta->dump_json;
const string fifo_to_read = dta->running_trailmix ? dta->fifo_A : dta->fifo_C;
#ifdef VERBOSE
    const string strn(fifo_to_read);
    const int idx = strn.find_last_of("/");
    const string to_print= strn.substr(idx + 1);
    if (!dta->quiet) {cerr<<"Reading GAM file "<<to_print<<endl;}
#endif

    auto normal_cin = std::cin.rdbuf();
    ifstream cin_(fifo_to_read);
    std::cin.rdbuf(cin_.rdbuf());

    unsigned int n_reads=0;
    vector<AlignmentInfo*> v;
    shared_ptr<vector<AlignmentInfo*>> read_vec = make_shared<vector<AlignmentInfo*>>(v);
    ofstream json_out(dta->jsonfilename);


    string emit_file = dta->running_trailmix ? dta->rpvg_gamfilename : "/dev/null";
    auto emitter = vg::io::VGAlignmentEmitter(emit_file, "GAM", dta->n_threads);

    vector<vg::Alignment> to_emit;
    if (dta->running_trailmix){
       cerr << "Emitting to: " << dta->rpvg_gamfilename << endl;
                             }

//ofstream fastq_out("aligned_reads.fastq");

function<void(vg::Alignment&)> lambda = [&n_reads, &read_vec, &dump_json, &json_out, &to_emit](vg::Alignment& a) {
    //static std::ofstream aligned_file("aligned.tsv");  // Open the file once and use it in the lambda
    if (dump_json) {
        json_out << pb2json(a) << "\n";
    }
    ++n_reads;
    if (n_reads % 20 == 0){
        //std::cerr << "ON READ: " << n_reads << std::endl;
                          }

    to_emit.emplace_back(a);
    AlignmentInfo* ai = new AlignmentInfo();
    ai->seq = a.sequence();
    ai->path = a.path();
    ai->name = a.name();
    ai->mapping_quality = a.mapping_quality();
    ai->quality_scores = a.quality();
    ai->is_paired = a.read_paired();
    ai->identity = a.identity();

    read_vec->emplace_back(ai);
};

    vg::get_input_file("-", [&](istream& in) {vg::io::for_each(in, lambda);});
if (dta->running_trailmix){
    emitter.emit_mapped_single(move(to_emit));
                          }

    std::cin.rdbuf(normal_cin);

    //if (std::filesystem::is_fifo(dta->fifo_A) && !dta->dump_json) {remove(dta->fifo_A);}
    if (read_vec->empty()){throw std::runtime_error("Error, no reads found in GAM file");}
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


