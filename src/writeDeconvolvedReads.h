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

using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;

void Trailmix::writeDeconvolvedReads(shared_ptr<Trailmix_struct> &dta, const string &input_file, const string &output_file, vector<int> reads_to_include){

    auto normal_cin = std::cin.rdbuf();
    ifstream cin_(input_file);
    std::cin.rdbuf(cin_.rdbuf());

    unsigned int read_counter=0;
    vector<AlignmentInfo*> v;
    shared_ptr<vector<AlignmentInfo*>> read_vec = make_shared<vector<AlignmentInfo*>>(v);

        auto emitter = vg::io::VGAlignmentEmitter(output_file, "GAM", dta->n_threads);
        cerr << "Emitting to: " << output_file << endl;
        vector<vg::Alignment> to_emit;

        function<void(vg::Alignment&)> lambda = [&read_counter, &to_emit, &reads_to_include](vg::Alignment& a) {
                                                     if (std::find(reads_to_include.begin(), reads_to_include.end(), read_counter) != reads_to_include.end()){
                                                         to_emit.emplace_back(a);
                                                                                                                                                             }
                                                     ++read_counter;
		    	            		    };

        vg::get_input_file("-", [&](istream& in) {vg::io::for_each(in, lambda);});
        emitter.emit_mapped_single(move(to_emit));

        std::cin.rdbuf(normal_cin);

    cerr << "RETURNING FROM WRITING DECONVOLVED READS" << endl;
    return;
}


