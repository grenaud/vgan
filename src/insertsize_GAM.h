#pragma once
#include "bdsg/odgi.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>

#include "vg/vg.pb.h"
#include "vg/io/stream.hpp"

#include "utility.hpp"
#include "alignment.hpp"
#include "AlignmentInfo.h"


using namespace std;
using namespace google::protobuf;


static vector<int> insize_histogram(vector<int> sizes, int max_len = 238);
static int insertsize(const auto &path, const string &algnseq);
static vector<int> insertsizeGAM(const string &gamfilename);

static vector<int> histogram(vector<int> sizes, int max_len = 238)
{
    vector<int> hist;
    hist.resize(max_len, 0);
    for (auto s : sizes)
    {
        hist[s]++;
    }
    return hist;
}

static vector<int> insertsizeGAM(const string &gamfilename)
{
    vector<int> sizes;
    function<void(vg::Alignment &)> lambda = [&sizes](vg::Alignment &a)
    {
        // if identity is 0 the read is unmapped.
        if (a.identity() != 0)
        {
            sizes.emplace_back(insertsize(a.path(), a.sequence()));
        }
    };

    vg::get_input_file(gamfilename, [&](istream &in)
                           { vg::io::for_each(in, lambda); });
    return sizes;
}

static int insertsize(const auto &path, const string &algnseq)
{
    int len_sequence = 0;
    auto mppgs = path.mapping();
    // Loop over mappings to nodes
    for (auto mppg : mppgs)
    {
        RepeatedPtrField<vg::Edit> ed = mppg.edit();
        // Loop over each edit in a node mapping
        for (const auto edit : ed)
        {
            if (vg::io::edit_is_match(edit) | vg::io::edit_is_sub(edit) | vg::io::edit_is_deletion(edit))
                len_sequence += edit.from_length();
           else if (vg::io::edit_is_insertion(edit))
                len_sequence += edit.to_length();
        }
    }
    return len_sequence;
}

