#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gzstream.h>
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
#include "libgab.h"
#include "Clade.h"
#include "Euka.h"
#include "Trailmix_struct.h"

class Gam2prof{
private:

    vector<Clade *> *  load_clade_info(const string &clade_info_path, const int lengthToProf);
    vector<vector<tuple<int, int, double, double > > > load_clade_chunks(const string &clade_chunk_path);

public:

    Gam2prof();
    Gam2prof(const Gam2prof & other);
    ~Gam2prof();
    Gam2prof & operator= (const Gam2prof & other);

    const string usage() const;
    const int run(int argc, char *argv[] , const string &cwdProg);

};