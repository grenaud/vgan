#pragma once

#include "insertsize_GAM.h"
#include "bdsg/odgi.hpp"
#include "subcommand/subcommand.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <gzstream.h>

class InSize
{
public:
    const vector<vector<bool>> load_path_supports_InSize(const string &pathsupportfile);

    InSize();
    InSize(const InSize &other);
    ~InSize();
    InSize &operator=(const InSize &other);

    const string usage() const;
    const int run(int argc, char *argv[], const string &cwdProg);
};
