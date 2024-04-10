#include <sys/wait.h>
#include "InSize.h"
#include "miscfunc.h"
#include "crash.hpp"
#include "config/allocator_config.hpp"
#include "preflight.hpp"
#include "io/register_libvg_io.hpp"

using namespace vg;
using namespace std;

InSize::InSize()
{
}

InSize::~InSize()
{
}

const string InSize::usage() const
{
    return string(
	string("")
        + "vgan insize <database prefix> <reads.gam> <options>\n"
        + "\n" + "Reads a GAM file and produces the insert sizes. Ignores unmapped reads."
        + "\n\n" + "Mandatory arguments:\n\n"
        + "<database prefix> is the prefix of the database to which the sequences were aligned\n"
        + "For instance, for:\n" + "\teuka_tetraanthro.vg  euka_tetraanthro.clade euka_tetraanthro.bins\n"
        + "in: /disk/software/euka/share/DB/\n"
        + "\t<database prefix> would be /disk/software/euka/share/DB/euka_tetraanthro\n"
        + "\n" + "<reads.gam> are the sorted gam alignments"
	+ "\n\nOptions:"
	+ " --histogram: return the number of reads of each size instead of returning sizes."""
    );
}

const int InSize::run(int argc, char *argv[], const string& cwdProg)
{
    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    temp_file::set_system_dir();

    int lastOpt = 1;
    string gamfilename = argv[1];
    string tmpdir = "/tmp/";

    bool make_hist = false;
    for(int i = 1; i < argc; i++)
    {
	if (string(argv[i]) == "--histogram"){
	    make_hist = true;
	}
    }

    string first_fifo = tmpdir + random_string(7);
    const char *fifo_A = first_fifo.c_str();
    mkfifo(fifo_A, 0666);
    pid_t wpid;
    int status = 0;
    pid_t pid1 = fork();
    if (pid1 == -1)
    {
        throw std::runtime_error("Error in fork");
    }

    const vg::subcommand::Subcommand *sc = NULL;
    if (pid1 == 0)
    {
        if (gamfilename != "")
        {
            ifstream src(gamfilename);
            ofstream dst(fifo_A);
            dst << src.rdbuf();
            exit(0);
        }
    }

    vector<int> sizes = insertsizeGAM(fifo_A);
    if (make_hist){
        vector<int> hist = histogram(sizes);
        for (const auto &h: hist){
            cerr << h << endl;
        }
    }

    else{
        for (const auto &s: sizes){
            cerr << s << endl;
        }
    }

    remove(fifo_A);
    return 0;
}
