#pragma once
#include "baseshift.h"
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <cerrno>

#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << v[i] << '\t';}cerr << endl << endl;

std::vector<std::string> buildRpvgArgumentsHaplotypeTranscripts(const std::shared_ptr<Trailmix_struct>& dta) {
    std::vector<std::string> rpvgarguments = {
        "rpvg", "-g", dta->graph_dir + dta->graph_prefix + ".xg", "-p", dta->graph_dir + dta->graph_prefix + ".gbwt",
        "-a", dta->rpvg_gamfilename, "-o", dta->tmpdir + "rpvg_ht", "-i", "haplotype-transcripts",
        "-f", dta->graph_dir + "/pantranscriptome.txt", "-u", "-s", "-t",
        std::to_string(dta->n_threads), "--path-node-cluster", "--use-hap-gibbs",
    };

    if (dta->rng_seed != "NONE") {
        rpvgarguments.emplace_back("-r");
        rpvgarguments.emplace_back(dta->rng_seed);
    }

    rpvgarguments.insert(rpvgarguments.end(), {
        "-y", std::to_string(dta->k), "-m", "60", "-d", "20", "--score-not-qual",
        "--min-noise-prob", "1e-10", "--vgan-temp-dir", dta->tmpdir, "-b",
            "--prob-precision", "1e-8", "--filt-best-score", "1e-4"
    });

    if (dta->strand_specific) {
        rpvgarguments.insert(rpvgarguments.end(), {"-e", "fr"});
    }

    return rpvgarguments;
}


std::vector<std::string> buildRpvgArguments(const std::shared_ptr<Trailmix_struct>& dta) {
    std::vector<std::string> rpvgarguments = {
        "rpvg", "-g", dta->graph_dir + dta->graph_prefix + ".xg", "-p", dta->graph_dir + dta->graph_prefix + ".gbwt",
        "-a", dta->rpvg_gamfilename, "-o", dta->tmpdir + "rpvg_hap", "-i", "haplotypes",
        "-u", "-s", "-t", std::to_string(dta->n_threads), "--filt-best-score", "1e-4",
        "--min-noise-prob", "1e-10", "-m", "60", "-d", "20", "--path-node-cluster", "--use-hap-gibbs",
    };

    if (dta->rng_seed != "NONE") {
        rpvgarguments.emplace_back("-r");
        rpvgarguments.emplace_back(dta->rng_seed);
    }

    rpvgarguments.insert(rpvgarguments.end(), {
        "-y", std::to_string(dta->k), "--score-not-qual",
        "--vgan-temp-dir", dta->tmpdir,
        "--prob-precision", "1e-10",
    });

    if (dta->strand_specific) {
        rpvgarguments.insert(rpvgarguments.end(), {"-e", "fr"});
    }

    return rpvgarguments;
}

// Function to find the paths that go through a given node
vector<string> paths_through_node(const bdsg::HashGraph& graph, const bdsg::handle_t& node) {
    vector<string> paths;
    graph.for_each_step_on_handle(node, [&](const bdsg::step_handle_t& step) {
        // Get the path associated with this step
        bdsg::path_handle_t path_handle = graph.get_path_handle_of_step(step);
        // Get the path name and add it to the result vector
        paths.emplace_back(graph.get_path_name(path_handle));
    });
    return paths;
}

std::vector<std::vector<double>> convertMapsToVector(std::vector<AlignmentInfo*>* & gam) {
    std::vector<std::vector<double>> result;

    if (gam->empty()) {
        throw std::runtime_error("Error: GAM file is empty.");
    }

    for (const auto& map : *gam) {
        std::vector<double> row;
        for (const auto& pair : map->pathMap) {
            row.emplace_back(pair.second);
        }

        // Compute the sum of the elements in the row
        double sum = 0.0;
        for (int i = 0; i < row.size(); i++) {
            sum += row[i];
        }

        // Push the row into the result vector
        result.emplace_back(row);
    }

    return result;
}

void write_fq_read(auto & dummyFASTQFile, int offset, const int window_size, const string &fastaseq, char dummyqualscore) {
        string seq_to_write, qual_to_write;
        for (const char base : fastaseq.substr(min(offset, int(fastaseq.size())), window_size)) {
            if (base != 'N') {
                seq_to_write += base;
                qual_to_write += dummyqualscore;
                             }
            else {
                  seq_to_write += 'A';
                  qual_to_write += '!';
                 }
                                                               }

        dummyFASTQFile << '@' << random_string(7) << '\n';
        dummyFASTQFile << seq_to_write;
        dummyFASTQFile << "\n+\n";
        dummyFASTQFile << qual_to_write;
        dummyFASTQFile << '\n';
        offset += window_size;
                                                                                                                                   }


const char get_dummy_qual_score(const double &background_error_prob) {
    // Given a background error probability, return a dummy quality score for the artificial FASTQ reads
    string illumina_encodings = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI";
    const int Q = -10 * log10(background_error_prob);
    return illumina_encodings[Q];
                                                                           }


const string fa2fq(const string & fastaseq, const char & dummyqualscore, const string & tmpdir) {

    // Given a consensus sequence and a dummy quality score, convert to FASTQ

    // https://stackoverflow.com/questions/7560114/random-number-c-in-some-range

    const string prefix = tmpdir;
    const string tempfqfilename= prefix+random_string(7);
    ofstream dummyFASTQFile;
    dummyFASTQFile.open(tempfqfilename);
    unsigned int window_size = ceil(fastaseq.size()/100);
    unsigned int offset = 0;
    string seq_to_write, qual_to_write;

     for (int i = 0; i < 101; ++i) {
        write_fq_read(dummyFASTQFile, offset, window_size, fastaseq, dummyqualscore);
        offset += 100;
                                   }

     for (int i = 1; i < 101; ++i) {
        write_fq_read(dummyFASTQFile, offset, window_size, fastaseq, dummyqualscore);
        offset += 100;
                                   }

     dummyFASTQFile.close();
     return tempfqfilename;
                                                                                     }

// SAFARI (grenaud/SAFARI, "Sensitive Alignments from a RYmer Index") is a damage-aware fork
// of `vg giraffe`, forked from vg v1.44 ("safari" subcommand). Its core mapping algorithm
// changes are deep (gbwtgraph minimizer-key encoding, snarl seed clustering, minimizer_mapper
// internals) and were never ported forward to vg v1.75 -- doing so in-process would mean
// reconciling those changes against a substantially restructured (zipcode-based) giraffe.
// Instead, TrailMix invokes a standalone, prebuilt safari-capable vg v1.44 binary
// (dep/safari_vg/bin/vg) as a SUBPROCESS for just the mapping step. It writes a GAM file,
// which the rest of the (vg v1.75-based) pipeline reads back exactly as it already does for
// the in-process giraffe call fifo -- GAM/protobuf format is stable across versions (verified:
// vg.proto differs by exactly one new, purely-additive field between v1.44 and v1.75).
static void run_safari_subprocess(vector<string> &arguments, const string &output_fifo, shared_ptr<Trailmix_struct> &dta){
    const string safari_bin = getFullPath(dta->cwdProg + "../dep/safari_vg/bin/vg");

    pid_t pid = fork();
    if (pid == -1) {
        throw runtime_error("[TrailMix] fork() failed for safari subprocess");
    }

    if (pid == 0) {
        // Child: redirect stdout to the GAM fifo, then exec the standalone safari binary.
        int fd = open(output_fifo.c_str(), O_WRONLY);
        if (fd == -1) {
            string msg = "[TrailMix] Child: open() failed on fifo " + output_fifo + ": " + string(strerror(errno)) + "\n";
            write(STDERR_FILENO, msg.c_str(), msg.size());
            _exit(126);
        }
        dup2(fd, STDOUT_FILENO);
        close(fd);

        vector<char*> argv;
        argv.reserve(arguments.size() + 1);
        for (auto &a : arguments) { argv.emplace_back(const_cast<char*>(a.c_str())); }
        argv.emplace_back(nullptr);

        execv(safari_bin.c_str(), argv.data());
        // execv only returns on failure.
        string msg = "[TrailMix] Child: execv() failed on " + safari_bin + ": " + string(strerror(errno)) + "\n";
        write(STDERR_FILENO, msg.c_str(), msg.size());
        _exit(127);
    }

    // Parent: caller is responsible for reading the fifo; just wait for the child here.
    int status = 0;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
        std::cerr << "[TrailMix] safari subprocess exited with status " << WEXITSTATUS(status) << std::endl;
    }
}

void Trailmix::map_giraffe(const string &fastaseq, shared_ptr<Trailmix_struct> &dta){

    if (!dta->quiet) {cerr << "Mapping reads..." << endl;}

    vector<string> arguments;
    arguments.emplace_back("vg");
    arguments.emplace_back("safari");

    // Map with safari to generate a GAM file
    string minimizer_to_use=  dta->graph_dir + dta->graph_prefix + ".min";
    string rymer_to_use=  dta->graph_dir + dta->graph_prefix + ".ry";


if (dta->fastq1filename != "" && dta->fastq2filename != "")
   {
        arguments.emplace_back("-f");
        arguments.emplace_back(dta->fastq1filename);
        arguments.emplace_back("-f");
        arguments.emplace_back(dta->fastq2filename);
        arguments.emplace_back("-Z");
        arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".giraffe.gbz"));
        arguments.emplace_back("-d");
        arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".dist"));
        arguments.emplace_back("-m");
        arguments.emplace_back(minimizer_to_use);
        arguments.emplace_back("-q");
        arguments.emplace_back(rymer_to_use);
        arguments.emplace_back("--deam-3p");
        arguments.emplace_back(dta->deam3pfreqE);
        arguments.emplace_back("--deam-5p");
        arguments.emplace_back(dta->deam5pfreqE);
        arguments.emplace_back("-j");
        arguments.emplace_back("0.5");

        run_safari_subprocess(arguments, dta->fifo_A, dta);

    }

else if (dta->fastq1filename != "" && dta->fastq2filename == "")
    {
               arguments.emplace_back("-f");
               arguments.emplace_back(dta->fastq1filename);
               arguments.emplace_back("-Z");
               arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".giraffe.gbz"));
               arguments.emplace_back("-d");
               arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".dist"));
               arguments.emplace_back("-m");
               arguments.emplace_back(minimizer_to_use);
               arguments.emplace_back("-q");
               arguments.emplace_back(rymer_to_use);
               arguments.emplace_back("--deam-3p");
               arguments.emplace_back(dta->deam3pfreqE);
               arguments.emplace_back("--deam-5p");
               arguments.emplace_back(dta->deam5pfreqE);
               arguments.emplace_back("-j");
               arguments.emplace_back("0.5");

           if (dta->interleaved) {
               arguments.emplace_back("-i");
                            }

               run_safari_subprocess(arguments, dta->fifo_A, dta);

    }

else if (fastaseq != ""){
    string fastapath;
    const char dummyq = get_dummy_qual_score(dta->background_error_prob);

    string fasta_cmd;
    const string dummy_fastq_file = fa2fq(fastaseq, dummyq, dta->tmpdir);


    arguments.emplace_back("-f");
    arguments.emplace_back(dummy_fastq_file);
    arguments.emplace_back("-Z");
    arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".giraffe.gbz"));
    arguments.emplace_back("-d");
    arguments.emplace_back(getFullPath(dta->graph_dir + dta->graph_prefix + ".dist"));
    arguments.emplace_back("-m");
    arguments.emplace_back(minimizer_to_use);
    arguments.emplace_back("-q");
    arguments.emplace_back(rymer_to_use);
    arguments.emplace_back("--deam-3p");
    arguments.emplace_back(dta->deam3pfreqE);
    arguments.emplace_back("--deam-5p");
    arguments.emplace_back(dta->deam5pfreqE);
    arguments.emplace_back("-j");
    arguments.emplace_back("0.5");

    run_safari_subprocess(arguments, dta->fifo_A, dta);

    remove(dummy_fastq_file.c_str());
   }

    if (!dta->quiet) {std::cerr << "Reads mapped" << endl;}
}

char** convert_to_char_array(const vector<string>& arguments) {
    char** argvtopass = new char*[arguments.size()];
    for (size_t i = 0; i < arguments.size(); i++) {
        argvtopass[i] = const_cast<char*>(arguments[i].c_str());
    }
    return argvtopass;
}

void Trailmix::run_gam2prof(shared_ptr<Trailmix_struct>& dta) {
    //vector<string> g2parguments = {"vgan", "gam2prof", "--running-trailmix"};
    //char** g2pargvtopass = convert_to_char_array(g2parguments);

    //Gam2prof().run(g2parguments.size(), g2pargvtopass, dta->cwdProg, dta);

    //delete[] g2pargvtopass;
}

// Function to find the paths that go through a given node
vector<string> Trailmix::paths_through_node(const bdsg::HashGraph& graph, const bdsg::handle_t& node) {
    vector<string> paths;
    graph.for_each_step_on_handle(node, [&](const bdsg::step_handle_t& step) {
        // Get the path associated with this step
        bdsg::path_handle_t path_handle = graph.get_path_handle_of_step(step);
        // Get the path name and add it to the result vector
        paths.emplace_back(graph.get_path_name(path_handle));
    });
    return paths;
}

void Trailmix::initializeParams(RunTreeProportionParams &params, shared_ptr<Trailmix_struct>& dta){

       params.freqs['A'] = 0.313122;
       params.freqs['C'] = 0.255707;
       params.freqs['G'] = 0.153784;
       params.freqs['T'] = 0.277387;
       params.freqs['R'] = params.freqs['A'] + params.freqs['G'];
       params.freqs['Y'] = params.freqs['C'] + params.freqs['T'];
       params.freqs['M'] = 1/(2*(   ((22)*(params.freqs['A']*params.freqs['G'])) + ((22)*(params.freqs['C']*params.freqs['T'])) + \
       (params.freqs['A']*params.freqs['C'] + (params.freqs['A']*params.freqs['T']) \
                         + (params.freqs['G']*params.freqs['C'] + (params.freqs['G'] * params.freqs['T']) )) ) );
        params.tr = dta->tree;
        params.soibean = false;
        params.root = dta->tree->root;
        params.chains = dta->chains;
        params.maxIter = dta->iter;
        params.burn = dta->burnin;
}

// Helper function to trim whitespace from both ends of a string
std::string trim(const std::string& str) {
    std::string out = str;
    out.erase(out.begin(), std::find_if(out.begin(), out.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    out.erase(std::find_if(out.rbegin(), out.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), out.end());
    return out;
}

// Helper function to remove special characters from a string
std::string remove_special_characters(const std::string& str) {
    std::string out;
    std::copy_if(str.begin(), str.end(), std::back_inserter(out), [](unsigned char c) {
        return std::isalnum(c) || c == '_'; // Keep alphanumeric and underscores, remove other characters
    });
    return out;
}

const int find_index(const std::vector<std::pair<std::string, double>>& vec, const std::string& str) {

    auto it = std::find_if(vec.begin(), vec.end(), [&str](const auto& p) {
        return p.first == str;
    });

    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    } else {
        return -1; // indicate not found
    }
}

void normalize(vector<double>& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x;
    }
    for (double& x : v) {
        x /= sum;
    }
}

void Trailmix::get_seed_source_estimates(shared_ptr<Trailmix_struct> &dta, const vector<string> &sigpaths){

vector<double> props(sigpaths.size(), 0.0);

for (unsigned int source=0; source < sigpaths.size(); ++source){
    string path_name = sigpaths[source];
    modifyPathNameInPlace(dta, path_name);
    //path_name.erase(std::remove(path_name.begin(), path_name.end(), ' '), path_name.end());
          const int find_idx = find_index(dta->tpms, path_name);
          if(find_idx == -1){
              props[source] = 0.0;
                            }
          else{
              props[source] = dta->tpms[find_idx].second;
              }
                                                            }
if (isnan(props[0])){
    for (auto &el : props){el = (1.0 / (dta->k));}
}

normalize(props);
if (dta->k > 1){
cerr << "Initial source proportion estimate: " << endl;
PRINTVEC(props)
               }
assert(props.size() == dta->k);

// Check for nan and replace with 1/k if needed
bool has_nan = false;
for (double val : props) {
    if (isnan(val)) {
        has_nan = true;
        break;
    }
}

if (has_nan) {
    double value = 1.0 / dta->k;
    for (double& val : props) {
        val = value;
    }
}


dta->seed = props;
return;
}
