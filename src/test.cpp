#define BOOST_TEST_MODULE haplocart_test
#define BOOST_TEST_MAIN
#include "Dup_Remover.h"
#include "HaploCart.h"
#include "readGAM.h"
#include "Euka.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "vgan_utils.h"
#include "bdsg/odgi.hpp"
#define PRINTVEC(v) for (int i=0; i<20; ++i){cerr << v[i] << '\t';}cerr << endl << endl;
using namespace vg;

bool is_convertible_to_double(const std::string& str) {
    try {
        std::stod(str);
        return true;
    } catch (const std::invalid_argument& ia) {
        // the string is not convertible to a double because it does not represent a valid number
        return false;
    } catch (const std::out_of_range& oor) {
        // the string is not convertible to a double because it is outside the range of representable values for a double
        return false;
    }
}

// Gam2Prof
std::pair<double, double> load_gam2prof(const std::string& filename)
{
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> data;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string field;
        std::vector<double> row;

        while (std::getline(iss, field, '\t')) {
            if (!is_convertible_to_double(field)){continue;}
            row.push_back(std::stod(field));
        }

        data.push_back(row);
    }

    double gam2prof_2nd_row_6th_col = data[1][5];
    double gam2prof_22nd_row_7th_col = data[21][6];

    //std::cerr << "gam2prof_2nd_row_6th_col: " << gam2prof_2nd_row_6th_col << std::endl;
    //std::cerr << "gam2prof_22nd_row_7th_col: " << gam2prof_22nd_row_7th_col << std::endl;

    return std::make_pair(gam2prof_2nd_row_6th_col, gam2prof_22nd_row_7th_col);
}

struct Fixture {
  Fixture() {
    boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_messages);
  }
};


// Euka
pair<vector<string>, vector<vector<double>>> load_detected_taxa_file(const string &detected_file_path){
    vector<string> taxa_names;
    vector<vector<double>> abundance_estimates;
    ifstream tsv_file(detected_file_path);
    string line;
    getline(tsv_file, line); // skip header line

    while (getline(tsv_file, line)) {
        stringstream ss(line);
        string taxa_name;
        string detected;
        int number_of_reads;
        double proportion_estimate;
        double ci85_lower;
        double ci85_upper;
        double ci95_lower;
        double ci95_upper;

        ss >> taxa_name >> detected >> number_of_reads >> proportion_estimate >> ci85_lower >> ci85_upper >> ci95_lower >> ci95_upper;
        if (taxa_name == ""){continue;}
        taxa_names.push_back(taxa_name);
        abundance_estimates.push_back(vector<double>{proportion_estimate, ci85_lower, ci85_upper, ci95_lower, ci95_upper});
    }
 return make_pair(taxa_names, abundance_estimates);

                                                         }

// HaploCart
const string get_haplocart_pred(const string &output_path) {
igzstream myfile;
myfile.open(output_path.c_str(), ios::in);
string line;
vector<string> tokens;
while (getline(myfile, line))
{
    if (line.size() == 0){continue;}
    getline(myfile, line);
    tokens= allTokensWhiteSpaces(line);
}
return tokens[1];
                                                           }


vector<string> get_haplocart_preds(const string& output_path) {
    vector<string> haplogroups;
    ifstream input_file(output_path);
    string line;

    // Skip the header line
    getline(input_file, line);

    // Read in the haplogroups from the second column of each subsequent line
    while (getline(input_file, line)) {
        stringstream ss(line);
        string field;

        // Skip the first field
        getline(ss, field, '\t');

        // Read the haplogroup from the second field
        getline(ss, field, '\t');
        haplogroups.push_back(field);
    }

    return haplogroups;
}


string get_sample_name(const string &output_path)   {
igzstream myfile;
myfile.open(output_path.c_str(), ios::in);
string line;
vector<string> tokens;
while (getline(myfile, line)) {
    if (line.size() == 0){continue;}
    getline(myfile, line);
    tokens= allTokensWhiteSpaces(line);
                              }
return tokens[0];
                                                    }


void hc_run_thread_test(const string & output_path, Haplocart * hc, const int n_threads) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back("../test/input_files/rCRS.fa");
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back(to_string(n_threads));
char** argvtopass = new char*[arguments.size()];
for (int i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }

hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                                         }



void hc_run_custom_posterior_output(const string & posterior_output_path, Haplocart * hc) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back("../test/input_files/rCRS.fq");
arguments.emplace_back("-o");
arguments.emplace_back("/dev/null");
arguments.emplace_back("-t");
arguments.emplace_back("20");
//arguments.emplace_back("-p");
arguments.emplace_back("-pf");
arguments.emplace_back(posterior_output_path);
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");

                                                                                           }

void hc_run_fasta(string input_path, string output_path, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                         }

void hc_run_fasta_bep(string input_path, string output_path, const double background_error_prob, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
arguments.emplace_back("-e");
arguments.emplace_back(to_string(background_error_prob));
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                         }

void hc_run_fq_single(string input_path, string output_path, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                             }

void hc_run_interleaved(string input_path, string output_path, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back(input_path);
arguments.emplace_back("-i");
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                             }


void hc_run_fq_paired(string input_path_1, string input_path_2, string output_path, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back(input_path_1);
arguments.emplace_back("-fq2");
arguments.emplace_back(input_path_2);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                                                    }

void hc_run_gam(string input_path_1, string output_path, Haplocart * hc, bool quiet) {
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-g");
arguments.emplace_back(input_path_1);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("20");
if(quiet){arguments.emplace_back("-q");}
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, getCWD(".")+"bin/");
                                                                         }

BOOST_AUTO_TEST_SUITE(haplocart)

BOOST_AUTO_TEST_CASE(fq_single_zipped)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/Q1_1.fq.gz";
  const string output_path = cwdProg + "test/output_files/Q1.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  BOOST_ASSERT(get_haplocart_pred(output_path) == "Q1");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "Q1_1.fq.gz");
}


BOOST_AUTO_TEST_CASE(multifasta)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/multifasta.fa";
  const string output_path = cwdProg + "test/output_files/multifasta.txt";
  vector<string> truth{"J2a1a1", "Z", "H2a2a1g"};
  hc_run_fasta(input_path, output_path, &hc, false);
  vector<string> preds = get_haplocart_preds(output_path);
  BOOST_CHECK_EQUAL_COLLECTIONS(preds.begin(), preds.end(), truth.begin(), truth.end());
}

BOOST_AUTO_TEST_CASE(load)
{
  Haplocart hc;
  const vector<vector<bool>> path_supports = hc.load_path_supports(getCWD(".")+"/share/vgan/hcfiles/");
  BOOST_CHECK_EQUAL(path_supports.size(), 11825);
  const vector<double> mappabilities = hc.load_mappabilities(getCWD(".")+"/share/vgan/hcfiles/");
  for (double x : mappabilities) {BOOST_CHECK_EQUAL((x <= 1 && x >= 0), true);}
  const vector<string> paths = hc.load_paths(getCWD(".")+"/share/vgan/hcfiles/");
  BOOST_CHECK_EQUAL(paths.size(), 5179);
  const map<string, vector<string>> parents = hc.load_parents(getCWD(".")+"/share/vgan/hcfiles/");
  BOOST_CHECK_EQUAL(parents.size(), 5437);
  const map<string, vector<string>> children = hc.load_children(getCWD(".")+"/share/vgan/hcfiles/");
  BOOST_CHECK_EQUAL(children.size(), 5438);
}

BOOST_AUTO_TEST_CASE(check_nodevec)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD(".")+"bin/");
  const auto [nodevec, minid, graph] = hc.readPathHandleGraph(cwdProg+"../share/vgan/hcfiles/graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
  BOOST_ASSERT(!nodevec.empty());
  for (int i = 0; i < nodevec.size(); ++i){
      BOOST_ASSERT(nodevec.at(i)->seq.size() <= 8);
                                          }
}


BOOST_AUTO_TEST_CASE(check_graph)
{
 Haplocart hc;
 const vector<string> paths = hc.load_paths(getCWD(".")+"/share/vgan/hcfiles/");
 bdsg::ODGI graph;
 const string cwdProg = getFullPath(getCWD(".")+"bin/");
 graph.deserialize(cwdProg+"../share/vgan/hcfiles/graph.og");
 const int minid = graph.min_node_id();
 const int maxid = graph.max_node_id();
 const int nodecount = graph.get_node_count();
 const int edgecount = graph.get_node_count();
 BOOST_CHECK_EQUAL(nodecount, 11821);
 BOOST_CHECK_EQUAL(minid, 1);
 BOOST_CHECK_EQUAL(maxid, 11821);
 for (const string & path : paths) {
     const string path_ = path.size() > 1 ? path : path + "_";
     BOOST_CHECK_EQUAL(graph.has_path(path_), true);
     BOOST_CHECK_EQUAL(graph.is_empty(graph.get_path_handle(path_)), false);
                                   }

     for (int i = minid; i<maxid; ++i) {
         BOOST_ASSERT(graph.get_degree(graph.get_handle(i), false) + graph.get_degree(graph.get_handle(i), true) > 0);
                                       }
}

BOOST_AUTO_TEST_CASE(invalid_bep1)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fa";
  BOOST_CHECK_THROW(hc_run_fasta_bep(input_path, "/dev/null", 2, &hc, true), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(invalid_bep2)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fa";
  BOOST_CHECK_THROW(hc_run_fasta_bep(input_path, "/dev/null", -2, &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(valid_bep)
{
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fa";
  Haplocart hc;
  BOOST_CHECK_NO_THROW(hc_run_fasta_bep(input_path, "/dev/null", 0.42, &hc, false));
}

BOOST_AUTO_TEST_CASE(fq_single_rcrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fq";
  const string output_path = cwdProg + "test/output_files/rCRS_fq.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "rCRS.fq");
}

BOOST_AUTO_TEST_CASE(fq_single_rsrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/RSRS.fq";
  const string output_path = cwdProg + "test/output_files/RSRS.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  const string out = get_haplocart_pred(output_path);
  BOOST_ASSERT(out == "L1'2'3'4'5'6" || out == "mt-MRCA");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "RSRS.fq");
}

BOOST_AUTO_TEST_CASE(consensus_rcrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fa";
  const string output_path = cwdProg + "test/output_files/rCRS_consensus.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(another_consensus)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/H2a2a1g.fa";
  const string output_path = cwdProg + "test/output_files/H2a2a1g.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1g");
}


BOOST_AUTO_TEST_CASE(zipped_consensus)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/Z.fa.gz";
  const string output_path = cwdProg + "test/output_files/Z.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Z");
}

BOOST_AUTO_TEST_CASE(consensus_wrong_format)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/rCRS.fq";
  const string output_path = cwdProg + "test/output_files/rCRS_bad.txt";
  BOOST_CHECK_THROW(hc_run_fasta(input_path, output_path, & hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(zipped_paired_fastq)

{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input1_path = cwdProg + "test/input_files/Q1_1.fq.gz";
  const string input2_path = cwdProg + "test/input_files/Q1_2.fq.gz";
  const string output_path = cwdProg + "test/output_files/Q1_paired.txt";
  hc_run_fq_paired(input1_path, input2_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Q1");
}

BOOST_AUTO_TEST_CASE(gam)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/alignments/J2a1a1a1.gam";
  const string output_path = cwdProg + "test/output_files/J2a1a1a1.txt";
  hc_run_gam(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "J2a1a1a1");
}

BOOST_AUTO_TEST_CASE(check_thread_minus_one)
{
  Haplocart hc;
  hc_run_thread_test("../test/output_files/thread_minus_one.txt", &hc, -1);
  BOOST_CHECK_EQUAL(get_haplocart_pred("../test/output_files/thread_minus_one.txt"), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(check_thread_zero)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string output_path = cwdProg + "test/output_files/thread_zero.txt";
  BOOST_CHECK_THROW(hc_run_thread_test(output_path, &hc, 0), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(check_thread_too_many)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string output_path = cwdProg + "test/output_files/thread_too_many.txt";
  hc_run_thread_test(output_path, &hc, 424242);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(custom_posterior_output)
{
  Haplocart hc;
  hc_run_custom_posterior_output("../test/output_files/rCRS_posterior.txt", & hc);
  BOOST_CHECK_EQUAL(std::filesystem::exists("../test/output_files/rCRS_posterior.txt"), true);
}

BOOST_AUTO_TEST_CASE(missing_input_consensus)
{
  Haplocart hc;
  BOOST_CHECK_THROW(hc_run_fasta("not_a_real_file.fa", "/dev/null", &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(missing_input_fq1)
{
  Haplocart hc;
  BOOST_CHECK_THROW(hc_run_fq_single("not_a_real_file.fq", "/dev/null", &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(missing_input_fq2)
{
  Haplocart hc;
  BOOST_CHECK_THROW(hc_run_fq_paired("not_a_real_file.fq", "also_not_a_real_file.fq", "/dev/null", &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(missing_input_gam)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  BOOST_CHECK_THROW(hc_run_gam("not_a_real_file.gam", "/dev/null", &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(multifasta_zipped)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/multifasta.fa.gz";
  const string output_path = cwdProg + "test/output_files/multifasta_zipped.txt";
  hc_run_fasta(input_path, output_path, &hc, false);
  vector<string> truth {"J2a1a1", "Z", "H2a2a1g"};
  vector<string> preds = get_haplocart_preds(output_path);
  BOOST_CHECK_EQUAL_COLLECTIONS(preds.begin(), preds.end(), truth.begin(), truth.end());
}

BOOST_AUTO_TEST_CASE(interleaved)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/Q1_interleaved.fq";
  const string output_path = cwdProg + "test/output_files/Q1_interleaved.txt";
  hc_run_interleaved(input_path, output_path, &hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Q1");
}

BOOST_AUTO_TEST_SUITE_END()


/// BEGIN TESTING GRAPH RECONSTRUCTION ///

BOOST_AUTO_TEST_SUITE(reconstruction)

/// BEGIN TESTING GRAPH RECONSTRUCTION ///

BOOST_AUTO_TEST_CASE(plus_strand_perfect_match)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(0)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(0)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_mismatch)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(1)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(1)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
     BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCGTGAGTAGGGTCCACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_insert_in_read)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(2)->path;
    cerr << "algnvector size: " << read_info->size() << endl;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(2)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
     BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCCCCGTGAGTAGGGTCGACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATA---CCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_deletion_in_read)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(3)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(3)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TGGGTGGAGCGCGCCCCATTTATACCGTGAGTAGGGTCGACCAAGAACCGCAAGA");
     BOOST_CHECK_EQUAL(read_seq, "TGGGTGGAGCGCGCCCCAT--------TGAGTAGGGTCGACCAAGAACCGCAAGA");
}

BOOST_AUTO_TEST_CASE(plus_strand_softclip)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(4)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(4)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "SSSSSSSSSSSSSSSSSSSSSSSSCGGATATAAACGCCAGGTTGAATCCGCATTT");
     BOOST_CHECK_EQUAL(read_seq, "CGGCTGTCAGCTGCCGTCTGCGTACGGATATAAACGCCAGGTTGAATCCGCATTT");
}

BOOST_AUTO_TEST_CASE(minus_strand_perfect_match)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(5)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(5)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_mismatch)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(6)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(6)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGCGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_insert_in_read)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(7)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(7)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTC------------GACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCCAGTCAGTCAGTGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_deletion_in_read)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(8)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(8)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTA---------TAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_softclip)
{
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    auto read_info = readGAM((cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str(), false, false, "/dev/null", "");
    const auto path = read_info->at(9)->path;
    Haplocart hc;
    const auto pathhandlegraphtuple = hc.readPathHandleGraph(cwdProg+"../test/reconstructInputSeq/target_graph.og", 1, cwdProg + "../share/vgan/hcfiles/");
    auto [nodevector, minid, graph] = pathhandlegraphtuple;
    auto ret_tuple = reconstruct_graph_sequence(graph, path, read_info->at(9)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "SSSSSSSSSSSSSSSSSSSSSSSSSSCACCGTAATCCATGCTTGATTGAGACCGCC");
     BOOST_CHECK_EQUAL(read_seq, "CTAGCTGCAGTCGCGCTCGTCATGCACACCGTAATCCATGCTTGATTGAGACCGCC");
}


/// END TESTING GRAPH RECONSTRUCTION ///
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////// TEST EUKA //////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE(euka)


                 ////////////// Output sanity checks //////////////////////

void run_detect_correct_taxa_no_damage(Euka* ek){
vector<string> euka_argvec;
euka_argvec.emplace_back("vgan");
euka_argvec.emplace_back("euka");
euka_argvec.emplace_back("-fq1");
euka_argvec.emplace_back(getCWD(".")+"bin/" + "../test/input_files/euka/three_none_100.fq.gz");
euka_argvec.emplace_back("-t");
euka_argvec.emplace_back("-1");
euka_argvec.emplace_back("-o");
euka_argvec.emplace_back(getCWD(".")+"bin/" + "../test/output_files/euka/three_none_100");

char** argvtopass = new char*[euka_argvec.size()];
for (int i=0;i<euka_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(euka_argvec[i].c_str());
                                       }

ek->run(euka_argvec.size(), argvtopass, getCWD(".")+"bin/");
                                                  }

BOOST_AUTO_TEST_CASE(detect_correct_taxa_no_damage){

Euka ek;
run_detect_correct_taxa_no_damage(&ek);

auto [taxa_names, abundance_estimates] = load_detected_taxa_file(getCWD(".")+"bin/" + "../test/output_files/euka/three_none_100_detected.tsv");

assert(taxa_names[0] == "Bovidae");
assert(taxa_names[1] == "Myotis");
assert(taxa_names[2] == "Ursidae");
assert(abundance_estimates[0][0] >= 0.01 && abundance_estimates[0][0] <= 0.05);
assert(abundance_estimates[1][0] >= 0.2 && abundance_estimates[1][0] <= 0.3);
assert(abundance_estimates[2][0] >= 0.68 && abundance_estimates[2][0] <= 0.78);
                                                   }


BOOST_AUTO_TEST_CASE(detect_correct_taxa_mid_damage){
Euka ek;
vector<string> euka_argvec;
euka_argvec.emplace_back("vgan");
euka_argvec.emplace_back("euka");
euka_argvec.emplace_back("-fq1");
euka_argvec.emplace_back(getCWD(".") +"bin/" + "../test/input_files/euka/three_dmid_100.fq.gz");
euka_argvec.emplace_back("-t");
euka_argvec.emplace_back("-1");
euka_argvec.emplace_back("-o");
euka_argvec.emplace_back(getCWD(".") +"bin/" + "../test/output_files/euka/three_dmid_100");

char** argvtopass = new char*[euka_argvec.size()];
for (int i=0;i<euka_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(euka_argvec[i].c_str());
                                       }
ek.run(euka_argvec.size(), argvtopass, getCWD(".")+"bin/");

auto [taxa_names, abundance_estimates] = load_detected_taxa_file(getCWD(".")+"bin/" + "../test/output_files/euka/three_dmid_100_detected.tsv");

BOOST_ASSERT(taxa_names[0] == "Bovidae");
BOOST_ASSERT(taxa_names[1] == "Myotis");
BOOST_ASSERT(taxa_names[2] == "Ursidae");

BOOST_ASSERT(abundance_estimates[0][0] >= 0.01 && abundance_estimates[0][0] <= 0.05);
BOOST_ASSERT(abundance_estimates[1][0] >= 0.2 && abundance_estimates[1][0] <= 0.3);
BOOST_ASSERT(abundance_estimates[2][0] >= 0.68 && abundance_estimates[2][0] <= 0.78);

                                                   }

BOOST_AUTO_TEST_CASE(detect_correct_taxa_high_damage){
Euka ek;
vector<string> euka_argvec;
euka_argvec.emplace_back("vgan");
euka_argvec.emplace_back("euka");
euka_argvec.emplace_back("-fq1");
euka_argvec.emplace_back(getCWD(".") + "bin/" +"../test/input_files/euka/three_dhigh_100.fq.gz");
euka_argvec.emplace_back("-t");
euka_argvec.emplace_back("-1");
euka_argvec.emplace_back("-o");
euka_argvec.emplace_back(getCWD(".")+"bin/" + "../test/output_files/euka/three_dhigh_100");

char** argvtopass = new char*[euka_argvec.size()];
for (int i=0;i<euka_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(euka_argvec[i].c_str());
                                       }
ek.run(euka_argvec.size(), argvtopass, getCWD(".")+"bin/");

auto [taxa_names, abundance_estimates] = load_detected_taxa_file(getCWD(".") +"bin/" + "../test/output_files/euka/three_dhigh_100_detected.tsv");

BOOST_ASSERT(taxa_names[0] == "Bovidae");
BOOST_ASSERT(taxa_names[1] == "Myotis");
BOOST_ASSERT(taxa_names[2] == "Ursidae");

BOOST_ASSERT(abundance_estimates[0][0] >= 0.01 && abundance_estimates[0][0] <= 0.05);
BOOST_ASSERT(abundance_estimates[1][0] >= 0.2 && abundance_estimates[1][0] <= 0.3);
BOOST_ASSERT(abundance_estimates[2][0] >= 0.68 && abundance_estimates[2][0] <= 0.78);
                                                   }

BOOST_AUTO_TEST_CASE(detect_formicidae){
Euka ek;
vector<string> euka_argvec;
euka_argvec.emplace_back("vgan");
euka_argvec.emplace_back("euka");
euka_argvec.emplace_back("-fq1");
euka_argvec.emplace_back(getCWD(".") + "test/input_files/euka/formicidae.fq.gz");
euka_argvec.emplace_back("-t");
euka_argvec.emplace_back("-1");
euka_argvec.emplace_back("-o");
euka_argvec.emplace_back(getCWD(".") + "test/output_files/euka/formicidae");
euka_argvec.emplace_back("--minBins");
euka_argvec.emplace_back("2");

char** argvtopass = new char*[euka_argvec.size()];
for (int i=0;i<euka_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(euka_argvec[i].c_str());
                                       }
ek.run(euka_argvec.size(), argvtopass, getCWD(".")+"bin/");

auto [taxa_names, abundance_estimates] = load_detected_taxa_file(getCWD(".") + "bin/" +"../test/output_files/euka/formicidae_detected.tsv");
cerr << abundance_estimates[0][0] << " " << taxa_names[0] << endl;
BOOST_ASSERT(abundance_estimates[0][0] == 1);
BOOST_ASSERT(taxa_names[0] == "Formicidae");

                                    }

BOOST_AUTO_TEST_CASE(detect_formicidae2){
Euka ek;
vector<string> euka_argvec;
euka_argvec.emplace_back("vgan");
euka_argvec.emplace_back("euka");
euka_argvec.emplace_back("-fq1");
euka_argvec.emplace_back(getCWD(".") + "test/input_files/euka/formicidae.fq.gz");
euka_argvec.emplace_back("-t");
euka_argvec.emplace_back("-1");
euka_argvec.emplace_back("-o");
euka_argvec.emplace_back(getCWD(".") + "test/output_files/euka/formicidae2");
euka_argvec.emplace_back("--entropy");
euka_argvec.emplace_back("0.9");

char** argvtopass = new char*[euka_argvec.size()];
for (int i=0;i<euka_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(euka_argvec[i].c_str());
                                       }
ek.run(euka_argvec.size(), argvtopass, getCWD(".")+"bin/");

auto [taxa_names, abundance_estimates] = load_detected_taxa_file(getCWD(".")+"bin/" + "../test/output_files/euka/formicidae2_detected.tsv");

BOOST_ASSERT(abundance_estimates[0][0] == 1);
BOOST_ASSERT(taxa_names[0] == "Formicidae");

                                    }

               ////////////// End output sanity checks //////////////////////


BOOST_AUTO_TEST_CASE(load) {
    const string path_support_path = "../share/vgan/euka_dir/euka_db_graph_path_supports";
    const string clade_info_path = "../share/vgan/euka_dir/euka_db.clade";
    const string clade_chunk_path = "../share/vgan/euka_dir/euka_db.bins";

    Euka ek;
    auto clade_info = ek.load_clade_info(clade_info_path, 5);
    auto path_supports = ek.load_path_supports_Euka(path_support_path);
    auto bins = ek.load_clade_chunks(clade_chunk_path);
    BOOST_ASSERT(clade_info->size() > 0);
    BOOST_ASSERT(path_supports.size() > 0);
    BOOST_ASSERT(bins.size() > 0);
                           }

BOOST_AUTO_TEST_SUITE_END()
//////////////////////////////////////////// END TEST EUKA ///////////////////////////////////////////////////////

