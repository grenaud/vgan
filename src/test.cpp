#define BOOST_TEST_MODULE vgan_test
#define BOOST_TEST_MAIN
#include "Dup_Remover.h"
#include "HaploCart.h"
#include "readGAM.h"
#include "Euka.h"
#include "soibean.h"
#include "TrailMix.h"
#include "MCMC.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "vgan_utils.h"
#include "bdsg/odgi.hpp"
#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << setprecision(10) << v[i] << '\t';}cerr << endl << endl;
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




/////////////////////////////////////////// BEGIN TEST SOIBEAN ///////////////////////////////////////////////////
struct BranchPlacementRecord {
    std::string source;
    int chain;
    double meanBranchPosition;
    double meanBranchPositionCI;
    double medianBranchPosition;
    double medianBranchPositionCI;
    double effectiveSampleSize;
    double autocorrelationVariance;
};

struct ProportionEstimateRecord {
    std::string source;
    int chain;
    double meanProportionEstimate;
    double meanProportionEstimateCI5;
    double medianProportionEstimate;
    double medianProportionEstimateCI95;
    double effectiveSampleSize;
    double autocorrelationVariance;
};

// Struct for Chain Diagnostics data
struct ChainDiagnosticsRecord {
    std::string source;
    double highestLogLikelihood;
    int rhatProportionEstimate;
    double rhatBranchPositionEstimate;
};


// Function to split a string by whitespace
std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> result;
    std::istringstream iss(s);
    for (std::string token; iss >> token; )
        result.push_back(token);
    return result;
}

auto load_result_file(const string &file_path){

                                                                                      }

auto load_trace_file(const string &file_path){

                                                                                      }

std::vector<ChainDiagnosticsRecord> loadChainDiagnosticsData(const std::string& filename) {
    std::vector<ChainDiagnosticsRecord> records;
    std::ifstream file(filename);
    std::string line;

    // Read and discard the header line
    getline(file, line);

    while (getline(file, line)) {
        auto tokens = split(line);
        if (tokens.size() >= 5) {
            ChainDiagnosticsRecord record;
            record.source = tokens[0];
            record.highestLogLikelihood = std::stod(tokens[1]);
            record.rhatProportionEstimate = std::stod(tokens[2]);
            record.rhatBranchPositionEstimate = std::stod(tokens[3]);
            records.push_back(record);
        }
    }
    return records;
}

std::vector<ProportionEstimateRecord> load_prop_diagnostics_file(const string &file_path) {
    std::vector<ProportionEstimateRecord> records;
    std::ifstream file(file_path);
    std::string line;

    // Read and discard the header line
    getline(file, line);

    while (getline(file, line)) {
        auto tokens = split(line);

        if (tokens.size() >= 9) {
            ProportionEstimateRecord record;
            record.source = tokens[0];
            record.chain = std::stoi(tokens[1]);
            record.meanProportionEstimate = std::stod(tokens[2]);
            record.meanProportionEstimateCI5 = std::stod(tokens[3]);
            record.medianProportionEstimate = std::stod(tokens[4]);
            record.medianProportionEstimateCI95 = std::stod(tokens[5]);

            try {
                record.effectiveSampleSize = std::stod(tokens[6]);
                record.autocorrelationVariance = std::stod(tokens[7]);
            } catch (const std::invalid_argument&) {
                record.effectiveSampleSize = -1;
                record.autocorrelationVariance = -1;
            }
            records.push_back(record);
        }
    }
    return records;
}


std::vector<BranchPlacementRecord> load_branch_placement_diagnostics_file(const std::string &file_path) {
    std::vector<BranchPlacementRecord> records;
    std::ifstream file(file_path);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file: " << file_path << std::endl;
        return records; // Return empty if file can't be opened
    }

    int line_number = 0;
    while (std::getline(file, line)) {
        line_number++;

        // Skip headers and empty lines
        if (line.empty() || line.rfind("Source", 0) == 0) {
            std::cerr << "Skipping header or empty line." << std::endl;
            continue;
        }

        std::istringstream iss(line);
        BranchPlacementRecord r;

        if (iss >> r.source >> r.chain >> r.meanBranchPosition >> r.meanBranchPositionCI
                >> r.medianBranchPosition >> r.medianBranchPositionCI
                >> r.effectiveSampleSize >> r.autocorrelationVariance) {
            //std::cerr << "Parsed record: " << std::endl;
            //std::cerr << "\tSource: " << r.source << std::endl;
            //std::cerr << "\tChain: " << r.chain << std::endl;
            //std::cerr << "\tMean Branch Position: " << r.meanBranchPosition << std::endl;
            //std::cerr << "\tMean Branch Position CI: " << r.meanBranchPositionCI << std::endl;
            //std::cerr << "\tMedian Branch Position: " << r.medianBranchPosition << std::endl;
            //std::cerr << "\tMedian Branch Position CI: " << r.medianBranchPositionCI << std::endl;
            //std::cerr << "\tEffective Sample Size: " << r.effectiveSampleSize << std::endl;
            //std::cerr << "\tAutocorrelation Variance: " << r.autocorrelationVariance << std::endl;

            records.emplace_back(r);
            //std::cerr << "Record added to records vector." << std::endl;
        } else {
            //std::cerr << "Warning: Failed to parse line " << line_number << ": " << line << std::endl;
        }
    }

    //std::cerr << "Finished reading file. Total records loaded: " << records.size() << std::endl;
    return records;
}



///////////////////////////////////////////// TEST SOIBEAN //////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE(Soibean)


void run_k1(soibean * sb){


vector<string> soibean_argvec;
soibean_argvec.emplace_back("vgan");
soibean_argvec.emplace_back("soibean");
soibean_argvec.emplace_back("-fq1");
soibean_argvec.emplace_back(getCWD(".")+"bin/" + "../test/input_files/soibean/k1.fq.gz");
soibean_argvec.emplace_back("-t");
soibean_argvec.emplace_back("50");
soibean_argvec.emplace_back("-o");
soibean_argvec.emplace_back(getCWD(".")+"bin/" + "../test/output_files/soibean/k1");
soibean_argvec.emplace_back("--dbprefix");
soibean_argvec.emplace_back("Ursidae");
soibean_argvec.emplace_back("--iter");
soibean_argvec.emplace_back("1000");
soibean_argvec.emplace_back("--burnin");
soibean_argvec.emplace_back("150");
soibean_argvec.emplace_back("--chains");
soibean_argvec.emplace_back("1");


char** argvtopass = new char*[soibean_argvec.size()];
for (int i=0;i<soibean_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(soibean_argvec[i].c_str());
                                       }

sb->run(soibean_argvec.size(), argvtopass, getCWD(".")+"bin/");


                                                  }

void run_k2(soibean * sb){


vector<string> soibean_argvec;
soibean_argvec.emplace_back("vgan");
soibean_argvec.emplace_back("soibean");
soibean_argvec.emplace_back("-fq1");
soibean_argvec.emplace_back(getCWD(".")+"bin/" + "../test/input_files/soibean/k2.fq.gz");
soibean_argvec.emplace_back("-t");
soibean_argvec.emplace_back("50");
soibean_argvec.emplace_back("-o");
soibean_argvec.emplace_back(getCWD(".")+"bin/" + "../test/output_files/soibean/k2");
soibean_argvec.emplace_back("--dbprefix");
soibean_argvec.emplace_back("Ursidae");
soibean_argvec.emplace_back("--iter");
soibean_argvec.emplace_back("1000");
soibean_argvec.emplace_back("--burnin");
soibean_argvec.emplace_back("150");
soibean_argvec.emplace_back("--chains");
soibean_argvec.emplace_back("1");


char** argvtopass = new char*[soibean_argvec.size()];
for (int i=0;i<soibean_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(soibean_argvec[i].c_str());
                                       }

sb->run(soibean_argvec.size(), argvtopass, getCWD(".")+"bin/");


                                                  }
BOOST_AUTO_TEST_CASE(k2)
{
    soibean sb;

    run_k2(&sb);
    auto branchRecords = load_branch_placement_diagnostics_file(getCWD(".")+"bin/" + "../test/output_files/soibean/k2BranchEstimate2.txt");


    auto propRecords = load_prop_diagnostics_file(getCWD(".")+"bin/" + "../test/output_files/soibean/k2ProportionEstimates2.txt");
    double expectedValue = 0.5; // set the expected value
    double tolerancePercent = 1.0; // set the tolerance as a percentage (e.g., 0.01 for 1%)
    BOOST_CHECK_CLOSE(propRecords[0].meanProportionEstimate, expectedValue, tolerancePercent);
    BOOST_CHECK_CLOSE(propRecords[1].meanProportionEstimate, expectedValue, tolerancePercent);

    BOOST_ASSERT(branchRecords[0].meanBranchPosition > 0.9);
    BOOST_ASSERT(branchRecords[1].meanBranchPosition > 0.9);
    //BOOST_ASSERT(chainRecords[0].rhatBranchPositionEstimate < 10);
    BOOST_ASSERT(branchRecords[0].source == "NC_003426.1_Ursus_americanus_mitochondrion__complete_genome" || branchRecords[0].source == "NC_009971.1_Ursus_thibetanus_mitochondrion__complete_genome");
    BOOST_ASSERT(branchRecords[1].source == "NC_009971.1_Ursus_thibetanus_mitochondrion__complete_genome" || branchRecords[1].source == "NC_003426.1_Ursus_americanus_mitochondrion__complete_genome");

}

BOOST_AUTO_TEST_CASE(k1)
{
    soibean sb;
    run_k1(&sb);
    auto branchRecords = load_branch_placement_diagnostics_file(getCWD(".")+"bin/" + "../test/output_files/soibean/k1BranchEstimate1.txt");


    auto propRecords = load_prop_diagnostics_file(getCWD(".")+"bin/" + "../test/output_files/soibean/k1ProportionEstimates1.txt");

    BOOST_ASSERT(propRecords[0].meanProportionEstimate == 1.0);

    BOOST_ASSERT(branchRecords[0].meanBranchPosition > 0.9);
    //BOOST_ASSERT(chainRecords[0].rhatBranchPositionEstimate < 10);
    BOOST_ASSERT(branchRecords[0].source == "NC_003426.1_Ursus_americanus_mitochondrion__complete_genome");

}



BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Output for k 1: right branch ? prop = 1 and pos > 0.9
/// Output for k 2: right branches ? prop between 0.4 and 0.6 and sum to 1 and pos > 0.9 


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


// Utility to tokenize a string based on any whitespace (spaces or tabs)
std::vector<std::string> tokenize(const std::string &str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string get_haplocart_pred(const std::string &output_path) {
    std::ifstream myfile(output_path);
    if (!myfile.is_open()) {
        std::cerr << "Failed to open file: " << output_path << std::endl;
        return "";
    }

    std::string line;

    // 1. Skip any initial empty lines
    while (std::getline(myfile, line)) {
        if (line.empty()) {
            continue;
        } else {
            break; // Exit the loop after finding the first non-empty line
        }
    }

    // 3. Process the data lines
    while (std::getline(myfile, line)) {

        // Skip empty lines within data
        if (line.empty()) {
            continue;
        }

        // Remove potential carriage return (for Windows-formatted files)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Tokenize the line
        std::vector<std::string> tokens = tokenize(line);
        std::cerr << std::endl;

        // Assuming the predicted haplotype is the second token (index 1)
        if (tokens.size() > 1) {
            return tokens[1];  // Return the predicted haplotype
        }
    }

    // If no valid data lines are found
    std::cerr << "No valid data lines found." << std::endl;
    return "";
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
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back("../test/input_files/haplocart/rCRS.fa");
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back(to_string(n_threads));
char** argvtopass = new char*[arguments.size()];
for (int i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }

hc->run(arguments.size(), argvtopass, dta);
                                                                                         }



void hc_run_custom_posterior_output(const string & posterior_output_path, Haplocart * hc) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back("../test/input_files/haplocart/rCRS.fq");
arguments.emplace_back("-o");
arguments.emplace_back("/dev/null");
arguments.emplace_back("-t");
arguments.emplace_back("50");
arguments.emplace_back("-pf");
arguments.emplace_back(posterior_output_path);
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, dta);

                                                                                           }

void hc_run_fasta(string input_path, string output_path, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("50");
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, dta);
                                                                         }

void hc_run_fasta_bep(string input_path, string output_path, const double background_error_prob, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-f");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("50");
arguments.emplace_back("-e");
arguments.emplace_back(to_string(background_error_prob));
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, dta);
                                                                         }

void hc_run_fq_single(string input_path, string output_path, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back(input_path);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("50");
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, dta);
                                                                             }

void hc_run_interleaved(string input_path, string output_path, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-fq1");
arguments.emplace_back(input_path);
arguments.emplace_back("-i");
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("50");
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
 hc->run(arguments.size(), argvtopass, dta);
                                                                             }


void hc_run_fq_paired(string input_path_1, string input_path_2, string output_path, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
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
arguments.emplace_back("50");
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, dta);
                                                                                                    }

void hc_run_gam(string input_path_1, string output_path, Haplocart * hc, bool quiet) {
shared_ptr dta = make_unique<Trailmix_struct>();
vector<string> arguments;
arguments.emplace_back("vgan");
arguments.emplace_back("haplocart");
arguments.emplace_back("-g");
arguments.emplace_back(input_path_1);
arguments.emplace_back("-o");
arguments.emplace_back(output_path);
arguments.emplace_back("-t");
arguments.emplace_back("50");
if(quiet){arguments.emplace_back("-q");}
arguments.emplace_back("-z");
arguments.emplace_back("tempdir");
char** argvtopass = new char*[arguments.size()];
for (size_t i=0;i<arguments.size();i++) {
                   argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                     }
hc->run(arguments.size(), argvtopass, dta);
                                                                         }

BOOST_AUTO_TEST_SUITE(haplocart)

BOOST_AUTO_TEST_CASE(fq_single_zipped)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/Q1_1.fq.gz";
  const string output_path = cwdProg + "test/output_files/haplocart/Q1.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  BOOST_ASSERT(get_haplocart_pred(output_path) == "Q1");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "Q1_1.fq.gz");
}


BOOST_AUTO_TEST_CASE(multifasta)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/multifasta.fa";
  const string output_path = cwdProg + "test/output_files/haplocart/multifasta.txt";
  vector<string> truth{"J2a1a1", "Z", "H2a2a1g"};
  hc_run_fasta(input_path, output_path, &hc, false);
  vector<string> preds = get_haplocart_preds(output_path);
  BOOST_CHECK_EQUAL_COLLECTIONS(preds.begin(), preds.end(), truth.begin(), truth.end());
}

BOOST_AUTO_TEST_CASE(load)
{
  Haplocart hc;
  shared_ptr dta = make_unique<Trailmix_struct>();
  hc.load_path_supports(dta);
  BOOST_CHECK_EQUAL(dta->path_supports.size(), 15725);
  hc.load_mappabilities(dta);
  for (double x : dta->mappabilities) {BOOST_CHECK_EQUAL((x <= 1 && x >= 0), true);}
  hc.load_path_names(dta);
  BOOST_CHECK_EQUAL(dta->path_names.size(), 5179);
  hc.load_parents(dta);
  BOOST_CHECK_EQUAL(dta->parents.size(), 5437);
  hc.load_children(dta);
  BOOST_CHECK_EQUAL(dta->children.size(), 5438);
}

BOOST_AUTO_TEST_CASE(check_nodevec)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD(".")+"bin/");
  shared_ptr dta = make_unique<Trailmix_struct>();
  hc.readPathHandleGraph(dta);
  BOOST_ASSERT(!dta->nodevector.empty());
  for (int i = 0; i < dta->nodevector.size(); ++i){
      BOOST_ASSERT(dta->nodevector.at(i)->seq.size() <= 8);
                                          }
}


BOOST_AUTO_TEST_CASE(check_graph)
{
 Haplocart hc;
 shared_ptr dta = make_unique<Trailmix_struct>();
 dta->running_trailmix = false;
 hc.read_PHG(dta);
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
 for (const string & path : dta->path_names) {
     const string path_ = path.size() > 1 ? path : path + "_";
     BOOST_CHECK_EQUAL(graph.has_path(path), true);
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
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fa";
  BOOST_CHECK_THROW(hc_run_fasta_bep(input_path, "/dev/null", 2, &hc, true), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(invalid_bep2)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fa";
  BOOST_CHECK_THROW(hc_run_fasta_bep(input_path, "/dev/null", -2, &hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(valid_bep)
{
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fa";
  Haplocart hc;
  BOOST_CHECK_NO_THROW(hc_run_fasta_bep(input_path, "/dev/null", 0.42, &hc, false));
}

BOOST_AUTO_TEST_CASE(fq_single_rcrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fq";
  const string output_path = cwdProg + "test/output_files/haplocart/rCRS_fq.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "rCRS.fq");
}

BOOST_AUTO_TEST_CASE(fq_single_rsrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/RSRS.fq";
  const string output_path = cwdProg + "test/output_files/haplocart/RSRS.txt";
  hc_run_fq_single(input_path, output_path, & hc, false);
  const string out = get_haplocart_pred(output_path);
  BOOST_ASSERT(out == "L1'2'3'4'5'6" || out == "mt-MRCA");
  BOOST_CHECK_EQUAL(get_sample_name(output_path), "RSRS.fq");
}

BOOST_AUTO_TEST_CASE(consensus_rcrs)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fa";
  const string output_path = cwdProg + "test/output_files/rCRS_consensus.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(another_consensus)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/H2a2a1g.fa";
  const string output_path = cwdProg + "test/output_files/H2a2a1g.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1g");
}


BOOST_AUTO_TEST_CASE(zipped_consensus)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/Z.fa.gz";
  const string output_path = cwdProg + "test/output_files/Z.txt";
  hc_run_fasta(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Z");
}

BOOST_AUTO_TEST_CASE(consensus_wrong_format)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/rCRS.fq";
  const string output_path = cwdProg + "test/output_files/haplocart/rCRS_bad.txt";
  BOOST_CHECK_THROW(hc_run_fasta(input_path, output_path, & hc, true), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(zipped_paired_fastq)

{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input1_path = cwdProg + "test/input_files/haplocart/Q1_1.fq.gz";
  const string input2_path = cwdProg + "test/input_files/haplocart/Q1_2.fq.gz";
  const string output_path = cwdProg + "test/output_files/haplocart/Q1_paired.txt";
  hc_run_fq_paired(input1_path, input2_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Q1");
}

BOOST_AUTO_TEST_CASE(gam)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/alignments/J2a1a1a1.gam";
  const string output_path = cwdProg + "test/output_files/haplocart/J2a1a1a1.txt";
  hc_run_gam(input_path, output_path, & hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "J2a1a1a1");
}

BOOST_AUTO_TEST_CASE(check_thread_minus_one)
{
  Haplocart hc;
  hc_run_thread_test("../test/output_files/haplocart/thread_minus_one.txt", &hc, -1);
  BOOST_CHECK_EQUAL(get_haplocart_pred("../test/output_files/haplocart/thread_minus_one.txt"), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(check_thread_zero)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string output_path = cwdProg + "test/output_files/haplocart/thread_zero.txt";
  BOOST_CHECK_THROW(hc_run_thread_test(output_path, &hc, 0), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(check_thread_too_many)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string output_path = cwdProg + "test/output_files/haplocart/thread_too_many.txt";
  hc_run_thread_test(output_path, &hc, 424242);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "H2a2a1");
}

BOOST_AUTO_TEST_CASE(custom_posterior_output)
{
  Haplocart hc;
  hc_run_custom_posterior_output("../test/output_files/haplocart/rCRS_posterior.txt", & hc);
  BOOST_CHECK_EQUAL(std::filesystem::exists("../test/output_files/haplocart/rCRS_posterior.txt"), true);
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
  const string input_path = cwdProg + "test/input_files/haplocart/multifasta.fa.gz";
  const string output_path = cwdProg + "test/output_files/haplocart/multifasta_zipped.txt";
  hc_run_fasta(input_path, output_path, &hc, false);
  vector<string> truth {"J2a1a1", "Z", "H2a2a1g"};
  vector<string> preds = get_haplocart_preds(output_path);
  BOOST_CHECK_EQUAL_COLLECTIONS(preds.begin(), preds.end(), truth.begin(), truth.end());
}

BOOST_AUTO_TEST_CASE(interleaved)
{
  Haplocart hc;
  const string cwdProg = getFullPath(getCWD("."));
  const string input_path = cwdProg + "test/input_files/haplocart/Q1_interleaved.fq";
  const string output_path = cwdProg + "test/output_files/haplocart/Q1_interleaved.txt";
  hc_run_interleaved(input_path, output_path, &hc, false);
  BOOST_CHECK_EQUAL(get_haplocart_pred(output_path), "Q1");
}

BOOST_AUTO_TEST_SUITE_END()


/// BEGIN TESTING GRAPH RECONSTRUCTION ///

BOOST_AUTO_TEST_SUITE(reconstruction)

/// BEGIN TESTING GRAPH RECONSTRUCTION ///

BOOST_AUTO_TEST_CASE(plus_strand_perfect_match)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(0)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(0)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_mismatch)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(1)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(1)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
     BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCGTGAGTAGGGTCCACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATACCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_insert_in_read)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(2)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(2)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
     BOOST_CHECK_EQUAL(read_seq, "CCCCATTTATACCCCCGTGAGTAGGGTCGACCAAGAAC");
    BOOST_CHECK_EQUAL(graph_seq, "CCCCATTTATA---CCGTGAGTAGGGTCGACCAAGAAC");
}

BOOST_AUTO_TEST_CASE(plus_strand_deletion_in_read)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(3)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(3)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TGGGTGGAGCGCGCCCCATTTATACCGTGAGTAGGGTCGACCAAGAACCGCAAGA");
     BOOST_CHECK_EQUAL(read_seq, "TGGGTGGAGCGCGCCCCAT--------TGAGTAGGGTCGACCAAGAACCGCAAGA");
}

BOOST_AUTO_TEST_CASE(plus_strand_softclip)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(4)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(3)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "SSSSSSSSSSSSSSSSSSSSSSSSCGGATATAAACGCCAGGTTGAATCCGCATTT");
     BOOST_CHECK_EQUAL(read_seq, "CGGCTGTCAGCTGCCGTCTGCGTACGGATATAAACGCCAGGTTGAATCCGCATTT");
}

BOOST_AUTO_TEST_CASE(minus_strand_perfect_match)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(5)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(5)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_mismatch)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(6)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(6)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGCGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_insert_in_read)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(7)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(7)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTC------------GACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCCAGTCAGTCAGTGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_deletion_in_read)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(8)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(8)->seq);
    auto [graph_seq, read_seq, mppg_sizes] = ret_tuple;
    BOOST_CHECK_EQUAL(graph_seq, "TCTTGCGGTTCTTGGTCGACCCTACTCACGGTATAAATGGGGCGCGCTCCAT");
     BOOST_CHECK_EQUAL(read_seq, "TCTTGCGGTTCTTGGTCGACCCTA---------TAAATGGGGCGCGCTCCAT");
}

BOOST_AUTO_TEST_CASE(minus_strand_softclip)
{
    shared_ptr dta = make_unique<Trailmix_struct>();
    const string cwdProg = getFullPath(getCWD(".")+"bin/");
    dta->fifo_C = (cwdProg + "../test/reconstructInputSeq/test_reads.gam").c_str();
    auto read_info = readGAM(dta);
    const auto path = read_info->at(9)->path;
    Haplocart hc;
    dta->graph.deserialize(cwdProg+"../test/reconstructInputSeq/target_graph.og");
    auto ret_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->at(9)->seq);
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
euka_argvec.emplace_back("4");
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
euka_argvec.emplace_back("4");
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
euka_argvec.emplace_back("4");
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
euka_argvec.emplace_back("4");
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
euka_argvec.emplace_back("4");
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

/////////////////////////////////////////// BEGIN TEST TRAILMIX /////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(TrailMix)

void run_k1_HV4b(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-fq1");
trailmix_argvec.emplace_back("../test/input_files/trailmix/HV4b.fq.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("50");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/k1_HV4b");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");
trailmix_argvec.emplace_back("--tm-files");
trailmix_argvec.emplace_back("../share/vgan/smaller_tmfiles/");
trailmix_argvec.emplace_back("--dbprefix");
trailmix_argvec.emplace_back("pub.graph");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");

                         }


void run_k1_G1a1a4_FASTA(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-f");
trailmix_argvec.emplace_back("../test/input_files/trailmix/G1a1a4.fa.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("50");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/G1a1a4");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");
trailmix_argvec.emplace_back("--tm-files");
trailmix_argvec.emplace_back("../share/vgan/smaller_tmfiles/");
trailmix_argvec.emplace_back("--dbprefix");
trailmix_argvec.emplace_back("pub.graph");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");
                                }

void run_k2_HV4b_S3(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-fq1");
trailmix_argvec.emplace_back("../test/input_files/trailmix/HV4b_S3.fq.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("50");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("2");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/k2_HV4b_S3");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");
trailmix_argvec.emplace_back("--tm-files");
trailmix_argvec.emplace_back("../share/vgan/smaller_tmfiles/");
trailmix_argvec.emplace_back("--dbprefix");
trailmix_argvec.emplace_back("pub.graph");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");
                         }


void run_k2_HV4b_S3_randStart(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-fq1");
trailmix_argvec.emplace_back("../test/input_files/trailmix/HV4b_S3.fq.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("50");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("2");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/k2_HV4b_S3");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");
trailmix_argvec.emplace_back("--tm-files");
trailmix_argvec.emplace_back("../share/vgan/smaller_tmfiles/");
trailmix_argvec.emplace_back("--dbprefix");
trailmix_argvec.emplace_back("pub.graph");
trailmix_argvec.emplace_back("--randStart");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");
                         }

void run_k2_wrong_matrix(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-fq1");
trailmix_argvec.emplace_back("../test/input_files/trailmix/HV4b_S3.fq.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("50");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("2");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/k2_wrong_matrix");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");
trailmix_argvec.emplace_back("--tm-files");
trailmix_argvec.emplace_back("../share/vgan/smaller_tmfiles/");
trailmix_argvec.emplace_back("--dbprefix");
trailmix_argvec.emplace_back("pub.graph");
trailmix_argvec.emplace_back("--deam3p");
trailmix_argvec.emplace_back("../share/vgan/damageProfiles/dhigh3p.prof");
trailmix_argvec.emplace_back("--deam5p");
trailmix_argvec.emplace_back("../share/vgan/damageProfiles/dhigh5p.prof");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");
                         }


void run_k2_full_graph(Trailmix * tm){
vector<string> trailmix_argvec;
trailmix_argvec.emplace_back("vgan");
trailmix_argvec.emplace_back("trailmix");
trailmix_argvec.emplace_back("-fq1");
trailmix_argvec.emplace_back("../test/input_files/trailmix/full_graph.fq.gz");
trailmix_argvec.emplace_back("-t");
trailmix_argvec.emplace_back("80");
trailmix_argvec.emplace_back("-k");
trailmix_argvec.emplace_back("2");
trailmix_argvec.emplace_back("-o");
trailmix_argvec.emplace_back("../test/output_files/trailmix/k2_full_graph");
trailmix_argvec.emplace_back("--iter");
trailmix_argvec.emplace_back("5000");
trailmix_argvec.emplace_back("--burnin");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("--chains");
trailmix_argvec.emplace_back("1");
trailmix_argvec.emplace_back("-z");
trailmix_argvec.emplace_back("tempdir");

char** argvtopass = new char*[trailmix_argvec.size()];
for (int i=0;i<trailmix_argvec.size();i++) {
                   argvtopass[i] = const_cast<char*>(trailmix_argvec[i].c_str());
                                       }

tm->run(trailmix_argvec.size(), argvtopass, getCWD(".")+"bin/");
                         }


BOOST_AUTO_TEST_CASE(k1_HV4b)
{
    Trailmix tm;
    run_k1_HV4b(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/k1_HV4bBranchEstimate.txt");
    BOOST_ASSERT(branchRecords[0].source == "HV4b");
}

BOOST_AUTO_TEST_CASE(k1_fasta)
{
    Trailmix tm;
    run_k1_G1a1a4_FASTA(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/G1a1a4BranchEstimate.txt");
    BOOST_ASSERT(branchRecords[0].source == "G1a1a4");
}

BOOST_AUTO_TEST_CASE(k2_HV4b_S3)
{
    Trailmix tm;
    run_k2_HV4b_S3(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/k2_HV4b_S3BranchEstimate.txt");
    BOOST_ASSERT((branchRecords[0].source == "HV4b" && branchRecords[1].source == "S3") || (branchRecords[0].source == "S3" && branchRecords[1].source == "HV4b"));
}

BOOST_AUTO_TEST_CASE(k2_HV4b_S3_randStart)
{
    Trailmix tm;
    run_k2_HV4b_S3_randStart(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/k2_HV4b_S3BranchEstimate.txt");
    BOOST_ASSERT((branchRecords[0].source == "HV4b" && branchRecords[1].source == "S3") || (branchRecords[0].source == "S3" && branchRecords[1].source == "HV4b"));
}

BOOST_AUTO_TEST_CASE(wrong_matrix)
{
    Trailmix tm;
    run_k2_wrong_matrix(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/k2_wrong_matrixBranchEstimate.txt");
    BOOST_ASSERT((branchRecords[0].source == "HV4b" && branchRecords[1].source == "S3") || (branchRecords[0].source == "S3" && branchRecords[1].source == "HV4b"));
}

BOOST_AUTO_TEST_CASE(k2_full_graph)
{
    Trailmix tm;
    run_k2_full_graph(&tm);
    auto branchRecords = load_branch_placement_diagnostics_file("../test/output_files/trailmix/k2_full_graphBranchEstimate.txt");
    BOOST_ASSERT((branchRecords[0].source == "HV4b" && branchRecords[1].source == "S3") || (branchRecords[0].source == "S3" && branchRecords[1].source == "HV4b"));
}


BOOST_AUTO_TEST_SUITE_END()

/////////////////////////////////////////// END TEST TRAILMIX /////////////////////////////////////////////////
