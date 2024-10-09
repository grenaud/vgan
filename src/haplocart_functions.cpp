#include "HaploCart.h"
#define probably_true(x) __builtin_expect(!!(x), 1)

using namespace vg;

int main_filter(int argc, char** argv);
int main_gamsort(int argc, char** argv);

bool Haplocart::contains_no_inf(const std::vector<double>& v) {
    return std::all_of(v.begin(), v.end(), [](double x) {
        return x != std::numeric_limits<double>::infinity() &&
               x != -std::numeric_limits<double>::infinity();
    });
}

void Haplocart::filter(const int n_threads, const bool interleaved, char const * fifo_A, char const * fifo_B){

    auto normal_cout = cout.rdbuf();
    auto normal_cin = cin.rdbuf();

    ifstream cin(fifo_A);
    ofstream cout(fifo_B);
    std::cin.rdbuf(cin.rdbuf());
    std::cout.rdbuf(cout.rdbuf());

    // Filter for only mapped reads and run gamsort
    vector<string> arguments{"vg", "filter"};
    arguments.push_back("-r");
    arguments.push_back("1");
    arguments.push_back("-o");
    arguments.push_back("1");
    arguments.push_back("-t");
    arguments.push_back(to_string(n_threads));
    if (interleaved) {
        arguments.push_back("-i");
                     }
    arguments.push_back("-");

    char** argvtopass = new char*[arguments.size()];
    for (int i=0;i<arguments.size();i++) {
        argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                         }

    main_filter(arguments.size(), argvtopass);

    std::cout.rdbuf(normal_cout);
    std::cin.rdbuf(normal_cin);

    remove(fifo_A);
    delete[] argvtopass;
}

void Haplocart::gamsort(const int n_threads, const bool interleaved, char const * fifo_B, char const * fifo_C, const string & tmpdir){

    auto normal_cout = cout.rdbuf();
    auto normal_cin = cin.rdbuf();

    ifstream cin(fifo_B);
    ofstream cout(fifo_C);

    std::cin.rdbuf(cin.rdbuf());
    std::cout.rdbuf(cout.rdbuf());

    vg::temp_file::set_dir(tmpdir);

    // Filter for only mapped reads and run gamsort
    vector<string> arguments{"vg", "gamsort", "-t"};
    arguments.push_back(to_string(n_threads));
    arguments.push_back("-");

    char** argvtopass = new char*[arguments.size()];
    for (size_t i=0;i<arguments.size();i++) {
        argvtopass[i] = const_cast<char*>(arguments[i].c_str());
                                         }
    main_gamsort(arguments.size(), argvtopass);

    std::cout.rdbuf(normal_cout);
    std::cin.rdbuf(normal_cin);

    cout.close();
    cin.close();

    remove(fifo_B);
    delete[] argvtopass;
}




const double Haplocart::get_background_freq(const char &base)
     {

    //Background proportions of nucleotides in the hominin mitogenome (taken from frequencies within the graph)

    switch(base)
    {
    case 'A':
        return 0.27532;
    case 'C':
        return 0.30044;
    case 'G':
        return 0.16644;
    case 'T':
        return 0.25780;
    }
    return 0.25;
     }


constexpr double Haplocart::get_p_incorrectly_mapped(const int &Q)
{return pow(10, ((-1 * Q)*0.1));}


void Haplocart::precompute_incorrect_mapping_probs(shared_ptr<Trailmix_struct> &dta){
for (int Q=0; Q!=100; ++Q) {
    dta->incorrect_mapping_vec.emplace_back(Haplocart::get_p_incorrectly_mapped(Q));
                           }
return;
}




