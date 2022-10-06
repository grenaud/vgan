#include "HaploCart.h"
#include <algorithm>
#include <thread>
#include <chrono>
#include <omp.h>
#include <numeric>
#include <sys/wait.h>
#include "crash.hpp"
#include "preflight.hpp"
#include "Dup_Remover.h"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "miscfunc.h"

using namespace google::protobuf;
using namespace vg;

Haplocart::Haplocart(){

}

Haplocart::~Haplocart(){

}

const string Haplocart::usage() const{

    return string(string("") + "vgan haplocart [options]"+
                  "\n\nHaplocart predicts the mitochondrial haplogroup for reads originating from an uncontaminated modern human sample."+
		  "\n\nExamples:\n"+
		  "\n\tvgan haplocart --hc-files /home/username/share/hcfiles/ -f myfasta.fa\n"+
		  "\n\tvgan haplocart --hc-files /home/username/share/hcfiles/ -fq1 seqreads_fwd.fq.gz -fq2 seqreads_rev.fq.gz\n"+
                  "\n\n"+
                  "Options:\n\n"+
                  "  Algorithm parameters\n"+
                  "  \t"+"-e [FLOAT]" + "\t\t\t"+ "Background error probability for FASTA input (default 0.0001)\n" +
                  "  Input/Output\n"+
                  "  \t"+"--hc-files [STR]" + "\t\t" + "HaploCart graph directory location (default: \"../share/hcfiles/\")\n" +
                  "  \t"+"-f [STR]" + "\t\t\t" + "FASTA consensus input file\n" +
                  "  \t"+"-fq1 [STR]" + "\t\t\t" + "FASTQ input file\n" +
                  "  \t"+"-fq2 [STR]" + "\t\t\t" + "FASTQ second input file (for paired-end reads)\n" +
                  "  \t"+"-g [STR]" + "\t\t\t" + "GAM input file\n" +
                  "  \t"+"-i " + "\t\t\t\t" + "Input FASTQ (-fq1) is interleaved\n" +
                  "  \t"+"-jf " + "\t\t\t\t" + "JSON output file (must be used with -j)\n" +
                  "  \t"+"-o [STR]" + "\t\t\t" + "Output file (default: stdout)\n" +
                  "  \t"+"-s [STR]" + "\t\t\t" + "Sample name\n" +
                  "  \t" + "-pf" + "\t\t\t\t" + "[STR] Posterior output file (default: stdout)\n" +
                  "  \t"+"-z"+ "\t\t\t\t" + "Temporary directory (default: /tmp/)\n" +
                  "  Non-algorithm parameters\n"+
                  "  \t"+"-j" + "\t\t\t\t"+ "Output JSON file of alignments \n" +
                  "  \t"+"-np" + "\t\t\t\t"+ "Do not compute clade-level posterior probabilities \n" +
                  "  \t"+"-t" + "\t\t\t\t"+ "Number of threads (-t -1 for all available)\n" +
                  "  \t"+"-q" + "\t\t\t\t"+ "Quiet mode\n"
		  );

}

const int Haplocart::run(int argc, char *argv[], const string &cwdProg){

    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    temp_file::set_system_dir();

    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        exit(1);
                                      }

    bool debug=false;
    bool quiet = false;
    bool interleaved = false;
    bool dump_json = false;
    bool webapp = false;
    bool compute_posteriors = true;
    bool rmdup = true;
    bool invoked_samplename = false;
    string posteriorfilename = "/dev/stdout";
    bool   hcfiledirspecified = false;
    string hcfiledir = "../share/hcfiles/";
    double background_error_prob = 0.0001;
    string gamfilename, fastafilename, fastq1filename, fastq2filename, samplename, jsonfilename;
    string outputfilename = "/dev/stdout";
    string tmpdir = "/tmp/";
    int n_threads = 1;

    for(int i=1;i<(argc);++i){

    if(string(argv[i]) == "-") {
        cerr << Haplocart::usage() << endl;
        exit(0);
                               }

    if(string(argv[i]) == "-d"){
            debug=true;
            continue;
                               }

    if(string(argv[i]) == "--hc-files"){
            hcfiledir=argv[i+1];
	    hcfiledirspecified=true;
            if (hcfiledir.back() != '/'){hcfiledir += '/';} 
            continue;
                               }


    if(string(argv[i]) == "-e") {
            background_error_prob = stod(argv[i+1]);
            if (background_error_prob < 0 || background_error_prob > 1) {
                       throw std::runtime_error("[HaploCart] Error, option -e is not a valid probability.");
                                                                        }
            continue;
                                }

    if(string(argv[i]) == "-g"){
            gamfilename = argv[i+1];
            samplename = gamfilename;
            continue;
                                }

    if(string(argv[i]) == "-f"){
            fastafilename = argv[i+1];
            samplename=fastafilename;
            continue;
                                }


    if(string(argv[i]) == "-fq1"){
            fastq1filename = argv[i+1];
            const int idx = fastq1filename.find_last_of("/");
            samplename=fastq1filename.substr(idx + 1);
            samplename=fastq1filename;
            continue;
                                }

    if(string(argv[i]) == "-fq2"){
            fastq2filename = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "-i"){
            interleaved = true;
            continue;
                               }

    if(string(argv[i]) == "-j"){
            dump_json = true;
            continue;
                               }

    if(string(argv[i]) == "-jf"){
            jsonfilename = argv[i+1];
            continue;
                                }

    if(string(argv[i]) == "-o"){
            outputfilename = argv[i+1];
            continue;
                               }

    if(string(argv[i]) == "-np"){
            compute_posteriors=false;
            continue;
                                }


    if(string(argv[i]) == "-pf"){
            posteriorfilename = argv[i+1];
            continue;
                                }

    if(string(argv[i]) == "-q"){
            quiet = true;
            continue;
                               }

    if(string(argv[i]) == "-s"){
            samplename = argv[i+1];
            invoked_samplename = true;
            continue;
                                 }

    if(string(argv[i]) == "-t"){
            if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw std::runtime_error("[HaploCart] Error, invalid number of threads");}
            if (stoi(argv[i+1]) == -1) {n_threads = std::thread::hardware_concurrency();}
            else if (stoi(argv[i+1]) <= std::thread::hardware_concurrency()) {
                n_threads = stoi(argv[i+1]);
                                                                             }
            else {
                cerr << "[HaploCart] Warning, specified number of threads is greater than the number available. Using " << n_threads << " threads\n";
                n_threads = std::thread::hardware_concurrency();
                 }
            continue;
                                 }

    if(string(argv[i]) == "-w"){
            webapp = true;
            std::cerr.rdbuf(NULL);
            n_threads = 62;
            continue;
                               }

    if(string(argv[i]) == "-z"){
            tmpdir = argv[i+1];
            if (tmpdir.back() != '/') {tmpdir += '/';}
            continue;
                               }

                              }

    ofstream outputFile(outputfilename, ios::app);

    if (!webapp && !quiet)
     cerr <<
     "╦ ╦┌─┐┌─┐┬  ┌─┐╔═╗┌─┐┬─┐┌┬┐\n"
     "╠═╣├─┤├─┘│  │ │║  ├─┤├┬┘ │\n"
     "╩ ╩┴ ┴┴  ┴─┘└─┘╚═╝┴ ┴┴└─ ┴\n"
     << "\n" << endl;


    if (!fs::is_directory(tmpdir) || !fs::exists(tmpdir)){
        std::filesystem::create_directory(tmpdir);
                                                         }

    const int idx_ = samplename.find_last_of("/");
    samplename=samplename.substr(idx_ + 1);


    // Handle erroneous input
    if (fastq1filename == "" && fastq2filename != "") {throw std::runtime_error("[HaploCart] Error, cannot invoke -fq2 without -fq1");}
    if (dump_json && jsonfilename == "") {throw(std::runtime_error("[HaploCart] Error, cannot invoke -j without -jf"));}
    if (!compute_posteriors && posteriorfilename != "/dev/stdout"){throw std::runtime_error("[HaploCart] Error, cannot invoke -pf without -p");}
    if (fastafilename.ends_with(".fq") && (!quiet)){cerr << "[HaploCart] Warning, input file is named like a FASTQ but the FASTA flag (-f) was invoked" << endl;}

    if(!hcfiledirspecified){
	hcfiledir  =  cwdProg +  hcfiledir;
    }
    
    // Check that input files exist
    if (!(filesystem::exists(hcfiledir)))
        {throw std::runtime_error("[HaploCart] Error, hc file directory (tried "+hcfiledir+" ) not found. This may be because you have downloaded the requisite files. \n Please see README.md for further instructions. ");}
    
    if (gamfilename != "") {if(!(filesystem::exists(gamfilename)))
                      {throw std::runtime_error("[HaploCart] Error, GAM input file " + gamfilename + " does not exist");}}
    if (fastafilename != "") {if(!(filesystem::exists(fastafilename)))
                      {throw std::runtime_error("[HaploCart] Error, consensus FASTA input file " + fastafilename + " does not exist");}}
    if (fastq1filename != "") {if(!(filesystem::exists(fastq1filename)))
                      {throw std::runtime_error("[HaploCart] Error, FASTQ1 input file " + fastq1filename + " does not exist");}}
    if (fastq2filename != "") {if(!(filesystem::exists(fastq2filename)))
                      {throw std::runtime_error("[HaploCart] Error, FASTQ2 input file " + fastq2filename + " does not exist");}}


    int n_input_files = 0;
    if (gamfilename != "") {n_input_files+=1;}
    if (fastafilename != "") {n_input_files+=1;}
    if (fastq1filename != "") {n_input_files+=1;}
    if (fastq2filename != "") {n_input_files+=1;}
    if (n_input_files == 0) {cerr << "[HaploCart] Error, no input file given" << '\n'; return 1;}
    if (n_input_files > 1 && !(fastq1filename != "" && fastq2filename != "" && gamfilename == "" && fastafilename == "")) {
        cerr << "[HaploCart] Error, cannot accept multiple input files." << '\n'; return 1;
                                                                                                                          }
    if (quiet == false) {cerr << "Predicting sample: " << samplename << '\n';}
    if (quiet == false) {cerr << "Using " << n_threads << " threads" << '\n';}
    const string & graphfilename = getFullPath(hcfiledir + "graph.og");

    int n_samples = 1;
    vector<string> fasta_seqs{""};
    vector<string> fasta_ids{""};
    if (fastafilename != "") {
        std::tie(fasta_seqs, fasta_ids) = Haplocart::read_fasta(fastafilename);
        if (fasta_seqs.size() > 1 && invoked_samplename) {throw std::runtime_error("[HaploCart] Error, cannot invoke -s on multifasta input.");}
        n_samples = fasta_seqs.size();
        if (dump_json && n_samples > 1) {throw std::runtime_error("[HaploCart] Error, cannot invoke -j for multifasta input");}
        if (!quiet) cerr << "Found " << n_samples << " sequences in FASTA file\n";
                             }


    // Load a bunch of stuff
    const map<const string, int> pangenome_map = load_pangenome_map(hcfiledir);
    const tuple<vector<NodeInfo *>, const int, bdsg::ODGI> pathhandlegraphtuple = Haplocart::readPathHandleGraph(graphfilename, min(n_threads, 7), hcfiledir);
    const auto [nodevector, minid, graph] = pathhandlegraphtuple;
    const vector<double> incorrect_mapping_vec = Haplocart::precompute_incorrect_mapping_probs();
    const vector<string> path_names = Haplocart::load_paths(hcfiledir);
    const int nbpaths = path_names.size();
    const vector<double> mappabilities = load_mappabilities(hcfiledir);
    const vector<double> qscore_vec = get_qscore_vec();

    const vg::subcommand::Subcommand* sc = NULL;
    int sample;
    pid_t wpid, wpid2, wpid3;
    int status = 0;
    int status2 = 0;
    int status3 = 0;
    pid_t pid1 = 0;
    pid_t pid2 = 0;
    pid_t pid3 = 0;
    vector<long double> log_likelihood_vec(5179, 0); // Vector of log likelihood
    vector<long double> empty_vec(5179, 0); // Vector of log likelihood
    vector<long double> final_vec(5179, 0); // Vector of log likelihood
    vector<AlignmentInfo *> * algnvector;
    #pragma omp wait_policy dynamic
    #pragma omp parallel for num_threads(1) private(log_likelihood_vec, sample, pid1, pid2, pid3) schedule(static) shared(sc)
    for (sample=0; sample<n_samples; ++sample) {
    if(!quiet){cerr << "Processing sample " << sample+1 << " of " << n_samples << endl;}
    if (fasta_ids.size() > 0 && fastafilename != "" && invoked_samplename == false) {samplename = fasta_ids[sample];}

   if (webapp && (sample > 0)){outputFile << "--------------------------------------------------------------------------------------------------------------------------" << endl;}

    /////////////// GIRAFFE ///////////////////////////

    // If we are given somethng other than GAM we need to map it first

    string first_fifo = tmpdir+random_string(7);
    string second_fifo = tmpdir+random_string(7);
    string third_fifo = tmpdir+random_string(7);
    const char * fifo_A = first_fifo.c_str();
    const char * fifo_B = second_fifo.c_str();
    const char * fifo_C = third_fifo.c_str();

    mkfifo(fifo_A, 0666);
    mkfifo(fifo_B, 0666);
    mkfifo(fifo_C, 0666);

    pid1 = fork();
    if (pid1 == -1) {
        throw std::runtime_error("Error in fork");
                    }


    if(pid1 == 0) {
        // Child
        if (gamfilename == "") {
            Haplocart::map_giraffe(fasta_seqs[sample], fastq1filename, fastq2filename, n_threads,
                                    interleaved, background_error_prob, samplename, fifo_A, sc, tmpdir, hcfiledir, quiet);
           while ((wpid = wait(&status)) > 0);
           exit(0);
                               }

        else {
               // Redirect buffer in case of GAM input
               ifstream src(gamfilename);
               ofstream dst(fifo_A);
               dst << src.rdbuf();
               exit(0);
             }

                  }

    else {

         ////////////////////////// FILTER ////////////////////

           pid2 = fork();
           if (pid2 == -1) {
                throw std::runtime_error("Error in fork");
                           }

           if (pid2 == 0) {
               Haplocart::filter(n_threads, interleaved, fifo_A, fifo_B);
               while ((wpid2 = wait(&status2)) > 0);
               exit(0);
                          }

           else {

        ///////////////////////// GAMSORT ////////////////////

             pid3 = fork();
             if (pid3 == -1) {
                 throw std::runtime_error("Error in fork");
                             }

             if (pid3 == 0) {
             Haplocart::gamsort(n_threads, interleaved, fifo_B, fifo_C, tmpdir);
             while ((wpid3 = wait(&status3)) > 0);
             exit(0);
                            }

        else {

        //////////////////////// NOW BACK TO YOUR REGULARLY SCHEDULED PROGRAMMING ////////////////////////

        algnvector = readGAM(fifo_C, quiet, dump_json, jsonfilename, fastafilename);


        int n_reads = algnvector->size();
        if (n_reads == 0)  {throw std::runtime_error("[HaploCart] Error, no reads mapped");}

        if (!quiet && fastafilename=="") {cerr << "Removing PCR duplicates ..." << '\n';}
        const auto dup_pair = Dup_Remover().remove_duplicates_internal(algnvector, n_threads, quiet);
        algnvector = dup_pair.first;
        n_reads = algnvector->size();

        bool use_background_error_prob = false;
        bool is_consensus_fasta = false;
        if (fastafilename != ""){
            if (quiet == false) {cerr << "Using background error probability of " << background_error_prob << '\n';}
                use_background_error_prob = true;
                is_consensus_fasta = true;
                                }

        const map<string, vector<string>> parents = load_parents(hcfiledir);
        const map<string, vector<string>> children = load_children(hcfiledir);

        if (!quiet && fastafilename=="") {cerr << "Computing haplogroup likelihoods from " << n_reads << " reads." << '\n';}

        int i;
        #pragma omp parallel for num_threads(n_threads) private(log_likelihood_vec, i)
        for(long unsigned int i=0;i!=n_reads;++i){
            // Discard unmapped reads
            if (algnvector->at(i) -> identity < 1e-10) {continue;}

            // Loop through mapped reads and update log likelihood vector accordingly
            log_likelihood_vec = update(i, pangenome_map, nodevector, algnvector->at(i), empty_vec, qscore_vec, mappabilities,
                                        nbpaths, quiet, use_background_error_prob, background_error_prob, incorrect_mapping_vec,
                                        n_reads, minid, is_consensus_fasta, n_threads, graph
                                       );

        #pragma omp critical
        for (int j=0;j<nbpaths;++j) {final_vec[j] += (log_likelihood_vec[j]);}
                                                 }

        const int maxh_index = std::max_element(final_vec.begin(),final_vec.end()) - final_vec.begin();
        const string predicted_haplotype = path_names[maxh_index];

        // Write output to file
        if (fastafilename != ""){n_reads = 1;}

        ofstream outputFile(outputfilename, ios::app);
        replace(samplename.begin(), samplename.end(), ' ', '_');

        if (webapp == false)            {
            if ((outputfilename == "/dev/stdout" && !quiet))             {
                 outputFile <<"\n\n";
                                                                         }

                if (sample==0 || quiet) {outputFile << "#sample\tpredicted haplogroup\treads" << endl;}
                outputFile << samplename << '\t' << predicted_haplotype << '\t' << n_reads << endl;
                                        }

        else {
            outputFile <<"\n\n" << endl;
            outputFile <<"<table>" << endl;
            if (fastafilename != ""){
                outputFile << "<tr><td>" << "#sample" << "</td><td>&emsp;" << "Haplogroup" << endl;
                                    }
            else {
                outputFile << "<tr><td>" << "Sample Number" << "</td><td>&emsp;" << "Haplogroup" << endl;
                 }
            if (fastafilename != "") {
                outputFile <<"<tr><td>" << samplename << "</td><td>&emsp;" << "<strong>"+predicted_haplotype+"</strong>" << endl;
                                     }
            else {
                outputFile <<"<tr><td>" << sample+1 << "</td><td>&emsp;" << "<strong>"+predicted_haplotype+"</strong>" << endl;
                 }
            outputFile <<"</table><br><br>" << endl;
             }

        if (compute_posteriors) {
            Haplocart::get_posterior(final_vec, path_names, parents, children, samplename, predicted_haplotype, \
                                                 posteriorfilename, webapp, sample);
                                }

        if (debug){
            cerr << "Writing log likelihoods to disk" << endl;
            std::vector<int> indices(final_vec.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
            [&](int A, int B) -> bool {
                return final_vec[A] > final_vec[B];
            });
            cerr << setprecision(10) << endl;

	    //@josh is this for you only?
	    ogzstream debug_out((cwdProg + "../data/Debug_log_likelihoods.txt.gz").c_str());
            for (int k=0; k!=final_vec.size(); ++k)       {
                debug_out << path_names[indices[k]] << '\t' << final_vec[indices[k]] << '\n';
                                                          }
                  }
           algnvector->clear();
           fill(begin(final_vec), end(final_vec), 0.0);
           fill(begin(empty_vec), end(empty_vec), 0.0);
                         }
                        }
                     }
                  }
            return 0;
}
