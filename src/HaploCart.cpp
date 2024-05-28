#include "HaploCart.h"
#include <algorithm>
#include <thread>
#include <chrono>
#include <omp.h>
#include <numeric>
#include <sys/wait.h>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "Dup_Remover.h"

using namespace google::protobuf;
using namespace vg;

Haplocart::Haplocart(){

}

Haplocart::~Haplocart(){

}

const string Haplocart::usage() {

    return string(string("") + "vgan Haplocart [options]"+
                  "\n\nHaplocart predicts the mitochondrial haplogroup for reads originating from an uncontaminated modern human sample."+
                  "\n\n"+
                  "Options:\n\n"+
                  "  Algorithm parameters\n"+
                  "  \t"+"-e [FLOAT]" + "\t\t\t"+ "Background error probability for FASTA input (default 0.0001)\n" +
                  "  Input/Output\n"+
                  "  \t"+"-f [STR]" + "\t\t\t" + "FASTA consensus input file\n" +
                  "  \t"+"-fq1 [STR]" + "\t\t\t" + "FASTQ input file\n" +
                  "  \t"+"-fq2 [STR]" + "\t\t\t" + "FASTQ second input file (for paired-end reads)\n" +
                  "  \t"+"-g [STR]" + "\t\t\t" + "GAM input file\n" +
                  "  \t"+"-i " + "\t\t\t\t" + "Input FASTQ (-fq1) is interleaved\n" +
                  "  \t"+"-o [STR]" + "\t\t\t" + "Output file (default: stdout)\n" +
                  "  \t"+"-s [STR]" + "\t\t\t" + "Sample name\n" +
                  "  \t" + "-pf" + "\t\t\t\t" + "[STR] Posterior output file (default: stdout)\n" +
                  "  \t"+"-z"+ "\t\t\t\t" + "Temporary directory (default: /tmp/)\n"
                  "  Non-algorithm parameters\n"+
                  "  \t"+"-t" + "\t\t\t\t"+ "Number of threads (-t -1 for all available)\n" +
                  "  \t"+"-q" + "\t\t\t\t"+ "quiet mode\n"
                  "  \t"+"-p" + "\t\t\t\t" + "Compute posterior probabilities for the prediction and its ancestral clades\n"
		  );

}

void Haplocart::run(int argc, char *argv[], shared_ptr<Trailmix_struct> &dta){

/*
    std::ios_base::sync_with_stdio(false);
    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    temp_file::set_system_dir();

    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        exit(1);
                                      }
*/

    dta->background_error_prob=0.0001;
    int i=0;
    int n_input_files=0;
    int idx_ = 0;
    vector<string> fasta_seqs{""};
    vector<string> fasta_ids{""};
    auto start = std::chrono::system_clock::now();

      // Load a bunch of stuff
        if (!dta->running_trailmix) {
        Haplocart::load_pangenome_map(dta);
        Haplocart::load_path_names(dta);
        Haplocart::readPathHandleGraph(dta);
        Haplocart::precompute_incorrect_mapping_probs(dta);
        dta->nbpaths = dta->path_names.size();
        load_mappabilities(dta);
        dta->qscore_vec = get_qscore_vec();
                                    }

    vector<double> log_likelihood_vec;
    log_likelihood_vec.resize(dta->nbpaths, 0.0);

    if(dta->reads_already_processed){goto infer;}

    for(unsigned int i=1;i<(argc);++i){

    if(string(argv[i]) == "-") {
        cerr << Haplocart::usage() << endl;
        exit(0);
                               }

    if(string(argv[i]) == "-e") {
            dta->background_error_prob = stod(argv[i+1]);
            if (dta->background_error_prob < 0 || dta->background_error_prob > 1) {
                       throw std::runtime_error("[HaploCart] Error, option -e is not a valid probability.");
                                                                        }
            continue;
                                }

    if(string(argv[i]) == "-g"){
            dta->gamfilename = argv[i+1];
            dta->samplename = dta->gamfilename;
            continue;
                                }

    if(string(argv[i]) == "-f"){
            dta->fastafilename = argv[i+1];
            dta->samplename=dta->fastafilename;
            continue;
                                }


    if(string(argv[i]) == "-fq1"){
            dta->fastq1filename = argv[i+1];
            const int idx = dta->fastq1filename.find_last_of("/");
            dta->samplename=dta->fastq1filename.substr(idx + 1);
            dta->samplename=dta->fastq1filename;
            continue;
                                }

    if(string(argv[i]) == "-fq2"){
            dta->fastq2filename = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "--hc-files"){
            dta->graph_prefix=argv[i+1];
	    dta->graph_dir_specified=true;
            if (dta->graph_prefix.back() != '/'){dta->graph_prefix += '/';}
            continue;
                                       }

    if(string(argv[i]) == "-i"){
            dta->interleaved = true;
            continue;
                                 }

    if(string(argv[i]) == "-o"){
            dta->outputfilename = argv[i+1];
            continue;
                                 }

    if(string(argv[i]) == "-p"){
            dta->compute_posteriors = true;
            continue;
                               }

    if(string(argv[i]) == "-pf"){
            dta->posteriorfilename = argv[i+1];
            continue;
                                }

    if(string(argv[i]) == "-q"){
            dta->quiet = true;
            continue;
                               }

    if(string(argv[i]) == "-s"){
            dta->samplename = argv[i+1];
            dta->invoked_samplename = true;
            continue;
                                 }

    if(string(argv[i]) == "-t"){
            if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw std::runtime_error("[HaploCart] Error, invalid number of threads");}
            if (stoi(argv[i+1]) == -1) {dta->n_threads = std::thread::hardware_concurrency();}
            else if (stoi(argv[i+1]) <= std::thread::hardware_concurrency()) {
                dta->n_threads = stoi(argv[i+1]);
                                                                             }
            else {
                cerr << "[HaploCart] Warning, specified number of threads is greater than the number available. Using " << dta->n_threads << " threads\n";
                dta->n_threads = std::thread::hardware_concurrency();
                 }
            continue;
                                 }

    if(string(argv[i]) == "-w"){
            dta->webapp = true;
            std::cerr.rdbuf(NULL);
            continue;
                               }

    if(string(argv[i]) == "-z"){
            dta->tmpdir = argv[i+1];
            if (dta->tmpdir.back() != '/') {dta->tmpdir += '/';}
            continue;
                               }

                              }


    if (!dta->webapp && !dta->quiet && !dta->running_trailmix)
     cerr <<
     "╦ ╦┌─┐┌─┐┬  ┌─┐╔═╗┌─┐┬─┐┌┬┐\n"
     "╠═╣├─┤├─┘│  │ │║  ├─┤├┬┘ │\n"
     "╩ ╩┴ ┴┴  ┴─┘└─┘╚═╝┴ ┴┴└─ ┴\n"
     << "\n\n\n" << endl;


    if (!fs::is_directory(dta->tmpdir) || !fs::exists(dta->tmpdir)){
        std::filesystem::create_directory(dta->tmpdir);
                                                         }

    idx_ = dta->samplename.find_last_of("/");
    dta->samplename=dta->samplename.substr(idx_ + 1);

    // Handle erroneous input
    //if (dta->fastq1filename == "" && dta->fastq2filename != "") {cerr << "[HaploCart] Error, cannot invoke -fq2 without -fq1" << '\n'; return 1;}

    // Check that input files exist

    if (!(filesystem::exists(dta->graph_dir))) {
        throw std::runtime_error("[HaploCart] Error, hc file directory (tried "+dta->graph_dir+" ) not found. This may be because you have downloaded the requisite files. \n Please see README.md for further instructions. ");
                                               }

    if (dta->gamfilename != "") {if(!(filesystem::exists(dta->gamfilename)))
                      {throw std::runtime_error("[HaploCart] Error, GAM input file " + dta->gamfilename + " does not exist");}}
    if (dta->fastafilename != "") {if(!(filesystem::exists(dta->fastafilename)))
                      {throw std::runtime_error("[HaploCart] Error, consensus FASTA input file " + dta->fastafilename + " does not exist");}}
    if (dta->fastq1filename != "") {if(!(filesystem::exists(dta->fastq1filename)))
                      {throw std::runtime_error("[HaploCart] Error, FASTQ1 input file " + dta->fastq1filename + " does not exist");}}
    if (dta->fastq2filename != "") {if(!(filesystem::exists(dta->fastq2filename)))
                      {throw std::runtime_error("[HaploCart] Error, FASTQ2 input file " + dta->fastq2filename + " does not exist");}}

    if (dta->gamfilename != "") {n_input_files+=1;}
    if (dta->fastafilename != "") {n_input_files+=1;}
    if (dta->fastq1filename != "") {n_input_files+=1;}
    if (dta->fastq2filename != "") {n_input_files+=1;}
    if (n_input_files == 0) {throw std::runtime_error("[HaploCart] Error, no input file given");}
    if (n_input_files > 1 && !(dta->fastq1filename != "" && dta->fastq2filename != "" && dta->gamfilename == "" && dta->fastafilename == "")) {
        throw std::runtime_error("[HaploCart] Error, cannot accept multiple input files.");
                                                                                                                                              }
    if (dta->quiet == false) {cerr << "Predicting sample: " << dta->samplename << '\n';}
    if (dta->quiet == false) {cerr << "Using " << dta->n_threads << " threads" << '\n';}

    if (dta->fastafilename != "") {
        std::tie(fasta_seqs, fasta_ids) = Haplocart::read_fasta(dta->fastafilename);
        if (fasta_seqs.size() > 1 && dta->invoked_samplename) {throw std::runtime_error("[HaploCart] Error, cannot invoke -s on multifasta input.");}
        dta->n_samples = fasta_seqs.size();
        if (!dta->quiet) cerr << "Found " << dta->n_samples << " sequences in FASTA file\n";
                                   }

    for (i=0; i<dta->n_samples; ++i) {
        if(!dta->quiet){cerr << "Processing sample " << i+1 << " of " << dta->n_samples << endl;}
        if (fasta_ids.size() > 0 && dta->fastafilename != "") {dta->samplename = fasta_ids[i];}


    /////////////// GIRAFFE ///////////////////////////

    // If we are given something other than GAM we need to map it first

    dta->first_fifo = dta->tmpdir+random_string(7);
    dta->second_fifo = dta->tmpdir+random_string(7);
    dta->third_fifo = dta->tmpdir+random_string(7);
    dta->fifo_A = dta->first_fifo.c_str();
    dta->fifo_B = dta->second_fifo.c_str();
    dta->fifo_C = dta->third_fifo.c_str();

    mkfifo(dta->fifo_A, 0666);
    mkfifo(dta->fifo_B, 0666);
    mkfifo(dta->fifo_C, 0666);

    dta->pid1 = fork();
    if (dta->pid1 == -1) {
        throw std::runtime_error("Error in fork");
                    }


    if(dta->pid1 == 0) {
        // Child
        if (dta->gamfilename == "") {
            if (!dta->reads_already_processed) {
                Haplocart::map_giraffe(fasta_seqs[i], dta->fastq1filename, dta->fastq2filename, dta->n_threads,
                                    dta->interleaved, dta->background_error_prob, dta->samplename, dta->fifo_A, dta->sc, dta->tmpdir, \
                                    dta->graph_dir, dta->quiet);
                                               }
           while ((dta->wpid = wait(&dta->status)) > 0);
           exit(0);
                               }

        else {
               // Redirect buffer in case of GAM input
               ifstream src(dta->gamfilename);
               ofstream dst(dta->fifo_A);
               dst << src.rdbuf();
               exit(0);
             }

                  }

    else {

         ////////////////////////// FILTER ////////////////////

           dta->pid2 = fork();
           if (dta->pid2 == -1) {
                throw std::runtime_error("Error in fork");
                           }

           if (dta->pid2 == 0) {
               Haplocart::filter(dta->n_threads, dta->interleaved, dta->fifo_A, dta->fifo_B);
               while ((dta->wpid2 = wait(&dta->status2)) > 0);
               exit(0);
                          }

           else {

        ///////////////////////// GAMSORT ////////////////////

             dta->pid3 = fork();
             if (dta->pid3 == -1) {
                 throw std::runtime_error("Error in fork");
                             }

             if (dta->pid3 == 0) {
             Haplocart::gamsort(dta->n_threads, dta->interleaved, dta->fifo_B, dta->fifo_C, dta->tmpdir);
             while ((dta->wpid3 = wait(&dta->status3)) > 0);
            exit(0);
                            }

           else {

        //////////////////////// NOW BACK TO YOUR REGULARLY SCHEDULED PROGRAMMING ////////////////////////

        //assert(!dta->reads_already_processed);
        dta->algnvector = readGAM(dta);
        assert(!dta->algnvector->empty());
        if(dta->running_trailmix){
            dta->reads_already_processed=true;return;
                                 }

        if (dta->algnvector->size() == 0) {throw std::runtime_error("[HaploCart] Error, no reads mapped");}

        if (dta->quiet == false && dta->fastafilename=="") {cerr << "Removing PCR duplicates ..." << '\n';}

//shared_ptr<vector<bool>> thing = Dup_Remover().remove_duplicates_internal(dta->algnvector, dta->n_threads, dta->quiet);

        if (dta->fastafilename != ""){
            if (dta->quiet == false) {cerr << "Using background error probability of " << dta->background_error_prob << '\n';}
                dta->use_background_error_prob = true;
                dta->is_consensus_fasta = true;
                                     }
        //if (dta->quiet == false) {cerr << "Computing haplogroup likelihoods from " << dta->algnvector->size() << " reads." << '\n';}

infer:

    dta->reads_already_processed=true;

        unsigned int i;
        vector<double> empty_vec(log_likelihood_vec.size(), 0.0);
        vector<double> final_vec(log_likelihood_vec.size(), 0.0);
        #pragma omp parallel for num_threads(dta->n_threads) private(i) private(log_likelihood_vec) schedule(dynamic)
        for(size_t i=0;i!=dta->algnvector->size();++i){
            //cerr << "On read: " << i << endl;
            // Discard unmapped reads
            if (dta->algnvector->at(i) -> identity < 1e-10) {continue;}

            // Loop through mapped reads and update log likelihood vector accordingly
            log_likelihood_vec = update(dta, i, empty_vec);

            #pragma omp critical
            for (size_t j=0;j<log_likelihood_vec.size();++j) {final_vec[j] += log_likelihood_vec[j];}
                log_likelihood_vec.clear();
                                                 }

        assert(contains_no_inf(final_vec));



         if (true){
            cerr << "Writing log likelihoods to disk" << endl;
            std::vector<int> indices(final_vec.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
            [&](int A, int B) -> bool {
                return final_vec[A] > final_vec[B];
            });
            cerr << setprecision(10) << endl;
            ofstream debug_out("debug.txt");
            for (int k=0; k!=final_vec.size(); ++k)       {
                debug_out << dta->path_names[indices[k]] << '\t' << final_vec[indices[k]] << '\n';
                                                          }
                  }


        const unsigned int maxh_index = max_element(final_vec.begin(), final_vec.end()) - final_vec.begin();
        const string proposed_haplotype = dta->path_names[maxh_index];

        // Write output to file

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        ofstream outputFile(dta->outputfilename, ios::app);

        replace(dta->samplename.begin(), dta->samplename.end(), ' ', '_');


        const size_t n_reads = dta->algnvector->size();
        if (!dta->running_trailmix) {
        if (dta->webapp == false) {
            if (dta->outputfilename == "") {
                cout << "\n\nSample\tPredicted haplotype\t# threads\ttime(ms)\t# reads" << endl;
                cout << dta->samplename << '\t' << proposed_haplotype << '\t' << dta->n_threads << '\t' << elapsed.count() << '\t' << n_reads << "\n\n" << endl;
                                      }
            else {
                outputFile << "\n\nSample\tPredicted haplotype\t# threads\ttime(ms)\t# reads" << endl;
                outputFile << dta->samplename << '\t' << proposed_haplotype << '\t' << dta->n_threads << '\t' << elapsed.count() << '\t' << n_reads << "\n\n" << endl;
                 }
                             }

        else {
            cout << "\n\n";
            cout << "<table>" << '\n';
            cout << "<tr><td>" << dta->samplename << "</td><td>" << proposed_haplotype << "</td><td>"  << "</td></tr>" << '\n';
            cout << "</table>" << '\n';
            for (unsigned int i=0;i<8;++i){cout << '\n';}
             }

        if (dta->compute_posteriors) {
            Haplocart::get_posterior(final_vec, dta->path_names, dta->parents, dta->children, dta->samplename, proposed_haplotype, dta->posteriorfilename, dta->webapp);
                                     }

                        }

                                 }
                         }

                        }
                   }
                                 }
