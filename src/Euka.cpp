#include "Euka.h"
#include "HaploCart.h"
#include <sys/wait.h>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "miscfunc.h"
#include "readOG_Euka.h"
#include "MCMC.h"
#include "damage.h"
#include "baseshift.cpp" // Mikkel code

using namespace google::protobuf;
using namespace vg;

#define VERBOSE_S
//#define DEBUGDEAM

Euka::Euka(){

}

Euka::~Euka(){

}

// functions: 

/**
 * get_qscore_vec gives a vector of
 *
 * This function creates a vector of probabilities for each phred score. It is used to look up the probability of a sequenceing error for each base
 *
 * @return vector of probabilities for phred scores 0-70
 *
 */
vector<double> Euka::get_qscore_vec() {

    vector<double> qscore_vec;
    for (int Q=0; Q<100;++Q) {
        if (Q >= 2) {
           qscore_vec.emplace_back(pow(10, ((-1 * Q)*0.1)));
                              }
        else{
           qscore_vec.emplace_back(0.25);
            }

    }
    return qscore_vec;
}

// function gets the average for per clade damage profils (specific for either c->t or g->a as well as 5' or 3' end). 
vector <double> Euka::get_avg(vector<double> dam, int l, int no_clades){

    vector <double> sum(l, 0.0);
    vector<double> av; 
    
    for (int i = 0; i<dam.size(); i++){

        sum[i%l] += dam.at(i); 
    }

    for (int j = 0; j<sum.size(); j++){
        av.emplace_back(sum[j]/no_clades);
    }
    
    return av;
}


const string Euka::usage() const{

    return string(string("")+ 

		   "\n"+
		    "              __        \n"+
		    "  ___  __  __/ /______ _\n"+
		    " / _ \\/ / / / /_/ __ `/\n"+
		    "/  __/ /_/ / ,< / /_/ / \n"+
		    "\\___/\\__,_/_/|_|\\__,_/  \n"+

                  "\neuka performs abundance estimation of eukaryotic taxa from an environmental DNA sample.\n"+
		  "\n"+

          "vgan euka [options] \n"+
          "\n\nNo damage example:\n"+
                  "\n\tvgan euka -fq1 seqreads.fq.gz\n"+
          "\n\nDamage example:\n"+
                  "\n\tvgan euka -fq1 seqreads.fq.gz --deam5p ../share/damageProfiles/dhigh5.prof --deam3p ../share/damageProfiles/dhigh3.prof\n"+
          "\n\nUser specific MCMC example:\n"+
                  "\n\tvgan euka -fq1 seqreads.fq.gz -iter 100000 -burnin 1000"

          

                  "\n\n"+
                  "Input options:\n"+
                  "\t\t"+""  +"" +"--euka_dir [STR]" + "\t" + "euka database location (default: current wroking directory)\n"+  
                  "\t\t"+""  +"" +"--dbprefix [STR]" + "\t" + "database prefix name (defualt: euka_db)\n"+
                  "\t\t"+""  +"" +"-fq1 [STR]"   +"\t\t" + "Input FASTQ file (for merged and single-end reads)"+"\n"+
                  "\t\t"+""  +"" +"-fq2 [STR]"   +"\t\t" + "Second input FASTQ file (for paired-end reads)"+"\n"+
                  "\t\t"+""  +"" +"-i"   +"\t\t\t" + "Paired-end reads are interleaved (default: false)"+"\n"+
                  "\t\t"+""  +"" +"-g [STR]"   +"\t\t" + "GAM input file"+"\n"+ 
                  "\t\t"+""  +"" +"-M [STR]"   +"\t\t" + "Alternative minimizer prefix input (defualt: euka_db)"+"\n"+
                  "\t\t"+""  +"" +"-o [STR]"   +"\t\t" + "Output file prefix (default: euka_output) "+"\n"+
                  "\t\t"+""  +"" +"-t"   +"\t\t\t" + "Number of threads (-1 for all available)"+"\n"+
                  "\t\t"+""  +"" +"-Z"   +"\t\t\t" + "Temporary directory (default: /tmp)"+"\n"+
                  "\n"+
                  "Filter options:\n"+
                  "\t\t"+""  +"" +"--minMQ [INT]"    +"\t\t"      +"Set the mapping quality minimum for a fragment (default: 29)"+"\n"+
                  "\t\t"+""  +"" +"--minFrag [INT]"    +"\t\t"      +"Minimum amount of fragments that need to map to a group (default: 10)"+"\n"+
                  "\t\t"+""  +"" +"--entropy [double]"    +"\t\t"      +"Minimum entropy score for a bin to be considered (default: 1.17)"+"\n"+
                  "\t\t"+""  +"" +"--minBins [INT]"    +"\t\t"      +"Minimum number of bins with an entropy higher than the thresold (default: 6)"+"\n"+
                  "\t\t"+""  +"" +"--maxBins [INT]"    +"\t\t"      +"Maximum number of bins without coverage (default: 0)"+"\n"+
                  "\n"+
                  "Damage options:"+"\n"+
		          "\t\t"+""  +"" +"--deam5p"   +"\t\t"      +"[.prof]"  +"\t"+"5p deamination frequency for eukaryotic species (default: no damage)"+"\n"+
		          "\t\t"+""  +"" +"--deam3p"   +"\t\t"      +"[.prof]"  +"\t"+"3p deamination frequency for eukaryotic species (default: no damage)"+"\n"+
                  "\t\t"+""  +"" +"-l [INT]"   +"\t\t" + "Set length for substitution matrix (default: 5)"+"\n"+ // Mikkels argument
                  "\t\t"+""  +"" +"--out_dir [STR]"   +"\t\t" + "Path for output prof-file (default: current wroking directory)"+"\n"+ // Mikkels argument
                  "\n"+
                  "Markov chain Monte Carlo options:"+"\n"+
                  "\t\t"+""  +"" +"--no-mcmc"     +"\t\t"   +"The MCMC does not run (default: false)"+"\n"+
                  "\t\t"+""  +"" +"--iter [INT]"      +"\t\t"   +"Define the number of iterartions for the MCMC (default: 10000)"+"\n"+
                  "\t\t"+""  +"" +"--burnin [INT]"      +"\t\t"   +"Define the burnin period for the MCMC (default: 100)"+"\n"+
                  "\n"+
                  "Additional output option:\n"+
                  "\t\t"+""  +"" +"--outFrag"     +"\t\t"   +"Outputs a file with all read names per taxonomic group (default: false)"+"\n"
                  "\t\t"+""  +"" +"--outGroup [string]"     +"\t\t"   +"Outputs all information about a taxonmic group of interest (default: empty)"+"\n"

                  //"\t"+"-b" + "\t\t\t"+"Produce binary compressed output       (default: "+booleanAsString(uncompressed)+")\n"+
		  ""
		  );

}



const int Euka::run(int argc, char *argv[], const string cwdProg){

    preflight_check();
    configure_memory_allocator();
    enable_crash_handling();
    temp_file::set_system_dir();

    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        exit(1);
                                      }
    int lastOpt=1;

    string deam5pfreqE;
    string deam3pfreqE;

    ifstream file(cwdProg + "../share/vgan/damageProfiles/none.prof");

    if (file) {
        // Load deamination profiles
        deam5pfreqE = getFullPath(cwdProg + "../share/vgan/damageProfiles/none.prof");
        deam3pfreqE = getFullPath(cwdProg + "../share/vgan/damageProfiles/none.prof");
    } else {
        // Handle error when getFullPath fails
        cerr << "Warning: deamination profile not found. Proceeding without the assumption of ancient damage." << endl;
        deam5pfreqE = ""; // Set to empty to handle it later
        deam3pfreqE = "";
    }

    bool specifiedDeam=false;
    bool interleaved=false;
    bool run_mcmc=true;
    int n_threads = 1;
    int iter=10000;
    int burnin=100;
    string fastq1filename, fastq2filename, gamfilename, samplename;
    string tmpdir = "/tmp/";
    string outputfilename = "euka_output";
    bool euka_dirspecified = false;
    string euka_dir = "../share/vgan/euka_dir/";
    string dbprefixS = "euka_db";
    int lengthToProf = 5;
    string prof_out_file_path = getFullPath(cwdProg+"../");
    unsigned int MINNUMOFBINS = 6; 
    unsigned int MINNUMOFREADS = 10; 
    unsigned int MINIMUMMQ = 29;
    int MAXIMUMOFBINS = 0;
    double ENTROPY_SCORE_THRESHOLD = 1.17;
    bool outFrag = false;
    string outGroup = "";
    string altMin = dbprefixS;
    bool safari = false;


    for(int i=1;i<(argc);i++){

        if(string(argv[i]) == "--euka_dir"){
            euka_dir = argv[i+1];
            euka_dirspecified=true; 
            if(euka_dir.back() != '/'){euka_dir += '/';}
            continue;
        }

        if(string(argv[i]) == "--dbprefix"){
            dbprefixS = argv[i+1];
            continue;
        }

        if(string(argv[i]) == "-fq1"){
            fastq1filename = argv[i+1];
            const int idx = fastq1filename.find_last_of("/");
            samplename=fastq1filename.substr(idx + 1);
            samplename=fastq1filename;
            if (fastq1filename.ends_with(".fa") || fastq1filename.ends_with(".fasta") || fastq1filename.ends_with(".fa.gz") || fastq1filename.ends_with(".fasta.gz"))
                {throw runtime_error("[euka] Input file must be FASTQ, not FASTA");}
            continue;
                                }

        if(string(argv[i]) == "-fq2"){
            fastq2filename = argv[i+1];
            if (fastq2filename.ends_with(".fa") || fastq2filename.ends_with(".fasta") || fastq2filename.ends_with(".fa.gz") || fastq2filename.ends_with(".fasta.gz"))
                {throw runtime_error("[euka] Input file must be FASTQ, not FASTA");}
            continue;
        }

        if(string(argv[i]) == "-i"){
            interleaved = true;
            if (fastq2filename != ""){throw runtime_error("[euka] If interleaved option chosen, Euka expects only one FASTQ file");}
            continue;
                               }

        if(string(argv[i]) == "-g"){
            gamfilename = argv[i+1];
            samplename = gamfilename;
            continue;
                                }
        if(string(argv[i]) == "-M"){
            altMin = argv[i+1];
            continue;
                                }
        if(string(argv[i]) == "--deam5p"  ){
            deam5pfreqE=string(argv[i+1]);
	    specifiedDeam=true;
            continue;
        }

        if(string(argv[i]) == "--deam3p"  ){
            deam3pfreqE=string(argv[i+1]);
	    specifiedDeam=true;
            continue;
        }
        if(string(argv[i]) == "--no-mcmc"  ){
            run_mcmc=false;
            continue;
        }
        if(string(argv[i]) == "--iter"){
            iter=stoi(argv[i+1]);
            assert(iter >= 0);
            continue;
        }
        if(string(argv[i]) == "--burnin"){
            burnin=stoi(argv[i+1]);
            assert(burnin >= 0);
            continue;
        }
        if(string(argv[i]) == "--entropy"){
            ENTROPY_SCORE_THRESHOLD=stod(argv[i+1]);
            assert(ENTROPY_SCORE_THRESHOLD >= 0);
            if (ENTROPY_SCORE_THRESHOLD > 5.0){throw runtime_error("[euka] Error, entropy thresold is too stringent");}
            continue;
        }
        if(string(argv[i]) == "--minBins"){
            MINNUMOFBINS=stoi(argv[i+1]);
            assert(MINNUMOFBINS >= 0);
            if (MINNUMOFBINS > 20){throw runtime_error("[euka] Error, minimum number of bins exceeds the total number of bins");}
            continue;
        }
        if(string(argv[i]) == "--maxBins"){
            MAXIMUMOFBINS=stoi(argv[i+1]);
            assert(MAXIMUMOFBINS >= 0);
            if (MAXIMUMOFBINS > 20){throw runtime_error("[euka] Error, maximum number of bins exceeds the total number of bins");}
            continue;
        }
        if(string(argv[i]) == "--minMQ"){
            MINIMUMMQ=stoi(argv[i+1]);
            assert(MINIMUMMQ >= 0 && MINIMUMMQ <= 60);
            continue;
        }
        if(string(argv[i]) == "--minFrag"){
            MINNUMOFREADS=stoi(argv[i+1]);
            assert(MINNUMOFREADS >= 0);
            continue;
        }
        if(string(argv[i]) == "--outFrag"  ){
            outFrag=true;
            continue;
        }
        if(string(argv[i]) == "--outGroup"){
            outGroup = argv[i+1];
            continue;
        }
        if(string(argv[i]) == "-S"  ){
            safari=true;
            continue;
        }


        if(string(argv[i]) == "-t"){
        if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw std::runtime_error("[euka] Error, invalid number of threads");}
        if (stoi(argv[i+1]) == -1) {n_threads = std::thread::hardware_concurrency();}
        else if (stoi(argv[i+1]) <= std::thread::hardware_concurrency()) {
                n_threads = stoi(argv[i+1]);
                                                                         }
        else {
               cerr << "[euka] Warning, specified number of threads is greater than the number available. Using " << n_threads << " threads\n";
               n_threads = std::thread::hardware_concurrency();
             }
            continue;
                                   }

        if(string(argv[i]) == "-o"){
            outputfilename = argv[i+1];
            continue;
                               }

        if(string(argv[i]) == "-z"){
            tmpdir = argv[i+1];
            if (tmpdir.back() != '/') {tmpdir += '/';}
            continue;
                               }

        // // // Mikkel code begin // // //
        
        if(string(argv[i]) == "-l"){
            lengthToProf = stoi(argv[i+1]);
            continue;
        }
        
        if(string(argv[i]) == "--out_dir"){
            prof_out_file_path = string(argv[i+1]);
            continue;
        }

        // // // Mikkel code ends // // //


    }
        cerr <<
        "              __        \n"
        "  ___  __  __/ /______ _\n"
        " / _ \\/ / / / /_/ __ `/\n"
        "/  __/ /_/ / ,< / /_/ / \n"
        "\\___/\\__,_/_/|_|\\__,_/  \n"
        << '\n' << endl;



    bool populateAlignmentVector=false;

    if(!euka_dirspecified){
        euka_dir = cwdProg + euka_dir;
    }


    string dbprefix              = euka_dir + dbprefixS;
    string altMinFull            = euka_dir + altMin;


    string ogfilename     = dbprefix+".og";
    string vgfilename     = dbprefix+".vg";
    string cladefilename  = dbprefix+".clade";
    string binsfilename   = dbprefix+".bins";
    //string pathsupportfile = dbprefix+"_graph_path_supports";
    string gbtwfilename = dbprefix+".gbwt";

    // Checking if all mandatory graph input files exist
    if (isFile(ogfilename) == false){
        throw(std::runtime_error(ogfilename + " does not exist."));

    } else if (isFile(cladefilename) == false){
        throw(std::runtime_error(cladefilename + " does not exist."));
    } else if (isFile(binsfilename) == false){
        throw(std::runtime_error(binsfilename + " does not exist."));
    } else if (isFile(gbtwfilename) == false){
    	throw(std::runtime_error(gbtwfilename + " does not exist."));
     } 


    // Check that user input files exist
    if (gamfilename != "") {if(!(filesystem::exists(gamfilename)))
                      {throw std::runtime_error("[euka] Error, GAM input file " + gamfilename + " does not exist");}}
    if (fastq1filename != "") {if(!(filesystem::exists(fastq1filename)))
                      {throw std::runtime_error("[euka] Error, FASTQ1 input file " + fastq1filename + " does not exist");}}
    if (fastq2filename != "") {if(!(filesystem::exists(fastq2filename)))
                      {throw std::runtime_error("[euka] Error, FASTQ2 input file " + fastq2filename + " does not exist");}}


    Damage dmg;
    dmg.initDeamProbabilities(deam5pfreqE,deam3pfreqE);


    /////// WRITE PATH SUPPORTS FROM VG INPUT (ONLY NEED TO DO THIS ONCE) /////////
    //auto readvgtuple = readVG(vgfilename);
    //vector<NodeInfo *> * nodevector=std::get<0>(readvgtuple);
    //write_path_supports("euka_db_graph_path_supports", nodevector);
    /////////// DONE WRITING PATH SUPPORTS ///////////////////////////////////////

    //read information about different clades
    cerr << "Reading in taxa information ..." << endl;
    cerr << "\t-------------------------------" << endl;
    vector<Clade *> * clade_vec = load_clade_info(cladefilename, lengthToProf);
    cerr << "... done!" << endl;
    if (clade_vec->empty()) {
    throw std::runtime_error("Error: The clade vector is empty. Unable to proceed. Check if the soibean.clade file is not empty.");
    }
        cerr << "-------------------------------" << endl;


    cerr << "Reading in variation graph ..." << endl;

    auto [nodevector, minid, graph, node_path_matrix, path_names] = Euka().readPathHandleGraph(ogfilename, 1, gbtwfilename, dbprefixS, clade_vec);


    

    // Check outgroup makes sense
    if (outGroup != ""){
        bool outgroup_found = false;
        for (unsigned int clade_idx = 0; clade_idx < clade_vec->size(); ++clade_idx){
            const string clade_name = clade_vec->at(clade_idx)->name;
            if (outGroup == clade_name){outgroup_found = true;}
                                                                                    }
   if (!outgroup_found){throw runtime_error("[euka] Outgroup not found in reference graph");}
                       }

    //read information about bin complexity
    cerr << "Reading in bin information ..." << endl;
    cerr << "\t-------------------------------" << endl;
    vector<vector<tuple<int, int, double, double >>>  chunks = load_clade_chunks(binsfilename);
    cerr << "... done!" << endl;
    cerr << "\t-------------------------------" << endl;
    if(chunks.empty()){throw runtime_error("Bins file is empty unable to proceed");}



    const vector<double> qscore_vec = get_qscore_vec();
    // array initialised in the Euka.h by setting it to the highest bit (T = 84 + 1 for offset). All bit characters are given the base freq value. Is used in readGAM to index the model2 log-likelihood.
    base_freq['A'] = log(0.362815);
    base_freq['C'] = log(0.207743);
    base_freq['G'] = log(0.116809);
    base_freq['N'] = log(0.25);
    base_freq['T'] = log(0.312435);

    // array to store the transversion/transition rates for the mitogenome. Used by calculating the propability of a mismatch due to a mutation
    t_T_ratio['A']['A'] = 1.0;
    t_T_ratio['C']['C'] = 1.0;
    t_T_ratio['G']['G'] = 1.0;
    t_T_ratio['T']['T'] = 1.0;
    t_T_ratio['A']['C'] = 0.02381;
    t_T_ratio['A']['G'] = 0.95238;
    t_T_ratio['A']['T'] = 0.02381;
    t_T_ratio['C']['A'] = 0.02381;
    t_T_ratio['C']['G'] = 0.02381;
    t_T_ratio['C']['T'] = 0.95238;
    t_T_ratio['G']['A'] = 0.95238;
    t_T_ratio['G']['C'] = 0.02381;
    t_T_ratio['G']['T'] = 0.02381;
    t_T_ratio['T']['A'] = 0.02381;
    t_T_ratio['T']['C'] = 0.95238;
    t_T_ratio['T']['G'] = 0.02381;


    // array to store bools for rare base call: WSMKRYBDHVN
    rare_bases['W'] = true;
    rare_bases['M'] = true;
    rare_bases['K'] = true;
    rare_bases['R'] = true;
    rare_bases['Y'] = true;
    rare_bases['B'] = true;
    rare_bases['D'] = true;
    rare_bases['H'] = true;
    rare_bases['V'] = true;
    rare_bases['T'] = false;
    rare_bases['A'] = false;
    rare_bases['C'] = false;
    rare_bases['G'] = false;
    rare_bases['S'] = false;
    rare_bases['N'] = false; 

    ////////////////////////////////////  GIRAFFE //////////////////////////////////////


    string first_fifo = tmpdir + random_string(7);
    const char * fifo_A = first_fifo.c_str();
    mkfifo(fifo_A, 0666);
    pid_t wpid;
    int status = 0;
    pid_t pid1 = fork();
    if (pid1 == -1) {
        throw std::runtime_error("Error in fork");
                    }

    const vg::subcommand::Subcommand* sc = NULL;
    if(pid1 == 0) {
    // Child process
        if (gamfilename == "") {

        cerr << "Mapping reads..." << endl;
       
        Euka::map_giraffe(fastq1filename, fastq2filename, n_threads, interleaved,

                      fifo_A, sc, tmpdir, euka_dir, dbprefix, altMinFull);
    

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


#ifdef VERBOSE_S
    //cerr<<"found "<<nodevector.size()<<" nodes "<<endl;
    //cerr<<"Input file "<<samplename<<endl;
#endif
    cerr << "Estimating clades: Please be patient! Depending on the size of your input file, this process can take some time." << endl; 
    cerr << "\t-------------------------------" << endl;
    auto [readgam, clade_id_list] = readGAM3(graph,fifo_A,populateAlignmentVector,clade_vec,\
                                             nodevector,qscore_vec, base_freq, t_T_ratio, rare_bases,\
                                             chunks, dmg.subDeamDiNuc, lengthToProf, prof_out_file_path, MINIMUMMQ,\
                                             MINNUMOFREADS, MINNUMOFBINS, ENTROPY_SCORE_THRESHOLD, MAXIMUMOFBINS);
    remove(fifo_A);
    cerr << " .. done!" << endl;
    cerr << "\t-------------------------------" << endl;

//////// Write file with all fragment names if specified //////////

    if (outFrag == true && !clade_id_list.empty() && outGroup == ""){


        ofstream outnameStorage((outputfilename + "_FragNames.tsv").c_str(), ios::trunc);
        for (int i = 0; i < clade_id_list.size(); ++i){
            outnameStorage << clade_vec->at(clade_id_list.at(i)*6+1)->name << '\t';
            for (int s = 1; s<clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.size(); ++s){
                    outnameStorage << clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.at(s);
                    if (s != clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.size() -1){
                        outnameStorage << '\t';
                    }
                }
                outnameStorage << endl;


        }
    }else if (outFrag == true && !clade_id_list.empty() && outGroup != ""){

        int extra_id = -1;
        for (int j; j<chunks.size(); j++){
            if (clade_vec->at(j*6+1)->name == outGroup){
                extra_id =clade_vec->at(j*6+1)->id;
            }
        }

        clade_id_list.emplace_back(extra_id);
        ofstream outnameStorage((outputfilename + "_FragNames.tsv").c_str(), ios::trunc);
        for (int i = 0; i < clade_id_list.size(); ++i){
            outnameStorage << clade_vec->at(clade_id_list.at(i)*6+1)->name << '\t';
            for (int s = 1; s<clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.size(); ++s){
                    outnameStorage << clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.at(s);
                    if (s != clade_vec->at(clade_id_list.at(i)*6+1)->nameStorage.size() -1){
                        outnameStorage << '\t';
                    }
                }
                outnameStorage << endl;

        }
    }


///////// RUN MCMC ////////


    if (clade_id_list.size() < 2 || run_mcmc == false){


        const vector<long double> init_vec = Euka::compute_init_vec(clade_vec, clade_id_list);

        ofstream outbin((outputfilename +"_coverage.tsv").c_str(), ios::trunc);
        ofstream out((outputfilename + "_abundance.tsv").c_str(), ios::trunc);
        ofstream outsurv((outputfilename + "_detected.tsv").c_str(), ios::trunc); 
        ofstream outinSize((outputfilename + "_inSize.tsv").c_str(), ios::trunc);

        
        out << "#Taxa" << '\t' << "detected"<< '\t'<< "Number_of_reads" << '\t' << "proportion_estimate" << '\n';
        outsurv << "#Taxa" << '\t' << "detected"<< '\t'<< "Number_of_reads" << '\t' << "proportion_estimate" << '\n';
        outbin << "#Taxa" << '\t';
            for (int bins=0; bins<21; bins++){
                outbin << "bin" << bins << '\t' <<"entropy"; 
                if (bins != 20){
                    outbin << '\t';
                }
            }
            outbin << endl;


        vector<int> clade_bins;
        vector<int> clade_list_id;
        vector<int> clade_list_count;
        int extra_id = -1;

        for (int i = 0; i<chunks.size(); i++)
        {

            //save id for specified outGroup if specified:
            if (outGroup == clade_vec->at(i*6+1)->name){
                extra_id =clade_vec->at(i*6+1)->id;
            }


            vector<int> check_for_zero;

            for (int k = 0; k<chunks[i].size()-1; k++){
                if (get<2>(chunks.at(i).at(k)) > ENTROPY_SCORE_THRESHOLD){
                    check_for_zero.emplace_back(get<3>(chunks[i].at(k)));
                }
            }
            // a group is disregarded if they have an empty bin, their minimum number of bins above the entropy threshold is too low or their total amount of reads is to low.
            if (std::count(check_for_zero.begin(), check_for_zero.end(), 0.0) || check_for_zero.size() < MINNUMOFBINS  || clade_vec->at(i*6+1)->count < MINNUMOFREADS)
            {

                // in case a group is defined to be analysed no matter what this chunk of code will produce the coverage file
                out << clade_vec->at(i*6+1)->name << '\t' << "no" << '\t' << clade_vec->at(i*6+1)->count << '\t' << 0 << endl;
                if (outGroup == clade_vec->at(i*6+1)->name)
                {
                                

                    outbin << clade_vec->at(i*6+1)->name << '\t';

                    for (size_t j = 0; j<chunks[i].size()-1; j++){
                        outbin << std::fixed << std::setprecision(5);
                        outbin << get<3>(chunks[i].at(j)) << '\t' << get<2>(chunks[i].at(j));
                        if (j != chunks[i].size() -2){
                            outbin << '\t';
                        }

                    }

                    outbin << endl;

                    outinSize << clade_vec->at(i*6+1)->name << '\t';
                    for (int s = 1; s<clade_vec->at(i*6+1)->inSize.size(); s++){
                        outinSize << clade_vec->at(i*6+1)->inSize.at(s);
                        if (s != clade_vec->at(i*6+1)->inSize.size() -1){
                            outinSize << '\t';
                        }
                    }
                    outinSize << endl;
                }

            }
            // a group is detected: 
            else
            {

                // creating a list of ids 
                clade_list_id.emplace_back(clade_vec->at(i*6+1)->id);
                // creating a list of all fragments counts per group. We only added a count to the clade if it passed the likelihood test and the mapping quality filter
                // this is defined in readGAM_Euka.h line 515. We are not accessing any fragments that do not pass the first set of filters. 
                clade_list_count.emplace_back(clade_vec->at(i*6+1)->count);
                
                outbin << clade_vec->at(i*6+1)->name << '\t';

                for (int j = 0; j<chunks[i].size()-1; j++){
                    outbin << std::fixed << std::setprecision(5);
                    outbin << get<3>(chunks[i].at(j)) << '\t' << get<2>(chunks[i].at(j));
                    if (j != chunks[i].size() -1){
                        outbin << '\t';
                    }

                }
                outbin << '\n';

                out << clade_vec->at(i*6+1)->name << '\t' << "yes" << '\t' << clade_vec->at(i*6+1)->count << '\t';
                outsurv << clade_vec->at(i*6+1)->name << '\t' << "yes" << '\t' << clade_vec->at(i*6+1)->count << '\t';
                outinSize << clade_vec->at(i*6+1)->name << '\t';


                for (int s = 1; s<clade_vec->at(i*6+1)->inSize.size(); s++){
                    outinSize << clade_vec->at(i*6+1)->inSize.at(s);
                    if (s != clade_vec->at(i*6+1)->inSize.size() -1){
                        outinSize << '\t';
                    }
                }
                outinSize << endl; 

                for (int index = 0; index<clade_list_id.size();index++) {
                    out << init_vec.at(index);
                    outsurv << init_vec.at(index);
                    if (index != clade_list_id.size()-1){
                        out << '\t';
                        outsurv <<'\t';
                    }
                }
                
            
                out << endl;
                outsurv << endl;
                clade_bins.clear();
            }

            
        } // end of loop through all the groups. We created the abundance file, the detected file, the bins file, the inSize file. 


        /// Mikkel code begin ///
        // Print the clades that have been found
        // check if the program should write output files or to stdout
        string prof_out_file = outputfilename;
        string count_out_file = outputfilename;
        if (prof_out_file_path != outputfilename) {
            // Check for correct end of file-path
            if (prof_out_file_path.back() != '/') {
                prof_out_file_path = prof_out_file_path + "/";
            }
            // check if out-folder exits, if not then create folder
            if (!fs::is_directory(prof_out_file_path) || !fs::exists(prof_out_file_path)) { // Check if src folder exists
                fs::create_directory(prof_out_file_path); // create src folder
            }
        }
            
        vector<double>end5ct;
        vector<double>end3ct;
        vector<double>end5ga;
        vector<double>end3ga;

        // Initiate the 
        // Print the clades that have been found
        for (size_t i = 0; i < clade_list_id.size(); i++){

            // Find loaction of baseshift data array in clade_vec
	    
            unsigned int** baseshift_data_array_location = clade_vec->at(clade_list_id[i]*6+1)->baseshift_clade_array; 
            
            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);
            
            // make filename for clades
            prof_out_file = outputfilename + "_" + clade_vec->at(clade_list_id[i]*6+1)->name + ".prof";
            
            vector <vector <double >> dam;
            // Display substitution matrix for cladeqq
            dam = baseshift_data_array.display_prof(prof_out_file);

            

            int half = dam.at(0).size()/2;
            for (size_t a = 0; a<dam.at(0).size();a++){

                if (a < half){
                    end5ct.emplace_back(dam.at(0).at(a));
                }
                else{
                    end3ct.emplace_back(dam.at(0).at(a));
                }
             
            }
            int half1 = dam.at(1).size()/2;
            for (size_t a = 0; a<dam.at(1).size();a++){

                if (a < half1){
                    end5ga.emplace_back(dam.at(1).at(a));
                }
                else{
                    end3ga.emplace_back(dam.at(1).at(a));
                }
             
            }
                
            
            // Display entire baseshift array for clade
            cerr << endl;
            //baseshift_data_array.print_counts(count_out_file);

        }

        if (extra_id != -1){
            //cout << extra_id << endl; 
            unsigned int** baseshift_data_array_location = clade_vec->at(extra_id*6+1)->baseshift_clade_array; 
            
            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);
            
            // make filename for clades
            prof_out_file = outputfilename + "_" + clade_vec->at(extra_id*6+1)->name + ".prof";
            
            vector <vector <double >> rest;
            // Display substitution matrix for cladeqq
            rest = baseshift_data_array.display_prof(prof_out_file);


        }
        
        // // Mikkel code ends // //

        vector<double> end5ct_av = get_avg(end5ct, lengthToProf, clade_id_list.size());
        vector<double> end3ct_av = get_avg(end3ct, lengthToProf, clade_id_list.size());
        vector<double> end5ga_av = get_avg(end5ga, lengthToProf, clade_id_list.size());
        vector<double> end3ga_av = get_avg(end3ga, lengthToProf, clade_id_list.size());

        /// output combined damagefile for secondary euka run:

        ofstream outdamall((outputfilename +"_5p.prof").c_str(), ios::trunc);
        outdamall <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << endl;
        for (int pas = 0; pas<end5ct_av.size(); pas++){

            //cout << " is it here " << endl; 

            for (int col = 0; col<12; col++){

                if (col == 5){
                    outdamall << end5ct_av.at(pas) << '\t';
                }
                else if(col == 11){
                    outdamall << 0 << endl;
                }
                else {
                    outdamall << 0 << '\t';
                }

            }

        }

        ofstream outdamall2((outputfilename +"_3p.prof").c_str(), ios::trunc);
        outdamall2 <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << endl;

        reverse(end3ga_av.begin(), end3ga_av.end());

        for (int pas = 0; pas<end3ga_av.size(); pas++){

            //cout << " is it here " << endl; 

            for (int col = 0; col<12; col++){

                if (col == 6){
                    outdamall2 << end3ga_av.at(pas) << '\t';
                }
                else if(col == 11){
                    outdamall2 << 0 << endl;
                }
                else {
                    outdamall2 << 0 << '\t';
                }

            }

        }


        cerr << "Abundance estimation completed! " << endl; 
        cerr << '\n';
        cerr << "No MCMC was computed. It was either specified by the user or less than 2 groups were present in the sample." << endl; 
        cerr << "You can find all four output files ("+outputfilename + "_abundance.tsv, "+outputfilename +"_detected.tsv, and "+outputfilename+"_coverage.tsv, and the damage profiles) in your current working directory!" << endl;
        cerr << '\n';
        cerr << '\n';
        cerr << "You have two options to visualize euka’s output. All scripts necessary for the visualization can be found in vgan/share/vgan/plottingScripts/." << endl; 
        cerr << '\n';
        cerr << '\t' << '\t' << "1) ./visualize_detected_taxa.sh "+ outputfilename << endl;
        cerr << '\n';
        cerr << "This will provide you with a summary plot of each detected taxa, including their damage profiles, estimated coverage across the pangenome graph and fragment length distribution." << endl;
        cerr << '\n';
        cerr << '\t' << '\t' << "2) python vgan/share/vgan/plottingScripts/make_tree_from_output.py " +outputfilename+"_abundance.tsv or python vgan/share/vgan/plottingScripts/make_tree_from_output.py " +outputfilename+"_detected.tsv" << endl; 
        cerr << '\n';
        cerr << "These two commands will plot a taxonomic tree with all ("+outputfilename+"_abundance.tsv) taxa or only the detected (" +outputfilename+"_detected.tsv) taxa." << endl; 

    }
    // there are more than 1 group detected or the mcmc was not turned off. 
    else
    {

        const vector<long double> init_vec = Euka::compute_init_vec(clade_vec, clade_id_list);
        vector<long double > clade_res = MCMC().run(iter, burnin, 0.01, init_vec, clade_vec, clade_id_list);

    /////// CREATE OUTPUT FILES ///////
        
        ofstream outbin((outputfilename +"_coverage.tsv").c_str(), ios::trunc);
        ofstream out((outputfilename + "_abundance.tsv").c_str(), ios::trunc);
        ofstream outsurv((outputfilename + "_detected.tsv").c_str(), ios::trunc); 
        ofstream outinSize((outputfilename + "_inSize.tsv").c_str(), ios::trunc);

        
        out << "#Taxa" << '\t' << "detected"<< '\t'<< "Number_of_reads" << '\t' << "proportion_estimate" << '\t' << "85%_confidence_interval_lower_bound" << '\t' << "85%_confidence_interval_higher_bound" << '\t' << "95%_confidence_interval_lower_bound" << '\t' << "95%_confidence_interval_higher_bound" << '\n';
        outsurv << "#Taxa" << '\t' << "detected"<< '\t'<< "Number_of_reads" << '\t' << "proportion_estimate" << '\t' << "85%_confidence_interval_lower_bound" << '\t' << "85%_confidence_interval_higher_bound" << '\t' << "95%_confidence_interval_lower_bound" << '\t' << "95%_confidence_interval_higher_bound" << '\n';
        outbin << "#Taxa" << '\t';
            for (int bins=0; bins<21; bins++){
                outbin << "bin" << bins << '\t' <<"entropy"; 
                if (bins != 20){
                    outbin << '\t';
                }
            }
            outbin << endl;

        vector<int> clade_bins;
        vector<int> clade_list_id;
        vector<int> clade_list_count;
        int extra_id = -1;

        for (int i = 0; i<chunks.size(); i++) {

            //save id for specified outGroup if specified:
            if (outGroup == clade_vec->at(i*6+1)->name){
                extra_id =clade_vec->at(i*6+1)->id;
            }


            vector<int> check_for_zero;
            // -1, because we do not count the last bin do to the fact there is a chance it is not 1500 bases. 
            for (int k = 0; k<chunks[i].size()-1; k++){
                if (get<2>(chunks.at(i).at(k)) > ENTROPY_SCORE_THRESHOLD){
                    check_for_zero.emplace_back(get<3>(chunks[i].at(k)));

                }
            }
            // disregarded groups except its defined as outgroup
            if (std::count(check_for_zero.begin(), check_for_zero.end(), 0.0) || check_for_zero.size() < MINNUMOFBINS || clade_vec->at(i*6+1)->count < MINNUMOFREADS)
            {

                
                out << clade_vec->at(i*6+1)->name << '\t' << "no" << '\t' << clade_vec->at(i*6+1)->count << '\t' << 0 << '\t' << 0 << '\t'<< 0 << '\t'<< 0 << '\t'<< 0 << endl;
                if (outGroup == clade_vec->at(i*6+1)->name){

                    //cout << outGroup << endl; 

                    outbin << clade_vec->at(i*6+1)->name << '\t';

                    for (size_t j = 0; j<chunks[i].size()-1; j++){
                        outbin << std::fixed << std::setprecision(5);
                        outbin << get<3>(chunks[i].at(j)) << '\t' << get<2>(chunks[i].at(j));
                        if (j != chunks[i].size() -2){
                            outbin << '\t';
                        }

                    }

                    outbin << endl;

                    outinSize << clade_vec->at(i*6+1)->name << '\t';
                    for (int s = 1; s<clade_vec->at(i*6+1)->inSize.size(); s++){
                        outinSize << clade_vec->at(i*6+1)->inSize.at(s);
                        if (s != clade_vec->at(i*6+1)->inSize.size() -1){
                            outinSize << '\t';
                        }
                    }
                    outinSize << endl;
                }


                      
            } else {

                // creating a list of ids 
                clade_list_id.emplace_back(clade_vec->at(i*6+1)->id);
                // creating a list of all fragments counts per clade
                clade_list_count.emplace_back(clade_vec->at(i*6+1)->count);
                
                outbin << clade_vec->at(i*6+1)->name << '\t';

                for (size_t j = 0; j<chunks[i].size()-1; j++){
                    outbin << std::fixed << std::setprecision(5);
                    outbin << get<3>(chunks[i].at(j)) << '\t' << get<2>(chunks[i].at(j));
                    if (j != chunks[i].size() -2){
                        outbin << '\t';
                    }

                }

                outbin << endl;

                out << clade_vec->at(i*6+1)->name << '\t' << "yes" << '\t' << clade_vec->at(i*6+1)->count << '\t';
                outsurv << clade_vec->at(i*6+1)->name << '\t' << "yes" << '\t' << clade_vec->at(i*6+1)->count << '\t';
                outinSize << clade_vec->at(i*6+1)->name << '\t';
                for (int s = 1; s<clade_vec->at(i*6+1)->inSize.size(); s++){
                    outinSize << clade_vec->at(i*6+1)->inSize.at(s);
                    if (s != clade_vec->at(i*6+1)->inSize.size() -1){
                        outinSize << '\t';
                    }
                }
                outinSize << endl; 

                for (int index = (clade_list_id.size()-1)*5; index<clade_list_id.size()*5;index++) {
                    out << clade_res.at(index);
                    outsurv << clade_res.at(index);
                    if (index != clade_res.size()-1){
                        out << '\t';
                        outsurv << '\t';
                    }
                }

                out << endl;
                outsurv << endl;
                clade_bins.clear();
            }

        }

        cerr << "Abundance estimation completed! " << endl; 
        cerr << '\n'; 
        // // Mikkel code begin // //
        
        // Print the clades that have been found
        // check if the program should write output files or to stdout
        string prof_out_file = euka_dir;
        string count_out_file = euka_dir;
        if (prof_out_file_path != outputfilename) {
            // Check for correct end of file-path
            if (prof_out_file_path.back() != '/') {
                prof_out_file_path = prof_out_file_path + "/";
            }
            // check if out-folder exits, if not then create folder
            if (!fs::is_directory(prof_out_file_path) || !fs::exists(prof_out_file_path)) { // Check if src folder exists
                fs::create_directory(prof_out_file_path); // create src folder
            }
        }
            
        vector<double>end5ct;
        vector<double>end3ct;
        vector<double>end5ga;
        vector<double>end3ga;

        // Initiate the 
        // Print the clades that have been found
        for (size_t i = 0; i < clade_list_id.size(); i++){

            // Find loaction of baseshift data array in clade_vec
            unsigned int** baseshift_data_array_location = clade_vec->at(clade_list_id[i]*6+1)->baseshift_clade_array; 
            
            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);
            
            // make filename for clades
            prof_out_file = outputfilename + "_" + clade_vec->at(clade_list_id[i]*6+1)->name + ".prof";
            //count_out_file = outputfilename + clade_vec->at(clade_list_id[i]*6+1)->name + ".count";
            
            vector <vector <double >> dam;
            // Display substitution matrix for cladeqq
            dam = baseshift_data_array.display_prof(prof_out_file);
            //baseshift_data_array.print_counts(count_out_file);

            

            int half = dam.at(0).size()/2;
            for (size_t a = 0; a<dam.at(0).size();a++){

                if (a < half){
                    end5ct.emplace_back(dam.at(0).at(a));
                }
                else{
                    end3ct.emplace_back(dam.at(0).at(a));
                }
             
            }
            int half1 = dam.at(1).size()/2;
            for (size_t a = 0; a<dam.at(1).size();a++){

                if (a < half1){
                    end5ga.emplace_back(dam.at(1).at(a));
                }
                else{
                    end3ga.emplace_back(dam.at(1).at(a));
                }
             
            }
                
            

            // Display entire baseshift array for clade
            cerr << endl;
            //baseshift_data_array.print_counts(count_out_file);

        }
        if (extra_id != -1){
            unsigned int** baseshift_data_array_location = clade_vec->at(extra_id*6+1)->baseshift_clade_array; 
            
            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);
            
            // make filename for clades
            prof_out_file = outputfilename + "_" + clade_vec->at(extra_id*6+1)->name + ".prof";
            
            vector <vector <double >> rest;
            // Display substitution matrix for cladeqq
            rest = baseshift_data_array.display_prof(prof_out_file);


        }
        
        // // Mikkel code ends // //

            vector<double> end5ct_av = get_avg(end5ct, lengthToProf, clade_id_list.size());
            vector<double> end3ct_av = get_avg(end3ct, lengthToProf, clade_id_list.size());
            vector<double> end5ga_av = get_avg(end5ga, lengthToProf, clade_id_list.size());
            vector<double> end3ga_av = get_avg(end3ga, lengthToProf, clade_id_list.size());

            /// output combined damagefile for secondary euka run:

            ofstream outdamall((outputfilename +"_5p.prof").c_str(), ios::trunc);
            outdamall <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << endl;
            for (size_t pas = 0; pas<end5ct_av.size(); pas++){

                //cerr << " is it here " << endl; 

                for (int col = 0; col<12; col++){

                    if (col == 5){
                        outdamall << end5ct_av.at(pas) << '\t';
                    }
                    else if(col == 11){
                        outdamall << 0 << endl;
                    }
                    else {
                        outdamall << 0 << '\t';
                    }

                }

            }

            ofstream outdamall2((outputfilename +"_3p.prof").c_str(), ios::trunc);
            outdamall2 <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << endl;

            reverse(end3ga_av.begin(), end3ga_av.end());

            for (size_t pas = 0; pas<end3ga_av.size(); pas++){

                //cerr << " is it here " << endl; 

                for (int col = 0; col<12; col++){

                    if (col == 6){
                        outdamall2 << end3ga_av.at(pas) << '\t';
                    }
                    else if(col == 11){
                        outdamall2 << 0 << endl;
                    }
                    else {
                        outdamall2 << 0 << '\t';
                    }

                }

            }



        std::size_t pos = outputfilename.find_last_of("/");
        string file_in_dir = outputfilename.substr(pos+1);

        cerr << "You can find all four output files ("+file_in_dir + "_abundance.tsv, "+file_in_dir +"_detected.tsv, "+file_in_dir+"_coverage.tsv, and the damage profils) in your current working directory!" << endl;
        cerr << '\n';
        cerr << '\n';
        cerr << "You have two options to visualize euka’s output. All scripts necessary for the visualization can be found in vgan/share/vgan/plottingScripts/." << endl; 
        cerr << '\n';
        cerr << '\t' << '\t' << "1) ./visualize_detected_taxa.sh "+ file_in_dir << endl;
        cerr << '\n';
        cerr << "This will provide you with a summary plot of each detected taxa, including their damage profiles, estimated coverage across the pangenome graph and fragment length distribution." << endl;
        cerr << '\n';
        cerr << '\t' << '\t' << "2) python vgan/share/vgan/plottingScripts/make_tree_from_output.py " +file_in_dir+"_abundance.tsv or python vgan/share/vgan/plottingScripts/make_tree_from_output.py " +file_in_dir+"_detected.tsv" << endl; 
        cerr << '\n';
        cerr << "These two commands will plot a taxonomic tree with all ("+file_in_dir+"_abundance.tsv) taxa or only the detected (" +file_in_dir+"_detected.tsv) taxa." << endl; 


  
    }


    
    return 0;

}

