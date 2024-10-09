#include "soibean.h"
#include <omp.h>
#include <sys/wait.h>
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/vg.pb.h>
#include <vg.hpp>
#include <boost/algorithm/string.hpp>
#include <handlegraph/handle_graph.hpp>
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"
#include "miscfunc.h"
#include "vgan_utils.h"
#include "Clade.h"
#include "MCMC.h"
#include "damage.h"


//#define DEBUGLCA
//#define PRINTLCA
//#define SANITYCHECK

using namespace vg;

soibean::soibean(){

}

soibean::~soibean(){

}

// Function to print the header
void soibeanSetup() {

 preflight_check();
  configure_memory_allocator();
  enable_crash_handling();
  temp_file::set_system_dir();

  if (!vg::io::register_libvg_io()) {
      cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
      exit(1);
                                    }

    std::cerr << "                 _ \n"
                 "              _ ( )\n"
                 "  ___    _   (_)| |_      __     _ _   ___  \n"
                 "/',__) /'_`\\ | || '_`\\  /'__`\\ /'_` )/' _ `\\ \n"
                 "\\__, \\( (_) )| || |_) )(  ___/( (_| || ( ) | \n"
                 "(____/\\___/'(_)(_,__/'`\\____)`\\__,_)(_) (_) \n"
              << std::endl;
}

bool validateGraphInputFiles(const string& ogfilename, const string& cladefilename, const string& binsfilename, const string& gbwtfilename, \
                             const string& gamfilename, const string& fq1filename, const string& fq2filename) {
    if (!isFile(ogfilename)) {
        throw runtime_error(ogfilename + " does not exist.");
    }
    if (!isFile(cladefilename)) {
        throw runtime_error(cladefilename + " does not exist.");
    }
    if (!isFile(binsfilename)) {
        throw runtime_error(binsfilename + " does not exist.");
    }
    if (!isFile(gbwtfilename)) {
        throw runtime_error(gbwtfilename + " does not exist.");
    }
      if (!isFile(gamfilename)) {
        throw runtime_error(cladefilename + " does not exist.");
    }
    if (!isFile(fq1filename)) {
        throw runtime_error(binsfilename + " does not exist.");
    }
    if (!isFile(fq2filename)) {
        throw runtime_error(gbwtfilename + " does not exist.");
    }

    return true;
}


string soibean::usage() const {
    return string("") +
                 "                 _ \n"
                 "              _ ( )\n"
                 "  ___    _   (_)| |_      __     _ _   ___  \n"
                 "/',__) /'_`\\ | || '_`\\  /'__`\\ /'_` )/' _ `\\ \n"
                 "\\__, \\( (_) )| || |_) )(  ___/( (_| || ( ) | \n"
                 "(____/\\___/'(_)(_,__/'`\\____)`\\__,_)(_) (_) \n"
              +
        "\n"+
        "\n"+
        "Usage: soibean [options]\n" +
        "Options:\n" +
        "  --soibean_dir <directory>   Specify the directory containing the soibean files\n" +
        "  --dbprefix <prefix>         Specify the prefix for the database files\n" +
        "  -o [STR]                    Output file prefix (default: beanOut) "+"\n"+
        "  -fq1 <filename>             Specify the input FASTQ file (single-end or first pair)\n" +
        "  -fq2 <filename>             Specify the input FASTQ file (second pair)\n" +
        "  -g <filename>               Specify the input GAM file\n" +
        "  -t <threads>                Specify the number of threads to use (default: 1)\n" +
        "  -z <directory>              Specify the temporary directory for intermediate files (default: /tmp/)\n" +
        "  -l <length>                 Specify the length threshold for profiling (default: 5)\n" +
        "  -i                          Enable interleaved input mode\n" +
        "Markov chain Monte Carlo options:\n" +
        "  --no-mcmc                   The MCMC does not run (default: false)\n" +
        "  --chains [INT]              Define the number of chains for the MCMC (default: 4)\n" +
        "  --iter [INT]                Define the number of iterations for the MCMC (default: 1.000.000)\n" +
        "  --burnin [INT]              Define the burn-in period for the MCMC (default: 100.000)\n";
}

///////////// Begin of function storage - temporary ////////////////////////////

std::vector<std::vector<double>> soibean::convertMapsToVector(std::vector<AlignmentInfo*>* & gam) {
    std::vector<std::vector<double>> result;

    if (gam->empty()) {
        throw std::runtime_error("Error: GAM file is empty.");
    }

    for (const auto& map : *gam) {
        std::vector<double> row;
        for (const auto& pair : map->pathMap) {
            row.push_back(pair.second); 
        }

        // Compute the sum of the elements in the row
        double sum = 0.0;
        for (int i = 0; i < row.size(); i++) {
            sum += row[i];
        }

        // Push the row into the result vector
        result.push_back(row);
    }

    return result;
}

const int soibean::run(int argc, char *argv[], const string & cwdProg)
{

  soibeanSetup();

  int lastOpt=1;
  bool run_mcmc=true;
  int n_threads = 1;
  bool interleaved=false;
  string fastq1filename, fastq2filename, gamfilename, samplename;
  string tmpdir = "/tmp/";
  bool sbdirspecified = false;
  string sbdir = "../share/soibean_dir/";
  string dbprefixS = "soibean_dir";
  string outputfilename = "beanOut";
  int lengthToProf = 5;
  unsigned int iter=50000;
  unsigned int burnin=5000;
  unsigned int chains = 1;
  bool specifiedDeam=false;

  string deam5pfreqE  = getFullPath(cwdProg+"../share/damageProfiles/none.prof");
  string deam3pfreqE  =  getFullPath(cwdProg+"../share/damageProfiles/none.prof");


  for(int i=1;i<(argc);i++)
  {

    if(string(argv[i]) == "--soibean_dir"){
        sbdir = argv[i+1];
        sbdirspecified=true; 
        if(sbdir.back() != '/'){sbdir += '/';}
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
            {throw runtime_error("[soibean] Input file must be FASTQ, not FASTA");}
        continue;
                            }

    if(string(argv[i]) == "-fq2"){
        fastq2filename = argv[i+1];
        if (fastq2filename.ends_with(".fa") || fastq2filename.ends_with(".fasta") || fastq2filename.ends_with(".fa.gz") || fastq2filename.ends_with(".fasta.gz"))
            {throw runtime_error("[soibean] Input file must be FASTQ, not FASTA");}
        continue;
    }

    if(string(argv[i]) == "-g"){
        gamfilename = argv[i+1];
        samplename = gamfilename;
        continue;
                            }

    if(string(argv[i]) == "-t"){
      if (stoi(argv[i+1]) < -1 || stoi(argv[i+1]) == 0) {throw runtime_error("[soibean] Error, invalid number of threads");}
      if (stoi(argv[i+1]) == -1) {n_threads = thread::hardware_concurrency();}
      else if (stoi(argv[i+1]) <= thread::hardware_concurrency()) {
              n_threads = stoi(argv[i+1]);
                                                                       }
    else {
           cerr << "[soibean] Warning, specified number of threads is greater than the number available. Using " << n_threads << " threads\n";
           n_threads = thread::hardware_concurrency();
         }
        continue;
                                   }

      if(string(argv[i]) == "-z"){
          tmpdir = argv[i+1];
          if (tmpdir.back() != '/') {tmpdir += '/';}
          continue;
                             }
      if(string(argv[i]) == "-l"){
        lengthToProf = stoi(argv[i+1]);
        continue;
      }
      if(string(argv[i]) == "-i"){
            interleaved = true;
            if (fastq2filename != ""){throw runtime_error("[soibean] If interleaved option chosen, soibean expects only one FASTQ file");}
            continue;
                               }
       if(string(argv[i]) == "--no-mcmc"  ){
            run_mcmc=false;
            continue;
        }
        if(string(argv[i]) == "--iter" || string(argv[i]) == "--iterations"){
            iter=stoi(argv[i+1]);
            assert(iter >= 0);
            continue;
        }
        if(string(argv[i]) == "--burnin"){
            burnin=stoi(argv[i+1]);
            assert(burnin >= 0);
            continue;
        }
        if(string(argv[i]) == "--chains"){
            chains=stoi(argv[i+1]);
            assert(burnin >= 0);
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
        if(string(argv[i]) == "-o"){
            outputfilename = argv[i+1];
            continue;
                               }

    }


    if(!sbdirspecified){
        sbdir = cwdProg + sbdir;
    }

    Euka ek;

    string dbprefix              = sbdir + dbprefixS;

    string ogfilename     = dbprefix+".og";
    string vgfilename     = dbprefix+".vg";
    string cladefilename  = sbdir + "soibean_db.clade";
    string binsfilename   = sbdir + "soibean_db.bins";
    string gbwtfilename = dbprefix+".gbwt";

    Damage dmg;
    dmg.initDeamProbabilities(deam5pfreqE,deam3pfreqE);

    vector<Clade *> * clade_vec = ek.load_clade_info(cladefilename, lengthToProf);
   
    cerr << "Reading in variation graph ..." << endl;

  
    auto [nodevector, minid, graph, node_path_matrix, path_names] = ek.readPathHandleGraph(ogfilename, 1, gbwtfilename, dbprefixS, clade_vec);
    if (nodevector.empty()) {
	throw std::runtime_error("Error: The nodevector is empty. Unable to proceed.");
    }
    if (path_names.empty()) {
	throw std::runtime_error("Error: The path_names vector is empty. Unable to proceed.");
    }
    if (node_path_matrix.empty()) {
	throw std::runtime_error("Error: The node_path_matrix is empty. Unable to proceed.");
    }
    
#ifdef PRINTLCA
    cerr << "TOTAL NUMBER OF PATHS: " << path_names.size() << endl;
    for (int i = 0; i < path_names.size(); ++i){
        cerr << i << '\t' << path_names[i] << endl;
    }
#endif


    //read information about different taxa

    vector<vector<tuple<int, int, double, double >>> chunks =  ek.load_clade_chunks(binsfilename);

    const vector<double> qscore_vec = ek.get_qscore_vec();
    

    cerr << "Collecting path information ..." << endl;

    vector<vector<string>> nodepaths;

    for (int i = minid; i < minid+nodevector.size(); ++i){

        //printprogressBarCerr( float(i)/float(nodevector.size()) );

        vector<string> paths;
        //uint64_t node = reinterpret_cast<uint64_t>(nodevector[i]);
        uint64_t node = i;
        // Get the handle for the node
        bdsg::handle_t node_handle = graph.get_handle(node);
        // Find the paths that go through the given node
        paths = soibean::paths_through_node(graph, node_handle);

        nodepaths.push_back(paths);
    }

    //cout << nodepaths.size()<< endl;
    //cout << nodepaths[1].size() << endl;


    cerr << " ... done!" << endl; 

#ifdef DEBUGLCA
    cerr << nodepaths.size() <<" " << nodepaths.at(1).size() << endl; 
    for (int i = 0; i < nodepaths.at(1).size(); ++i){
        cerr << nodepaths.at(6752352)[i] << endl;
    }
#endif


    ////////////////////////////////////  GIRAFFE //////////////////////////////////////

    string first_fifo = tmpdir + random_string(9);
    const char * fifo_A = first_fifo.c_str();
    mkfifo(fifo_A, 0666);
    pid_t wpid;
    const vg::subcommand::Subcommand* sc = NULL;

    pid_t pid1 = fork();
    
    if (pid1 == -1) {
    throw runtime_error("Error in fork");
    }
    
    if(pid1 == 0) {    // Child process
	// Code for writing to the FIFO
	if (gamfilename == "") {
	    cerr << "mapping reads..." << endl;
	    ek.map_giraffe(fastq1filename, fastq2filename, n_threads, interleaved, fifo_A, sc, tmpdir, sbdir, dbprefix);
	    cerr << "and done" << endl;
	    exit(0);
	} else {
	    // Redirect buffer in case of GAM input
	    ifstream src(gamfilename);
        ofstream dst(fifo_A);
        dst << src.rdbuf();
        exit(0);
	}
    }

    int status;
    // Block the parent process until the child process finishes writing
    //while ((wpid = wait(&status)) > 0);

    cerr << "done mapping reads..." << endl;

    /////////////////////////////////////////////////////////////////////////////////

    bool entire_graph = false;
    double taxaMu = 0.0;
    if(dbprefixS == "soibean_db"){
        entire_graph = true;
    }else{
        for (int i = 1; i<clade_vec->size(); i+=6){
            if(clade_vec->at(i)->name == dbprefixS){
                taxaMu = clade_vec->at(i)->dist;
            }
        }
    }

    shared_ptr dta = make_unique<Trailmix_struct>();
    auto gam = precompute_GAM(graph, fifo_A, clade_vec, nodepaths, path_names, qscore_vec, false, minid, false, dmg.subDeamDiNuc, dta);

    //turn all alignment likelihoods for each path into a vector<vector<int>> for the MCMC input
    vector<vector<double>> probMatrix = convertMapsToVector(gam);

    cerr << "Number of paths: " << probMatrix[0].size() << endl;
    cerr << "Number of reads: " << gam->size() << endl;

    cerr << "Loading tree ... " << endl;
    spidir::Tree taxatree = NULL;  // Initialize taxatree pointer to nullptr
    string treename = sbdir+"/tree_dir/"+dbprefixS+".new.dnd";
    spidir::readNewickTree(treename.c_str(), &taxatree);

    if (taxatree.nodes.size() == 0){
        throw runtime_error("The tree is empty");
    }
    unordered_map <string, int> pathToTreeNode;
    unordered_map <string , string > parentMap;
    for (int j = 0; j< taxatree.nodes.size(); ++j) {
	pathToTreeNode[taxatree.nodes[j]->longname] = taxatree.nodes[j]->name;
	
        if (taxatree.nodes[j]->parent == NULL){
            parentMap[taxatree.nodes[j]->longname] = taxatree.nodes[j]->longname;
        }else{
            parentMap[taxatree.nodes[j]->longname] = taxatree.nodes[j]->parent->longname;
        }   
    }

    cerr << " ... done!" << endl;


    //get the signature node frequencies for the k estimnation
    cerr << "Finding the initial estimate ..." << endl;
    unordered_map <string, int> frequencies;
    for (const auto &read :*(gam))
    {

        //cout << read->pathMap["N11"] << endl;

        if (read->mostProbPath.size() == 1)
        {
            //cout << "There is a sig node" << endl;
            frequencies[read->mostProbPath[0]]++;
        }
    }
    

    MCMC mcmc;

    std::string freqFilename = sbdir+"soibean_db.baseFreq"; 
    std::ifstream freqFile(freqFilename);

    if (!freqFile.is_open()) {
        std::cerr << "Failed to open the base frequency file." << std::endl;
        return 1;
    }

    std::string searchName = dbprefixS;
    std::array<double, 256> freqs = {}; // Initialize all values to 0
    
    std::string line;
    while (std::getline(freqFile, line)) {
        std::istringstream iss(line);
        std::string name;
        double a, c, g, t;

        iss >> name >> a >> c >> g >> t;

        if (name == searchName) {
            freqs['A'] = a;
            freqs['C'] = c;
            freqs['G'] = g;
            freqs['T'] = t;
            break; // Exit the loop once the match is found
        }
    }

    

    // Assign the frequencies to the map this is a test set for the bears only! 
    // freqs['A'] = 0.313122;
    // freqs['C'] = 0.255707;
    // freqs['G'] = 0.153784;
    // freqs['T'] = 0.277387;
    freqs['R'] = freqs['A'] + freqs['G'];
    freqs['Y'] = freqs['C'] + freqs['T'];
    freqs['M'] = 1/(2*(   ((1/22)*(freqs['A']*freqs['G'])) + ((1/22)*(freqs['C']*freqs['T'])) + (freqs['A']*freqs['C'] + (freqs['A']*freqs['T']) \
                                                           + (freqs['G']*freqs['C'] + (freqs['G'] * freqs['T']) )) ) );


    vector<pair<string, int>> freqVec(frequencies.begin(), frequencies.end());
    sort(freqVec.begin(),freqVec.end(), [](const auto & a, const auto & b){
        return a.second > b.second;
    });
    
    vector<int> sigNodes;
    vector<string> sigPaths;
    double thres = gam->size()*0.01;
    cout << thres << endl;
    for (const auto & pair : freqVec)
    {
        if(pair.second >= thres) {
            sigNodes.emplace_back(pathToTreeNode[pair.first]);
            sigPaths.emplace_back(pair.first);

            cerr << "Identified signature paths: "<< pair.first << " with tree node: " << pathToTreeNode[pair.first] << endl;
        }

    }

    cerr << "... done! The Identified signature paths are used as input for the MCMC. " << endl;

    //for (int num = 0; num < 28; ++num){
    //    sigNodes = {num};
    
    //sigNodes = {10,21};
    string num = outputfilename;
    
    for (size_t i = 0; i < sigNodes.size(); ++i) {

        
        std::vector<int> subVector(sigNodes.begin(), sigNodes.begin() + i + 1);

        double logLike = 0.0L;
        double freq = log(1.0/sigNodes.size());
        if (sigPaths.size() == 1){


            for (const auto &read :*(gam))
            {
                logLike += read->pathMap[sigPaths[i]];
            }
            cerr << "PathMap loglikelihood: " << logLike << endl;

        }else{

            for (const auto &read :*(gam))
            {

                double inter = freq + read->pathMap[sigPaths[0]];
                for (int j = 1; j<subVector.size(); ++j)
                {
                    inter = oplusInitnatl(inter, (freq + read->pathMap[sigPaths[j]]));
                
                }
                logLike += inter;
            }

        }


        spidir::Tree* taxatreePtr = &taxatree;
        RunTreeProportionParams params(taxatreePtr);
        params.probMatrix = probMatrix;
        params.sources = vector<unsigned int>{16};//subVector;
        params.align = gam;
        params.burn = burnin;
        params.maxIter = iter;
        params.chains = chains;
        params.logLike = logLike;
        params.freqs = freqs;

        // Define the number of chains.
        std::vector<std::vector<MCMCiteration>> MCMCiterationsVec(chains);

        if (run_mcmc){
            
            unsigned int chainIndex = 0;
            for (auto& chainVec : MCMCiterationsVec) {
                std::cerr << "Running chain number: " << chainIndex << std::endl;
                //mcmc.run_tree_proportion(params, chainVec, graph, nodepaths, num, dta);
                chainIndex++;
            }

             cerr << "DONE WITH MCMC" << endl;
          
        //mcmc.processMCMCiterations(MCMCiterationsVec, i);
        //}
    }
}

    remove(fifo_A);

return 0;
}


