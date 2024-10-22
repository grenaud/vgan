
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
#include "MCMC_sb.h"
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

    std::cerr << "                  _ \n"
                 "               _ ( )\n"
                 "   ___    _   (_)| | _     __     _ _   ___  \n"
                 " /',__) /'_`\\ | || '_`\\  /'__`\\ /'_` )/' _ `\\ \n"
                 " \\__, \\( (_) )| || |_) )(  ___/( (_| || ( ) | \n"
                 " (____/ \\___/'(_)(_,__/'`\\____)`\\__,_)(_) (_) \n"
              << std::endl;
}


const string soibean::usage(const std::string& cwdProg) const {
    // Assuming getFullPath(cwdProg) returns the full directory path of the executable
    std::string execDir = getFullPath(cwdProg);

    // Construct the paths by appending the specific paths to the executable directory path
    std::string makeGraphFilesCmd = execDir + "../share/vgan/soibean_dir/make_graph_files.sh";


    return string("") +
                 "                  _ \n"
                 "               _ ( )\n"
                 "   ___    _   (_)| | _     __     _ _   ___  \n"
                 " /',__) /'_`\\ | || '_`\\  /'__`\\ /'_` )/' _ `\\ \n"
                 " \\__, \\( (_) )| || |_) )(  ___/( (_| || ( ) | \n"
                 " (____/ \\___/'(_)(_,__/'`\\____)`\\__,_)(_) (_) \n"
              +
        "\n"+
        "\n"+
        "Usage: soibean [options]\n" +
        "\n"+
        "First, to create a taxon of interest for the --dbprefix option please use:\n"+
        "      " + makeGraphFilesCmd + " [taxon name]\n" +
        "The taxon name must be from our database. To get an overview of the available taxa use:\n"
        "      " + makeGraphFilesCmd + " taxa\n" +
        "\n"+
        "\n"+
        "No damage example:\n"+
        "\tvgan soibean -fq1 [input.fq.gz] --dbprefix [taxon name] -o [output prefix]\n"
        "\n"+
        "Damage example:\n"+
        "\tvgan soibean -fq1 [input.fq.gz] --dbprefix [taxon name] -o [output prefix] --deam5p ../share/damageProfiles/dhigh5.prof --deam3p ../share/damageProfiles/dhigh3.prof\n"
        "\n"+
        "User defined MCMC:\n"+
        "\tvgan soibean -fq1 [input.fq.gz] --dbprefix [taxon name] -o [output prefix] --iter 1000000 --burnin 100000 -k 3\n"
        "\n"+

        "Options:\n" +
        "  --soibean_dir <directory>   Specify the directory containing the soibean files\n" +
        "  --tree_dir <directory>      Specify the directory containing the HKY trees\n" +
        "  --dbprefix <prefix>         Specify the prefix for the database files\n" +
        "  -o [STR]                    Output file prefix (default: beanOut) "+"\n"+
        "  -fq1 <filename>             Specify the input FASTQ file (single-end or first pair)\n" +
        "  -fq2 <filename>             Specify the input FASTQ file (second pair)\n" +
        "  -t <threads>                Specify the number of threads to use (default: 1)\n" +
        "  -z <directory>              Specify the temporary directory for intermediate files (default: /tmp/)\n" +
        "  -i                          Enable interleaved input mode\n" +

        "Damage options:"+"\n"+
        "  --deam5p [.prof]            5p deamination frequency for eukaryotic species (default: no damage)"+"\n"+
        "  --deam3p [.prof]            3p deamination frequency for eukaryotic species (default: no damage)"+"\n"+


        "Markov chain Monte Carlo options:\n" +
        "  --no-mcmc                   The MCMC does not run (default: false)\n" +
        "  --chains [INT]              Define the number of chains for the MCMC (default: 4)\n" +
        "  --iter [INT]                Define the number of iterations for the MCMC (default: 500.000)\n" +
        "  --burnin [INT]              Define the burn-in period for the MCMC (default: 75.000)\n"+
        "  --randStart [bool]          Set to get random starting nodes in the tree instead of the signature nodes (default: false)\n"+
        "  -k [INT]                    User defined value of k (k = number of expected sources) (default: not defined)\n";
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

std::vector<unsigned int> soibean::generateRandomNumbers(const int maxNum, const int k) {
    std::vector<unsigned int> sigNodes;
    std::random_device rd;
    std::mt19937 gen(rd());
    unsigned int min_value = 0;
    unsigned int max_value = static_cast<unsigned int>(maxNum) - 1;  // Adjusted to prevent out-of-range
    std::uniform_int_distribution<> distrib(min_value, max_value);

    unsigned int n = static_cast<unsigned int>(k);  // Number of random numbers to generate
    for(unsigned int i = 0; i < n; ++i) {
        unsigned int random_number = static_cast<unsigned int>(distrib(gen));
        sigNodes.emplace_back(random_number);
    }

    return sigNodes;
}

double soibean::calculateRhat(const std::vector<double>& means, const std::vector<double>& variances, int chainLength) {
    int numChains = means.size();
    if (numChains < 2) {
        std::cerr << "Warning: The Rhat could not be computed, due to insufficant number of chains. " << std::endl;
        return -1;
    }

    if (means.size() != variances.size()) {
        std::cerr << "Error: Mismatched sizes of means and variances vectors." << std::endl;
        return -1;
    }

    // Compute within-chain variance W
    double W = std::accumulate(variances.begin(), variances.end(), 0.0) / numChains;

    // Compute between-chain variance B
    double grandMean = std::accumulate(means.begin(), means.end(), 0.0) / numChains;
    double B = 0.0;
    for (int i = 0; i < numChains; ++i) {
        B += std::pow(means[i] - grandMean, 2);
    }
    B *= chainLength / (numChains - 1);

    // Compute R-hat
    double varEstimate = ((chainLength - 1.0) * W + B) / chainLength;
    double Rhat = std::sqrt(varEstimate / W);

    return Rhat;
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
  bool treedirspecified = false;
  string sbdir = getFullPath(cwdProg+"../share/vgan/soibean_dir/");
  string treedir = getFullPath(cwdProg+"../share/vgan/soibean_dir/tree_dir/");
  string dbprefixS = "soibean_db";
  string outputfilename = "beanOut";
  int lengthToProf = 5;
  unsigned int iter=100000;
  unsigned int burnin=15000;
  unsigned int chains = 4;
  bool dbprefixFound=false;
  bool specifiedDeam=false;
  bool randStart = false;
  bool specifiedk = false;
  double con = 0.0;
  int k = 1;
  size_t cutk = 0;

  string deam5pfreqE  = getFullPath(cwdProg+"../share/vgan/damageProfiles/none.prof");
  string deam3pfreqE  =  getFullPath(cwdProg+"../share/vgan/damageProfiles/none.prof");


  for(int i=1;i<(argc);i++)
  {

    if(string(argv[i]) == "--soibean_dir" || string(argv[i]) == "--soibean-dir"){
        sbdir = argv[i+1];
        sbdirspecified=true; 
        if(sbdir.back() != '/'){sbdir += '/';}
        continue;
    }
    if(string(argv[i]) == "--tree_dir" || string(argv[i]) == "--tree-dir"){
        treedir = argv[i+1];
        treedirspecified=true; 
        if(treedir.back() != '/'){treedir += '/';}
        continue;
    }

    if(string(argv[i]) == "--dbprefix"){
        dbprefixS = argv[i+1];
        dbprefixFound=true;
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
            assert(chains >= 0);
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
        if(string(argv[i]) == "--randStart"  || string(argv[i]) == "--randstart"){
            randStart=true;
            continue;
        }
        if(string(argv[i]) == "-k"){
            k=stoi(argv[i+1]);
            assert(k >= 0);
            specifiedk = true;
            continue;
        }
        if(string(argv[i]) == "--con"){
            con=stod(argv[i+1]);
            //assert(con >= 0);
            continue;
        }
        if(string(argv[i]) == "-cutk"){
            cutk=stoi(argv[i+1]);
            assert(cutk >= 0);
            continue;
        }


    }



    Euka ek;

    string dbprefix              = sbdir + dbprefixS;

    string ogfilename     = dbprefix+".og";
    string vgfilename     = dbprefix+".vg";
    string cladefilename  = sbdir + "soibean_db.clade";
    string binsfilename   = sbdir + "soibean_db.bins";
    string gbwtfilename = dbprefix+".gbwt";

    if (specifiedDeam){
        if (deam5pfreqE.size() == 0 || deam3pfreqE.size() == 0){
            throw runtime_error("Error the damage profiles do not exist. Unable to proceed.");
        }
    }
    Damage dmg;
    dmg.initDeamProbabilities(deam5pfreqE,deam3pfreqE);
    
    if (!dbprefixFound){
        throw runtime_error("No database specified. Please choose a taxon of interest and specify it with the --dbprefix option. You can create a graph for you taxon of interested by using the make_graph_file.sh script.");
    }


    vector<Clade *> * clade_vec = ek.load_clade_info(cladefilename, lengthToProf);
    if (clade_vec->empty()) {
    throw std::runtime_error("Error: The clade vector is empty. Unable to proceed. Check if the soibean.clade file is not empty.");
    }

    if (iter < burnin) {
    throw std::runtime_error("The number of iterations must be higher than the burn-in period. Unable to proceed.");
    }

   
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

    //vector<vector<tuple<int, int, double, double >>> chunks =  ek.load_clade_chunks(binsfilename);

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
	    //cerr << "Mapping reads..." << endl;
	    ek.map_giraffe(fastq1filename, fastq2filename, n_threads, interleaved, fifo_A, sc, tmpdir, sbdir, dbprefix);
	    //cerr << "and done" << endl;
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
    
    cerr << "...done!" << endl;
    
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

    auto gam = analyse_GAM(graph,fifo_A,clade_vec,nodevector, nodepaths, path_names, qscore_vec, false, minid, false, dmg.subDeamDiNuc);
    
    //turn all alignment likelihoods for each path into a vector<vector<int>> for the MCMC input
    vector<vector<double>> probMatrix = convertMapsToVector(gam);
    
    cerr << "Number of paths: " << probMatrix[0].size() << endl;
    cerr << "Number of reads: " << gam->size() << endl;
    
    
    cerr << "Loading tree ... " << endl;
    spidir::Tree taxatree = NULL;  // Initialize taxatree pointer to nullptr
    string treename = treedir+dbprefixS+".new.dnd";
    spidir::readNewickTree(treename.c_str(), &taxatree);
    
    if (taxatree.nodes.size() == 0){
        throw runtime_error("The tree is empty");
    }
    unordered_map <string, int> pathToTreeNode;
    unordered_map <string , string > parentMap;
    int leafcounter = 0;
    double shortestBranch = taxatree.nodes[0]->dist;
    for (int j = 0; j< taxatree.nodes.size(); ++j) {
        if (taxatree.nodes[j]->isLeaf()){
            leafcounter++;
        }
        //cerr << setprecision(14)<< taxatree.nodes[j]->longname << " " << taxatree.nodes[j]->name << " " << taxatree.nodes[j]->dist << endl;
        //cerr << setprecision(14)<< "Parent " << taxatree.nodes[j]->parent->longname << " " << taxatree.nodes[j]->parent->name << " " << taxatree.nodes[j]->parent->dist << endl;
    	pathToTreeNode[taxatree.nodes[j]->longname] = taxatree.nodes[j]->name;
        if (taxatree.nodes[j]->dist < shortestBranch && taxatree.nodes[j]->dist != 0.0){shortestBranch = taxatree.nodes[j]->dist;}

	
        if (taxatree.nodes[j]->parent == NULL){
            parentMap[taxatree.nodes[j]->longname] = taxatree.nodes[j]->longname;
        }else{
            parentMap[taxatree.nodes[j]->longname] = taxatree.nodes[j]->parent->longname;
        }   
    }
    if (shortestBranch != 0 && shortestBranch < 1){
        con = shortestBranch;
    }else{
        con = 0.01;
    }

    cerr << " ... done!" << endl;
    cerr << "Number of tree nodes " << taxatree.nodes.size() << endl;
    if(taxatree.nodes.size() != probMatrix[0].size()){
        throw runtime_error("The number of tree nodes and paths in the graph is unequal. Unable to proceed. Exiting...");
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

    
    freqs['R'] = freqs['A'] + freqs['G'];
    freqs['Y'] = freqs['C'] + freqs['T'];
    freqs['M'] = 1/(2*(   ((22)*(freqs['A']*freqs['G'])) + ((22)*(freqs['C']*freqs['T'])) + (freqs['A']*freqs['C'] + (freqs['A']*freqs['T']) + (freqs['G']*freqs['C'] + (freqs['G'] * freqs['T']) )) ) );

    vector<unsigned int> sigNodes;
    vector<string> sigPaths;

    if(specifiedk){
        if (k > probMatrix[0].size()){
            cerr << "Number for k cannot be larger than the number of tree nodes present. The tree has " << taxatree.nodes.size() << " nodes. Please adjust k accordingly. Exiting..." << endl;
            throw runtime_error("Invalid number of k.");
        }
        cerr << "User specified number of sources was set to k = " << k << ". A random start is being initiated."  << endl;
        sigNodes = soibean::generateRandomNumbers(probMatrix[0].size(), k);
        cerr << "Random starting nodes: ";
        
        for(int i = 0; i < sigNodes.size(); ++i) {
            cerr << sigNodes[i] << " ";
            sigPaths.emplace_back(taxatree.nodes[sigNodes[i]]->longname);
           
        }
        cerr << endl;
    }else{
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
    
        //cout << freqs['M'] << endl;
        vector<pair<string, int>> freqVec(frequencies.begin(), frequencies.end());
        sort(freqVec.begin(),freqVec.end(), [](const auto & a, const auto & b){
            return a.second > b.second;
        });
        
       
        double thres = gam->size()*0.01;
        //cout << thres << endl;
        for (const auto & pair : freqVec)
        {
            if(pair.second >= thres) {
                sigNodes.emplace_back(pathToTreeNode[pair.first]);
                sigPaths.emplace_back(pair.first);
            }

        }

        if(cutk > 0){
            sigNodes.resize(cutk);
            sigPaths.resize(cutk);
        }

        if(sigNodes.empty()){
            for (const auto & pair : freqVec)
            {
                cerr << "No signature node-sets could be identified. Rerunning with no minimum threshold..." << endl;
                sigNodes.emplace_back(pathToTreeNode[pair.first]);
                sigPaths.emplace_back(pair.first);
                
            }
            if (sigNodes.empty()){
                k = 3;
                specifiedk = true;
                cerr << "Still, no signature node-sets could be identified. Initiating the MCMC with k = 3 and random starting nodes." << endl;
            }
        }else{
            for (int sig = 0; sig < sigPaths.size(); ++sig){
                cerr << "Identified signature paths: "<< sigPaths[sig] << " with tree node: " << sigNodes[sig] << endl;
            }
        }

            cerr << "... done! The Identified signature paths are used as input for the MCMC. " << endl;

    
        if (randStart){

            sigNodes = soibean::generateRandomNumbers(probMatrix[0].size(), sigPaths.size());
            cerr << "Random starting nodes: ";
            
            for(int i = 0; i < sigNodes.size(); ++i) {
                cerr << sigNodes[i] << " ";
               
            }
            cerr << endl;

        }
    }


    string num = outputfilename;
    
    for (size_t i = 0; i < sigNodes.size(); ++i) {

        
        std::vector<unsigned int> subVector(sigNodes.begin(), sigNodes.begin() + i + 1);

        double logLike = 0.0L;
        double freq = log(1.0/sigNodes.size());
        if (sigPaths.size() == 1){


            for (const auto &read :*(gam))
            {
                logLike += read->pathMap[sigPaths[i]];
            }
            cerr << "Initial log-likelihood: " << logLike << endl;

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
            cerr << "Initial log-likelihood: " << logLike << endl;

        }
        

        spidir::Tree* taxatreePtr = &taxatree;
        RunTreeProportionParams params(taxatreePtr);
        params.probMatrix = probMatrix;
        params.sources = subVector;
        params.align = gam;
        params.burn = burnin;
        params.maxIter = iter;
        params.chains = chains;
        params.logLike = logLike;
        params.freqs = freqs;

        // Define the number of chains.
        std::vector<std::vector<MCMCiteration>> MCMCiterationsVec(chains);
        
        if (run_mcmc)
        {

            std::ofstream diagnostics;
            // This map will store the vectors of statistics for all chains, with branch names as keys
            unordered_map<string, vector<vector<vector<double>>>> branchStatsMap;
            unsigned int chainIndex = 0;
            vector<double> chainLogLikes;
            // Open diagnostics file for this chain
            diagnostics.open(num + "Diagnostics" + to_string(subVector.size()) + to_string(chainIndex) + ".txt");
            diagnostics << "Source\tHighest log-likelihood\tfor chain\tRhat for the proportion estimate\tRhat for the branch position estimate" << endl;
            // Iterate over each chain
            for (auto& chainVec : MCMCiterationsVec) {

                // Set random sources if not the first chain
                if (chainIndex != 0){
                    params.sources = soibean::generateRandomNumbers(probMatrix[0].size(), subVector.size());
                }

                std::cerr << "Running chain number: " << chainIndex << std::endl;
                // Run the MCMC process for this chain

                std::vector<MCMCiteration> chainiter = mcmc.run_tree_proportion_sb(params, chainVec, graph, nodepaths, num, \
                                                                                   n_threads, path_names.size(), chainIndex, con);

                shared_ptr dta = make_unique<Trailmix_struct>();

                // Process the MCMC iterations to get the statistics map for this chain
                pair<unordered_map<string, vector<vector<double>>>, double> intermStatsMapPair = mcmc.processMCMCiterations_sb(chainiter, subVector.size(), num, \
                                                                                                                            chainIndex, taxatreePtr, leafcounter);

                unordered_map<string, vector<vector<double>>> intermStatsMap = intermStatsMapPair.first;
                chainLogLikes.emplace_back(intermStatsMapPair.second);
                // For each branch name, add the current chain's statistics vector to the main map
                for (const auto& branchStat : intermStatsMap) {
                    const auto& branchName = branchStat.first;
                    const auto& statsForCurrentChain = branchStat.second;

                    // Ensure that we have a vector to store the stats for this branch
                    if (branchStatsMap.find(branchName) == branchStatsMap.end()) {
                        branchStatsMap[branchName] = vector<vector<vector<double>>>();
                    }

                    // Add the stats for the current chain to the corresponding branch in the map
                    branchStatsMap[branchName].emplace_back(statsForCurrentChain);
                    //cerr << "Chain " << chainIndex << " " << branchName << endl;
                }

                chainIndex++;
            }


            // Now, calculate R-hat statistics for each branch across all chains
            int numChains = branchStatsMap.begin()->second.size(); // Assuming all branches have stats for the same number of chains
            int chainLength = iter - burnin; // Assume all chains have the same length

            // Iterate over each branch
            int source_idx = 0;
            for (const auto& branchStat : branchStatsMap) {

                //cerr << "Looking at branch " << branchStat.first << endl;
                const auto& branchName = branchStat.first;
                const auto& allChainStats = branchStat.second;
                if(allChainStats.empty()){throw runtime_error("all Chains vector is empty");}
                std::vector<double> Propmeans(numChains);
                std::vector<double> Propvariances(numChains);
                std::vector<double> Posmeans(numChains);
                std::vector<double> Posvariances(numChains);

                // Collect the statistics for each chain for this branch the branch can have only one name for now

                for (int chain = 0; chain < numChains; ++chain) {
                    if(!allChainStats[chain].empty()){
                    Propmeans[chain] = allChainStats[chain][0][0];
                    Propvariances[chain] = allChainStats[chain][0][1]; 
                    Posmeans[chain] = allChainStats[chain][0][2]; 
                    Posvariances[chain] = allChainStats[chain][0][3]; }
                }
                double maxLogLike = chainLogLikes[0];

                int maxIndex = 0;
                for (int h = 0; h <chainLogLikes.size(); ++h){
                    if(chainLogLikes[h] > maxLogLike){
                        maxLogLike = chainLogLikes[h];
                        maxIndex = h;
                    }
                }

                // Calculate R-hat for proportions
                double PropRhat = soibean::calculateRhat(Propmeans, Propvariances, chainLength);
                if(isnan(PropRhat)){
                    cerr << "Warning: R-hat for the proportion estimate is not computed because we are handling a single source." << endl;
                }
                else if (PropRhat == -1){
                    cerr << "Warning: The R-hat cannot be computed for a single chain." << endl;
                }
                else if(PropRhat > 1.05){
                    std::cerr << "Warning: R-hat for proportion of branch " << branchName << " is above 1.05, indicating that the chains have not converged for the parameter." << std::endl;
                }

                // Calculate R-hat for positions
                double PosRhat = soibean::calculateRhat(Posmeans, Posvariances, chainLength);
                if(isnan(PosRhat)){
                    cerr << "Warning: R-hat cannot be computed. The value for the parameter is identical in every iteration." << endl;
                }
                else if (PosRhat == -1){
                    cerr << "Warning: The R-hat cannot be computed for a single chain." << endl;
                }
                else if(PosRhat > 1.05){
                    std::cerr << "Warning: R-hat for position of branch " << branchName << " is above 1.05, indicating that the chains have not converged for the parameter." << std::endl;
                }

                // Write the R-hat statistics for this branch to the diagnostics file
                diagnostics << branchName << '\t' << maxLogLike << '\t' << maxIndex << '\t' << PropRhat << '\t' << PosRhat << endl;
                source_idx++;
            }
        }

    }

    remove(fifo_A);

return 0;
}


