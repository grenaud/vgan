#include "gam2prof.h"
#include "bdsg/odgi.hpp"


Gam2prof::Gam2prof(){

}

Gam2prof::~Gam2prof(){

}


const string Gam2prof::usage() const{

    return string(string("") +

        "\nvgan gam2prof [options] <reads.gam>\n"+
        "\nGam2prof is a tool that reads a GAM file and produces a deamination profile for the 5' and 3' ends.\n"+

        "\n\n"+
        "Input options:\n"+
        "\t"+""  +"" +"--euka_dir [STR]" + "\t" + "Euka database location (default: \"../share/euke_dir/\")\n"+  
        "\t"+""  +"" +"--dbprefix [STR]" + "\t" + "Database prefix name (defualt: all)\n"+
        "\t"+""  +"" +"-l [INT]"          +"\t\t" + "Set length for substitution matrix (default: 5)"+"\n"+ 
        "\t"+""  +"" +"-out_dir [STR]"    +"\t\t" + "Path for output prof-file(s) (default: current location, '/dev/stdout' for stdout)"+"\n"+ 
        "\t"+""  +"" +"-prof [STR]"       +"\t\t" + "Set which determination profiles you want: (default: both)"+"\n"+
                                                    "\t\t\t\t\t{both, all, 5p3p, 5p, 3p}"+"\n"+ 
                                                    "\t\t\t\t\tboth: print everything in one file"+"\n"+ 
                                                    "\t\t\t\t\tall:  print everything in one file, plus 5p in one file and 3p in one file"+"\n"+
                                                    "\t\t\t\t\t5p3p: print 5p in one file and 3p in one file"+"\n"+ 
                                                    "\t\t\t\t\t5p:   print 5p in one file"+"\n"+  
                                                    "\t\t\t\t\t3p:   print 3p in one file"+"\n"+ 
        "\t"+""  +"" +"-count "         +"\t\t" + "If stated: print basecounts for each position to out_dir"+"\n"+ 
    ""
    );
}

const int Gam2prof::run(int argc, char *argv[], const string cwdProg){


    // Check option
    bool euka_dirspecified = false; 
    string euka_dir = "../share/euka_dir/"; 
    string dbprefix = "all";
    int lengthToProf = 5;
    string prof_out_file_path = cwdProg;
    string print_ends = "both";
    vector<string> options = {"both", "all", "5p3p", "5p", "3p"};
    bool print_counts_to_user = false;

    for(int i=1;i<(argc);i++){

        if(string(argv[i]) == "--euka_dir"){
            euka_dir = argv[i+1];
            euka_dirspecified=true; 
            if(euka_dir.back() != '/'){euka_dir += '/';}
            continue;
        }

        if(string(argv[i]) == "--dbprefix"){
            dbprefix = argv[i+1];
            continue;
        }

        if(string(argv[i]) == "-l"){
            lengthToProf = stoi(argv[i+1]);
            continue;
        }
        
        if(string(argv[i]) == "-out_dir"){
            prof_out_file_path = string(argv[i+1]);
            continue;
        }

        if(string(argv[i]) == "-prof"){
            
            if (std::find(options.begin(), options.end(), string(argv[i+1])) != options.end() ) {
                print_ends = string(argv[i+1]);
            } else {
                cout << "Cannot recognize '" << string(argv[i+1]) << "' as input." << endl;
                cout << "Use either 'both', 'all' '5p3p', '5p' or '3p'." << endl;
                return 1;
            }
            continue;
        }

        if(string(argv[i]) == "-count"){
            print_counts_to_user = true;
            continue;
        }

        
    }

    string gamfilename;
    // Get gamfile name
    if  (string(argv[argc-1])[0] == '/') {
        gamfilename  = string( argv[argc-1]);
    } else  {
        gamfilename  = cwdProg + string( argv[argc-1]);
    }
    

    if(!euka_dirspecified){
        euka_dir = cwdProg + euka_dir;
    }    

    // Get database names
    dbprefix              = euka_dir + dbprefix;
    string ogfilename     = dbprefix+".og";
    string cladefilename  = dbprefix+".clade";
    string binsfilename   = dbprefix+".bins";

    // Check if files exits
    if (isFile(gamfilename) == false){
        throw(std::runtime_error(gamfilename + " does not exist."));
    } else if (isFile(ogfilename) == false){
        throw(std::runtime_error(ogfilename + " does not exist."));
    } else if (isFile(cladefilename) == false){
        throw(std::runtime_error(cladefilename + " does not exist."));
    } else if (isFile(binsfilename) == false){
        throw(std::runtime_error(binsfilename + " does not exist."));
    } 

    // Handle path graph
    bdsg::ODGI graph;
    graph.deserialize(ogfilename);

    //read information about different clades
    vector<Clade *> * clade_vec = load_clade_info(cladefilename, lengthToProf);

    //read information about bin complexity
    vector<vector<tuple<int, int, double, double >>>  chunks = load_clade_chunks(binsfilename);



    
    int n_reads=0;
    int n_map_reads=0;

    
    function<void(vg::Alignment&)> lambda = [&graph, &n_reads, &n_map_reads, &clade_vec,  \
         &chunks,  &lengthToProf, &prof_out_file_path](vg::Alignment& a) { 

        ++n_reads;


        // if Identity is 0 the read is unmapped.
        if (a.identity() != 0){
            
            ++n_map_reads; // we count both the mapped and the unmapped reads

            // get the first node id for the mapped read
            int n_index = a.path().mapping()[0].position().node_id();

            // the node ID for the first node of the mapping is saved as n_index
            // the n_index is used to find the clade ID by looping through the saved bin indexes. Every clade has their own bins defind by the node IDs.
            // As soon as the fist node ID of the mapping is found within the bins the number of the vector is returned, which corresponds to the clade ID.
            // The clade ID (c_n) is then used to get the clade name, clade pairwise average distance, etc.
            int c_n = 0;

            for (long unsigned int i=0; i<chunks.size(); i++){

                for (long unsigned int j=0; j<chunks.at(i).size();++j){

                    //adds a count to the bin the read falls into for the specific clade.
                    if (n_index == std::clamp(n_index, get<0>(chunks.at(i).at(j)), get<1>(chunks.at(i).at(j)))) {

                    //cerr << "This is the clade id: "<< i << endl;
                        c_n = i;

                    }else{
                        continue;
                    }
                }
            }


            // with the clade id access the pairwise average distance computed for each clade
            double pair_dist = clade_vec->at(c_n*3 +1)->dist;  // *3+1 is to index the vector of pointers (every clade has 3 entries (*3) but the counting starts at 0)

            // reconstructing the graph and the read sequence
            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());


            // Init baseshift data array in clade_vec
            int** baseshift_data_array_location = clade_vec->at(c_n*3+1)->baseshift_clade_array; 
        
            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);

            // Calculate baseshifts
            baseshift_data_array.baseshift_calc(graph_seq, read_seq);
            

            // going through the mapping again to count all node ids that the read hits.
            for (long unsigned int i=0; i<a.path().mapping().size(); ++i) {
                // node id to look for c_n is the clade id to get the right bins
                int n_id = a.path().mapping()[i].position().node_id();

                for (long unsigned int j=0; j<chunks.at(c_n).size();++j){
                    //adds a count to the bin the read falls into for the specific clade.
                    if (n_id == std::clamp(n_id, get<0>(chunks.at(c_n).at(j)), get<1>(chunks.at(c_n).at(j)))) {
                    // we are adding a count for each node id that a mapping hits, there can be multiple mappings per read
                    // therefore multiple nodes can be hit and more than one count can be issues per read
                        std::cerr.precision(20);
                        get<3>(chunks.at(c_n).at(j)) += 1;
                    }
                }
            }
             
        
        } //END if identify not zero
    
    }; //END lambda function

    // run lamba function on gamfile
    vg::get_input_file(gamfilename, [&](istream& in) { vg::io::for_each(in, lambda);    });


    // check if clade could be detected or not
    vector<int> clade_list_id;
    int k1 =0; // test code
    int k2 = 0; // test code
    int k3 = 0; // test code
    for (int i = 0; i< chunks.size(); i++) {

        vector<int> check_for_zero;
        for (int k = 0; k<chunks[i].size(); k++){
            check_for_zero.emplace_back(get<3>(chunks[i].at(k)));
            //cout << get<3>(chunks[i].at(k)) << endl; // test code
            k1++; // test code
        }
        if (std::count(check_for_zero.begin(), check_for_zero.end(), 0.0)){
            k2++; // test code
            continue;
        } else {
            // creating a list of ids for detected clades
            clade_list_id.emplace_back(clade_vec->at(i*3+1)->id);
            k3++; // test code
        }
    }
    //cerr << "clade_list_id size: " << clade_list_id.size() << endl; // test code
    //cerr << "k1: " << k1 << endl; // test code
    //cerr << "k2: " << k2 << endl; // test code
    //cerr << "k3: " << k3 << endl; // test code
    
    
    // check if the program should write output files or to stdout
    string prof_out_file = cwdProg;
    string count_out_file = cwdProg;
    if (prof_out_file_path != cwdProg && prof_out_file_path != "/dev/stdout") {
        // Check for correct end of file-path
        if (prof_out_file_path.back() != '/') {
            prof_out_file_path = prof_out_file_path + "/";
        }
        // check if out-folder exits, if not then create folder
        if (!fs::is_directory(prof_out_file_path) || !fs::exists(prof_out_file_path)) { // Check if src folder exists
            fs::create_directory(prof_out_file_path); // create src folder
        }
    }
        
    // Print the clades that have been found
    for (int i = 0; i < clade_list_id.size(); i++){

        // Find loaction of baseshift data array in clade_vec
        int** baseshift_data_array_location = clade_vec->at(clade_list_id[i]*3+1)->baseshift_clade_array; 
        
        // init baseshift data array in baseshift class
        Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);
        
        // make filename for clades
        if (prof_out_file_path != "/dev/stdout"){
            prof_out_file = prof_out_file_path + clade_vec->at(clade_list_id[i]*3+1)->name + ".prof";
            count_out_file = prof_out_file_path + clade_vec->at(clade_list_id[i]*3+1)->name + ".count";
        } else {
            prof_out_file = "/dev/stdout";
            count_out_file = "/dev/stdout";
        }
        


        // Display substitution matrix for clade

        if (print_ends == "both") {
            baseshift_data_array.print_prof(prof_out_file);
        }
        else if (print_ends == "all") {
            baseshift_data_array.print_prof(prof_out_file);
            baseshift_data_array.print_prof(prof_out_file, "5p");
            baseshift_data_array.print_prof(prof_out_file, "3p");
        }
        else if (print_ends == "5p3p") {
            baseshift_data_array.print_prof(prof_out_file, "5p");
            baseshift_data_array.print_prof(prof_out_file, "3p");
        }
        else if (print_ends == "5p") {
            baseshift_data_array.print_prof(prof_out_file, "5p");
        }
        else if (print_ends == "3p") {
            baseshift_data_array.print_prof(prof_out_file, "3p");
        }

        if (print_counts_to_user) {
            baseshift_data_array.print_counts(count_out_file);
        }
        

    }

    cout << "Outfiles can be found in " << prof_out_file_path << endl;

    return 0;

}

