#include "gam2prof.h"
#include "bdsg/odgi.hpp"
#include "Trailmix_struct.h"

static shared_ptr<Trailmix_struct> NULLPTR_TRAILMIX_STRUCT = make_unique<Trailmix_struct>();

Gam2prof::Gam2prof(){

}

Gam2prof::~Gam2prof(){

}


const string Gam2prof::usage() const{

    return string(string("") +"\nvgan gam2prof [options] <database prefix> <reads.gam>\n"+
           "\n"+
                  "\n\n"+
                  //"Options:\n"+
                  "\t\t"+""  +"" +"-l"   +"\t\t" + "Set length for substitution matrix (default: 5)"+"\n"+
                  "\t\t"+""  +"" +"-p"   +"\t\t" + "Path for output prof-file (default: stdout)"+"\n"+
          ""
          );
}

const int Gam2prof::run(int argc, char *argv[], const string &cwdProg){

    // Check options
    int lengthToProf = 5;
    string prof_out_file_path = "/dev/stdout";

    for(int i=1;i<(argc);i++){

        if(string(argv[i]) == "-l"){
            lengthToProf = stoi(argv[i+1]);
            continue;
        }

        if(string(argv[i]) == "-p"){
            prof_out_file_path = string(argv[i+1]);
            continue;
        }

    }

    // Get file names
    const string gamfilename  = cwdProg + string(argv[argc-1]);
    const string dbprefix    = cwdProg + string(argv[argc-2]);
    const string ogfilename     = dbprefix+".og";
    const string cladefilename  = dbprefix+".clade";
    const string binsfilename   = dbprefix+".bins";

    // Check if files exits
    if (isFile(gamfilename) == false){throw(std::runtime_error(gamfilename + " does not exist."));}
    if (isFile(ogfilename) == false){throw(std::runtime_error(ogfilename + " does not exist."));}
                       


    if (isFile(cladefilename) == false){throw(std::runtime_error(cladefilename + " does not exist."));}
    if (isFile(binsfilename) == false){throw(std::runtime_error(binsfilename + " does not exist."));}
                      

    // Deserialize handlegraph
    bdsg::ODGI graph;


    graph.deserialize(ogfilename);
    

    Euka ek;

    //read information about different clades
    vector<Clade *> * clade_vec = ek.load_clade_info(cladefilename, lengthToProf);

    //read information about bin complexity
    vector<vector<tuple<int,int,double,double>>> chunks = load_clade_chunks(binsfilename);


    int n_reads=0;
    int n_map_reads=0;


    function<void(vg::Alignment&)> lambda = [&graph, &n_reads, &n_map_reads, &clade_vec,  \
         &chunks,  &lengthToProf, &prof_out_file_path](vg::Alignment& a) {

        ++n_reads;


        // if Identity is 0 the read is unmapped.
        if (a.identity() != 0){
            ++n_map_reads; // we count both the mapped and the unmapped reads

            // get the first node id for the mapped read
            const int n_index = a.path().mapping()[0].position().node_id();

          // reconstructing the graph and the read sequence
          auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());

            // the node ID for the first node of the mapping is saved as n_index
            // the n_index is used to find the clade ID by looping through the saved bin indexes. Every clade has their own bins defind by the node IDs.
            // As soon as the fist node ID of the mapping is found within the bins the number of the vector is returned, which corresponds to the clade ID.
            // The clade ID (c_n) is then used to get the clade name, clade pairwise average distance, etc.
            int c_n = 0;

            for (size_t i=0; i<chunks.size(); ++i){
                for (size_t j=0; j<chunks.at(i).size();++j){
                    //adds a count to the bin the read falls into for the specific clade.
                    if (n_index == std::clamp(n_index, get<0>(chunks.at(i).at(j)), get<1>(chunks.at(i).at(j)))) {c_n = i;}
                    else{continue;}
                                                           }
                                                  }


            // with the clade id access the pairwise average distance computed for each clade
            const double pair_dist = clade_vec->at(c_n*6 +1)->dist;  // *3+1 is to index the vector of pointers (every clade has 3 entries (*3)
                                                                     //  but the counting starts at 0)

            // Init baseshift data array in clade_vec
            unsigned int** baseshift_data_array_location = clade_vec->at(c_n*6+1)->baseshift_clade_array;

            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);

            // Calculate baseshifts
            baseshift_data_array.baseshift_calc(graph_seq, read_seq);

            // going through the mapping again to count all node ids that the read hits.
            for (size_t i=0; i<a.path().mapping().size(); ++i) {
                // node id to look for c_n is the clade id to get the right bins
                const int n_id = a.path().mapping()[i].position().node_id();

                for (size_t j=0; j<chunks.at(c_n).size();++j){
                    //adds a count to the bin the read falls into for the specific clade.
                    if (n_id == std::clamp(n_id, get<0>(chunks.at(c_n).at(j)), get<1>(chunks.at(c_n).at(j)))) {
                    // we are adding a count for each node id that a mapping hits, there can be multiple mappings per read
                    // therefore multiple nodes can be hit and more than one count can be issues per read
                        std::cerr.precision(20);
                        get<3>(chunks.at(c_n).at(j)) += 1;
                                                                                                              }
                                                             }
               }
             // END EUKA-SPECIFIC CODE

        } //END if identify not zero
    };



        vg::get_input_file(gamfilename,[&](istream& in) {vg::io::for_each(in, lambda);});
                           
                              

    // check if the program should write output files or to stdout
    string prof_out_file = "/dev/stdout";
    string all_out_file = "/dev/stdout";
    if (prof_out_file_path != "/dev/stdout") {
        // Check for correct end of file-path
        if (prof_out_file_path.back() != '/') {
            prof_out_file_path = prof_out_file_path + "/";
        }
        // check if out-folder exits, if not then create folder
        if (!fs::is_directory(prof_out_file_path) || !fs::exists(prof_out_file_path)) { // Check if src folder exists
            fs::create_directory(prof_out_file_path); // create src folder
                                                                                      }
                                              }


/////////////////////////////////////////////////// FOR EUKA ////////////////////////////////////////////////////////////


    // check if clade could be detected or not
    vector<int> clade_list_id;
    for (size_t i=0; i<chunks.size(); i++) {
        vector<int> check_for_zero;
        for (size_t k = 0; k<chunks[i].size(); k++){
            check_for_zero.emplace_back(get<3>(chunks[i].at(k)));
                                                   }
        if (std::count(check_for_zero.begin(), check_for_zero.end(), 0.0)){continue;}
        else {
            // creating a list of ids for detected clades
            clade_list_id.emplace_back(clade_vec->at(i*6+1)->id);
             }
                                              }

    // Print the clades that have been found
    for (size_t i=0; i<clade_list_id.size(); ++i){

        // Find loaction of baseshift data array in clade_vec
        unsigned int** baseshift_data_array_location = clade_vec->at(clade_list_id[i]*6+1)->baseshift_clade_array;

        // init baseshift data array in baseshift class
        Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);

        // make filename for clades of outfiles if chosen
        if (prof_out_file_path != "/dev/stdout") {
            // Make filenames
            prof_out_file = prof_out_file_path + clade_vec->at(clade_list_id[i]*6+1)->name + ".prof";
            all_out_file = prof_out_file_path + clade_vec->at(clade_list_id[i]*6+1)->name + ".all";
                                                 }

        // Display substitution matrix for clade
        baseshift_data_array.display_prof(prof_out_file);

        // Display entire baseshift array for clade
        cerr << '\n';
        baseshift_data_array.display_counts(all_out_file);

                                                }


//////////////////////////////////////////////// END FOR EUKA ////////////////////////////////////////


    if (prof_out_file_path != "/dev/stdout") {
        cerr << "Outfiles can be found in " << prof_out_file_path << endl;
                                             }

    return 0;

}
