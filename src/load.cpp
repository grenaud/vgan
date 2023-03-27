#include "HaploCart.h"
#include "Euka.h"
#include "gam2prof.h" // Mikkel code


const vector<double> Haplocart::load_mappabilities(const string &hcfiledir) {
vector<double> mappabilities;
const string mappability_path = getFullPath(hcfiledir+"mappability.tsv");
igzstream myfile;
myfile.open(mappability_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
    {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     const string start_pos = tokens[1];
     const string end_pos = tokens[2];
     const string mappability = tokens[3];
     for (int i=stoi(start_pos); i < stoi(end_pos); ++i) {
         mappabilities.emplace_back(stod(mappability));
     }
    }

return mappabilities;
}


const map<const string, int> Haplocart::load_pangenome_map(const string &hcfiledir){
map<const string, int> pangenome_map;
const string pangenome_map_path = getFullPath(hcfiledir+"parsed_pangenome_mapping");
igzstream myfile;
myfile.open(pangenome_map_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
    {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     const string key = tokens[0];
     const int val = stoi(tokens[1])+1;
     pangenome_map.insert(make_pair(key, val));
    }
return pangenome_map;
}


const vector<string> Haplocart::load_paths(const string &hcfiledir){
vector<string> path_names;
string path_name_path = getFullPath(hcfiledir+"graph_paths");
igzstream myfile;
myfile.open(path_name_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
  {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     string path = tokens[0].substr(0, path.find("."));
     path_names.emplace_back(path);
  }

return path_names;
}


/*
 * Load bin information from file clade_chunk_input.csv
 *
 * This function loads information about the different bins for each of the clades from a file called clade_chunk_input.csv
 * It then populated a vector of vectors of tuples with the bin ranges in node IDs, their entropy score and a counter.
 *
 * @param clade_chunk_path The full path to the file containing the node ID range and entropy scores as well as the name of the clade
 * @return a vector of vectors of tuples
 *
 */
vector<vector<tuple<int, int, double, double > > > Euka::load_clade_chunks(string clade_chunk_path){
    vector<vector<tuple<int, int, double, double > > > chunks;
    vector<tuple<int, int, double, double > > bin_range;
    igzstream myfile;
    myfile.open(clade_chunk_path.c_str(), ios::in);
    string line;
    while ( getline (myfile,line) )
    {
        vector<string> tokens= allTokensWhiteSpaces(line);

        for(int j=1;j<tokens.size();j+=3) {
            tuple<int, int, double, double > minmax = make_tuple(0, 0, 0.0, 0.0);
            get<0>(minmax) = stoi(tokens[j]);
            get<1>(minmax) = stoi(tokens[j+1]);
            get<2>(minmax) = stod(tokens[j+2]);
            get<3>(minmax) = 0.0;
            bin_range.emplace_back(minmax);
        }

        chunks.emplace_back(bin_range);
        bin_range.clear();

    }
    return chunks;
}

/*
 * Load clade information from file clade_info.csv
 *
 * This function loads clade information from a file called clade_info.csv
 * It then populated a vector of Clade pointers which contain the ID+distance+name
 *
 * @param clade_info_path The full path to the file containing the ID, distance and names of the clades
 * @return a pointer to a vector of clade pointers
*/ 
 

vector<Clade *> * Euka::load_clade_info(const string clade_info_path, int lengthToProf){ // Mikkel last argument
    vector<Clade *> * clade_vec = new vector<Clade *>();
    igzstream myfile;
    myfile.open(clade_info_path.c_str(), ios::in);
    string line;
    int** init_array; // Mikkel code
    while ( getline (myfile,line) )
        {
            vector<string> tokens= allTokensWhiteSpaces(line);
            assert(tokens.size() == 3);
            for(unsigned int i=0;i<tokens.size();i++){
                Clade * clade_info = new Clade(0, "", 0.0, 0, {0.0}, {0.0}, {0}, init_array, {""}); //Mikkel code last argument
                //cout<<"tokens["<<i<<"]= "<<tokens[i]<<endl;
                clade_info->id = stoi(tokens[0]);
                clade_info->name = tokens[1];
                clade_info->dist = stod(tokens[2]);
                clade_info->count = 0;
                clade_info->clade_like = {0.0};
                clade_info->inSize = {0};
                clade_info->nameStorage = {""};


                // Mikkel code begins

                // Construct baseshift clade_Array
                int **baseshift_clade_array;
                baseshift_clade_array = new int*[lengthToProf*2];

                for (int p = 0; p < lengthToProf*2; p++){
                    baseshift_clade_array[p] = new int[16];

                    for (int i = 0; i < 16; i++){
                        baseshift_clade_array[p][i] = 0;
                    }
                }; 
                
                clade_info->baseshift_clade_array = baseshift_clade_array;

                // Mikkel code ends

                clade_vec->emplace_back(clade_info);
            }

        }

    return clade_vec;
}

// // // Mikkel code begin - This should be done better

vector<Clade *> * Gam2prof::load_clade_info(const string clade_info_path, int lengthToProf){ // Mikkel last argument
    vector<Clade *> * clade_vec = new vector<Clade *>();
    igzstream myfile;
    myfile.open(clade_info_path.c_str(), ios::in);
    string line;
    int** init_array; // Mikkel code
    while ( getline (myfile,line) )
        {
            vector<string> tokens= allTokensWhiteSpaces(line);
            for(unsigned int i=0;i<tokens.size();i++){
                Clade * clade_info = new Clade(0, "", 0.0, 0, {0.0}, {0.0}, {0}, init_array, {""}); //Mikkel code last argument
                //cout<<"tokens["<<i<<"]= "<<tokens[i]<<endl;
                clade_info->id = stoi(tokens[0]);
                clade_info->name = tokens[1];
                clade_info->dist = stod(tokens[2]);
                clade_info->count = 0;
                clade_info->clade_like = {0.0};
                clade_info->inSize = {0};
                clade_info->nameStorage = {""};


                // Mikkel code begins

                // Construct baseshift clade_Array
                int **baseshift_clade_array;
                baseshift_clade_array = new int*[lengthToProf*2];

                for (int p = 0; p < lengthToProf*2; p++){
                    baseshift_clade_array[p] = new int[16];

                    for (int i = 0; i < 16; i++){
                        baseshift_clade_array[p][i] = 0;
                    }
                }; 
                
                clade_info->baseshift_clade_array = baseshift_clade_array;

                // Mikkel code ends

                clade_vec->emplace_back(clade_info);
            }

        }

    return clade_vec;
}

vector<vector<tuple<int, int, double, double > > > Gam2prof::load_clade_chunks(string clade_chunk_path){
    vector<vector<tuple<int, int, double, double > > > chunks;
    vector<tuple<int, int, double, double > > bin_range;
    igzstream myfile;
    myfile.open(clade_chunk_path.c_str(), ios::in);
    string line;
    while ( getline (myfile,line) )
    {
        vector<string> tokens= allTokensWhiteSpaces(line);

        for(int j=1;j<tokens.size();j+=3) {
            tuple<int, int, double, double > minmax = make_tuple(0, 0, 0.0, 0.0);
            get<0>(minmax) = stoi(tokens[j]);
            get<1>(minmax) = stoi(tokens[j+1]);
            get<2>(minmax) = stod(tokens[j+2]);
            get<3>(minmax) = 0.0;
            bin_range.emplace_back(minmax);
        }

        chunks.emplace_back(bin_range);
        bin_range.clear();

    }
    return chunks;
}

// // Mikkel code end



const vector<vector<bool>> Euka::load_path_supports_Euka(const string & pathsupportfile) {
    vector<vector<bool>> path_supports(6925366);
    const string path_support_path = pathsupportfile;
    igzstream myfile;
    myfile.open(path_support_path.c_str(), ios::in);
    string line;
    int index = 0;
    while (getline(myfile, line))
    {
        for (int j=0; j<line.size(); ++j) {
        char supported = line[j];
        if (supported == '1') {path_supports[index].emplace_back(true);}
        else {path_supports[index].emplace_back(false);}
        }
        index += 1;
    }
    return path_supports;
}


const vector<vector<bool>> Haplocart::load_path_supports(const string &hcfiledir) {
vector<vector<bool>> path_supports(11825);
const string path_support_path = getFullPath(hcfiledir+"path_supports");
igzstream myfile;
myfile.open(path_support_path.c_str(), ios::in);
string line;
int index = 0;
while (getline(myfile, line))
  {
     for (int j=0; j<5179; ++j) {
         char supported = line[j];
         if (supported == '1') {path_supports[index].emplace_back(true);}
         else {path_supports[index].emplace_back(false);}
                                 }
   index += 1;
  }
return path_supports;
}


const map<string, vector<string>> Haplocart::load_parents(const string &hcfiledir) {
map<string, vector<string>> parents;
igzstream myfile;
myfile.open(getFullPath(hcfiledir+"parents.txt").c_str(), ios::in);
string line;
while (getline(myfile, line)) {
    const vector<string> tokens= allTokensWhiteSpaces(line);
    if (tokens.size() == 0) {continue;}
    string node = tokens[0];
    vector<string> parent_vec;
    for (int j=1; j<tokens.size(); ++j) {
        if (tokens[j].find('[') == std::string::npos) {
            parent_vec.emplace_back(tokens[j]);
                                                      }
                                        }
    pair<string, vector<string>> parent_pair = make_pair(node, parent_vec);
    parents.insert(parent_pair);
                              }

return parents;
}

const map<string, vector<string>> Haplocart::load_children(const string &hcfiledir) {
map<string, vector<string>> children;
igzstream myfile;
myfile.open(getFullPath(hcfiledir+"children.txt").c_str(), ios::in);
string line;
while (getline(myfile, line)) {
    const vector<string> tokens= allTokensWhiteSpaces(line);
    if (tokens.size() == 0) {continue;}
    string node = tokens[0];
    vector<string> child_vec;
    for (int j=1; j<tokens.size(); ++j) {
        if (tokens[j].find('[') == std::string::npos) {
            child_vec.emplace_back(tokens[j]);
                                                      }
                                        }
    const pair<string, vector<string>> child_pair = make_pair(node, child_vec);
    children.insert(child_pair);
                              }

return children;
}
