#include "TrailMix.h"
#include "HaploCart.h"
#include "Euka.h"
#include "gam2prof.h" // Mikkel code

vector<string> split_string(const string &input, char delimiter) {
    stringstream ss(input);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delimiter)) {
        tokens.emplace_back(item);
    }
    return tokens;
}

string decompress_gz_file(const string &filepath) {
    gzFile file = gzopen(filepath.c_str(), "rb");
    if (!file) {
        throw runtime_error("Failed to open the file for decompression.");
    }

    stringstream decompressed_data;
    char buffer[4096];
    int bytes_read;
    while ((bytes_read = gzread(file, buffer, sizeof(buffer))) > 0) {
        decompressed_data.write(buffer, bytes_read);
    }

    gzclose(file);
    return decompressed_data.str();
}

bool is_number(const std::string& s) {
    return !s.empty() && std::find_if(s.begin(),
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

void Trailmix::load_read_probs(std::shared_ptr<Trailmix_struct>& dta) {
    std::string read_prob_file = dta->tmpdir + "rpvg_ht_probs.txt.gz";
    
    if (!filesystem::exists(read_prob_file)) {
        throw std::runtime_error("[TrailMix] Error: RPVG read-level output file not found.");
    }

    std::vector<double> probs;
    std::vector<std::vector<unsigned int>> idxs;

    std::string line;
    std::string decompressed_data = decompress_gz_file(read_prob_file);
    std::stringstream ss(decompressed_data);
    
    if (!std::getline(ss, line)) {
        throw std::runtime_error("[TrailMix] Error: Failed to read header from RPVG file.");
    }
    
    if (!std::getline(ss, line)) {
        throw std::runtime_error("[TrailMix] Error: Failed to read data from RPVG file.");
    }

    dta->RPVG_hap_names = split(line, ' ');
    
    for (auto& hn : dta->RPVG_hap_names) {
        size_t pos = hn.find('_');
        if (pos != std::string::npos) {
            hn = hn.substr(0, pos);
        }
    }

    int counter = 0;
    while (std::getline(ss, line)) {
        ++counter;
        std::vector<unsigned int> inner_idxs;
        std::vector<std::string> tokens = split_string(line, ' ');

        if (tokens.size() < 3 || !is_number(tokens[0])) {
            continue;
        }

        try {
            unsigned int multiplicity = std::stoul(tokens[0]);
            for (size_t i = 2; i < tokens.size(); ++i) {
                size_t colon_pos = tokens[i].find(':');
                if (colon_pos == std::string::npos) {
                    continue;
                }

                double prob = std::stod(tokens[i].substr(0, colon_pos)) * multiplicity;
                std::vector<std::string> inner_tokens = split(tokens[i].substr(colon_pos + 1), ',');
                for (const std::string& token : inner_tokens) {
                    inner_idxs.emplace_back(std::stoul(token));
                }
                probs.emplace_back(prob);
            }
            idxs.emplace_back(inner_idxs);
        } catch (const std::invalid_argument& ia) {
            throw std::runtime_error("[TrailMix] Error: Invalid argument - " + std::string(ia.what()));
        } catch (const std::out_of_range& oor) {
            throw std::runtime_error("[TrailMix] Error: Out of range - " + std::string(oor.what()));
        } catch (...) {
            throw std::runtime_error("[TrailMix] Error: Unknown exception occurred.");
        }
    }

    dta->read_probs = std::make_pair(probs, idxs);
}

void Trailmix::load_hap_combos(shared_ptr<Trailmix_struct> &dta) {

//cerr << "temp dir (load hap combos): " << dta->tmpdir << endl;

if (!std::filesystem::exists(dta->tmpdir+"rpvg_hap.txt")){
    throw runtime_error("[TrailMix] Error, no RPVG output was produced. This is probably due to very low confidence in placing the data. It may be worthwhile to consider alternative values for the number of contributing sources k.\n");
                                                         }

    dta->hap_combos.clear();
    const string hap_combo_path = getFullPath(dta->tmpdir + "rpvg_hap.txt");
    igzstream myfile;
    myfile.open(hap_combo_path.c_str(), ios::in);
    string line;
    getline(myfile, line);
    while (getline(myfile, line)) {
        vector<string> tokens = allTokensWhiteSpaces(line);
        vector<string> newtokens;
        for (const auto &token : tokens) {
            // Handling specific token formats
            if (token.starts_with("_gbwt_ref_")) {
                // Ensure that this manipulation doesn't alter formats like "MT576650.1"
                // Modify the logic here based on the expected format of your tokens
                // Example: Extract substring but preserve specific formats
                string modifiedToken = token.substr(10);
                size_t pos = modifiedToken.find("_0_0");
                if (pos != string::npos) {
                    modifiedToken = modifiedToken.substr(0, pos);
                }
                newtokens.emplace_back(modifiedToken);
            } else if (token != tokens.back()) {
                // Check and preserve specific formats like "MT576650.1"
                size_t dotPos = token.find('.');
                if (dotPos != string::npos) {
                    newtokens.emplace_back(token.substr(0, dotPos + 2)); // Preserving the format up to ".1"
                } else {
                    newtokens.emplace_back(token);
                }
            } else {
                newtokens.emplace_back(token);
            }
        }
        dta->hap_combos.emplace_back(newtokens);
    }

    return;
}

void Trailmix::load_tpms(shared_ptr<Trailmix_struct> &dta) {
    if (!std::filesystem::exists(dta->tmpdir+"rpvg_ht.txt")){
        throw runtime_error("[TrailMix] Error, no RPVG output was produced.");
                                                            }
    // Open the file for reading
    std::ifstream file(dta->tmpdir + "rpvg_ht.txt");
    // Check if the file was successfully opened
    if (!file.is_open()) {
        throw std::runtime_error("[TrailMix] Error opening RPVG file");
    }
    // Read the file line by line
    std::string line;
    getline(file, line);
    while (std::getline(file, line)) {
    // Split the line into tokens
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, '\t')) {
        tokens.emplace_back(token);
    }
    if (tokens[0] == "Unknown"){break;}
    // Parse the tokens and add them to dta->tpms
    if (tokens.size() == 7) {
        std::string name = tokens[0].substr(10).substr(0, tokens[0].substr(10).find("_0_0")).substr(0, tokens[0].substr(10).find('.'));
        Trailmix::modifyPathNameInPlace(dta, name);
        const double tpm = std::stod(tokens[6]);
        dta->tpms.emplace_back(std::make_pair(name, tpm));
    }
     }
    // Close the file
    file.close();
    return;
                                                                        }


void Haplocart::load_path_names(shared_ptr<Trailmix_struct> &dta){
string path_name_path = getFullPath(dta->graph_dir+"graph_paths");
igzstream myfile;
myfile.open(path_name_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
  {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     string path = tokens[0];
     dta->path_names.emplace_back(path.substr(0, path.find('.')));
  }

if (dta->path_names.empty()){throw runtime_error("[HaploCart] No path name file found");}
return;
}

void Haplocart::load_mappabilities(std::shared_ptr<Trailmix_struct> &dta) {
    std::string mappability_path = getFullPath(dta->graph_dir + "mappability.tsv");
    igzstream myfile;
    myfile.open(mappability_path.c_str(), ios::in);
    std::string line;
    while (std::getline(myfile, line)) {
        std::vector<std::string> tokens = allTokensWhiteSpaces(line);
        if (tokens.size() != 4) {
            throw std::runtime_error("Invalid line in mappability file: " + line);
        }
        int start_pos = std::stoi(tokens[1]);
        int end_pos = std::stoi(tokens[2]);
        const double mappability = max(0.1, std::stod(tokens[3]));
        if (start_pos >= end_pos) {
            throw std::runtime_error("Invalid start/end positions in mappability file: " + line);
        }
        dta->mappabilities.resize(end_pos, mappability);
        std::fill(dta->mappabilities.begin() + start_pos, dta->mappabilities.begin() + end_pos, mappability);
        assert(!dta->mappabilities.empty());
    }
    return;
}


void Haplocart::load_pangenome_map(shared_ptr<Trailmix_struct> &dta){
const string pangenome_map_path = getFullPath(dta->graph_dir+"parsed_pangenome_mapping");
igzstream myfile;
myfile.open(pangenome_map_path.c_str(), ios::in);
string line;
while (getline(myfile, line))
    {
     const vector<string> tokens= allTokensWhiteSpaces(line);
     const string key = tokens[0];
     const int val = stoi(tokens[1])+1;
     dta->pangenome_map.insert(make_pair(key, val));
    }
assert(!dta->pangenome_map.empty());
return;
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
    unsigned int** init_array; // Mikkel code
    while ( getline (myfile,line) )
        {
            vector<string> tokens= allTokensWhiteSpaces(line);
            assert(tokens.size() == 6);
            for(unsigned int i=0;i<tokens.size();i++){
                Clade * clade_info = new Clade(0, "", 0.0, 0, 0, 0, 0, {0.0}, {0.0}, {0}, init_array, {""}); //Mikkel code last argument
                //cout<<"tokens["<<i<<"]= "<<tokens[i]<<endl;
                clade_info->id = stoi(tokens[0]);
                clade_info->name = tokens[1];
                clade_info->dist = stod(tokens[2]);
                clade_info->noPaths = stoi(tokens[3]);
                clade_info->snode = stoi(tokens[4]);
                clade_info->enode = stoi(tokens[5]);
                clade_info->count = 0;
                clade_info->clade_like = {0.0};
                clade_info->inSize = {0};
                clade_info->nameStorage = {""};


                // Mikkel code begins

                // Construct baseshift clade_Array
                unsigned int **baseshift_clade_array;
                baseshift_clade_array = new unsigned int*[lengthToProf*2];

                for (int p = 0; p < lengthToProf*2; p++){
                    baseshift_clade_array[p] = new unsigned int[16];

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

// Overload for boost test

//vector<Clade *> * Euka::load_clade_info(const string &clade_info_path, const int lengthToProf, const bool){
//    return Euka::load_clade_info(clade_info_path,lengthToProf);
//                                                                                                          }

// // // Mikkel code begin - This should be done better
// vector<Clade *> * Gam2prof::load_clade_info(const string &clade_info_path, const int lengthToProf){ // Mikkel last argument
//     vector<Clade *> * clade_vec = new vector<Clade *>();
//     igzstream myfile;
//     myfile.open(clade_info_path.c_str(), ios::in);
//     string line;
//     unsigned int** init_array; // Mikkel code
//     while ( getline (myfile,line) )
//         {
//             vector<string> tokens= allTokensWhiteSpaces(line);
//             for(unsigned int i=0;i<tokens.size();i++){
//                 Clade * clade_info = new Clade(0, "", 0.0, 0, {0.0}, {0.0}, {0}, init_array); //Mikkel code last argument
//                 //cout<<"tokens["<<i<<"]= "<<tokens[i]<<endl;
//                 clade_info->id = stoi(tokens[0]);
//                 clade_info->name = tokens[1];
//                 clade_info->dist = stod(tokens[2]);
//                 clade_info->count = 0;
//                 clade_info->clade_like = {0.0};
//                 clade_info->inSize = {0};

//                 // Mikkel code begins

//                 // Construct baseshift clade_Array
//                 unsigned int **baseshift_clade_array;
//                 baseshift_clade_array = new unsigned int*[lengthToProf*2];

//                 for (int p = 0; p < lengthToProf*2; p++){
//                     baseshift_clade_array[p] = new unsigned int[16];

//                     for (int i = 0; i < 16; i++){
//                         baseshift_clade_array[p][i] = 0;
//                     }
//                 };
//                 clade_info->baseshift_clade_array = baseshift_clade_array;

//                 // Mikkel code ends

//                 clade_vec->emplace_back(clade_info);
//             }

//         }

//     return clade_vec;
// }




// const vector<vector<bool>> Euka::load_path_supports_Euka(const string & pathsupportfile) {
//     vector<vector<bool>> path_supports(6271738);
//     const string path_support_path = pathsupportfile;
//     igzstream myfile;
//     myfile.open(path_support_path.c_str(), ios::in);
//     string line;
//     int index = 0;
//     while (getline(myfile, line))
//     {
//         for (int j=0; j<line.size(); ++j) {
//         char supported = line[j];
//         if (supported == '1') {path_supports[index].emplace_back(true);}
//         else {path_supports[index].emplace_back(false);}
//         }
//         index += 1;
//     }
//     return path_supports;
// }



vector<vector<tuple<int, int, double, double > > > Gam2prof::load_clade_chunks(const string &clade_chunk_path){
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


void Haplocart::load_path_supports(shared_ptr<Trailmix_struct> &dta) {
vector<vector<bool>> path_supports(15725);
const string path_support_path = getFullPath(dta->graph_dir+"path_supports");
igzstream myfile;
myfile.open(path_support_path.c_str(), ios::in);
string line;
int index = 0;
while (getline(myfile, line))
  {
     for (int j=0; j<12852; ++j) {
         char supported = line[j];
         if (supported == '1') {path_supports[index].emplace_back(true);}
         else {path_supports[index].emplace_back(false);}
                                 }
   index += 1;
  }
dta->path_supports = move(path_supports);
}

void Haplocart::load_parents(shared_ptr<Trailmix_struct> &dta) {
igzstream myfile;
myfile.open(getFullPath(dta->graph_dir+"parents.txt").c_str(), ios::in);
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
    dta->parents.insert(parent_pair);
                              }
assert(!dta->parents.empty());
return;
}

void Haplocart::load_children(shared_ptr<Trailmix_struct> &dta) {
igzstream myfile;
myfile.open(getFullPath(dta->graph_dir+"children.txt").c_str(), ios::in);
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
    dta->children.insert(child_pair);
                              }
assert(!dta->children.empty());
return;
}
