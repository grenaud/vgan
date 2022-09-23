#include "HaploCart.h"

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
