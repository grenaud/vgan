#include "HaploCart.h"
#include "TrailMix.h"
#include <fstream>
#include <iostream>

const bool contains_only_valid_bases(const string &str) {
    return str.find_first_not_of("actgwsmkrybdhvnACTGWSMKRYBDHVN") ==
        std::string::npos;
}

const char get_invalid_base(const string &str) {
    return str[str.find_first_not_of("actgwsmkrybdhvnACTGWSMKRYBDHVN")];
                                               }

const pair<vector<string>, vector<string>> Trailmix::read_fasta(const string & fastafilename){
    // From https://rosettacode.org/wiki/FASTA_format#C.2B.2B

    igzstream gzinput;
    gzinput.open(fastafilename.c_str());

    string line, name, content, ret;
    vector<string> vec;
    vector<string> ids;

    while(getline(gzinput,line).good()){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if (content != "") {
                vec.emplace_back(content);
                                                 }
            if( !name.empty() ){ // Print out what we read from the last entry
                ret += content;
                ids.emplace_back(name);
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                if (!(contains_only_valid_bases(line))) {
                    string invalid_base;
                    invalid_base += get_invalid_base(line);
                    throw std::runtime_error("[HaploCart] Error, invalid base " + invalid_base);
                                                        }
                content += line;
            }
        }
    }

    if( !name.empty() ){
        ret += content;
        vec.emplace_back(content);
        ids.emplace_back(name);
    }

    if (vec.empty()){throw std::runtime_error("[TrailMix] Error, no sequences found in FASTA input file");}
    vector<string> sorted_ids = ids; sort(sorted_ids.begin(), sorted_ids.end());
    for (size_t i=1; i<sorted_ids.size(); ++i) {
        if (sorted_ids[i] == sorted_ids[i-1]) {cerr << "[TrailMix] Warning: Duplicate id in multiFASTA file: " + sorted_ids[i] << '\n';}
                                                                         }
    return make_pair(vec, ids);
}

const pair<vector<string>, vector<string>> Haplocart::read_fasta(const string & fastafilename){
    // From https://rosettacode.org/wiki/FASTA_format#C.2B.2B

    igzstream gzinput;
    gzinput.open(fastafilename.c_str());

    string line, name, content, ret;
    vector<string> vec;
    vector<string> ids;

    while(getline(gzinput,line).good()){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if (content != "") {
                vec.emplace_back(content);
                                                 }
            if( !name.empty() ){ // Print out what we read from the last entry
                ret += content;
                ids.emplace_back(name);
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                if (!(contains_only_valid_bases(line))) {
                    string invalid_base;
                    invalid_base += get_invalid_base(line);
                    throw std::runtime_error("[HaploCart] Error, invalid base " + invalid_base);
                                                        }
                content += line;
            }
        }
    }

    if( !name.empty() ){
        ret += content;
        vec.emplace_back(content);
        ids.emplace_back(name);
    }

    if (vec.empty()){throw std::runtime_error("[HaploCart] Error, no sequences found in FASTA input file");}
    vector<string> sorted_ids = ids; sort(sorted_ids.begin(), sorted_ids.end());
    for (size_t i=1; i<sorted_ids.size(); ++i) {
        if (sorted_ids[i] == sorted_ids[i-1]) {cerr << "[HaploCart] Warning: Duplicate id in multiFASTA file: " + sorted_ids[i] << '\n';}
                                                                         }
    return make_pair(vec, ids);
}
