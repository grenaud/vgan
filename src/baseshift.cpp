#include "baseshift.h"
#include "Euka.h"


using namespace std;

Baseshift::Baseshift(int lengthMax, int** baseshift_data_array){

    // define constants
    lengthToProf = lengthMax;
    dna2int['A']=0;
    dna2int['C']=1;
    dna2int['G']=2;
    dna2int['T']=3;

    // save data array as private data array
    this->baseshift_data_array = baseshift_data_array;

    return;
}

void Baseshift::baseshift_calc(string graph_seq, string read_seq){

    char graph_base;
    char read_base;
    int pos = 0;

    for (int p = 0; p < lengthToProf*2; p++) {
        // find bases from seq begining
        if (p < lengthToProf){
            pos = p;
            graph_base = toupper(graph_seq[pos]);
            read_base = toupper(read_seq[pos]);
        }
        // find bases from seq end
        if (p >= lengthToProf) {
            pos = -(lengthToProf*2)+p;
            graph_base = toupper(graph_seq[graph_seq.length()+pos]);
            //cout << "graph_base " << graph_base << endl;
            read_base = toupper(read_seq[read_seq.length()+pos]);
            //cout << "read_base " << read_base << endl; 
        }

        // Ignore softclips , indels and gaps
        if (graph_base == 'S' or read_base == 'S') continue;
        if (graph_base == 'I' or read_base == 'I') continue;
        if (graph_base == '-' or read_base == '-') continue;
        if (graph_base == 'N' or read_base == 'N') continue;
        //cout <<"graph base " << graph_base <<" "<< dna2int[graph_base] << " " <<"read base " << read_base <<" " << dna2int[read_base] << endl;
        baseshift_data_array[p][toIndex[dna2int[graph_base]][dna2int[read_base]]] ++;
    }

    return;
}

void Baseshift::print_counts(string prof_out_file) const{

    ofstream prof_out;
    prof_out.open(prof_out_file.c_str());

    prof_out <<"A>A\tA>C\tA>G\tA>T\tC>A\tC>C\tC>G\tC>T\tG>A\tG>C\tG>G\tG>T\tT>A\tT>C\tT>G\tT>T\tPosition"<<endl;

    for (int p = 0; p < lengthToProf*2; p++){
        for (int i = 0; i < 16; i++){
            prof_out << baseshift_data_array[p][i] << "\t";
        }

        if (p<lengthToProf) prof_out << "\t" << p << endl;
        if (p>=lengthToProf) prof_out << "\t" << -(lengthToProf*2)+p << endl;
        if (p==lengthToProf-1) prof_out << endl;
    }
    prof_out << endl;
    return;
}


vector<vector<double>> Baseshift::print_prof(string prof_out_file, string ends) const{

    if (ends != "both" && prof_out_file != "/dev/stdout") {
       prof_out_file.insert(prof_out_file.length()-5, "_"+ends);
    }

    ofstream prof_out;
    prof_out.open(prof_out_file.c_str());

    prof_out <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\tPosition"<<endl;

    int col = 0;
    int col_div = 0;

    vector <double> CtoT;
    vector <double> GtoA;

    for (int p = 0; p < lengthToProf*2; p++){
        // skip 5p if 3p is chosen
        if (ends == "3p") {
            if (p < lengthToProf) continue;
        }

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){

                if (i != j) {



                    prof_out << setprecision(4) << (double)baseshift_data_array[p][col] / 
                        (   (double)baseshift_data_array[p][0 + col_div] +
                            (double)baseshift_data_array[p][1+col_div] +
                            (double)baseshift_data_array[p][2+col_div] +
                            (double)baseshift_data_array[p][3+col_div] 
                            ) << "\t";

                        if (i == 1 && j == 3){
                            CtoT.emplace_back((double)baseshift_data_array[p][col] / 
                        (   (double)baseshift_data_array[p][0+col_div] +
                            (double)baseshift_data_array[p][1+col_div] +
                            (double)baseshift_data_array[p][2+col_div] +
                            (double)baseshift_data_array[p][3+col_div] 
                            ));
                        }
                        if (i == 2 && j == 0){
                            //cout << "pos " << p << endl;
                            //cout << "col " << col << endl; 
                            //cout << "i "<< i << endl; 
                            //cout << "j " << j << endl; 
                            GtoA.emplace_back((double)baseshift_data_array[p][col] / 
                        (   (double)baseshift_data_array[p][0+col_div] +
                            (double)baseshift_data_array[p][1+col_div] +
                            (double)baseshift_data_array[p][2+col_div] +
                            (double)baseshift_data_array[p][3+col_div] 
                            ));
                        }
                }
                col++;
            }
            col_div += 4;
        }


        if (p<lengthToProf) prof_out << p << endl; 
        if (p>=lengthToProf) prof_out << -(lengthToProf*2)+p << endl;
        if (p==lengthToProf-1) {
            if (ends == "5p") break; // Break loop if 5p is chosen
            prof_out <<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\tPosition"<<endl;
        } 
        col = 0;
        col_div = 0;


    }
    //cout << "size comparison "<< CtoT.size() <<" " << GtoA.size() << endl; 
    

    vector<vector<double>> damage_subs;
    damage_subs.emplace_back(CtoT);
    damage_subs.emplace_back(GtoA);

    // for (const auto el : GtoA){
    //     cout << el << endl; 
    // }

    return damage_subs;
}

