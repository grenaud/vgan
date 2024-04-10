#include "Euka.h"
#include "MCMC.h"

using namespace std;


//#define DEBUGINIT

const vector<long double> Euka::compute_init_vec(vector<Clade *> * clade_vec, vector<int> &clade_list_id){

// summarising the total number of fragment counts to receive the distribution between the clade that have been detected.
    long double total_frag_count = 0;

    for (int i=0; i<clade_list_id.size(); i++){

    	total_frag_count += clade_vec->at(clade_list_id.at(i)*6+1)->count;

    }

#ifdef DEBUGINIT
    cerr << "total fragment count "<< total_frag_count << endl;
#endif

    vector<long double> total_f_list;
    for (int n=0; n<clade_list_id.size(); n++){

        const long double frag_count = clade_vec->at(clade_list_id.at(n)*6+1)->count;

        const long double dis = frag_count/total_frag_count;

#ifdef DEBUGINIT
        cerr << "Clade id: " << clade_list_id.at(n) <<endl;
        cerr << "frag " << frag_count << " total count " << total_frag_count << endl;
        cerr << "distribution: " << dis <<endl;
#endif
        // creating a vector with the inital abundances of the different clades.

        //total_f_list.emplace_back(clade_list_id.at(n));
        total_f_list.emplace_back(dis);

    }
//     vector <long double> new_vec;
//     for (int j = 0; j<total_f_list.size(); j+=2){
//         // Make sure we have an even number of elements
//     	assert(total_f_list.size() % 2 == 0);

//     	int vec_len = total_f_list.size()/2;
//         long double frac = total_f_list.at(j+1);
//         long double total_frac = 0.0;
// #ifdef DEBUGINIT
//         cerr << "vec_len " << vec_len << endl; 
//         cerr << "frac " << frac << endl; 
// #endif
//     	for (int k = 1; k<clade_vec->at(total_f_list.at(j)*3+1)->clade_like.size(); k++){

//             cerr << "clade like at position 0: " <<clade_vec->at(total_f_list.at(j)*3+1)->clade_like.at(0) << endl; 

//             cerr << "clade like at postition " << k << ": " << clade_vec->at(total_f_list.at(j)*3+1)->clade_like.at(k) << endl; 

//             //double inv = (frac * clade_vec->at(total_f_list.at(j)*3+1)->clade_like.at(k)) + (clade_vec->at(total_f_list.at(j)*3+1)->clade_not_like.at(k)/ vec_len);
//     		total_frac += (frac * clade_vec->at(total_f_list.at(j)*3+1)->clade_like.at(k)) + (clade_vec->at(total_f_list.at(j)*3+1)->clade_not_like.at(k)/ vec_len);


//     	}
// #ifdef DEBUGINIT
//     	cerr << total_frac << endl;
//     	cerr << clade_vec->at(total_f_list.at(j)*3+1)->clade_like.size() << endl;
// #endif

//     	long double new_frac = (total_frac/clade_vec->at(total_f_list.at(j)*3+1)->clade_like.size());


//     	new_vec.emplace_back(new_frac);
// #ifdef DEBUGINIT
//     	cerr << "New fraction for clade " << clade_vec->at(total_f_list.at(j)*3+1)->name << " is " << new_frac << endl;
// #endif


//     }

    return total_f_list;


}
