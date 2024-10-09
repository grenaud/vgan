#include "Euka.h"
#include "MCMC.h"

using namespace std;


//#define DEBUGINIT

const vector<double> Euka::compute_init_vec(vector<Clade *> * clade_vec, vector<int> &clade_list_id){

// summarising the total number of fragment counts to receive the distribution between the clade that have been detected.
    double total_frag_count = 0;

    for (int i=0; i<clade_list_id.size(); i++){

    	total_frag_count += clade_vec->at(clade_list_id.at(i)*3+1)->count;

    }

#ifdef DEBUGINIT
    cerr << "total fragment count "<< total_frag_count << endl;
#endif

    vector<double> total_f_list;
    for (int n=0; n<clade_list_id.size(); n++){

        const double frag_count = clade_vec->at(clade_list_id.at(n)*3+1)->count;

        const double dis = frag_count/total_frag_count;

#ifdef DEBUGINIT
        cerr << "Clade id: " << clade_list_id.at(n) <<endl;
        cerr << "frag " << frag_count << " total count " << total_frag_count << endl;
        cerr << "distribution: " << dis <<endl;
#endif
        // creating a vector with the inital abundances of the different clades.

        //total_f_list.emplace_back(clade_list_id.at(n));
        total_f_list.emplace_back(dis);

    }
    return total_f_list;
}
