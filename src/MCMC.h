#pragma once

#include "Clade.h"
#include <vector>
//
using namespace std;

class MCMC {

	private:

	public:
		const vector<long double> generate_proposal(vector<long double> &current_vec, const double &alpha);
		const vector<long double> softmax(vector<long double> &log_vec);
		double get_average(vector<Clade *> * clade_vec, int id);

        vector<long double > run(int iter, int burnin, double tol, const vector<long double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id);
        long double get_proposal_likelihood(const vector <long double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id, const vector<long double> &init_vec);


	    MCMC();
	    MCMC(const MCMC & other);
	    ~MCMC();
	    MCMC & operator= (const MCMC & other);
};
