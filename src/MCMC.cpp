#include "MCMC.h"
#include "Euka.h"
#include <algorithm>


//#define DEBUGGENERATEVEC
//#define VERBOSE_MCMC


MCMC::MCMC(){

}

MCMC::~MCMC(){

}


using namespace std;


// computes the softmax of the proposal vector.
// normalises so each proposal vector sums up to 1
const vector<long double> MCMC::softmax(vector<long double> &log_vec){
	long double K = 0.0;

	for (int i = 0; i<log_vec.size(); i++){

		K += exp(log_vec.at(i));
#ifdef DEBUGGENERATEVEC
		cerr << "sum to k " << K << endl;
#endif

	}
#ifdef DEBUGGENERATEVEC
	cerr << "K " << K << endl;
#endif

	vector <long double> transformed_log_vec;
	for(int i = 0; i < log_vec.size(); i++){
		transformed_log_vec.emplace_back(exp(log_vec.at(i))/K);
	}

	return transformed_log_vec;
}


// function to generate the proposal vector for the MCMC runs
const vector<long double> MCMC::generate_proposal(vector<long double> &current_vec, const double &alpha){
	// checking that vector sums up to 1

	long double check = 0.0;
	for (int p = 0; p<current_vec.size(); p++){
#ifdef DEBUGGENERATEVEC
		cerr << current_vec.at(p) << endl; 
#endif
		check += current_vec.at(p);
	}
#ifdef DEBUGGENERATEVEC
	cerr << "check "<< check << endl; 
#endif
	long double lower = 0.99;
	long double upper = 1.01;


#ifdef DEBUGGENERATEVEC
	for (auto element:current_vec){
        cerr << "Current vec before log: "<< element << endl;
    }
#endif

    assert(check > lower && check < upper);
    
	vector <long double> current_vec_log;
	// transform vector into log space.
	for (int q = 0; q<current_vec.size(); q++){
		current_vec_log.emplace_back(log(current_vec.at(q)));

	}
#ifdef DEBUGGENERATEVEC
	for (auto elemet:current_vec_log){
        cerr << "Current vec after log: "<< elemet << endl;

    }
#endif
    // generate random number generator
	std::random_device rd;
	std::mt19937 g(rd());
	// looping through the vector of log transformed distributions.
	//For each element we will sample from a normal distribution to get a new draw
	vector<long double> proposal_vec;
	for (long double element:current_vec_log){

		std::normal_distribution<long double> d{element,alpha};

#ifdef DEBUGGENERATEVEC
		cerr << "exp proposal_vec " << exp(d(g)) << endl;
#endif
		proposal_vec.emplace_back(d(g));

	}

	// transforming the proposal_vec with the softmax function
	const vector<long double> projected_vec = MCMC::softmax(proposal_vec);
#ifdef DEBUGGENERATEVEC
	for (auto element:projected_vec){
        cerr << "exp projected_vec: " << element << endl;
    }
#endif
	check = accumulate(projected_vec.begin(), projected_vec.end(), 0.0);
	lower = 0.99;
	upper = 1.01;
	// Check that we have successfully projected back onto the unit simplex

	assert(check > lower && check < upper);

	return projected_vec;


}

double MCMC::get_average(vector<Clade *> * clade_vec, int id){

	long double sumC = 0.0;
	int n = 0;
	double av = 0.0;

	for (int i = 1; i<clade_vec->size(); i+=3){
		
		if (clade_vec->at(i)->id != id){
			
			for (int j = 0; j < clade_vec->at(i)->clade_not_like.size(); j++){

				sumC += clade_vec->at(i)->clade_not_like.at(j);
				n++; 
			}
			
		}
	

	}
	if (n > 0){
		av = sumC/n;
	}
	return av;
}




long double MCMC::get_proposal_likelihood(const vector <long double> &proposal_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id, const vector<long double> &init_vec){

	MCMC mcmc;

	long double proposal_log_likelihood = 0.0; 
	vector<double> total_f_list;
    for (int i=0; i<clade_list_id.size(); i++){

    	total_f_list.emplace_back(clade_list_id.at(i));
        total_f_list.emplace_back(proposal_vec.at(i));
        
    }

    for (int j = 0; j<total_f_list.size(); j+=2){
        // Make sure we have an even number of elements
    	assert(total_f_list.size() % 2 == 0);


    	
        double frac = total_f_list.at(j+1);
        double total_frac = 0.0;
        double av = 0.0;

        av = mcmc.get_average(clade_vec, total_f_list.at(j));

    	for (int k = 1; k<clade_vec->at(total_f_list.at(j)*3+1)->clade_like.size(); k++){

    		total_frac += log(((frac * clade_vec->at(total_f_list.at(j)*3+1)->clade_like.at(k)) + (clade_vec->at(total_f_list.at(j)*3+1)->clade_not_like.at(k)  * (1/335))));
    		
#ifdef DEBUGGENERATEVEC
    		cerr << "total_frac " << total_frac << endl;
#endif

    	}

    	proposal_log_likelihood += total_frac;

    	}

#ifdef DEBUGGENERATEVEC
    	cerr << "proposal_log_likelihood: " << proposal_log_likelihood << endl;

#endif

    return proposal_log_likelihood;

	
}





vector<long double > MCMC::run(int iter, int burnin, double tol, const vector<long double> &init_vec, vector<Clade *> * clade_vec, vector<int> &clade_list_id){

	MCMC mcmc;

	vector <long double> current_best = init_vec;
	//vector <long double> current_best = {0.0909090909,0.0909090909, 0.0909090909, 0.0909090909, 0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909,0.0909090909};
	vector <long double> proposal_vec;
    
    long double current_log_likelihood = -9999999;
    double acceptance_prob = 0;
    double u=0;
    std::random_device rd;
    std::mt19937 gen(rd());
    int n_reject = 0; 
    
    struct mcmc_moves{
    	long double log_like;
    	vector <long double> abund_vec;
    };
#ifdef DEBUGGENERATEVEC
    //ofstream outputFile("proposal_log_likelihood.txt", ios::trunc);
    //outputFile << "-1" << "\t" << current_log_likelihood << '\t' << current_best[0] << '\t' << current_best[1] << endl;
#endif
    mcmc_moves move[iter];
    int no = 0;

    cerr << "Computing MCMC:" << endl;
	for (int iteration = 0; iteration<iter; iteration++){

		printprogressBarCerr( float(iteration + 1)/float(iter) );

		proposal_vec = mcmc.generate_proposal(current_best, 0.1);
        long double proposal_log_likelihood = mcmc.get_proposal_likelihood(proposal_vec, clade_vec, clade_list_id, init_vec);

       // we are not adding proposal vectors to the struct before after the burin in period
		if (iteration > burnin){

        move[iteration].log_like = proposal_log_likelihood;
        move[iteration].abund_vec = proposal_vec;
         
    	}
        else{

        	continue; 
        }

        

#ifdef DEBUGGENERATEVEC
        //outputFile << iteration << '\t' << proposal_log_likelihood << '\t' << proposal_vec[0] << '\t' << proposal_vec[1] <<'\t' << proposal_vec[2] << '\t' << proposal_vec[3] <<'\t' << '\t' << proposal_vec[4] << '\t' << proposal_vec[5] <<'\t' << '\t' << proposal_vec[6] << '\t' << proposal_vec[7] <<'\t' << '\t' << proposal_vec[8] << '\t' << proposal_vec[9] <<'\t' << '\t' << proposal_vec[10] << '\t';

        cerr << "proposal_log_likelihood " << proposal_log_likelihood << endl; 
#endif
        acceptance_prob = min((long double)(1.0), expl(proposal_log_likelihood-current_log_likelihood));
#ifdef DEBUGGENERATEVEC
        cerr << "acceptance_prob " << acceptance_prob << endl; 
#endif
        std::uniform_real_distribution<> dis(0, 1);
        u = dis(gen);
#ifdef DEBUGGENERATEVEC
        cerr << u << endl;
        //proposal_log_likelihood > current_log_likelihood
#endif 
        if (u <= acceptance_prob || iteration == 0) {
#ifdef DEBUGGENERATEVEC
        	//outputFile << "accept" << endl;
#endif 
#ifdef VERBOSE_MCMC
            cerr << "ACCEPTING proposal." << endl;
            cerr << "Proposal log likelihood: " << proposal_log_likelihood << endl;
#endif
            current_log_likelihood = proposal_log_likelihood;
            current_best = proposal_vec;
                                                    }
        else {
#ifdef DEBUGGENERATEVEC
        	//outputFile << "reject" << endl; 
#endif
#ifdef VERBOSE_MCMC
            cerr << "REJECTING proposal" << endl;
#endif
            current_best = current_best;
            n_reject++;
             }
#ifdef VERBOSE_MCMC
        cerr << "\n\nCurrent log likelihood: " << current_log_likelihood << endl;
#endif
	}
	cerr<<endl;

#ifdef VERBOSE_MCMC
	cerr << "MCMC completed. The final log likelihood is: " << current_log_likelihood << endl; 
#endif


	
	int per85 = 0.85 * (iter - burnin); 
	int per95 = 0.95 * (iter - burnin);
	long double sums[proposal_vec.size()] = {0.0};
	vector<long double> posterior_estimate;
	vector<long double> sorted_clade; 

	// The posterior mean as well as the confidence intervalls will be calculated per clade (each clade is independent from each other)
	//We are lọoping through the abundance vector with j
	// 
	for (int j=0; j<proposal_vec.size(); j++){

		for (int i = burnin+1; i<iter; i++){
			
			sums[j] += move[i].abund_vec[j];
#ifdef DEBUGOUTPUT
			cerr << "before " <<move[i].abund_vec[j] <<endl;  
#endif		
			sorted_clade.emplace_back(move[i].abund_vec[j]);

		}
		// sorting the struct for highest posterior density intervall (HDI)

		sort(sorted_clade.begin(), sorted_clade.end());

		const int m = sorted_clade.size()/2;
		
		long double p = sorted_clade[m];


		// posterior point estimate for each clade (each fraction of the abdundance vector)
		//long double p = sums[j]/(iter-burnin);
#ifdef DEBUGOUTPUT
		cerr << "posterior_estimate " << p << endl; 
#endif 
		posterior_estimate.emplace_back(p);

		long double low_end_85 = quant(sorted_clade, 0.15);
		long double high_end_85 = quant(sorted_clade, 0.85);
		long double low_end_95 = quant(sorted_clade, 0.05);
		long double high_end_95 = quant(sorted_clade, 0.95);

		posterior_estimate.emplace_back(low_end_85);
		posterior_estimate.emplace_back(high_end_85);
		posterior_estimate.emplace_back(low_end_95); 
		posterior_estimate.emplace_back(high_end_95);

		sorted_clade.clear();


		
	}
	
	
	return posterior_estimate;

}
