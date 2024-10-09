#include "damage.h"
using namespace std;


Damage::Damage(){

}

Damage::~Damage(){

}

//#define DEBUGDEAM

//! A method to combined 5' deam rates and 3' deam rates
/*!
  This method is called by the initDeamProbabilities.
  It uses the "worse" deamination
*/
void Damage::combineDeamRates(double f1[4],double f2[4],double f[4],int b){
    double minFreq = MIN2(f1[b] , f2[b]);
    //cout<<b<<"\t"<<f1[b]<<"\t"<<f2[b]<<"\t"<<f[b]<<"\t"<<minFreq<<endl;

    if(f1[b] == minFreq){//use f1
        for(int i=0;i<4;i++)
            f[i] = f1[i];
                        }
    else{
        if(f2[b] == minFreq){//use f2
            for(int i=0;i<4;i++)
                f[i] = f2[i];
                            }
        else{
            cerr<<"ERROR in combineDeamRates(), wrong state"<<endl;
            exit(1);
            }
        }
                                                                                       }

void Damage::validateProbabilities() {
    // Helper function to check if a probability is valid
    auto isValidProb = [](double p) { return p >= 0.0 && p <= 1.0; };

    // Check each probability in your structures
    for(const auto& sub : {sub5p, sub3p}) {
        for(const auto& probSub : sub) {
            for(const double prob : probSub.s) {
                if (!isValidProb(prob)) {
                    cerr << "Invalid probability detected: " << prob << endl;
                    exit(1);
                }
            }
        }
    }

    for(const auto& subDiNuc : {sub5pDiNuc, sub3pDiNuc}) {
        for(const auto& diNuc : subDiNuc) {
            for(int nuc1 = 0; nuc1 < 4; nuc1++) {
                double sum = 0.0;
                for(int nuc2 = 0; nuc2 < 4; nuc2++) {
                    double prob = diNuc.p[nuc1][nuc2];
                    if (!isValidProb(prob)) {
                        cerr << "Invalid probability detected: " << prob << endl;
                        exit(1);
                    }
                    sum += prob;
                }
                // Allow a small margin for floating-point errors
                if (fabs(sum - 1.0) > 0.00001) {
                    cerr << "Probability sum error for di-nucleotide: " << sum << endl;
                    exit(1);
                }
            }
        }
    }
}


//! A method to initialize the deamination probabilities
/*!
  This method is called by the run/main function
*/
void Damage::initDeamProbabilities(const string & deam5pfreqE,const string & deam3pfreqE){

    vector<substitutionRates> sub5pT;
    vector<substitutionRates> sub3pT;

    readNucSubstitionRatesFreq(deam5pfreqE,sub5pT);
    readNucSubstitionRatesFreq(deam3pfreqE,sub3pT);


    //5'
    for(unsigned int i=0;i<sub5pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub5pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub5pT[i].s[ ind2 ];
	    }

	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;
	}
	sub5p.emplace_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub5p.size()-1);i<MAXLENGTHFRAGMENT;i++){
	sub5p.emplace_back( sub5p[sub5p.size()-1] );
    }

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub5p[i].s[nuc];
	    }
	}
	sub5pDiNuc.emplace_back(diNuctoadd);
    }


    //3'
    for(unsigned int i=0;i<sub3pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub3pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub3pT[i].s[ ind2 ];
	    }

	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;
	}
	sub3p.emplace_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub3p.size()-1);i<MAXLENGTHFRAGMENT;i++){
	sub3p.emplace_back( sub3p[sub3p.size()-1] );
    }

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub3p[i].s[nuc];
	    }
	}
	sub3pDiNuc.emplace_back(diNuctoadd);
    }

#ifdef DEBUGDEAM

    cerr<<endl<<"-- 5' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }

    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
    	}
    	cerr<<endl;
    }

    cerr<<"-- 3' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }


    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
    	}
    	cerr<<endl;
    }

#endif
    //dummy values
    for(unsigned int L=0;L<MINLENGTHFRAGMENT;L++){     //for each fragment length
	vector<probSubstition> subDeam_;
	vector<diNucleotideProb> subDeamDiNuc_;

	subDeam.emplace_back(      subDeam_      );
	subDeamDiNuc.emplace_back( subDeamDiNuc_ );
    }

    cerr<<"Computing substitutions due to ancient damage:"<<endl;
    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length
	//printprogressBarCerr( float(L-MINLENGTHFRAGMENT)/float(MAXLENGTHFRAGMENT-MINLENGTHFRAGMENT) );
	// if( (L%8)==0){
	//     cerr<<".";
	// }
	//vector<alleleFrequency > defaultDNA;
	vector<probSubstition> subDeam_;
	vector<diNucleotideProb> subDeamDiNuc_;
	for(unsigned int l=0;l<L;l++){     //position
	    probSubstition    subDeam__;
	    diNucleotideProb  subDeamDiNuc__;

	    for(int b1=0;b1<4;b1++){//original base
		//pick the highest deam rates for that position
		combineDeamRates(sub5pDiNuc[l].p[b1] ,  sub3pDiNuc[L-l-1].p[b1] , subDeamDiNuc__.p[b1], b1);

		for(int b2=0;b2<4;b2++){//post deam base
		    int b = b1*4+b2;
		    subDeam__.s[b]  = subDeamDiNuc__.p[b1][b2];
		}

	    }
	    subDeam_.emplace_back(      subDeam__      );
	    subDeamDiNuc_.emplace_back( subDeamDiNuc__ );
	}

	subDeam.emplace_back(      subDeam_      );
	subDeamDiNuc.emplace_back( subDeamDiNuc_ );
    }
    cerr<<endl; //flushing progress bar


#ifdef DEBUGDEAM

    cerr<<"-- per length deamination rates --"<<endl;

    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length

	cerr<<endl<<"L="<<L<<endl;
	for(unsigned int l=0;l<L;l++){     //position

	    cerr<<"l="<<l<<" - ";
	    for(int nuc1=0;nuc1<4;nuc1++){
		for(int nuc2=0;nuc2<4;nuc2++){
		    int nuc = nuc1*4+nuc2;
		    cerr<<subDeam[L][l].s[nuc]<<" ";
		}
		cerr<<" - ";
	    }
	    cerr<<endl;
	}
    }



    for(unsigned int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){     //for each fragment length

	cerr<<endl<<"L="<<L<<endl;
	for(unsigned int l=0;l<L;l++){     //position

	    cerr<<"l="<<l<<" - "<<endl;

	    for(int nuc1=0;nuc1<4;nuc1++){
		cerr<<"ACGT"[nuc1]<<"\t";
		for(int nuc2=0;nuc2<4;nuc2++){
		    cerr<<subDeamDiNuc[L][l].p[nuc1][nuc2]<<"\t";
		}
		cerr<<endl;
	    }
	    cerr<<endl;
	}
    }


#endif

    //exit(1);

    //if no deamination, cannot have a mismatch
    for(int b1=0;b1<4;b1++){
	for(int b2=0;b2<4;b2++){
	    int b = b1*4+b2;
	    if(b1==b2){
		defaultSubMatch.s[ b ]           = 1.0;
		defaultSubMatchMatrix.p[b1][b2]  = 1.0;
	    }else{
		defaultSubMatch.s[ b ]           = 0.0;
		defaultSubMatchMatrix.p[b1][b2]  = 0.0;
	    }
	}
    }
    validateProbabilities();
    cerr<<".. done"<<endl;
}//end initDeamProbabilities
