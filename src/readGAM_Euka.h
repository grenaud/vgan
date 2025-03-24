#pragma once

#include "bdsg/odgi.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <variant>
#include <numeric>
#include <functional>
//include
#include "vg/vg.pb.h"
#include "vg/io/basic_stream.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "vg/io/alignment_io.hpp"
#include "vg/io/json2pb.h"
#include "vg/io/stream.hpp"
//VG src/
#include "Euka.h"
#include "vgan_utils.h"
#include "utility.hpp"
#include "alignment.hpp"
#include "AlignmentInfo.h"
#include "libgab.h"
#include "miscfunc.h"
#include "baseshift.h" // Mikkel code
//#define VERBOSE_S
//#define DEBUGREADGAM
//#define DEBUGMUTDEAMERROR
//#define DEBUGDAMAGERATES
//#define DAMAGE_HIGH_LOW_TEST
//#define DEBUGDAMAGE2

using namespace std;
using namespace google::protobuf;
namespace fs = std::filesystem;




static pair <vector<AlignmentInfo *>*, vector<int>> readGAM3(const bdsg::ODGI &graph, const string &gamfilename, const bool populatevector,vector<Clade *> * clade_vec, \
 const vector<NodeInfo *> &nodevector, const vector<double> &qscore_vec, \
 const double * base_freq, const double t_T_ratio['T' + 1]['T' +1 ], \
 const bool * rare_bases, vector<vector<tuple<int, int, double, double >>> &chunks, \
 vector <vector<diNucleotideProb> > &subDeamDiNuc, int lengthToProf, string prof_out_file_path, const unsigned int &MINIMUMMQ, const unsigned int &MINNUMOFREADS, \
 const unsigned int &MINNUMOFBINS, const double &ENTROPY_SCORE_THRESHOLD, const int &MAXIMUMOFBINS)
{

    if (std::filesystem::exists(gamfilename) == false){
        throw(std::runtime_error(gamfilename + " does not exist. Aborting."));
    }

#ifdef VERBOSE_S
    cerr<<"reading GAM file "<<gamfilename<<endl;

#endif


    int n_reads=0;
    int n_map_reads=0;
    vector<AlignmentInfo *> * read_vec=NULL;
    if(populatevector){
        read_vec = new vector<AlignmentInfo *>();
    }


    function<void(vg::Alignment&)> lambda = [&graph, &n_reads,&read_vec,&populatevector, &n_map_reads, &clade_vec, &nodevector, &qscore_vec, \
       &base_freq, &t_T_ratio, &rare_bases, &chunks, &subDeamDiNuc, &lengthToProf, &prof_out_file_path, &MINIMUMMQ, &MINNUMOFREADS, &MINNUMOFBINS, \
       &ENTROPY_SCORE_THRESHOLD, &MAXIMUMOFBINS](vg::Alignment& a)   // Mikkel code two arguments lengthToProf & prof_out_file_path
    {

        ++n_reads;

    // if Identity is 0 the read is unmapped.
        if (a.identity() != 0){
        ++n_map_reads; // we count both the mapped and the unmapped reads


        if(populatevector){
            AlignmentInfo * ai = new AlignmentInfo();
            ai->seq = a.sequence();
            ai->name = a.name();
            ai->path = a.path();
            ai->mapping_quality = a.mapping_quality();
            ai->quality_scores = a.quality();
            read_vec->emplace_back(ai);

        }else
        {


            ///////////////////////
            //  Begin EUKA code  //
            ///////////////////////
            
#ifdef DEBUGREADGAM
            cerr<<a.name() << " "<<a.sequence()<<" "<<a.sequence().size()<<endl;
#endif



        // get the first node id for the mapped read
            int n_index = a.path().mapping()[0].position().node_id();


#ifdef DEBUGREADGAM

            cerr<<"n_index "<<n_index<<endl;

#endif
            // the node ID for the first node of the mapping is saved as n_index
            // the n_index is used to find the clade ID by looping through the saved bin indexes. Every clade has their own bins defind by the node IDs.
            // As soon as the fist node ID of the mapping is found within the bins the number of the vector is returned, which corresponds to the clade ID.
            // The clade ID (c_n) is then used to get the clade name, clade pairwise average distance, etc.
            int c_n = 0;

            double entropy_score = 0.0;


            for (long unsigned int i=0; i<chunks.size(); i++){

                for (long unsigned int j=0; j<chunks.at(i).size();++j){

                //adds a count to the bin the read falls into for the specific clade.
                    if (n_index == std::clamp(n_index, get<0>(chunks.at(i).at(j)), get<1>(chunks.at(i).at(j)))) {


                        entropy_score = get<2>(chunks.at(i).at(j));
  
                        c_n = i;
                    

                    }else{
                        continue;

                    }


                }
            }

            

#ifdef DEBUGREADGAM

            cerr<<"c_n "<<c_n<<endl;

#endif

            // with the clade id access the pairwise average distance computed for each clade

            double pair_dist = clade_vec->at(c_n*6 +1)->dist;  // *3+1 is to index the vector of pointers (every clade has 3 entries (*3) but the counting starts at 0)


            // reconstructing the graph and the read sequence
            auto [graph_seq, read_seq, mppg_sizes] = reconstruct_graph_sequence(graph, a.path(), a.sequence());


            //// MIKKEL CODE START HERE //////

            // Init baseshift data array in clade_vec

            unsigned int** baseshift_data_array_location = clade_vec->at(c_n*6+1)->baseshift_clade_array;

            // init baseshift data array in baseshift class
            Baseshift baseshift_data_array(lengthToProf, baseshift_data_array_location);

            //cerr << "Clade name: " << clade_vec->at(c_n*3+1)->name << endl;
            //cerr << "Id: "<< clade_vec->at(c_n*3+1)->id << endl;
            //cerr << "Clade baseshift_clade_array memory location: " << baseshift_data_array_location << endl;

            // Calculate baseshifts
            baseshift_data_array.baseshift_calc(graph_seq, read_seq);

            //// MIKKEL CODE END HERE //////



            double in_clade_lik     = 0.0;
            double not_in_clade_lik = 0.0;
            int edit_dist           = 0;
            int mis_dist            = 0;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // BEGIN CODE TO COMPUTE LIKELIHOOD OF CORRECT ALIGNMENT
            // We compute the likelihood of 2 models:
            //  model 1: the fragment came from some mitogenome in the clade, we believe the alignment to be genuine
            //  model 2: the fragment aligned at random, it came from bacteria or fungi or something else
            //
            //  We then compute the following likelihood ratio:
            //
            //          model 1 / model 2
            //
            //  if they are equal this means that we cannot trust the alignment
            //
            /////////////////////////////////////////////////////////////////////////////////////////////////////////

            double log_lik   = 0.0; //log likelihood of model 1: the fragment came from some mitogenome in the clade
            double log_lik_2 = 0.0; //log likelihood of model 2: the fragment aligned at random

            //temporary variables
            double first = 0.0;
            double second = 0;

            // initiating a softclip counter to penalise them correctly.
            int softclip_count = 0;

            unsigned int n=0; //this is the coordinate on the aligned fragment
            unsigned int Lseq = a.sequence().size(); //this is the length of the aligned fragment

            const bool isrev = a.path().mapping()[0].position().is_reverse(); //if the fragment is on the reverse strand
            //unsigned int incrementStrand=1;

            if (isrev == true){ // in case of a reversed fragment the counter starts at the last position of the fragment (5' end) and counts backwards
                n= a.sequence().size()-1;
            }


#ifdef DEBUGDAMAGERATES

            cerr<<a.name()<<'\t'<<a.sequence().size()<<endl;
            cerr<<"Original seq: "<<a.sequence()<<endl;
            cerr<<"G-seq: "<<graph_seq<<" R-seq: "<<read_seq<<endl;
#endif


           //for each base
            for ( unsigned int m=0; m<graph_seq.size();m++)
            {

#ifdef DEBUGREADGAM
                cerr<<m<<" G="<<graph_seq[m]<<" R="<<read_seq[m]<<" Q="<<int(a.quality()[m])<<" pos="<<n<<endl;
#endif

                    //This if either base is unresolved, we consider both to be
                if (graph_seq[m] == 'N' || read_seq[m]== 'N') {
                    // treated as sequencing error
                   log_lik   = base_freq[read_seq[m]];
                   log_lik_2 = base_freq[read_seq[m]];

               }

                    //gaps in either the fragment or graph
                else if (graph_seq[m] == '-' || read_seq[m] == '-'){
                    // propability of observing a gap in a real MSA.
                log_lik   = log(0.002);
                    // propability of observing a gap when we do random alignments against the graph
                log_lik_2 = log(0.2);
                }

                    //the UIPAC unresolved bases, unnecessary but just in case
                else if (rare_bases[graph_seq[m]] == true || rare_bases[read_seq[m]] == true) {
                    log_lik   = log((1-pair_dist)* 0.001);
                        //any rare bases should be treated as a sequencing error
                    log_lik_2 = log(0.001);

                }

                    // Softclips penalty
                    // Softclips will be penelised by assuming that 75% of the time they are actually a missmatch and 25% of the time they are a match.
                    // The same ratio will be applyed for model 2.

                else if (graph_seq[m] == 'S' || read_seq[m] == 'S') {
                    int base_quality = a.quality()[m];
                        //cerr << base_quality << endl;
                        //int base_quality = 10;

                    ++softclip_count;
                    //cerr << softclip_count << endl;
                    if (softclip_count%3==0) {
                        log_lik = log(1 - qscore_vec[base_quality]);

                    }
                    else {
                    log_lik = log(qscore_vec[base_quality]/3);
                    }

                log_lik_2 = log(0.25);

                }
                //This gets executed if you do not have a softclip or a gap or unresolved base
                //you either have a match or a mismatch
                else
                {
                    // access the base quality for every base of every read. Because it is a const it is initialised before.
                    int base_quality = a.quality()[m];
                    // match both the base from the fragment and the sequence in the graph agree


//////////////////////////  BEGIN CODE MODEL 1 //////////////////
                    double probBasePreDamage [4];

#ifdef DEBUGDAMAGERATES
                    if (graph_seq[m] == read_seq[m]){
                        cerr<<"match ";
                    }
                    else{
                        cerr<<"mismatch ";
                    }

                    cerr<<"graph="<<graph_seq[m]<<" read="<<read_seq[m]<<" q="<<qscore_vec[base_quality]<< " pos="<<n<<endl;
                    if(isrev == true){
                        cerr << "- strand" <<endl;
                    }
                    else {
                        cerr<<"+ strand"<<endl;
                    }
                    cerr<<"predamage:"<<endl;
#endif
                    // LEVEL 1 //
                    // filling a substitution matrix with the propbabilities of observing a base pre damage.
                    for(int bpo=0;bpo<4;bpo++){
                        if("ACGT"[bpo] == graph_seq[m]){// no mutation
                            probBasePreDamage[bpo]=(1 - pair_dist);
                        }else{ // mutation
                            probBasePreDamage[bpo]=pair_dist * t_T_ratio[graph_seq[m]]["ACGT"[bpo]];
                        }


#ifdef DEBUGDAMAGERATES
                        cerr<<bpo<<"\t"<<probBasePreDamage[bpo]<<endl;
#endif
                    }

                    // LEVEL 2 //
                    //initializing a post damage substitution matrix. First filled with 0.
                    double probBasePostDamage [4];

                    for(int bpd=0;bpd<4;bpd++){
                        probBasePostDamage[bpd]=0;
                    }

                    // filling post-damage substitution matrix with the new probabilities of observing a base. The probabilities are based on the damage profile for --deam5p and --deam3p

                    // If no damage profile is provided the pre-damage substitution matrix is multiplied by 0 and stays the same.

                    for(int bpd=0;bpd<4;bpd++){
                        //post damage = prob pre damage * damage rate from bpo to bpd
                        for(int bpo=0;bpo<4;bpo++){
                            probBasePostDamage[bpd]+=probBasePreDamage[bpo]*subDeamDiNuc[Lseq][n].p[bpo][bpd];
#ifdef DEBUGDAMAGE2
                            cerr << probBasePreDamage[bpo]<<"," << subDeamDiNuc[Lseq][n].p[bpo][bpd] << "," << probBasePostDamage[bpd]<< endl;
#endif

#ifdef DEBUGDAMAGERATES

                            cerr<<"probBasePreDamage: "<< probBasePreDamage[bpo] << " subDeamDiNuc: " << subDeamDiNuc[Lseq][n].p[bpo][bpd] << endl;

                            if (graph_seq[m] == 'C' && read_seq[m] == 'T'){

                                cerr << "DAMAGE" << endl;
                                cerr << Lseq << "," << n << endl;
                                cerr << subDeamDiNuc[Lseq][n].p[bpo][bpd] << endl;
                            }
                            else if (graph_seq[m] == 'G' && read_seq[m] == 'A'){
                                cerr << "DAMAGE" << endl;
                                cerr << Lseq << "," << n << endl;
                                cerr << subDeamDiNuc[Lseq][n].p[bpo][bpd] << endl;
                            }
                            else {

                                cerr << "NO DAMAGE" << endl; 

                            }
#endif

                        }
                    }

#ifdef DEBUGDAMAGERATES
                    cerr<<"postdamage:"<<endl;

                    for(int bpd=0;bpd<4;bpd++){
                        cerr<<bpd<<"\t"<<probBasePostDamage[bpd]<<endl;
                    }
#endif



                    double log_lik_marg   = 0.0; //sum of the likelihoods but in log space, this will be added to log_lik
                    // LEVEL 3 //
                    //as we do not know the base post damage, we marginalize over it
                    // we already know the probability of each base after damage. we marginalise over these probabilities and add the probability to see a seq error
                    // if the base in the read is still equal to the bpd (which are the base probility after damage from the graph) there also can't be a seq error (1 - e)
                    for(int bpd=0;bpd<4;bpd++){
                        if( "ACGT"[bpd]== read_seq[m]){// no sequencing error
                        //                                                    Prob of bpd                no seq error
                            log_lik_marg = oplusInitnatl( log_lik_marg , log(    probBasePostDamage[bpd]*(  1 - qscore_vec[base_quality]) ) );

                        }else{ // we are already marginalising over all possible bases, if the base is not matching it must be a seq error.
                        //                                                    Prob of bpd                seq error
                            log_lik_marg = oplusInitnatl( log_lik_marg , log(    probBasePostDamage[bpd]* (qscore_vec[base_quality]/3) ) ) ;
                        }
                    }

#ifdef DEBUGDAMAGERATES
                    cerr<<"log_lik_marg:"<<log_lik_marg<<" p="<<exp(log_lik_marg)<<endl;

#endif
        
                    log_lik = log_lik_marg;


#ifdef DAMAGE_HIGH_LOW_TEST

                    if (graph_seq[m] == 'C' && read_seq[m] == 'T'){

                        cout<<exp(log_lik)<<endl;
                    }
                    else if (graph_seq[m] == 'G' && read_seq[m] == 'A'){
                        cout<<exp(log_lik)<<endl;
                    }
                }
#endif

    ////////////////////  END CODE MODEL 1 ///////////////////

    /////////////////// BEGIN CODE MODEL 2 //////////////////

                    if (graph_seq[m] == read_seq[m]){
                                
                        
                        log_lik_2 = log(1-0.25536);
                    

                    }
                // mismatch the base from the fragment and the sequence in the graph do not agree
                    else {
                        // mismatch counter

                        mis_dist += 1;
                        //                         (1-u) (e/3)                       u * P (graph base to seq base) (1-e)
                        //first   = oplusnatl((log(1- pair_dist) + log(qscore_vec[base_quality]/3)), (log(pair_dist * t_T_ratio[graph_seq[m]][read_seq[m]]) + log(1 - qscore_vec[base_quality])));
                        //second  = oplusnatl(first,(log(pair_dist * qscore_vec[base_quality]/3)));
                        //log_lik = oplusnatl(second,(log(pair_dist * qscore_vec[base_quality]/3)));

        

                        // Model 2 is calculated based on the amount of mismatchtes that can be seen when random non-biological sequences get aligned to the graph.
                        // Every mismatch from every read that maps is counted (read lengths 20-150bp in 5 step scale 1000x) and the averge taken.
                        log_lik_2 = log(0.25536);

                    }
                


    ////////////////// END CODE MODEL 2 //////////////////////

                } // else statment ending of the S or rate bases case senario

                in_clade_lik     += log_lik;
                not_in_clade_lik += log_lik_2;
#ifdef DEBUGDAMAGE2
                cerr<<a.name()<<","<<(in_clade_lik - not_in_clade_lik)<< "," << clade_vec->at(c_n*6 +1)->name <<endl; 
#endif

                if(read_seq[m] != '-' && isrev == false){
                    n++;
                } else if (read_seq[m] != '-' && isrev == true){
                    n--;
                }


            } // END loop through each base


    //////////////////////////////// CLADE LIKELIHOOD ////////////////////////////////////////
    /// We are computing the clade likelihood with the following function: 
    /// 
    /// Likelihood function = (vi * (logit(M1/M2) * MQ) + (1 - (logit(M1/M2) * (MQ)) / len(v)))
    //
    // function to update the likelihood function.
    // The clade_like and clade_not_like vector are stored within the CLADE class and are estimated for each read that maps to a specific clade.
    // clade like = logit(M1/M2) * (MQ)
    // clade_not_like = 1 - (logit(M1/M2) * (MQ))
    // As soon as  a clade has passed all filters and therefore has been "detected", it is considered for the distribution vector (init_distr)
    // The init distr is calculated based on the number of total reads mapped to detected clades and than divided by the number of reads mapped to one clade
    // The number of fragments of the distribution = len(v)
    // The fragment of the distribution = vi
    /// the product of the mapping quality and the logit of like_ratio is stored in the clade object as a vector "clade_like". It is computed while going through the gam file. 
    /// The rest of the function will be computed after the lambda function went through the gam file. 
    /////////////////////////////////////////////////////////////////////////////////////////


            double ratio = (in_clade_lik - not_in_clade_lik); 
            

            double map_q = (1- get_p_incorrectly_mapped(a.mapping_quality())); //get_p_incorrectly_mapped is taken from Haplocart but also stored in miscfunc.h
            
            // saving the probability of each read coming from the clade and not coming from that clade.
            clade_vec->at(c_n*6+1)->clade_like.emplace_back(map_q * exp((in_clade_lik) - oplusInitnatl(in_clade_lik, not_in_clade_lik)));
            clade_vec->at(c_n*6+1)->clade_not_like.emplace_back(1 - (map_q *(exp((in_clade_lik) - oplusInitnatl(in_clade_lik, not_in_clade_lik)))));
    



#ifdef DEBUGREADGAM

            cerr << a.name() << "," << in_clade_lik << "," << not_in_clade_lik << "," << (in_clade_lik - not_in_clade_lik) << "," << clade_vec->at(c_n*6 +1)->name << "," << a.mapping_quality() << endl;
#endif

            //cout << a.name() << "," << in_clade_lik << "," << not_in_clade_lik << "," << (in_clade_lik - not_in_clade_lik) << "," << clade_vec->at(c_n*3 +1)->name << "," << a.mapping_quality() << endl;
            // the mapping (read and all nodes it mapped to) has a likelihood higher than 1 it will be counted and included in the coverage estimation
            if (
                ((in_clade_lik) - (not_in_clade_lik) > 1) && //if the likelihood of the in-clade model is greater than the random alignment one
                
                //cout << clade_vec->at(c_n*6 + 1)->name << '\t' << a.mapping_quality() << endl; 
                (a.mapping_quality() > MINIMUMMQ) //&&                   //AND the mapping quality is greater than 29
                //(entropy_score > ENTROPY_SCORE_THRESHOLD) 
                

                ){ // we only take detected reads and pass them on 
		    //cout << clade_vec->at(c_n*6 + 1)->name << '\t' << a.mapping_quality() << endl;
                    //cerr << a.name() << ",PASS"<< "\t" <<((in_clade_lik) - (not_in_clade_lik))<<endl;
                    clade_vec->at(c_n*6 +1)->count++;
                    clade_vec->at(c_n*6+1)->inSize.emplace_back(a.sequence().size());
                    clade_vec->at(c_n*6+1)->nameStorage.emplace_back(a.name());


#ifdef DEBUGREADGAM
                    cerr << a.name() << ",PASS"<< "\t" <<((in_clade_lik) - (not_in_clade_lik))<<endl;
#endif

                    // going through the mapping again to count all node ids that the read hits.
                     
                    for (long unsigned int i=0; i<a.path().mapping().size(); ++i) {

                        // node id to look for c_n is the clade id to get the right bins
                        int n_id = a.path().mapping()[i].position().node_id();

                        //cerr << chunks.at(c_n).size() << endl;
                        for (long unsigned int j=0; j<chunks.at(c_n).size();++j){

                            //adds a count to the bin the read falls into for the specific clade.
                            double newt = 0.0; 
                            if (n_id == std::clamp(n_id, get<0>(chunks.at(c_n).at(j)), get<1>(chunks.at(c_n).at(j)))) {
                            // we are adding a count for each node id that a mapping hits, there can be multiple mappings per read
                            // therefore multiple nodes can be hit and more than one count can be issues per read
                                std::cerr.precision(20);

                                newt = 1.0/a.path().mapping().size();
                                //newt = 1.0/; 
                                get<3>(chunks.at(c_n).at(j)) += newt;
                              
                                
                                
                            }
                        }
                    }




            }else{//end if a fragment has sufficient MQ and likelihood ratio

#ifdef DEBUGREADGAM
                cerr << a.name() << ",FAIL"<< endl;
#endif

                }



                
        

            } // END if not populating vector
        


        }//END if identify not zero
    




    }; //END lambda function
    
    

    vg::get_input_file(gamfilename, [&](istream& in) { vg::io::for_each(in, lambda);    });

    int sum_of_clades = 0;
    for (int i = 1;i<clade_vec->size();i+=6) {
        sum_of_clades += clade_vec->at(i)->count;
        //cout << "sum_of_clades" << sum_of_clades << endl;

       
    }

    vector<int> clade_bins;
    vector<int> clade_list_id;
    vector<int> clade_list_count;
    vector<int> extra;


    int max_allowed_zero_bins = MAXIMUMOFBINS;

    for (int i = 0; i < chunks.size(); i++) {

        vector<int> check_for_zero;
        for (int k = 0; k < chunks[i].size() - 1; k++) {
            if (get<2>(chunks.at(i).at(k)) > ENTROPY_SCORE_THRESHOLD) {
                // cerr <<  get<3>(chunks[i].at(k)) << endl;
                check_for_zero.emplace_back(get<3>(chunks[i].at(k)));
            }
        }

    // Count how many zeros are in the vector
    int num_zero_bins = std::count(check_for_zero.begin(), check_for_zero.end(), 0.0);

    // Modify the condition to allow a certain number of zero bins
    if (num_zero_bins > max_allowed_zero_bins || check_for_zero.size() < MINNUMOFBINS || clade_vec->at(i*6+1)->count < MINNUMOFREADS) {

#ifdef VERBOSE_S
            cerr << "The clade " << clade_vec->at(i*6+1)->name << " could not be detected." << endl;
            
#endif     
        }
        else {
            // creating a list of ids 
            clade_list_id.emplace_back(clade_vec->at(i*6+1)->id);



#ifdef VERBOSE_S
            cerr << "Clade " << clade_vec->at(i*6+1)->name << " is present in the sample." << endl;
#endif


            clade_bins.clear();
        }
    }

cerr << "Number of fragments in input file: " << n_reads << endl; 
cerr << "Number of mapped fragments: " << n_map_reads << endl;
cerr << "Number of fragments after filtering: " << sum_of_clades << endl; 

#ifdef VERBOSE_S

    cerr<<"done reading GAM file "<<gamfilename<<endl;
    cerr<<"total amount of reads "<<n_reads<<endl;
    cerr<<"amount of mapped reads "<<n_map_reads<<endl;
    cerr<<"A total of " <<n_map_reads - sum_of_clades<< " where excluded based on a low log-likelihood ratio\n" << endl;
#endif


    return make_pair(read_vec, clade_list_id);
}


