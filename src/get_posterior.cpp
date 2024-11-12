#include "HaploCart.h"
#include "TrailMix.h"


void Haplocart::write_posterior_log(const string &samplename, const string &posteriorfilename, vector<string> clade_vec, vector<double> confidence_vec,
                                    const bool &webapp) {

    ofstream posteriorFile(posteriorfilename, ios::app);

    if (webapp == false) {
        for (size_t i = 0; i < clade_vec.size(); ++i) {posteriorFile << samplename << '\t' << clade_vec[i] << '\t' << confidence_vec[i] << '\t' << i << '\n';}
                         }

    else {
        cout << "<table>" << '\n';
        for (int i = 0; i < clade_vec.size(); ++i) {
            cout << "<tr><td>" << samplename << "</td><td>" << clade_vec[i] << "</td><td>" << confidence_vec[i] << "</td><td>" << i << '\n';
                                                   }
        posteriorFile << "<table>" << '\n';
        for (int i = 0; i < clade_vec.size(); ++i) {
            posteriorFile << "<tr><td>" << samplename << "</td><td>" << clade_vec[i] << "</td><td>" << confidence_vec[i] << "</td><td>" << i << '\n';
                                                   }
        posteriorFile << "<table>" << '\n';
        cout << "<table>" << '\n';
         }
                                                                                                                         }


const set<string> Haplocart::get_children(set<string> preds, map<string, vector<string>> children) {

// Given a set of haplotypes, return all children
set<string> all_child_set;
for (const string & p : preds) {
    if (children.find(p)->second.size() > 0) {
        vector<string> child_vec = children.find(p)->second;
        for (const string & child : child_vec) {
            all_child_set.insert(child);
                                             }
                                             }
                                  }
return all_child_set;
                                                                                                         }

const vector<double> Haplocart::get_posterior_of_clade(vector<double> all_top, vector<double> final_vec, set<string> preds,
                                                             map<string, vector<string>> children, const double total_ll,
                                                             const vector<string> path_names, bool initial) {

    // Get sum of the log likelihoods for each haplotype subsumed within a given clade
    set<string> child_set;
    if (initial == false) {
        child_set = get_children(preds, children);

        int idx = 0;
        for (const string & path: path_names) {
                auto it = find(child_set.begin(), child_set.end(), path);
                if (it != child_set.end()) {
                    all_top.emplace_back(final_vec[idx]);
                                           }
                idx += 1;
                                              }

    if (child_set.size() > 0) {
        all_top = get_posterior_of_clade(all_top, final_vec, child_set, children, total_ll, path_names, initial);
                              }
                          }

    return all_top;

                                                                                              }

const double Haplocart::sum_log_likelihoods(vector<double> log_lik_vec) {

// Sum log likelihoods using Gabriel magic
double ret = log_lik_vec[0];
for (int i=1; i!=log_lik_vec.size(); ++i) {
    ret = oplusInitnatl(ret, log_lik_vec[i]);
                                          }
return ret;
}

void Haplocart::get_posterior(const vector<double> &final_vec, const vector<string> path_names,
                              map<string, vector<string>> parents, map<string, vector<string>> children,
                              const string & samplename, const string &predicted_haplotype, const string &posteriorfilename,
                              const bool webapp) {


   double total_ll = Haplocart::sum_log_likelihoods(final_vec);
   vector<string> parent_vec;
   if (parents.count(predicted_haplotype) > 0) {
       parent_vec = parents.find(predicted_haplotype)->second;
                                               }

   vector<string> clade_vec{predicted_haplotype};
   vector<double> confidence_vec;
   set<string> pred{predicted_haplotype};
   int predicted_haplotype_idx = find(path_names.begin(), path_names.end(), predicted_haplotype) - path_names.begin();
   vector<double> all_top{final_vec[predicted_haplotype_idx]};
   double considered_ll = Haplocart::sum_log_likelihoods(all_top);
   double top_ratio = exp(considered_ll - total_ll);
   confidence_vec.emplace_back(top_ratio);
   all_top.clear();
   pred.clear();
   for (int j=0; j<parent_vec.size(); ++j) {
       if (parent_vec[j] != parent_vec[j-1]) {
           clade_vec.emplace_back(parent_vec[j]);
                                             }
       pred.insert(parent_vec[j]);
       all_top = get_posterior_of_clade(all_top, final_vec, pred, children, total_ll, path_names, false);
       considered_ll = sum_log_likelihoods(all_top);
       top_ratio = exp(considered_ll - total_ll);
       if (parent_vec[j] != parent_vec[j-1]) {
           confidence_vec.emplace_back(top_ratio);
                                             }
       all_top.clear();
       pred.clear();

                                           }

   Haplocart::write_posterior_log(samplename, posteriorfilename, clade_vec, confidence_vec, webapp);

                                                            }
