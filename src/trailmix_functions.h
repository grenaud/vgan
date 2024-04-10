#pragma once
#include "TrailMix.h"
#include <unordered_set>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <tuple>
#include <cfloat>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>

#define PRINTVEC(v) for (int i=0; i<v.size(); ++i){cerr << v[i] << '\t';}cerr << endl << endl;

const double TRANSITION_RATE = 23.0L;
const double TRANSVERSION_RATE = 1.0L;

std::unordered_map<int, int> inverse_map(const std::unordered_map<int, int>& original) {
    std::unordered_map<int, int> inverse;
    for (const auto& [key, value] : original) {
        inverse.insert(std::make_pair(value, key));
    }
    return inverse;
}

const int get_hg_index(const vector<string> paths, const string & path){
    auto it = find(paths.begin(), paths.end(), path);
    if (it != paths.end()) {
       return std::distance(paths.begin(), it);
                           }
    else {
        return -1;
         }
                                                                       }

void Trailmix::create_path_node_map(shared_ptr<Trailmix_struct> &dta) {
    unordered_map<int, int> path_node_map;
    for (size_t i = 0; i < dta->tree->nodes.size(); i++) {
        const auto it = std::find_if(dta->path_names.begin(), dta->path_names.end(), [&](const string &path) {
            return path == dta->tree->nodes[i]->longname.substr(0, dta->tree->nodes[i]->longname.find('.'));
        });
        if (it != dta->path_names.end()) {
            const int found_idx = std::distance(dta->path_names.begin(), it);
            path_node_map[found_idx] = i;
        }
    }
    dta->path_node_map = move(path_node_map);
    dta->node_path_map = inverse_map(dta->path_node_map);
    //ofstream out("pnm");
    //for (const auto & el : dta->path_node_map){out << dta->path_names[el.first] << '\t' << el.second << '\n';}
    //ofstream out2("npm");
    //for (const auto & el2 : dta->node_path_map){out2 << el2.first << '\t' << dta->path_names[el2.second] << '\n';}
    //throw runtime_error("TESTING TREE");
}

const unsigned int branch_to_gbwt(shared_ptr<Trailmix_struct> &dta, unsigned int branch){
    return dta->path_to_gbwt[dta->node_path_map[branch]];
}

const double avg(const Eigen::Matrix4d &matrix) {
  double sum = 0;
  int count = 0;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (true) {
        sum += matrix(i, j);
        count++;
      }
    }
  }

  double average = sum / count;

  // Check for inf and nan
  if (std::isinf(average) || std::isnan(average)) {
    cerr << "AVERAGE IS NOT FINITE" << endl;
    // handle the error or return a default value
    return 0; // for example
  }

  return average;
}

const long double get_adjusted_p_mutation(const long double mutation_rate, const long double placement){
    return mutation_rate + placement;
}

bool has_duplicates(const vector<unsigned int>& combo) {
  set<unsigned int> unique_elements(combo.begin(), combo.end());
  return unique_elements.size() != combo.size();
}

const long double mutation_ratio(const Eigen::Matrix4d& count_matrix) {
    long double ct_mutations = count_matrix(1, 3);
    long double ga_mutations = count_matrix(2, 0);
    long double tc_mutations = count_matrix(3, 1);
    long double ag_mutations = count_matrix(0, 2);

    long double numerator = ct_mutations + ga_mutations;
    long double denominator = tc_mutations + ag_mutations;

    if (denominator == 0) {
        // Handle potential division by zero
        return std::numeric_limits<long double>::infinity();
    }

    return numerator / denominator;
}


Eigen::Matrix4d interpolateMatrices(const Eigen::Matrix4d& mat1, const Eigen::Matrix4d& mat2, double fraction) {
    // Clamp the fraction between 0 and 1
    fraction = std::max(0.0, std::min(1.0, fraction));

    // Perform the linear interpolation
    Eigen::Matrix4d result = (1.0 - fraction) * mat1 + fraction * mat2;
    return result;
}


double log_multivariate_normal_pdf(const Eigen::Matrix4d& observed, const Eigen::Matrix4d& expected, const Eigen::Matrix4d& covariance) {
    const int k = 4; // Dimension of the matrices (4)

    Eigen::Vector4d diff;
    for (int i = 0; i < 4; ++i) {
        diff(i) = (observed.row(i) - expected.row(i)).squaredNorm();
    }
    Eigen::Matrix4d covariance_inv = covariance.inverse();
    double log_det_covariance = covariance.determinant();
    double constant_term = k * std::log(2 * M_PI);

    double log_pdf = -0.5 * (diff.transpose() * covariance_inv * diff + std::log(std::abs(log_det_covariance)) + constant_term);
    return log_pdf;
}


void PRINTMATRIX(const Matrix4d& matrix) {
    assert(matrix.rows() == 4 && matrix.cols() == 4);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << matrix(i, j) << '\t';
        }
        cout << endl;
    }
}

void Trailmix::place_sources(shared_ptr<Trailmix_struct> &dta){
    Trailmix::enumerate_branch_combos(dta);
    assert(!dta->branch_combos.empty());
    Trailmix::infer(dta);
                                                              }

long double dirichlet_multinomial_likelihood(const Eigen::Matrix4d &observed, const Eigen::Matrix4d &alpha) {
    //cerr << "in dirichlet multinomial" << endl;
    long double likelihood = 0.0;

    for (int i = 0; i < 4; ++i) {
        long double row_likelihood = 1.0;
        long double alpha_sum = alpha.row(i).sum();
        long double observed_sum = observed.row(i).sum();

        row_likelihood *= boost::math::tgamma(alpha_sum) / boost::math::tgamma(observed_sum + alpha_sum);

        for (int j = 0; j < 4; ++j) {
            row_likelihood *= boost::math::tgamma(observed(i, j) + alpha(i, j)) / boost::math::tgamma(alpha(i, j));
        }


        if (row_likelihood < 0.0){return -DBL_MAX;}

        likelihood += log(row_likelihood);
    }

    if (isnan(likelihood) || isinf(likelihood)){return -DBL_MAX;}
    return likelihood;
}


bool is_positive(const Eigen::VectorXd& vec) {
    return (vec.array() >= 0).all();
}

double constraint_function_lower(unsigned n, const double *x, double *grad, void *data) {
    double sum = 0.0;
    for (unsigned i = 0; i < n/2; ++i) {
        sum += x[i];
    }
    if (grad) {
        for (unsigned i = 0; i < n/2; ++i) {
            grad[i] = 1.0;
        }
    }
    return sum - 1.0;  // The constraint is satisfied when the result is >= 0
}


/*
const std::tuple<size_t, size_t, vector<bool>> get_top_sources(vector<vector<long double>> &model_lls, const vector<vector<vector<bool>>> &all_sources) {

    long double max_likelihood = -std::numeric_limits<long double>::infinity();
    size_t max_i = 0;
    size_t max_j = 0;

    // Iterate through the 2D model_lls vector
    for (size_t i = 0; i < model_lls.size(); ++i) {
        for (size_t j = 0; j < model_lls[i].size(); ++j) {
            if (!isfinite(model_lls[i][j])) {
               model_lls[i][j] = -std::numeric_limits<long double>::max();
                                        }
            if (model_lls[i][j] > max_likelihood) {
                max_likelihood = model_lls[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    return std::make_tuple(max_i, max_j, all_sources[max_i][max_j]);
}
*/
std::tuple<size_t, size_t, vector<bool>> get_top_sources(vector<vector<long double>> &model_lls, const vector<vector<vector<bool>>> &all_sources) {

    long double max_likelihood = -std::numeric_limits<long double>::infinity();
    size_t max_i = 0;
    size_t max_j = 0;

    // Iterate through the 2D model_lls vector
    for (size_t i = 0; i < model_lls.size(); ++i) {
        for (size_t j = 0; j < model_lls[i].size(); ++j) {
            if (!isfinite(model_lls[i][j])) {
               model_lls[i][j] = -std::numeric_limits<long double>::max();
            }
            if (model_lls[i][j] > max_likelihood) {
                max_likelihood = model_lls[i][j];
                max_i = i;
                max_j = j;
            }
            std::cout << "At source: " << i << ", model: " << j << ", likelihood: " << model_lls[i][j] << ", current max_likelihood: " << setprecision(15) << max_likelihood << std::endl;
        }
    }

    std::cout << "Max likelihood found at source: " << max_i << ", model: " << max_j << ", likelihood: " << setprecision(15) << max_likelihood << std::endl;

    // Display the corresponding source assignment in all_sources
    std::cout << "Corresponding source assignment in all_sources: ";
    for (bool val : all_sources[max_i][max_j]) {
        std::cout << val << ' ';
    }
    std::cout << std::endl;

    return std::make_tuple(max_i, max_j, all_sources[max_i][max_j]);
}

void Trailmix::write_output(std::shared_ptr<Trailmix_struct> &dta) {
    size_t best_k;
    size_t best_index;
    vector<bool> best_sources;
    const string outfile = dta->outputfilename.empty() ? "/dev/stdout" : dta->outputfilename;
    ofstream tm_out(outfile, ios::app);
    // Write header to TSV
    tm_out << "source_number\tchild_node\tparent_node\tchild_name\tparent_name\tbranch_lengths\tbranch_placements\tsource_props\tassignment" << std::endl;
    if (dta->auto_mode) {
        std::tie(best_k, best_index, best_sources) = get_top_sources(dta->model_lls, dta->all_sources);

        std::cerr << "OPTIMAL K: " << best_k + 1 << std::endl;

        for (size_t n = 0; n <= best_k; ++n) {
            tm_out << n + 1 << "\t";
            tm_out << dta->tree->nodes[dta->best_branch_combos[best_k][best_index][n]]->name << '\t';
            tm_out << dta->tree->nodes[dta->best_branch_combos[best_k][best_index][n]]->parent->name << '\t';
            tm_out << dta->gbwt_paths[branch_to_gbwt(dta, dta->tree->nodes[dta->best_branch_combos[best_k][best_index][n]]->name)] << '\t';
            tm_out << dta->gbwt_paths[branch_to_gbwt(dta, dta->tree->nodes[dta->best_branch_combos[best_k][best_index][n]]->parent->name)] << '\t';

            // Branch lengths
            tm_out << dta->tree->nodes[dta->best_branch_combo[n]]->dist << '\t';

            // Branch placements
            tm_out << dta->all_branch_placements[dta->max_branch_combo_index][n] * dta->all_branch_placements[dta->max_branch_combo_index][n];

            // Source proportions
            tm_out << dta->optimized_thetas[dta->max_branch_combo_index][n] << '\t';

            // Assignment
            std::string assignment = best_sources[n] ? "Ancient" : "Modern";
            tm_out << assignment;
            tm_out << endl;
        }
    } else {
        if(dta->best_branch_combo.empty()){cerr << "No branch combos, not writing output." << endl; exit(0);};
        for (size_t n = 0; n < dta->k; ++n) {
            tm_out << n + 1 << "\t";
            tm_out << dta->tree->nodes[dta->best_branch_combo[n]]->name << '\t';
            tm_out << dta->tree->nodes[dta->best_branch_combo[n]]->parent->name << '\t';
            tm_out << dta->gbwt_paths[branch_to_gbwt(dta, dta->tree->nodes[dta->best_branch_combo[n]]->name)] << '\t';
            tm_out << dta->gbwt_paths[branch_to_gbwt(dta, dta->tree->nodes[dta->best_branch_combo[n]]->parent->name)] << '\t';

            dta->paths_to_surject.emplace_back(dta->gbwt_paths[branch_to_gbwt(dta, dta->tree->nodes[dta->best_branch_combo[n]]->name)]);

            // Branch lengths
            tm_out << dta->tree->nodes[dta->best_branch_combo[n]]->dist << '\t';

            // Branch placements
            tm_out << dta->all_branch_placements[dta->max_branch_combo_index][n] * dta->all_branch_placements[dta->max_branch_combo_index][n] << '\t';

            // Source proportions
            tm_out << dta->optimized_thetas[dta->max_branch_combo_index][n] << '\t';

            // Assignment
            tm_out << "NA";
            tm_out << endl;
        }
    }
}


double log_sum_exp_of_logs(double log_x, double log_y) {
    double max_val = std::max(log_x, log_y);
    return max_val + log(exp(log_x - max_val) + exp(log_y - max_val));
}

long double log_sum_exp(long double x, long double y) {
    if (x > y) {
        return x + std::log1p(std::exp(y - x));
    } else {
        return y + std::log1p(std::exp(x - y));
    }
}

Eigen::Matrix4d normalize_matrix(const Eigen::Matrix4d& mat) {
    double sum = mat.sum();
    return mat / sum;
}

long double log_weighted_average(const std::vector<long double>& log_likelihoods, const std::vector<long double>& weights) {
    if (log_likelihoods.size() != weights.size()) {
        throw std::invalid_argument("log_likelihoods and weights must have the same size");
    }

    std::vector<long double> log_weights(log_likelihoods.size());
    for (size_t i = 0; i < log_likelihoods.size(); ++i) {
        log_weights[i] = log_likelihoods[i] + std::log(weights[i]);
    }

    long double log_weighted_sum = log_weights[0];
    for (size_t i = 1; i < log_weights.size(); ++i) {
        log_weighted_sum = log_sum_exp_of_logs(log_weighted_sum, log_weights[i]);
    }

    return log_weighted_sum;
}

std::unordered_set<int> removeEntries(std::shared_ptr<Trailmix_struct>& dta, const std::vector<unsigned int>& inputs) {
    std::unordered_set<int> removed_indices;

    for (auto it = dta->mapping_supports.begin(); it != dta->mapping_supports.end();) {
        std::size_t count = 0;
        for (const auto& input : inputs) {
            if (std::find(it->second.begin(), it->second.end(), input) != it->second.end()) {
                ++count;
            }
            if (count > 1) {
                removed_indices.insert(it->first);
                it = dta->mapping_supports.erase(it);
                break;
            }
        }
        if (count <= 1) {
            ++it;
        }
    }
    //std::cerr << "UNIQUE MAPPINGS: " << dta->mapping_supports.size() << std::endl;
    return removed_indices;
}


double dirichlet_logpdf(const std::vector<long double>& x, const std::vector<long double>& alpha) {
    double sum_alpha = 0.0;
    double sum_x = 0.0;
    double logpdf = 0.0;
    int k = x.size();
    // Calculate the sum of the alpha values
    for (int i = 0; i < k; i++) {
        sum_alpha += alpha[i];
    }
    // Calculate the sum of the x values
    for (int i = 0; i < k; i++) {
        sum_x += x[i];
    }
    // Calculate the log of the beta function
    double log_beta = 0.0;
    for (int i = 0; i < k; i++) {
        log_beta += std::lgamma(alpha[i]);
    }
    log_beta -= std::lgamma(sum_alpha);
    // Calculate the log of the Dirichlet PDF
    for (int i = 0; i < k; i++) {
        logpdf += (alpha[i] - 1) * std::log(x[i]) - std::lgamma(alpha[i]);
    }
    logpdf += std::lgamma(sum_alpha);
    logpdf -= log_beta;
    // Calculate the log of the normalization constant
    double log_norm = -std::lgamma(sum_x);
 // Return the log of the Dirichlet PDF with normalization
    //assert(!isinf(logpdf + log_norm) && !isnan(logpdf + log_norm));
    if (isinf(logpdf + log_norm) || isnan(logpdf + log_norm)){return -DBL_MAX;}
    return logpdf + log_norm;
}


struct OptData {
    vector<reference_wrapper<Matrix4d>> observed_child;
    vector<reference_wrapper<Matrix4d>> observed_parent;
    vector<unsigned int> branch_combo;
    vector<double> branch_lengths;
    double read_coverage;
    unsigned int counter = 0;
    vector<bool> assignments;
    bool printed = false;
    shared_ptr<Trailmix_struct> dta;
};

void Trailmix::compute_read_coverage(shared_ptr<Trailmix_struct> &dta){
    long double n_aligned_bases = 0.0;
    for (unsigned int i=0; i<dta->algnvector->size(); ++i){
        n_aligned_bases += (long double)(dta->algnvector->at(i)->seq.size());
                                                         }

    dta->read_coverage = n_aligned_bases / 16565;
                                                                      }


const long double Trailmix::sum_log_likelihoods(vector<long double> &log_lik_vec) {

// Sum log likelihoods using Gabriel magic
long double ret = log_lik_vec[0];
for (int i=1; i!=log_lik_vec.size(); ++i) {
    ret = oplusInitnatl(ret, log_lik_vec[i]);
                                          }
return ret;
}


void rescale(std::vector<std::vector<long double>>& samples) {
  // Compute the mean and standard deviation of each sample
  std::vector<long double> means(samples.size());
  std::vector<long double> std_devs(samples.size());
  for (std::size_t i = 0; i < samples.size(); ++i) {
    long double sum = 0.0;
    for (auto x : samples[i]) {
      sum += x;
    }
    means[i] = sum / static_cast<long double>(samples[i].size());

    long double sum_sq_diff = 0.0;
    for (auto x : samples[i]) {
      sum_sq_diff += std::pow(x - means[i], 2.0);
    }
    std_devs[i] = std::sqrt(sum_sq_diff / static_cast<long double>(samples[i].size() - 1));
  }

  // Rescale each sample to have the same mean and standard deviation while preserving ordering
  for (std::size_t i = 0; i < samples.size(); ++i) {
    std::sort(samples[i].begin(), samples[i].end());
    long double med = samples[i][samples[i].size() / 2];
    long double std_dev = std_devs[i];
    long double scaling_factor = 100.0 / means[i]; // <-- adjust scaling factor to desired mean
    for (std::size_t j = 0; j < samples[i].size(); ++j) {
      samples[i][j] = scaling_factor * std::exp((samples[i][j] - med) / std_dev) * std::exp(med);
    }
  }
}


const int find_index(const std::vector<std::pair<std::string, double>>& vec, const std::string& str) {
    auto it = std::find_if(vec.begin(), vec.end(), [&str](const auto& p){ return p.first == str; });
    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    }
    else {
        return -1; // indicate not found
    }
}

void normalize(vector<long double>& v) {
    long double sum = 0.0;
    for (long double x : v) {
        sum += x;
    }
    for (long double& x : v) {
        x /= sum;
    }
}

vector<long double> normalize_copy(vector<long double>& v) {
    vector<long double> ret;
    long double sum = 0.0;
    for (long double x : v) {
        sum += x;
    }
    for (long double& x : v) {
        ret.emplace_back(min(max(0.0L, x / sum), 1.0L));
    }
    return ret;
}

void Trailmix::get_seed_source_estimates(shared_ptr<Trailmix_struct> &dta, const vector<unsigned int> &branch_combo){

vector<int> paths;
for (const auto & el : branch_combo)    {
    //cerr << "path: " << dta->path_names[dta->node_path_map[el]] << endl;
    //cerr << "path id: " << dta->node_path_map[el] << endl;
    paths.emplace_back(dta->node_path_map[el]);
                                        }


assert(paths.size() == dta->k);

vector<long double> props(paths.size(), 0.0);

for (unsigned int source=0; source < paths.size(); ++source){
    const string path_name = dta->path_names[paths[source]];
    for (unsigned int i = 0; i < dta->tpms.size(); ++i) {
          const int find_idx = find_index(dta->tpms, path_name);
          assert(find_idx != -1);
          props[source] = dta->tpms[find_idx].second;
                                                        }
                                                            }
cerr << "SEED: " << endl;
PRINTVEC(props)

if (isnan(props[0])){
    for (auto &el : props){el = (1.0 / (dta->k));}
}

normalize(props);
cerr << "SEED NORMALIZED: " << endl;
PRINTVEC(props)
assert(props.size() == dta->k);
dta->seed = props;
return;
                                                                          }


const double Trailmix::get_branch_combo_posterior(shared_ptr<Trailmix_struct> &dta, unsigned int branch_combo_idx){
    // Given index of branch combo, extract the corresponding hap combo posterior probability from RPVG
    //cerr << "n: " << n << endl;
    int current_idx = 0;
    for (const auto &v : dta->branch_combos){
        if (v.size() + current_idx >= branch_combo_idx){
            //cerr << "Returning index: " << current_idx << endl;
            return dta->hap_combo_posteriors[current_idx];
                                        }
        else{
            current_idx += v.size();
            }
                                            }

    throw std::runtime_error("We should have returned something but didn't");
                                                                                                                   }


void Trailmix::branch_combo_cartesian_product(vector<vector<unsigned int>>& v, shared_ptr<Trailmix_struct>& dta) {
  auto product = [](long long a, vector<unsigned int>& b) { return a * b.size(); };
  const long long N = accumulate(v.begin(), v.end(), 1LL, product);
  vector<unsigned int> u(v.size());

  for (long long n = 0; n < N; ++n) {
    lldiv_t q{n, 0};
    for (long long i = v.size() - 1; 0 <= i; --i) {
      q = div(q.quot, v[i].size());
      u[i] = v[i][q.rem];
    }
    if (!has_duplicates(u)) {
      assert(u.size() == dta->k);
      dta->branch_combos.emplace_back(u);
    }
  }
}

void Trailmix::enumerate_branch_combos(shared_ptr<Trailmix_struct> &dta){

    if (dta->debug) {cerr << "Node combo size: " << dta->node_combos.size() << endl;}
    for (const auto & combo : dta->node_combos){

        assert(combo.size() == dta->k);

        vector<vector<unsigned int>> options(combo.size(), vector<unsigned int>());
        for (size_t i=0; i<combo.size(); ++i){
            options[i].emplace_back(combo[i]->name);
            //cerr << "Number of children: " << combo[i]->nchildren << endl;
            for (unsigned int j=0; j<combo[i]->nchildren; ++j){
                options[i].emplace_back(combo[i]->children[j]->name);
                                                              }
                                             }

       // Exhaustively list all combinations

       assert(!options.empty());
       dta->branch_combos.clear();
       branch_combo_cartesian_product(options, dta);
       if (dta->debug){cerr << "Removing duplicate branch combos ..." << endl;}
       Trailmix::remove_duplicate_branch_combos(dta->branch_combos);
       if (dta->debug){cerr << "Duplicate branch combos removed" << endl;}
                                                }

                                                                         }

void Trailmix::get_dists_on_branch(const std::vector<unsigned int> &branch_combo, std::shared_ptr<Trailmix_struct> &dta) {
    for (const auto& source : branch_combo) {
        const long double branch_length = dta->tree->nodes[source]->dist;
        std::vector<long double> branch_dists;
        const unsigned int n_points = 10;
            for (double point = 0; point < n_points - 1; ++point) {
                const long double frac = (point+1)/n_points;
                branch_dists.emplace_back((branch_length * frac));
                                                                  }
        //cerr << "PLACING BACK DISTS: " << endl;
        //PRINTVEC_HIGH_PRECISION(branch_dists)
        dta->bc_dists.emplace_back(branch_dists);
                                            }
}

void Trailmix::remove_duplicate_branch_combos(vector<vector<unsigned int>> & branch_combos){
    sort(branch_combos.begin(), branch_combos.end());
    branch_combos.erase(std::unique(branch_combos.begin(), branch_combos.end()), branch_combos.end());
}

const long double log_sum_exp(const std::vector<long double>& log_vals) {
    assert(!log_vals.empty());
    if (log_vals.size() == 1) {
        return log_vals[0];
    }
    auto max_val_iter = std::max_element(log_vals.begin(), log_vals.end());
    long double max_val = *max_val_iter;
    long double sum = 0.0;
    for (const auto& val : log_vals) {
        if (val != max_val) {
            sum += exp(val - max_val);
        }
    }
    return max_val + log(1.0 + sum);
}


bool validate_data(shared_ptr<Trailmix_struct>& dta, int num_sources) {
    for (int k = 0; k < num_sources; ++k) {
        if (dta->model_lls[k].size() != dta->all_sources[k].size()) {
            cerr << "Mismatch at k = " << k << ": model_lls[k] size is " << dta->model_lls[k].size() << ", all_sources[k] size is " << dta->all_sources[k].size() << endl;
            return false;
        }
        for (size_t i = 0; i < dta->model_lls[k].size(); ++i) {
            if (dta->all_sources[k][i].size() != k + 1) {
                cerr << "Mismatch at k = " << k << ", i = " << i << ": all_sources[k][i] size is " << dta->all_sources[k][i].size() << ", expected size is " << k + 1 << endl;
                return false;
            }
        }
    }
    return true;
}


void housekeeping(shared_ptr<Trailmix_struct> &dta){
        for (Eigen::Matrix4d& matrix : dta->sub_vec) {
            matrix.setZero();
                                                     }
        dta->to_increment.clear();
        dta->to_increment_ancient.clear();
        dta->ancient_counts = 0;
        dta->modern_counts = 0;
}

void Trailmix::generate_combinations(shared_ptr<Trailmix_struct> &dta, int num_sources, int current_pos) {
    if (current_pos == num_sources) {
        cerr << "Current combination: ";
        for (bool source : dta->sources) {
            cerr << (source ? "\033[1;34m1 (ancient)\033[0m " : "\033[1;32m0 (modern)\033[0m ");
        }

        housekeeping(dta);
        run_trailmix(dta);
        dta->first_of_new_k=false;
        cerr << "MODEL_LIKELIHOOD: " << dta->cur_model_likelihood << endl;
        if (dta->cur_model_likelihood > dta->best_model_likelihood) {
            dta->best_model_likelihood = dta->cur_model_likelihood;
        }
        if (dta->k - 1>= dta->model_lls.size()) {
            std::vector<long double> innerVector;
            std::vector<std::vector<bool>> innerSources;
            dta->model_lls.emplace_back(innerVector);
            dta->all_sources.emplace_back(innerSources);
        }
        dta->model_lls[dta->k - 1].emplace_back(dta->cur_model_likelihood);
        dta->all_sources[dta->k - 1].emplace_back(dta->sources);
        assert(validate_data(dta, num_sources));
        cerr << "Storing likelihood for k: " << dta->k << ", likelihood: " << dta->cur_model_likelihood << endl;
        cerr << endl;
    } else {
        dta->sources[current_pos] = false;
        dta->current_source=0;
        generate_combinations(dta, num_sources, current_pos + 1);
        dta->sources[current_pos] = true;
        dta->current_source=1;
        generate_combinations(dta, num_sources, current_pos + 1);
    }
}

void Trailmix::mcmc_setup(shared_ptr<Trailmix_struct> &dta){
    cerr << "Setting up for MCMC..." << endl;
    Trailmix::enumerate_branch_combos(dta);
    assert(!dta->branch_combos.empty());
    // Iterate through the vectors and insert the index pairs into the map
    for (unsigned int i = 0; i < dta->gbwt_paths.size(); ++i) {
        for (unsigned int j = 0; j < dta->path_names.size(); ++j) {
            if (dta->gbwt_paths[i].substr(0, dta->gbwt_paths[i].find('.')) == dta->path_names[j]) {
                dta->path_to_gbwt[j] = i;
                break;
            }
        }
    }
    for (const auto &entry : dta->path_to_gbwt) {
        dta->gbwt_to_path[entry.second] = entry.first;
                                                }

    Trailmix::compute_read_coverage(dta);
    unsigned int branch_combo_idx=0;
    double best_branch_combo_ll = -DBL_MAX;
    double best_branch_dist = 0.0;

    cerr << endl << endl;
    dta->branch_combo_lls.clear();
    if (dta->branch_combos.empty()){
        cerr << "NO BRANCH COMBOS" << endl;
        dta->cur_model_likelihood = -DBL_MAX;
        return;
                                   }

        for (unsigned int branch_combo_idx=0; branch_combo_idx<dta->branch_combos.size(); ++branch_combo_idx){
           assert(dta->branch_combos[branch_combo_idx].size() == dta->k);
           cerr << "On branch combo " << branch_combo_idx+1 << " out of " << dta->branch_combos.size() << endl;
           Trailmix::get_seed_source_estimates(dta, dta->branch_combos[branch_combo_idx]);
           assert(!dta->branch_combos[branch_combo_idx].empty());
           vector<long double> branch_lengths;
           for (auto el : dta->branch_combos[branch_combo_idx]){branch_lengths.emplace_back(dta->tree->nodes[el]->dist);}
           //dta->mcmc_chain.branch_lengths = move(branch_lengths);
                                                                                                             }

         if (dta->k-1 >= dta->best_branch_combos.size()) {
            dta->best_branch_combos.resize(dta->k);
                                                         }

           cerr << "MCMC setup done" << endl;
}


void Trailmix::infer(shared_ptr<Trailmix_struct> &dta){

    // Iterate through the vectors and insert the index pairs into the map
    for (unsigned int i = 0; i < dta->gbwt_paths.size(); ++i) {
        for (unsigned int j = 0; j < dta->path_names.size(); ++j) {
            if (dta->gbwt_paths[i].substr(0, dta->gbwt_paths[i].find('.')) == dta->path_names[j]) {
                dta->path_to_gbwt[j] = i;
                break;
            }
        }
    }
    for (const auto &entry : dta->path_to_gbwt) {
        dta->gbwt_to_path[entry.second] = entry.first;
                                                }

    //dta->theta = vector<long double>{0.2, 0.8};
    Trailmix::compute_read_coverage(dta);
}
