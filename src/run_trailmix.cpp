#include "TrailMix.h"
#include "HaploCart.h"
#include "gam2prof.h"
#include "time.h"
#include <algorithm>
#include <functional>
#include "trailmix_functions.h"
#include "rpvg_main.hpp"
#include "damage.h"
#include "sys/wait.h"
#include "MCMC.h"

char** convert_to_char_array(const vector<string>& arguments) {
    char** argvtopass = new char*[arguments.size()];
    for (size_t i = 0; i < arguments.size(); i++) {
        argvtopass[i] = const_cast<char*>(arguments[i].c_str());
    }
    return argvtopass;
}

void Trailmix::run_trailmix(shared_ptr<Trailmix_struct>& dta) {
    MCMC mcmc;

    // Haplocart run
    //run_haplocart(dta);

    if (!dta->first_of_new_k) {
        goto place_sources;
    }

    // Run gam2prof
    //run_gam2prof(dta);

    // Run RPVG haplotypes
    mcmc.run_rpvg_haplotypes(dta);

    // Run RPVG haplotype-transcripts
    //mcmc.run_rpvg_haplotype_transcripts(dta);

    // Run Haplocart again
    //run_haplocart(dta);

    while ((dta->rpvg_pid = wait(&dta->rpvg_status)) > 0);
    while ((dta->rpvg_ht_pid = wait(&dta->rpvg_ht_status)) > 0);

    //exit(0);

    dta->rpvg_status = 0;
    dta->rpvg_ht_status = 0;

    dta->hap_combos.clear();
    dta->hap_combo_posteriors.clear();
    dta->node_combos.clear();

    if (dta->debug) {
        cerr << "Loading RPVG output" << endl;
    }

    //load_hap_combos(dta);
    //load_tpms(dta);
    //assert(!dta->tpms.empty());

    if (dta->debug) {
        cerr << "Number of hap combos: " << dta->hap_combos.size() << endl;
    }
    create_node_combos(dta);

place_sources:

    cerr << "LOADING READ PROBS" << endl;
    //load_read_probs(dta);

    cerr << "PLACING SOURCES..." << endl;
    // place_sources(dta);

    // Run MCMC
    //run_mcmc(dta);
}

void Trailmix::run_haplocart(shared_ptr<Trailmix_struct>& dta) {
    vector<string> hcarguments;

    if (!dta->fastafilename.empty()) {
        hcarguments = {"vgan", "haplocart", "-f", dta->fastafilename};
    } else if (!dta->fastq1filename.empty() && dta->fastq2filename.empty()) {
        hcarguments = {"vgan", "haplocart", "-fq1", dta->fastq1filename};
    } else if (!dta->gamfilename.empty()) {
        hcarguments = {"vgan", "haplocart", "-g", dta->gamfilename};
    }

    char** hcargvtopass = convert_to_char_array(hcarguments);

    cerr << "RUNNING HAPLOCART" << endl;
    Haplocart().run(hcarguments.size(), hcargvtopass, dta);
    hcarguments.clear();
    cerr << "\n" << endl;

    delete[] hcargvtopass;
}

void Trailmix::run_gam2prof(shared_ptr<Trailmix_struct>& dta) {
    vector<string> g2parguments = {"vgan", "gam2prof", "--running-trailmix"};
    char** g2pargvtopass = convert_to_char_array(g2parguments);

    Gam2prof().run(g2parguments.size(), g2pargvtopass, dta->cwdProg, dta);

    delete[] g2pargvtopass;
}

void MCMC::run_rpvg_haplotypes(shared_ptr<Trailmix_struct>& dta) {

    cerr << "RPVG GAMFILE: " << dta->rpvg_gamfilename << endl;

    assert(std::filesystem::exists(dta->rpvg_gamfilename));

    vector<string> rpvgarguments = {"rpvg", "-g", dta->graph_dir + dta->graph_prefix + ".og", "-p", dta->graph_dir + dta->graph_prefix + ".gbwt",
                                    "-a", dta->rpvg_gamfilename, "-o", dta->tmpdir + "rpvg_hap", "-i", "haplotypes",
                                    "--use-hap-gibbs", "-u", "-s", "-t", to_string(dta->n_threads), "--filt-best-score", "1e-30",
                                    "--min-noise-prob", "1e-30"};
/*
    if (dta->rng_seed != "NONE") {
        rpvgarguments.push_back("-r");
        rpvgarguments.push_back(dta->rng_seed);
    }
*/

    rpvgarguments.insert(rpvgarguments.end(), {"-y", to_string(dta->k), "-m", dta->mu, "-d", dta->sigma, "--score-not-qual",
                                               "-n", "1", "--vgan-temp-dir", dta->tmpdir, "--gibbs-thin-its", "1",
                                               "--min-noise-prob", "1e-30", "--prob-precision", "1e-30"});

/*
    if (dta->strand_specific) {
        rpvgarguments.insert(rpvgarguments.end(), {"-e", "fr"});
    }
*/

    char** rpvgargvtopass = convert_to_char_array(rpvgarguments);

    cerr << "RUNNING RPVG FOR K= " << dta->k << endl;
    int retcode1 = rpvg_main(rpvgarguments.size(), rpvgargvtopass);
    cerr << "DONE!" << endl;
    rpvgarguments.clear();

    delete[] rpvgargvtopass;
}

void MCMC::run_rpvg_haplotype_transcripts(shared_ptr<Trailmix_struct>& dta) {
    vector<string> rpvgarguments = {"rpvg", "-g", dta->graph_dir + dta->graph_prefix + ".xg", "-p", dta->graph_dir + dta->graph_prefix + ".gbwt",
                                    "-a", dta->rpvg_gamfilename, "-o", dta->tmpdir + "rpvg_ht", "-i", "haplotype-transcripts",
                                    "-f", dta->graph_prefix + "/pantranscriptome.txt", "--use-hap-gibbs", "-u", "-s", "-t",
                                    to_string(dta->n_threads)};

    if (dta->rng_seed != "NONE") {
        rpvgarguments.push_back("-r");
        rpvgarguments.push_back(dta->rng_seed);
    }

    rpvgarguments.insert(rpvgarguments.end(), {"-y", to_string(dta->k), "-m", dta->mu, "-d", dta->sigma, "--score-not-qual",
                                               "--vgan-temp-dir", dta->tmpdir, "-b", "--use-hap-gibbs", "--max-em-its", "1"});

    if (dta->strand_specific) {
        rpvgarguments.insert(rpvgarguments.end(), {"-e", "fr"});
    }

    char** rpvgargvtopass = convert_to_char_array(rpvgarguments);

    int retcode2 = rpvg_main(rpvgarguments.size(), rpvgargvtopass);
    rpvgarguments.clear();

    delete[] rpvgargvtopass;
}

// void Trailmix::run_mcmc(shared_ptr<Trailmix_struct>& dta) {
//     //RunTreeProportionParams params(dta);
//     //std::vector<std::vector<MCMCiteration>> MCMCiterationsVec(params.chains);

//     // Initialize the maximum posterior and corresponding node combo
//     double max_posterior = -1.0;
//     std::vector<int> max_node_combo_names;

//     // Find the node combo with the highest posterior
//     for (size_t i = 0; i < dta->node_combos.size(); ++i) {
//         double current_posterior = dta->hap_combo_posteriors[i];
//         if (current_posterior > max_posterior) {
//             max_posterior = current_posterior;

//             // Create a new vector to store the names of all nodes in the node_combo
//             max_node_combo_names.clear();
//             for (auto& node : dta->node_combos[i]) {
//                 max_node_combo_names.push_back(node->name);
//             }
//         }
//     }

    //params.chains = dta->chains;
    //MCMC mcmc;

//     for (unsigned int chain = 0; chain < params.chains; ++chain) {
//         cerr << "ON CHAIN: " << chain + 1 << "  OUT OF " << params.chains << endl;
//         params.tr = dta->tree;
//         params.soibean = false;
//         params.root = dta->tree->root;
//         params.chains = dta->chains;
//         params.maxIter = dta->iter;
//         params.burn = dta->burnin;
//         params.sources = vector<int>{42}; // max_node_combo_names;
//         params.align = dta->gam;
//         //mcmc.run_tree_proportion(params, MCMCiterationsVec[chain]);
//     }

//     if (dta->auto_mode) {
//         dta->rpvg_status = -1;
//         dta->rpvg_ht_status = -1;
//     }

//     //mcmc.processMCMCiterations(MCMCiterationsVec);
// }

void Trailmix::create_node_combos(shared_ptr<Trailmix_struct>& dta) {
    for (const auto& combo : dta->hap_combos) {
        vector<spidir::Node*> single_combo;
        single_combo.clear();
        for (size_t i = 0; i < combo.size(); ++i) {
            if (i < combo.size() - 2) {
                const int found_idx = find(dta->path_names.begin(), dta->path_names.end(), combo[i]) - dta->path_names.begin();
                assert(find(dta->path_names.begin(), dta->path_names.end(), combo[i]) != dta->path_names.end());
                single_combo.push_back(dta->tree->nodes[dta->path_node_map[found_idx]]);
            } else if (i == combo.size() - 1) {
                dta->hap_combo_posteriors.push_back(stod(combo[i]));
            }
        }
        assert(single_combo.size() == dta->k);
        dta->node_combos.push_back(single_combo);
    }
}



