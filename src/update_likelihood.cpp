#pragma once
#include "HaploCart.h"
#include "vgan_utils.h"
#include "process_mapping.cpp"
#include <gbwtgraph/gbwtgraph.h>
#include "path_string.hpp"
//#include <boost/range/adaptor/reversed.hpp>


const vector<double> Haplocart::update(shared_ptr<Trailmix_struct> &dta, const int i, vector<double> log_likelihood_vec) noexcept {
    const vector<double> ret = Haplocart::update_likelihood(dta, dta->algnvector->at(i), log_likelihood_vec);
    //if (i%100==0){cerr << "ON READ: " << i << endl;}
    return ret;
}


inline const vector<double> Haplocart::update_likelihood(shared_ptr<Trailmix_struct> &dta, const AlignmentInfo* read_info,
                             vector<double> log_likelihood_vec) noexcept
{
    const auto path = read_info->path;
    const tuple<string, string, vector<int>> seq_tuple = reconstruct_graph_sequence(dta->graph, path, read_info->seq);
    const vector<int> mppg_sizes = get<2>(seq_tuple);

    bool is_reverse = path.mapping(0).position().is_reverse();


      for (unsigned int path_idx = 0; path_idx < dta->path_names.size(); ++path_idx) {

        int position_in_read = 0;

        for (size_t i = 0; i < path.mapping().size(); ++i) {
            vector<int> quality_scores;
            quality_scores.reserve(read_info->seq.size());

            //cerr << "position in read: " << position_in_read << endl;
            size_t read_seq_len = (position_in_read + mppg_sizes[i] <= read_info->seq.size()) ? mppg_sizes[i] : read_info->seq.size() - position_in_read;
            //cerr << "read seq len: " << read_seq_len << endl;
            string read_seq = read_info->seq.substr(position_in_read, read_seq_len);

            string graph_seq = get<0>(seq_tuple).substr(position_in_read, mppg_sizes[i]);

            for (size_t j = position_in_read; j < position_in_read + read_seq.size(); ++j) {
                const int qscore = int(read_info->quality_scores[j]);
                if (qscore >= 90) {dta->use_background_error_prob = true;}
                quality_scores.emplace_back(qscore);
            }

            log_likelihood_vec = process_mapping(dta, read_info, log_likelihood_vec, path.mapping()[i], read_seq, quality_scores, graph_seq, path_idx);

            position_in_read += read_seq.size();
        }
    }

    return log_likelihood_vec;
}

