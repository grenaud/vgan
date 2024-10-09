#pragma once
#include "HaploCart.h"
#include "path_string.hpp"
#include <boost/range/adaptor/reversed.hpp>


inline const tuple<string, string, vector<int>> reconstruct_graph_sequence(const bdsg::ODGI &graph, const auto &path, const string &algnseq)

{
    string graph_seq = "";
    string read_seq = algnseq;
    vector<int> mppg_sizes;
    int mppg_counter=0;
    int edit_counter=0;

    string appended = "";

     // Loop over mappings to nodes
    string ps = vg::algorithms::path_string(graph, path);

      auto mppgs = path.mapping();
      int f=0;
      for (auto &mppg : mppgs) {
      string node_seq = graph.get_sequence(graph.get_handle(mppg.position().node_id(), mppg.position().is_reverse()));

      int aligned_length = 0; // Initialize aligned_length
      RepeatedPtrField<vg::Edit> ed = mppg.edit();
          edit_counter = 0;
          int offset = mppg.position().offset();
          // Loop over each edit in a node mapping
          for (const auto edit : ed) {
              const int32 to_length = edit.to_length();
              const int32 from_length = edit.from_length();
              const string edit_sequence = edit.sequence();
              

              // If we see an insertion in the read at the first edit of the first mapping, or last edit of the last mapping, it's a softclip.
              const bool softclip = (mppg_counter == 0 && offset == 0 && edit_counter == 0 && from_length == 0 && to_length > 0 && vg::io::edit_is_insertion(edit)) ||
                                (mppg_counter == mppgs.size()-1 && offset == 0 && edit_counter == ed.size() && from_length == 0 && to_length > 0 && vg::io::edit_is_insertion(edit));

          if (vg::io::edit_is_match(edit) || vg::io::edit_is_sub(edit)) {
              graph_seq += node_seq.substr(offset, from_length);
              aligned_length = node_seq.substr(offset, from_length).size();
             

                                                                        }

          else if (vg::io::edit_is_insertion(edit)) {
              if (softclip) {
                   // Annotate with an S for downstream use
                   graph_seq += string(to_length, 'S');
                   aligned_length = to_length;
                   

                            }
              else {
                   graph_seq += string(to_length, '-');
                   aligned_length = to_length;
    

                   }
                                                    }

          else if (vg::io::edit_is_deletion(edit)) {
              graph_seq += node_seq.substr(offset, from_length);
              aligned_length = node_seq.substr(offset, from_length).size();
              ps.insert(f, string(from_length, '-'));
                                                   }
          offset += from_length;
          f+=from_length;
          
          mppg_sizes.emplace_back(aligned_length);
                               }
                                          }


    const tuple<string, string, vector<int>> seq_tuple = make_tuple(graph_seq, ps, mppg_sizes);
    return seq_tuple;
}
