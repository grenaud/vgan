#pragma once
#include "path_string.hpp"
#include <boost/range/adaptor/reversed.hpp>

const tuple<string, string, vector<int>> reconstruct_graph_sequence(const bdsg::ODGI &graph, const auto &path, const string &algnseq)

{
    string graph_seq = "";
    string read_seq = "";
    vector<int> mppg_sizes;
    int mppg_counter=0;
    int edit_counter=0;
    bool is_reverse = false;

     // Loop over mappings to nodes

      auto mppgs = path.mapping();
      for (auto mppg : mppgs) {

      auto handle = graph.get_handle(mppg.position().node_id(), mppg.position().is_reverse());
      string node_seq = graph.get_sequence(handle);
      mppg_sizes.emplace_back(node_seq.size());
      RepeatedPtrField<vg::Edit> ed = mppg.edit();
          int offset = mppg.position().offset();
          // Loop over each edit in a node mapping
          for (const auto edit : ed) {
          const int32 to_length = edit.to_length();
          const int32 from_length = edit.from_length();
          const string edit_sequence = edit.sequence();

          // If we see an insertion in the read at the first edit of the first mapping, or last edit of the last mapping, it's a softclip.
          const bool softclip = (mppg_counter == 0 && edit_counter == 0 && from_length == 0 && to_length > 0) ||
                                (mppg_counter == mppgs.size() - 1 && edit_counter == ed.size() - 1 && from_length == 0 && to_length > 0);
          if (vg::io::edit_is_match(edit)) {
              graph_seq += node_seq.substr(offset, from_length);
              read_seq += node_seq.substr(offset, from_length);
                                                                        }
          else if (vg::io::edit_is_sub(edit)) {
              graph_seq += node_seq.substr(offset, from_length);
              read_seq += edit_sequence;
                                               }
          else if (vg::io::edit_is_insertion(edit)) {
              if (softclip) {
                   // Annotate with an S for downstream use
                   graph_seq += string(to_length, 'S');
                   read_seq += algnseq.substr(0, to_length);
                            }
              else {
                   graph_seq += string(to_length, '-');
                   read_seq += node_seq.substr(0, to_length);
                   }
                                                    }
          else if (vg::io::edit_is_deletion(edit)) {
              graph_seq += node_seq.substr(offset, from_length);
              read_seq += string(from_length, '-');
                                                   }
          offset += from_length;
                               }
                                          }

    const tuple<string, string, vector<int>> seq_tuple = make_tuple(graph_seq, read_seq, mppg_sizes);
    return seq_tuple;
}
