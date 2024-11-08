#include "Dup_Remover.h"
#include "../dep/vg/deps/progress_bar/progress_bar.hpp"
#include <vg/io/stream.hpp>
#include "readGAM.h"

Dup_Remover::Dup_Remover(){

}

Dup_Remover::~Dup_Remover(){

}

string_view Dup_Remover::usage(){
    return string("vgan rmdup [sorted_input.gam] > [output.gam]\n") +
           string("Remove duplicates from GAM file.\n") +
           string("Note that the GAM file must be sorted. If it is not please run vg gamsort first");
                                }


void Dup_Remover::update_for_single(shared_ptr<vector<bool>> &is_dup, const bool is_reverse, shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector, \
                                             long unsigned int i, const int n_threads) {

    const pair<long int, long int> start_coords = make_pair(sorted_algnvector->at(i)->path.mapping()[0].position().node_id(),
                                 sorted_algnvector->at(i)->path.mapping()[0].position().offset());

    pair<long int, long int> candidate_start_coords;

    for (size_t j=i+1; j!=sorted_algnvector->size(); ++j){

             candidate_start_coords = make_pair(sorted_algnvector->at(j)->path.mapping()[0].position().node_id(),
                                                sorted_algnvector->at(j)->path.mapping()[0].position().offset());

        if (i<j && candidate_start_coords == start_coords) {
            is_dup->at(j)=true;
                                                           }

                                                               }
}

void Dup_Remover::update_for_paired(shared_ptr<vector<bool>> &is_dup, const bool is_reverse, shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector,\
                                    const long unsigned int i, const int n_threads) {

    const pair<long int, long int> start_coords = make_pair(sorted_algnvector->at(i)->path.mapping()[0].position().node_id(),
                                 sorted_algnvector->at(i)->path.mapping()[0].position().offset());

    const int n_mappings = sorted_algnvector->at(i)->path.mapping().size();
    const pair<long int, long int> stop_coords =  make_pair(sorted_algnvector->at(i)->path.mapping()[n_mappings].position().node_id(),
                                  sorted_algnvector->at(i)->path.mapping()[n_mappings].position().offset());

        for (long unsigned int j=i+1; j!=sorted_algnvector->size(); ++j){
            const int candidate_n_mappings = sorted_algnvector->at(j)->path.mapping().size();
            const pair<long int, long int> candidate_start_coords = make_pair(sorted_algnvector->at(j)->path.mapping()[0].position().node_id(),
                                                    sorted_algnvector->at(j)->path.mapping()[0].position().offset());

            const pair<long int, long int> candidate_stop_coords = make_pair(sorted_algnvector->at(j)->path.mapping()[candidate_n_mappings].position().node_id(),
                                                    sorted_algnvector->at(j)->path.mapping()[candidate_n_mappings].position().offset());

            if (i<j && candidate_start_coords == start_coords && candidate_stop_coords == stop_coords) {
                is_dup->at(j)=true;
                                                                                                       }
                                                                         }
}

shared_ptr<vector<bool>> Dup_Remover::remove_duplicates_internal(shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector, \
                                             const int n_threads, const bool quiet){
    // Remove duplicates before processing GAM file
    int progress_step = 1000;
    static vector<AlignmentInfo*> sorted_unduplicated_algnvector;
    // Keep track of which read we mark as a duplicate
    shared_ptr<vector<bool>> is_dup = make_shared<vector<bool>>(sorted_algnvector->size(), false);
    bool is_paired;
    bool is_reverse;
    long unsigned int i;
    ProgressBar* progress = new ProgressBar(sorted_algnvector->size(), "Duplicate Removal");
    for (i=0; i < sorted_algnvector->size(); ++i){
        if (sorted_algnvector->size() > 1000 && i % progress_step == 0 && !quiet) {
            progress->Progressed(i);
                                                                                  }
        is_paired = sorted_algnvector->at(i)->is_paired;
        is_reverse = sorted_algnvector->at(i)->path.mapping()[0].position().is_reverse();
        if (!is_paired) {
            Dup_Remover::update_for_single(is_dup, is_reverse, sorted_algnvector, i, n_threads);
                        }

        else if (is_paired){
            Dup_Remover::update_for_paired(is_dup, is_reverse, sorted_algnvector, i, n_threads);
                           }

    for (long unsigned int i=0; i!=sorted_algnvector->size(); ++i) {
        if (is_dup->at(i) == false) {
            sorted_unduplicated_algnvector.emplace_back(sorted_algnvector->at(i));
                                    }
                                                                   }
                                                   }
//delete progress; progress = nullptr;
return is_dup;
                                                                                          }


void Dup_Remover::remove_duplicates(shared_ptr<Trailmix_struct> &dta, const char *gam_file){
  auto algninfo_gam = readGAM(dta);
  const auto algn_vg = readGAM(gam_file);
  shared_ptr<vector<bool>> is_dup = Dup_Remover::remove_duplicates_internal(algninfo_gam, 1, false);
  vector<vg::Alignment> to_write;
  for (size_t i=0; i<is_dup->size(); ++i) {
      if (!is_dup->at(i)){to_write.emplace_back(algn_vg.at(i));}
                                          }

  vg::io::write_buffered(cout,to_write,1);
                                                                                           }
