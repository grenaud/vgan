#pragma once
#include "Trailmix_struct.h"
#include "AlignmentInfo.h"

using namespace std;

class Dup_Remover{

public:
/**
 * @brief Print usage
 *
 */

string_view usage();

/**
 * @brief Remove PCR duplicates from a GAM file.
 *
 * Remove PCR duplicates from a GAM file.
 *
 */

void remove_duplicates(shared_ptr<Trailmix_struct>& dta, const char *gam_file);

shared_ptr<vector<bool>> remove_duplicates_internal(shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector,\
                                     const int n_threads, const bool quiet);

/**
 * @brief Helper function to update vector of booleans when encountering a paired read.
 *
 * Helper function to update vector of booleans when encountering a paired forward-end read.
 *
 */

void update_for_paired(shared_ptr<vector<bool>> &is_dup, const bool is_reverse, shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector,
const long unsigned int i, const int n_threads);

/**
 * @brief  Helper function to update vector of booleans when encountering a single-end read.
 *
 * Helper function to update vector of booleans when encountering a single-end read.
 *
 */

void update_for_single(shared_ptr<vector<bool>> &is_dup, const bool is_reverse, \
                       shared_ptr<vector<AlignmentInfo*>>& sorted_algnvector, long unsigned int i, const int n_threads);

Dup_Remover();
Dup_Remover(const Dup_Remover & other);
~Dup_Remover();
Dup_Remover & operator= (const Dup_Remover & other);

private:
                 };
