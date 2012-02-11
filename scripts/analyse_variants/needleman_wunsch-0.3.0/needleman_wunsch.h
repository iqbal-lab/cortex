/*
 needleman_wunsch.h
 project: NeedlemanWunsch
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 25-May-2011
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef NEEDLEMAN_WUNSCH_HEADER_SEEN
#define NEEDLEMAN_WUNSCH_HEADER_SEEN

#include "uthash.h"

typedef struct NW_SCORE NW_SCORE;
typedef struct NW_SCORING NW_SCORING;

struct NW_SCORE
{
  int id; // hash key
  int swap_score;
  UT_hash_handle hh; // makes this structure hashable
};

struct NW_SCORING
{
  int gap_open, gap_extend;

  // Turn these on to turn off penalties for gaps at the start/end of alignment
  char no_start_gap_penalty;
  char no_end_gap_penalty;

  // If swap_table != NULL, but char->char pair swap is not in the hashtable,
  // should we use match/mismatch values?
  char use_match_mismatch;
  int match, mismatch;

  char case_sensitive;

  NW_SCORE* swap_table;
};

/* Allocate memory for result */

// alloc memory for result (returns length of seq_a + seq_b)
int nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b);

// length is = length_a + length_b
int nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b);

/* Alignment */

int needleman_wunsch(char* seq_a, char* seq_b,
                     char* alignment_a, char* alignment_b,
                     NW_SCORING* scoring);

int score_alignment(char* alignment_a, char* alignment_b, NW_SCORING* scoring);

/* Scoring */

// Scoring set up
NW_SCORING* custom_scoring(int num_chars, char* chars, int* scores,
                           int gap_open, int gap_extend,
                           char no_start_gap_penalty, char no_end_gap_penalty,
                           char use_match_mismatch,
                           int match, int mismatch,
                           char case_sensitive);

NW_SCORING* simple_scoring(int match, int mismatch, int gap_open, int gap_extend,
                           char no_start_gap_penalty, char no_end_gap_penalty,
                           char case_sensitive);

// Some scoring systems
NW_SCORING* scoring_system_PAM30();
NW_SCORING* scoring_system_PAM70();
NW_SCORING* scoring_system_BLOSUM80();
NW_SCORING* scoring_system_BLOSUM62();
NW_SCORING* scoring_system_DNA_hybridization();

#endif
