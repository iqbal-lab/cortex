/*
 needleman_wunsch.c
 project: NeedlemanWunsch
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 25-May-2011
 
 To test/compile, see nw_test.c
 see README
 
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

// Turn on debugging output by defining DEBUG
//#define DEBUG

#define arr_lookup(arr,width,i,j) arr[((j)*(width)) + (i)]
#define MAX_3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : ((y) >= (z) ? (y) : (z)))
#define GENERATE_HASH(a,b) (int)(((int)(a) << 8) | (int)(b))

#define MATCH 0
#define GAP_A 1
#define GAP_B 2
#define MATRIX_NAME(x) ((x) == 0 ? "MATCH" : (\
                        (x) == 1 ? "GAP_A" : (\
                        (x) == 2 ? "GAP_B" : (\
                        "Unknown matrix"))))

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower

// DEBUG
#ifdef DEBUG
#define INT_MIN -99
#else
#include <limits.h>
#endif

#include "needleman_wunsch.h"

NW_SCORING* custom_scoring(int num_chars, char* chars, int* scores,
                           int gap_open, int gap_extend,
                           char no_start_gap_penalty, char no_end_gap_penalty,
                           char use_match_mismatch,
                           int match, int mismatch,
                           char case_sensitive)
{
  // Create hash table
  NW_SCORE* hashtable = NULL;

  int i, j;
  char a, b;
  for(i = 0; i < num_chars; i++)
  {
    a = case_sensitive ? chars[i] : tolower(chars[i]);

    for(j = 0; j < num_chars; j++)
    {
      b = case_sensitive ? chars[j] : tolower(chars[j]);

      NW_SCORE* new_entry = (NW_SCORE*) malloc(sizeof(NW_SCORE));
      new_entry->id = GENERATE_HASH(a, b);
      new_entry->swap_score = arr_lookup(scores,num_chars,i,j);

      HASH_ADD_INT(hashtable, id, new_entry);
    }
  }

  NW_SCORING* scoring = (NW_SCORING*) malloc(sizeof(NW_SCORING));

  scoring->swap_table = hashtable;

  // Gap of length 1 has penalty (gap_open+gap_extend)
  // of length N: (gap_open + gap_extend*N)
  scoring->gap_open = gap_open;
  scoring->gap_extend = gap_extend;

  scoring->no_start_gap_penalty = no_start_gap_penalty;
  scoring->no_end_gap_penalty = no_end_gap_penalty;

  scoring->use_match_mismatch = use_match_mismatch;
  scoring->match = match;
  scoring->mismatch = mismatch;
  
  scoring->case_sensitive = case_sensitive;

  return scoring;
}

NW_SCORING* simple_scoring(int match, int mismatch, int gap_open, int gap_extend,
                           char no_start_gap_penalty, char no_end_gap_penalty,
                           char case_sensitive)
{
  NW_SCORING* scoring = (NW_SCORING*) malloc(sizeof(NW_SCORING));

  scoring->swap_table = NULL;

  // Gap of length 1 has penalty (gap_open+gap_extend)
  // of length N: (gap_open + gap_extend*N)
  scoring->gap_open = gap_open;
  scoring->gap_extend = gap_extend;

  scoring->no_start_gap_penalty = no_start_gap_penalty;
  scoring->no_end_gap_penalty = no_end_gap_penalty;

  scoring->use_match_mismatch = 1;
  scoring->match = match;
  scoring->mismatch = mismatch;

  scoring->case_sensitive = case_sensitive;

  return scoring;
}

int score_lookup(NW_SCORING* scoring, char a, char b)
{
  if(!scoring->case_sensitive)
  {
    a = tolower(a);
    b = tolower(b);
  }

  if(scoring->swap_table == NULL)
  {
    return a == b ? scoring->match : scoring->mismatch;
  }

  // Look up in table
  int hash_key = GENERATE_HASH(a,b);
  NW_SCORE* result;

  HASH_FIND_INT(scoring->swap_table, &hash_key, result);

  if(result == NULL)
  {
    if(scoring->use_match_mismatch)
    {
      return a == b ? scoring->match : scoring->mismatch;
    }
    else
    {
      // Error
      fprintf(stderr, "Error: Unknown character pair (%c,%c) and "
                      "match/mismatch have not been set\n", a, b);
      exit(EXIT_FAILURE);
    }
  }
  
  return result->swap_score;
}

/*
 * - both strings should be of the same length and end with \0 char
 * - characters should be one of: aAcCgGtT
 * - gaps should be '-'
 */
int score_alignment(char* alignment_a, char* alignment_b, NW_SCORING* scoring)
{
  int score = 0;
  char in_gap_a = 0;
  char in_gap_b = 0;
  
  int start = 0;
  int end = strlen(alignment_a);
  
  if(scoring->no_start_gap_penalty)
  {
    // Move start
    if(alignment_a[end-1] == '-')
    {
      while(alignment_a[end-1] == '-') end--;
    }
    else if(alignment_b[end-1] == '-')
    {
      while(alignment_b[end-1] == '-') end--;
    }
  }
  
  if(scoring->no_end_gap_penalty)
  {
    // Move end
    if(alignment_a[start] == '-')
    {
      while(alignment_a[start] == '-') start++;
    }
    else if(alignment_b[start] == '-')
    {
      while(alignment_b[start] == '-') start++;
    }
  }
  
  int i;
  for(i = start; i < end; i++)
  {
    if(alignment_a[i] == '-')
    {
      if(in_gap_a) {
        score += scoring->gap_extend;
      }
      else {
        in_gap_a = 1;
        score += scoring->gap_open + scoring->gap_extend;
      }
      in_gap_b = 0;
    }
    else if(alignment_b[i] == '-')
    {
      if(in_gap_b) {
        score += scoring->gap_extend;
      }
      else {
        in_gap_b = 1;
        score += scoring->gap_open + scoring->gap_extend;
      }
      in_gap_a = 0;
    }
    else
    {
      in_gap_a = 0;
      in_gap_b = 0;

      score += score_lookup(scoring, alignment_a[i], alignment_b[i]);
    }
  }
  
  return score;
}


/* Allocate memory for alignment results */

int nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b)
{
  int length_a = strlen(seq_a);
  int length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  int longest_alignment = length_a + length_b;
  
  // longest_alignment + 1 to allow for \0
  *alignment_a = (char*) malloc((longest_alignment+1) * sizeof(char));
  *alignment_b = (char*) malloc((longest_alignment+1) * sizeof(char));

  return longest_alignment;
}

// length is = length_a + length_b
int nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b)
{
  // longest_alignment + 1 to allow for \0
  *alignment_a = realloc(*alignment_a, (length+1) * sizeof(char));
  *alignment_b = realloc(*alignment_b, (length+1) * sizeof(char));

  return length;
}

// Find backtrack start when scoring->no_end_gap_penalty is 1
char find_end_max(int *score_arr, int length_a, int length_b,
                  int *curr_score, int *seq_i, int *seq_j)
{
  int i, j, temp;
  char updated = 0;

  for(i = 1; i <= length_a; i++)
  {
    temp = arr_lookup(score_arr, length_a+1, i, length_b);
    if(temp > *curr_score)
    {
      *curr_score = temp;
      *seq_i = i-1;
      *seq_j = length_b-1;
      updated = 1;
    }
  }

  for(j = 1; j <= length_b; j++)
  {
    temp = arr_lookup(score_arr, length_a+1, length_a, j);
    if(temp > *curr_score)
    {
      *curr_score = temp;
      *seq_i = length_a-1;
      *seq_j = j-1;
      updated = 1;
    }
  }

  return updated;
}

#ifdef DEBUG
void print_scoring(NW_SCORING* scoring)
{
  printf("scoring:\n");
  printf("  match: %i; mismatch: %i; (use_match_mismatch: %i)\n",
         scoring->match, scoring->mismatch, scoring->use_match_mismatch);

  printf("  gap_open: %i; gap_extend: %i;\n",
         scoring->gap_open, scoring->gap_extend);

  printf("  no_start_gap_penalty: %i; no_end_gap_penalty: %i;\n",
         scoring->no_start_gap_penalty, scoring->no_end_gap_penalty);

  printf("  swap_table: %i\n", (scoring->swap_table == NULL));
}

void print_matrices(int* match_score, int* gap_a_score, int* gap_b_score,
                    int length_a, int length_b)
{
  int score_width = length_a+1;
  int i, j;

  printf("match_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(match_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_a_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(gap_a_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_b_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(gap_b_score, score_width, i, j));
    }
    printf("\n");
  }
}
#endif

/**
 * Align with gap start and continue penalties
 */
int needleman_wunsch(char* seq_a, char* seq_b,
                     char* alignment_a, char* alignment_b,
                     NW_SCORING* scoring)
{
#ifdef DEBUG
  print_scoring(scoring);
#endif

  int length_a = strlen(seq_a);
  int length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  int longest_alignment = length_a + length_b;
  
  int score_width = length_a+1;
  int score_height = length_b+1;
  
  int arr_size = score_width * score_height;
  
  // 2d array (length_a x length_b)
  // addressing [a][b]

  // Score having just matched
  int* match_score = (int*) malloc(arr_size * sizeof(int));
  // score having just deleted from seq_a
  int* gap_a_score = (int*) malloc(arr_size * sizeof(int));
  // score having just inserted into seq_a
  int* gap_b_score = (int*) malloc(arr_size * sizeof(int));
  
  int i, j, index;
  
  // [0][0]
  match_score[0] = 0;
  gap_a_score[0] = 0;
  gap_b_score[0] = 0;
  
  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_score[i] = INT_MIN;
    
    // Think carefully about which way round these are
    gap_a_score[i] = INT_MIN;
    gap_b_score[i] = scoring->no_start_gap_penalty ? 0
                     : scoring->gap_open + i * scoring->gap_extend;
  }
  
  // work down first column -> [0][j]
  for(j = 1; j < score_height; j++)
  {
    index = j*score_width;
    match_score[index] = INT_MIN;
    
    // Think carefully about which way round these are
    gap_a_score[index] = scoring->no_start_gap_penalty ? 0
                         : scoring->gap_open + j * scoring->gap_extend;
    gap_b_score[index] = INT_MIN;
  }
  
  // Update Dynamic Programming arrays
  int seq_i, seq_j, sub_penalty;
  int old_index, new_index;

  for (i = 1; i < score_width; i++)
  {
    for (j = 1; j < score_height; j++)
    {
      // It's an indexing thing...
      seq_i = i - 1;
      seq_j = j - 1;
      
      sub_penalty = score_lookup(scoring, seq_a[seq_i], seq_b[seq_j]);
      
      // Update match_score[i][j] with position [i-1][j-1]
      new_index = j*score_width + i;
      old_index = (j-1)*score_width + (i-1);
      
      // substitution
      match_score[new_index] = MAX_3(match_score[old_index], // continue alignment
                                     gap_a_score[old_index], // close gap in seq_a
                                     gap_b_score[old_index]) + sub_penalty;
                                     // ^ close gap in seq_b
      
      // Update gap_a_score[i][j] from position [i][j-1]
      old_index = (j-1)*score_width + i;
      
      // Long arithmetic since some INTs are set to INT_MIN and penalty is -ve
      // (adding as ints would cause an integer overflow)
      gap_a_score[new_index]
        = MAX_3((long)match_score[old_index] + scoring->gap_extend + scoring->gap_open,
                (long)gap_a_score[old_index] + scoring->gap_extend,
                (long)gap_b_score[old_index] + scoring->gap_extend + scoring->gap_open);
    
      // Update gap_b_score[i][j] from position [i-1][j]
      old_index = j*score_width + (i-1);
      
      gap_b_score[new_index]
        = MAX_3((long)match_score[old_index] + scoring->gap_extend + scoring->gap_open,
                (long)gap_a_score[old_index] + scoring->gap_extend + scoring->gap_open,
                (long)gap_b_score[old_index] + scoring->gap_extend);
    }
  }
   
  //
  // Trace back now (score matrices all calculated)
  //
  
  // work backwards re-tracing optimal alignment, then shift sequences into place
  
  char curr_matrix;
  int curr_score;
  
  // Position of next alignment character in buffer (working backwards)
  int next_char = longest_alignment-1;
  
  // Previous scores on each matrix
  int prev_match_score, prev_gap_a_score, prev_gap_b_score;
  
  // penalties if coming from each of the prev matrices
  int prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty;
  
  if(scoring->no_end_gap_penalty)
  {
    curr_matrix = MATCH;
    curr_score = arr_lookup(match_score, score_width, score_width-1, 0);
    
    find_end_max(match_score, length_a, length_b,
                 &curr_score, &seq_i, &seq_j);
    
    if(find_end_max(gap_a_score, length_a, length_b,
                    &curr_score, &seq_i, &seq_j))
    {
      curr_matrix = GAP_A;
    }
    
    if(find_end_max(gap_b_score, length_a, length_b,
                    &curr_score, &seq_i, &seq_j))
    {
      curr_matrix = GAP_B;
    }
    
    #ifdef DEBUG
    printf("no_end_gap_penalty: (matrix: %s, curr_score: %i, seq_i: %i, seq_j: %i)\n",
           MATRIX_NAME(curr_matrix), curr_score, seq_i, seq_j);
    #endif
    
    // Fill in last gap
    int i;
    for(i = length_a - 1; i > seq_i; i--, next_char--)
    {
      alignment_a[next_char] = seq_a[i];
      alignment_b[next_char] = '-';
    }

    int j;
    for(j = length_b - 1; j > seq_j; j--, next_char--)
    {
      alignment_a[next_char] = '-';
      alignment_b[next_char] = seq_b[j];
    }
  }
  else
  {
    // Get max score (and therefore current matrix)
    curr_matrix = MATCH;
    curr_score = match_score[arr_size-1];
    
    if(gap_a_score[arr_size-1] > curr_score)
    {
      curr_matrix = GAP_A;
      curr_score = gap_a_score[arr_size-1];
    }
    
    if(gap_b_score[arr_size-1] > curr_score)
    {
      curr_matrix = GAP_B;
      curr_score = gap_b_score[arr_size-1];
    }

    seq_i = length_a - 1;
    seq_j = length_b - 1;
  }

#ifdef DEBUG
  print_matrices(match_score, gap_a_score, gap_b_score, length_a, length_b);
#endif

  // Hold this value to return later
  long max_alignment_score = curr_score;
  
  // note: longest_alignment = strlen(seq_a) + strlen(seq_b)
  // seq_i is the index of the next char of seq_a to be added (working bckwrds)
  // seq_j is the index of the next char of seq_b to be added (working bckwrds)
  for(; seq_i >= 0 && seq_j >= 0; next_char--)
  {
    #ifdef DEBUG
    printf("matrix: %s (%i,%i) score: %i\n",
           MATRIX_NAME(curr_matrix), seq_i, seq_j, curr_score);
    #endif

    switch (curr_matrix)
    {
      case MATCH:
        alignment_a[next_char] = seq_a[seq_i];
        alignment_b[next_char] = seq_b[seq_j];
        
        // Match
        prev_match_penalty = score_lookup(scoring, seq_a[seq_i], seq_b[seq_j]);
        prev_gap_a_penalty = prev_match_penalty; // match
        prev_gap_b_penalty = prev_match_penalty; // match
        
        // Moving back on i and j
        seq_i--;
        seq_j--;
        break;

      case GAP_A:
        alignment_a[next_char] = '-';
        alignment_b[next_char] = seq_b[seq_j];
        
        prev_match_penalty = scoring->gap_extend + scoring->gap_open;
        prev_gap_a_penalty = scoring->gap_extend;
        prev_gap_b_penalty = scoring->gap_extend + scoring->gap_open;
        
        // Moving back on j
        seq_j--;
        break;

      case GAP_B:
        alignment_a[next_char] = seq_a[seq_i];
        alignment_b[next_char] = '-';
        
        prev_match_penalty = scoring->gap_extend + scoring->gap_open;
        prev_gap_a_penalty = scoring->gap_extend + scoring->gap_open;
        prev_gap_b_penalty = scoring->gap_extend;
        
        // Moving back on i
        seq_i--;
        break;

      default:
        fprintf(stderr, "Err: invalid matrix number\n");
        exit(EXIT_FAILURE);
    }
    
    // Current score matrix position is [seq_i+1][seq_j+1]
    int score_i = seq_i + 1;
    int score_j = seq_j + 1;

    // [score_i][score_j] is the next position in the score matrices
    
    prev_match_score = arr_lookup(match_score, score_width, score_i, score_j);
    prev_gap_a_score = arr_lookup(gap_a_score, score_width, score_i, score_j);
    prev_gap_b_score = arr_lookup(gap_b_score, score_width, score_i, score_j);
    
#ifdef DEBUG
    printf("  prev_match_score: %i; prev_gap_a_score: %i; prev_gap_b_score: %i\n",
           prev_match_score, prev_gap_a_score, prev_gap_b_score);
#endif
    
    // Now figure out which matrix we came from
    if((long)prev_match_score + prev_match_penalty == curr_score)
    {
      // Both sequences have a base
      curr_matrix = MATCH;
      curr_score = prev_match_score;
    }
    else if((long)prev_gap_a_score + prev_gap_a_penalty == curr_score)
    {
      // Gap in seq_a
      curr_matrix = GAP_A;
      curr_score = prev_gap_a_score;
    }
    else if((long)prev_gap_b_score + prev_gap_b_penalty == curr_score)
    {
      // Gap in seq_b
      curr_matrix = GAP_B;
      curr_score = prev_gap_b_score;
    }
    else
    {
      fprintf(stderr, "Fail\n");
      exit(EXIT_FAILURE);
    }
  }
  
  // Free memory
  free(match_score);
  free(gap_a_score);
  free(gap_b_score);
  
  // Gap in A
  while(seq_j >= 0)
  {
    alignment_a[next_char] = '-';
    alignment_b[next_char] = seq_b[seq_j];
    next_char--;
    seq_j--;
  }

  // Gap in B
  while(seq_i >= 0)
  {
    alignment_a[next_char] = seq_a[seq_i];
    alignment_b[next_char] = '-';
    next_char--;
    seq_i--;
  }

  // Shift alignment strings back into 0th position in char arrays
  int first_char = next_char+1;
  int alignment_len = longest_alignment - first_char;

#ifdef DEBUG
  printf("END: seq_i: %i; seq_j: %i; first_char %i; longest_alignment %i; length %i;\n",
         seq_i, seq_j, first_char, longest_alignment, alignment_len);
#endif

  int pos;
  for(pos = 0; pos < alignment_len; pos++)
  {
    alignment_a[pos] = alignment_a[pos+first_char];
    alignment_b[pos] = alignment_b[pos+first_char];
  }
  
  alignment_a[pos] = '\0';
  alignment_b[pos] = '\0';

  return max_alignment_score;
}

//
// Some scoring systems
//

// Scoring for protein comparisons of length <35bp
NW_SCORING* scoring_system_PAM30()
{
  int pam30[529] = {6,-7,-4,-3,-6,-4,-2,-2,-7,-5,-6,-7,-5,-8,-2,0,-1,-13,-8,-2,-3,
  -3,-3,-7,8,-6,-10,-8,-2,-9,-9,-2,-5,-8,0,-4,-9,-4,-3,-6,-2,-10,-8,-7,-4,-6,-4,
  -6,8,2,-11,-3,-2,-3,0,-5,-7,-1,-9,-9,-6,0,-2,-8,-4,-8,6,-3,-3,-3,-10,2,8,-14,
  -2,2,-3,-4,-7,-12,-4,-11,-15,-8,-4,-5,-15,-11,-8,6,1,-5,-6,-8,-11,-14,10,-14,
  -14,-9,-7,-6,-15,-14,-13,-13,-8,-3,-8,-15,-4,-6,-12,-14,-9,-4,-2,-3,-2,-14,8,
  1,-7,1,-8,-5,-3,-4,-13,-3,-5,-5,-13,-12,-7,-3,6,-5,-2,-9,-2,2,-14,1,8,-4,-5,
  -5,-9,-4,-7,-14,-5,-4,-6,-17,-8,-6,1,6,-5,-2,-9,-3,-3,-9,-7,-4,6,-9,-11,-10,
  -7,-8,-9,-6,-2,-6,-15,-14,-5,-3,-5,-5,-7,-2,0,-4,-7,1,-5,-9,9,-9,-6,-6,-10,-6,
  -4,-6,-7,-7,-3,-6,-1,-1,-5,-5,-5,-5,-7,-6,-8,-5,-11,-9,8,-1,-6,-1,-2,-8,-7,-2,
  -14,-6,2,-6,-6,-5,-6,-8,-7,-12,-15,-5,-9,-10,-6,-1,7,-8,1,-3,-7,-8,-7,-6,-7,
  -2,-9,-7,-6,-7,0,-1,-4,-14,-3,-4,-7,-6,-6,-8,7,-2,-14,-6,-4,-3,-12,-9,-9,-2,
  -4,-5,-5,-4,-9,-11,-13,-4,-7,-8,-10,-1,1,-2,11,-4,-8,-5,-4,-13,-11,-1,-10,-5,
  -5,-8,-9,-9,-15,-13,-13,-14,-9,-6,-2,-3,-14,-4,9,-10,-6,-9,-4,2,-8,-10,-13,-8,
  -2,-4,-6,-8,-8,-3,-5,-6,-4,-8,-7,-6,-8,-10,8,-2,-4,-14,-13,-6,-7,-4,-5,0,-3,
  0,-4,-3,-5,-4,-2,-6,-7,-8,-4,-5,-6,-2,6,0,-5,-7,-6,-1,-5,-3,-1,-6,-2,-5,-8,-5,
  -6,-6,-7,-2,-7,-3,-4,-9,-4,0,7,-13,-6,-3,-3,-6,-4,-13,-2,-8,-15,-15,-13,-17,
  -15,-7,-14,-6,-12,-13,-4,-14,-5,-13,13,-5,-15,-10,-14,-11,-8,-10,-4,-11,-4,
  -12,-8,-14,-3,-6,-7,-9,-11,2,-13,-7,-6,-5,10,-7,-6,-9,-7,-2,-8,-8,-8,-6,-7,
  -6,-5,-6,2,-2,-9,-1,-8,-6,-6,-3,-15,-7,7,-8,-6,-5,-3,-7,6,6,-12,-3,1,-3,-1,
  -6,-9,-2,-10,-10,-7,-1,-3,-10,-6,-8,6,0,-5,-3,-4,-3,1,-14,6,6,-5,-1,-6,-7,-4,
  -5,-13,-4,-5,-6,-14,-9,-6,0,6,-5,-3,-6,-3,-5,-9,-5,-5,-5,-5,-5,-6,-5,-5,-8,-5,
  -3,-4,-11,-7,-5,-5,-5,-5};

  char *bases = "ARNDCQEGHILKMFPSTWYVBZX";

  // Gap open -9, gap extend -1
  return custom_scoring(23, bases, pam30, -9, -1, 0,0,0,0,0,0);
}

// Scoring for protein comparisons of length 35-50
NW_SCORING* scoring_system_PAM70()
{
  int pam70[529] = {4,-2,-2,-2,-1,-1,-1,0,-2,-2,-2,-1,-1,-2,-1,1,0,-3,-2,0,-2,
  -1,-1,-2,6,-1,-2,-4,1,0,-3,0,-3,-3,2,-2,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-2,-1,6,
  1,-3,0,0,-1,0,-4,-4,0,-2,-3,-2,0,0,-4,-2,-3,3,0,-1,-2,-2,1,6,-4,-1,1,-2,-1,-4,
  -4,-1,-3,-4,-2,0,-1,-5,-4,-4,4,1,-2,-1,-4,-3,-4,9,-3,-4,-3,-4,-1,-2,-4,-2,-2,
  -3,-1,-1,-3,-3,-1,-4,-4,-2,-1,1,0,-1,-3,6,2,-2,1,-3,-2,1,0,-3,-2,0,-1,-2,-2,
  -2,0,3,-1,-1,0,0,1,-4,2,5,-2,0,-4,-3,1,-2,-4,-1,0,-1,-4,-3,-3,1,4,-1,0,-3,-1,
  -2,-3,-2,-2,6,-2,-4,-4,-2,-3,-4,-3,-1,-2,-3,-4,-4,-1,-2,-2,-2,0,0,-1,-4,1,0,
  -2,8,-4,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,-1,0,-1,-2,-3,-4,-4,-1,-3,-4,-4,-4,4,2,
  -3,1,0,-3,-3,-1,-3,-1,3,-4,-3,-1,-2,-3,-4,-4,-2,-2,-3,-4,-3,2,4,-3,2,0,-3,-3,
  -2,-2,-1,1,-4,-3,-1,-1,2,0,-1,-4,1,1,-2,-1,-3,-3,5,-2,-3,-1,0,-1,-3,-2,-3,-1,
  1,-1,-1,-2,-2,-3,-2,0,-2,-3,-2,1,2,-2,6,0,-3,-2,-1,-2,-1,1,-3,-2,-1,-2,-3,-3,
  -4,-2,-3,-4,-4,-1,0,0,-3,0,6,-4,-3,-2,1,3,-1,-4,-4,-2,-1,-2,-2,-2,-3,-2,-1,-3,
  -2,-3,-3,-1,-3,-4,8,-1,-1,-4,-3,-3,-2,-1,-2,1,-1,0,0,-1,0,0,-1,-1,-3,-3,0,-2,
  -3,-1,4,1,-3,-2,-2,0,0,-1,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-2,-1,-1,-2,-1,1,5,-3,
  -2,0,-1,-1,-1,-3,-3,-4,-5,-3,-2,-4,-3,-2,-3,-2,-3,-2,1,-4,-3,-3,11,2,-3,-4,-3,
  -3,-2,-2,-2,-4,-3,-2,-3,-4,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-2,-3,-2,-2,0,-3,-3,
  -4,-1,-2,-3,-4,-3,3,1,-3,1,-1,-3,-2,0,-3,-2,4,-3,-3,-1,-2,-1,3,4,-4,0,1,-1,-1,
  -4,-4,-1,-3,-4,-2,0,-1,-4,-3,-3,4,0,-1,-1,0,0,1,-4,3,4,-2,0,-3,-3,1,-2,-4,-1,
  0,-1,-3,-2,-3,0,4,-1,-1,-1,-1,-2,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,-1,-3,-2,
  -1,-1,-1,-1};
  
  char *bases = "ARNDCQEGHILKMFPSTWYVBZX";

  // Gap open -10, gap extend -1
  return custom_scoring(23, bases, pam70, -10, -1, 0,0,0,0,0,0);
}

// Scoring for protein comparisons of length 50-85
NW_SCORING* scoring_system_BLOSUM80()
{
  int blosum80[529] = {-2,0,-3,-3,-3,-1,-2,-4,-1,2,0,-5,-4,-1,-3,-2,-1,-3,9,-1,
  -3,-6,1,-1,-4,0,-5,-4,3,-3,-5,-3,-2,-2,-5,-4,-4,-2,0,-2,-3,-1,9,2,-5,0,-1,-1,
  1,-6,-6,0,-4,-6,-4,1,0,-7,-4,-5,5,-1,-2,-3,-3,2,10,-7,-1,2,-3,-2,-7,-7,-2,-6,
  -6,-3,-1,-2,-8,-6,-6,6,1,-3,-1,-6,-5,-7,13,-5,-7,-6,-7,-2,-3,-6,-3,-4,-6,-2,
  -2,-5,-5,-2,-6,-7,-4,-2,1,0,-1,-5,9,3,-4,1,-5,-4,2,-1,-5,-3,-1,-1,-4,-3,-4,
  -1,5,-2,-2,-1,-1,2,-7,3,8,-4,0,-6,-6,1,-4,-6,-2,-1,-2,-6,-5,-4,1,6,-2,0,-4,
  -1,-3,-6,-4,-4,9,-4,-7,-7,-3,-5,-6,-5,-1,-3,-6,-6,-6,-2,-4,-3,-3,0,1,-2,-7,
  1,0,-4,12,-6,-5,-1,-4,-2,-4,-2,-3,-4,3,-5,-1,0,-2,-3,-5,-6,-7,-2,-5,-6,-7,-6,
  7,2,-5,2,-1,-5,-4,-2,-5,-3,4,-6,-6,-2,-3,-4,-6,-7,-3,-4,-6,-7,-5,2,6,-4,3,0,
  -5,-4,-3,-4,-2,1,-7,-5,-2,-1,3,0,-2,-6,2,1,-3,-1,-5,-4,8,-3,-5,-2,-1,-1,-6,-4,
  -4,-1,1,-2,-2,-3,-4,-6,-3,-1,-4,-5,-4,2,3,-3,9,0,-4,-3,-1,-3,-3,1,-5,-3,-2,-4,
  -5,-6,-6,-4,-5,-6,-6,-2,-1,0,-5,0,10,-6,-4,-4,0,4,-2,-6,-6,-3,-1,-3,-4,-3,-6,
  -3,-2,-5,-4,-5,-5,-2,-4,-6,12,-2,-3,-7,-6,-4,-4,-2,-3,2,-2,1,-1,-2,-1,-1,-1,
  -2,-4,-4,-1,-3,-4,-2,7,2,-6,-3,-3,0,-1,-1,0,-2,0,-2,-2,-1,-2,-3,-3,-2,-3,-1,
  -1,-4,-3,2,8,-5,-3,0,-1,-2,-1,-5,-5,-7,-8,-5,-4,-6,-6,-4,-5,-4,-6,-3,0,-7,-6,
  -5,16,3,-5,-8,-5,-5,-4,-4,-4,-6,-5,-3,-5,-6,3,-3,-2,-4,-3,4,-6,-3,-3,3,11,-3,
  -5,-4,-3,-1,-4,-5,-6,-2,-4,-4,-6,-5,4,1,-4,1,-2,-4,-3,0,-5,-3,7,-6,-4,-2,-3,
  -2,5,6,-6,-1,1,-2,-1,-6,-7,-1,-5,-6,-4,0,-1,-8,-5,-6,6,0,-3,-2,0,-1,1,-7,5,6,
  -4,0,-6,-5,1,-3,-6,-2,-1,-2,-5,-4,-4,0,6,-1,-1,-2,-2,-3,-4,-2,-2,-3,-2,-2,-2,
  -2,-2,-3,-3,-1,-1,-5,-3,-2,-3,-1,-2};
  
  char *bases = "ARNDCQEGHILKMFPSTWYVBZX";

  // Gap open -10, gap extend -1
  return custom_scoring(23, bases, blosum80, -10, -1, 0,0,0,0,0,0);
}

// Scoring for protein comparisons of length >85
NW_SCORING* scoring_system_BLOSUM62()
{
  int blosum62[529] = {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,
  -1,0,-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-2,0,6,1,
  -3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-2,-2,1,6,-3,0,2,-1,-1,-3,-4,
  -1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,
  -1,-1,-2,-2,-1,-3,-3,-2,-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,
  3,-1,-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,0,-2,0,-1,-3,
  -2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-2,0,1,-1,-3,0,0,-2,8,-3,
  -3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,
  -2,-1,-3,-1,3,-3,-3,-1,-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,
  -4,-3,-1,-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-1,-1,
  -2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-2,-3,-3,-3,-2,-3,-3,
  -3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,
  -2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,
  -2,-2,0,0,0,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-3,
  -3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-2,-2,-2,-3,
  -2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,0,-3,-3,-3,-1,-2,-2,-3,
  -3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,
  -2,0,-1,-4,-3,-3,4,1,-1,-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,
  4,-1,0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1};
  
  char *bases = "ARNDCQEGHILKMFPSTWYVBZX";

  // Gap open -10, gap extend -1
  return custom_scoring(23, bases, blosum62, -10, -1, 0,0,0,0,0,0);
}

// Scoring system for predicting DNA hybridization
// "Optimization of the BLASTN substitution matrix for prediction of
//   non-specific DNA microarray hybridization" (2009)
// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2831327/
NW_SCORING* scoring_system_DNA_hybridization()
{
  int sub_matrix[64] = { 2, 2,-4,-4,-4,-4,-4,-4,
                         2, 2,-4,-4,-4,-4,-4,-4,
                        -4,-4, 5, 5,-4,-4,-4,-4,
                        -4,-4, 5, 5,-4,-4,-4,-4,
                        -4,-4,-4,-4, 5, 5,-4,-4,
                        -4,-4,-4,-4, 5, 5,-4,-4,
                        -4,-4,-4,-4,-4,-4, 2, 2,
                        -4,-4,-4,-4,-4,-4, 2, 2};

  char *bases = "AaCcGgTt";

  // Gap open -10, gap extend -10
  return custom_scoring(8, bases, sub_matrix, -10, -10, 0,0,0,0,0,0);
}

