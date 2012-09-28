/* This file is part of Repeat-finder

    Repeat-finder is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Repeat-finder is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Repeat-finder; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _MATCH_H_
#define _MATCH_H_

#include "fasta.h"

#define MAX_MATCH_LENGTH 1000
#define SECONDARY_MOTIF_INCREMENT 16
#define SUB_MATCH_INCREMENT 8

typedef struct {
	int index;
	int start;
	int end;
	int size;
	
	int sequence_index;
	
	int p_index;
	
	int substitutions;
	int deletions;
	int insertions;
	int errors; //Total of above three
	
	int truncated; // Indicates whether stats are reliable
	
	double repeat_count;
	double purity;

	int motif_index;
	int secondary_motifs[SECONDARY_MOTIF_INCREMENT];
	int num_secondary_motifs;
	int max_secondary_motifs;
	
	Boolean in_compound;
//	Boolean is_degenerative;

	struct Match** sub_matches;
	int num_sub_matches;
	int max_sub_matches;
	
	int total_gap_sizes;
	
	//motif position for *next* match in either direction.
	int l_motif_pos;
	int r_motif_pos;

	int hotspot_index;
	int hotspot_interval;
	int content[5];
	
	char sequence[MAX_MATCH_LENGTH];
	int left_sequence;
	int right_sequence;

} Match, *MatchPtr;

#define MATCH_INCREMENT 10240
// is mismatch symbol?
#define IS_MISMATCH(x) (x == '+' || x == '-' || x == '*')

MatchPtr findMatch(SeqStructPtr seqP, int pos, int motif_index);
void findMatches(SeqStructPtr seqP, char** repeats, int numRepeats);
MatchPtr copyMatch(MatchPtr match);
// Remove a match from the array
void remove_match(int index);
// Delete a match structure properly
void deleteMatch(MatchPtr m);
void remove_partial_ends(Match* m, int direction);
void extendMatches(SeqStructPtr seqP, MatchPtr* matches);
void shortenMatch(MatchPtr m, int amount, int direction);
void trim_match(MatchPtr m);

#define LEFT 0
#define RIGHT 1

int recursiveMatch(SeqStructPtr seqP, MatchPtr match, int direction,
int* error_distance);
void add_secondary_motif(MatchPtr m, int m_index);
Boolean isDegenerate(MatchPtr m, MatchPtr n);
void sort_matches();

void find_compound_repeats(FileStructPtr fsp, Boolean compound_repeat, Boolean degenerative, int compound_gap_size);

int remove_small_matches(int size);
void dump_match(int i);
void dump_matches();

Boolean are_end_bases_pure(MatchPtr m, int direction, int length);

#endif //_MATCH_H_
