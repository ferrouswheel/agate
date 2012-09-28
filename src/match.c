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

#include "match.h"
#include "motif.h"
#include "common.h"

extern char* database_motifs[MAX_MOTIFS];

int num_matches =  0;
int max_matches = MATCH_INCREMENT;
MatchPtr *matches;
int min_error_distance = 5;
int partial_ends = False;

int check_error_distance(MatchPtr m, int direction);

MatchPtr copyMatch(MatchPtr match)
{
	MatchPtr new_match;
	new_match = malloc(sizeof(Match));
	memcpy(new_match, match, sizeof(Match));
	return new_match;
}

void deleteMatch(MatchPtr m)
{
	if (m) {
		if (m->sub_matches)	free(m->sub_matches);
		//if (m->secondary_motifs) free(m->secondary_motifs);
		free(m);
	}
}

// finds pure matches of at least two repeats (at least three for mononucleotides).
MatchPtr findMatch(SeqStructPtr seqP, int pos, int motif_index) {
	int length;
	char* motif;
	
	MatchPtr new_match;
	int i;
	
	motif = database_motifs[motif_index];
	length = strlen(motif);

	new_match = malloc(sizeof(Match));
	memset(new_match, 0, sizeof(Match));
	new_match->motif_index = motif_index;
	new_match->hotspot_index = -1;
	new_match->sequence_index = seqP->index;
	
	i = pos;
	new_match->start = pos;	

	// At least two pure copies must be found...
	while (baseNCompare(seqP, i, motif, length) == 0)
		i += length;
		
	new_match->end = i - 1;
	new_match->size = i - new_match->start;
	new_match->repeat_count = (double) new_match->size / length;
	new_match->purity = 1.0;
	new_match->left_sequence = MAX_MATCH_LENGTH / 2;

	// ... and three pure copies if motif is only 1 base long
	if (length == 1 && new_match->repeat_count < 3) new_match->size = 0;
	else if (new_match->repeat_count < 2) new_match->size = 0;
	
	if (new_match->size == 0) {
		deleteMatch(new_match);
		return NULL;
	} else {
		baseNCopy(seqP, new_match->sequence + new_match->left_sequence, 
			new_match->start, new_match->size);
		new_match->right_sequence = new_match->left_sequence + new_match->size - 1;
		
		new_match->l_motif_pos = length - 1;
		new_match->r_motif_pos = 0;
		
		return new_match;
	}
}

// Macros shift motif position circularly left or right.
#define MOTIF_LEFT(x,length) (x) -= 1; if ((x) < 0) (x) = length - 1;
#define MOTIF_RIGHT(x,length) (x) += 1; if ((x) >= length) (x) = 0;

int recursiveMatch(SeqStructPtr seqP, MatchPtr match, int direction,
int* error_distance)
{
	int length;
	char* motif;
	int *motif_pos;
	int i;
	int start_i;
	// Initial position (start for left, end for right)
	// This is compared at end to check
	// a full match has been accomplished... otherwise
	// 0 is returned.
	int initial_position;

	motif = database_motifs[match->motif_index];
	length = strlen(motif);
		
	if (direction == LEFT) {
		i = match->start - 1;
		motif_pos = &(match->l_motif_pos);
		initial_position = match->start;
	} else {
		i = match->end+1;
		motif_pos = &(match->r_motif_pos);
		initial_position = match->end;

	}
	
	start_i = i;

#define ADD_LEFT(x,y) x->left_sequence--; x->sequence[x->left_sequence] = y;
#define ADD_RIGHT(x,y) x->right_sequence++; x->sequence[x->right_sequence] = y;
	
	while (i >= 0 && i < seqP->seqLen)
	{
		char base;
		base = getBaseAt(seqP, i);
		
		if (base == motif[*motif_pos]) {
			switch(direction){
			case LEFT:
				match->start = i;
				
				ADD_LEFT(match, base)
				
				i--;
				MOTIF_LEFT(*motif_pos, length)
				break;
			case RIGHT:
				match->end = i;

				ADD_RIGHT(match, base)							
				
				i++;
				MOTIF_RIGHT(*motif_pos, length)
				break;
			}
			match->size = (match->end - match->start) + 1;
			(*error_distance)++;
				
		} else {
			if (*error_distance < (min_error_distance-1)) {
				*error_distance = 0;
				//remove_partial_ends(match, motif_pos, direction);
				return abs(start_i - i);	
			}
			else {
				MatchPtr mismatches[3];
				int new_error_distances[3];
								
				int j;
				
				int distance = 0;
				
				int max_size = 0;
				int max_size_index = -1;
				
				
				
				*error_distance = 0;
				
				for (j=0; j < 3; j++) {
					mismatches[j] = copyMatch(match);
					new_error_distances[j] = *error_distance;
				}
			
				switch (direction)
				{
				case LEFT:
				
				    // Substitution moves both i and motifpos
					mismatches[0]->substitutions++;
				    mismatches[0]->start--;
				    MOTIF_LEFT(mismatches[0]->l_motif_pos, length)
				    if (mismatches[0]->start >= 0) {
						ADD_LEFT(mismatches[0], '*')

						distance = recursiveMatch(seqP, mismatches[0], direction,
						&new_error_distances[0]);

						if (distance == 0) {
							MOTIF_RIGHT(mismatches[0]->l_motif_pos, length)
							mismatches[0]->start++;
							mismatches[0]->left_sequence++;
							mismatches[0]->substitutions--;
						}
						
				    }

					//Insertion moves i
					mismatches[1]->insertions++;
					mismatches[1]->start--;
					if (mismatches[1]->start >= 0) {
						ADD_LEFT(mismatches[1], '+')
						
						distance = recursiveMatch(seqP, mismatches[1], direction,
						&new_error_distances[1]);
	
						if (distance == 0) {
							mismatches[1]->start++;
							mismatches[1]->left_sequence++;
							mismatches[1]->insertions--;
						}

					}
			    
			    	// Deletion moves only motifpos
					mismatches[2]->deletions++;
				    MOTIF_LEFT(mismatches[2]->l_motif_pos, length)
					if (mismatches[2]->start >= 0) {
						ADD_LEFT(mismatches[2], '-')
						
						distance = recursiveMatch(seqP, mismatches[2], direction,
						&new_error_distances[2]);
						
						if (distance == 0) {
							MOTIF_RIGHT(mismatches[2]->l_motif_pos, length)
							mismatches[2]->left_sequence++;
							mismatches[2]->deletions--;
						}
						
				    }				    	
			    			
					break;
				case RIGHT:
				    // Substitution moves both i and motifpos
					mismatches[0]->substitutions++;
					mismatches[0]->end++;
				    MOTIF_RIGHT(mismatches[0]->r_motif_pos, length)
				    if (mismatches[0]->end < seqP->seqLen) {
						ADD_RIGHT(mismatches[0], '*')
						
						distance = recursiveMatch(seqP, mismatches[0], direction,
						&new_error_distances[0]);
						
						if (distance == 0) {
							MOTIF_LEFT(mismatches[0]->r_motif_pos, length)
							mismatches[0]->end--;
							mismatches[0]->right_sequence--;
							mismatches[0]->substitutions--;
						}
				    }				    	

					//Insertion moves i
					mismatches[1]->insertions++;
					mismatches[1]->end++;
					if (mismatches[1]->end < seqP->seqLen) {
						ADD_RIGHT(mismatches[1], '+')
																		
						distance = recursiveMatch(seqP, mismatches[1], direction,
						&new_error_distances[1]);
						if (distance == 0) {
							mismatches[1]->end--;
							mismatches[1]->right_sequence--;							
							mismatches[1]->insertions--;
						}
					}
			    
			    	// Deletion moves only motifpos
					mismatches[2]->deletions++;			    	
					MOTIF_RIGHT(mismatches[2]->r_motif_pos, length)
					if (mismatches[2]->end < seqP->seqLen) {
						ADD_RIGHT(mismatches[2], '-')

						distance = recursiveMatch(seqP, mismatches[2], direction,
						&new_error_distances[2]);
						
						if (distance == 0) {
							MOTIF_LEFT(mismatches[2]->r_motif_pos, length)
							mismatches[2]->right_sequence--;
							mismatches[2]->deletions--;
						}
				    }    	
			    			
					break;
				}
				
				// Check which alternative got the furtherst
				if (strlen(motif) == 1)
				{
					for (j=0; j < 3; j++)
					{
						if (mismatches[j]->size > max_size) {
							max_size = mismatches[j]->size;
							max_size_index = j;
						}
					}
				}
				else
				{
					if (!partial_ends)
					{
						for (j=0; j < 3; j++)
						{
							remove_partial_ends(mismatches[j], direction);
						}
					}
					for (j=2; j >= 0; j--)
					{
						switch (direction)
						{
						case LEFT:
							if (mismatches[j]->start <= initial_position
								&& mismatches[j]->size > max_size) {
								max_size = mismatches[j]->size;
								max_size_index = j;
							}
							break;
						case RIGHT:
							if (mismatches[j]->end >= initial_position
								&& mismatches[j]->size > max_size) {
								max_size = mismatches[j]->size;
								max_size_index = j;
							}
							break;					
						}
					}	
				}
				
				// If no max size found it means that partial
				// match trimming has prevented mismatches from
				// getting longer than the initial match.
				if (max_size_index == -1) {
					for (j=0; j < 3; j++) free(mismatches[j]);
					return abs(start_i - i);
				} else {
					memcpy(match, mismatches[max_size_index], sizeof(Match));
					for (j=0; j < 3; j++) free(mismatches[j]);
				}
				
			}
			
			break;

		}
	}
	return abs(start_i - i);
}

void remove_partial_ends(Match* m, int direction)
{
	int length;
	char* motif;
	
	motif = database_motifs[m->motif_index];
	length = strlen(motif);

	switch (direction){
	case LEFT:
		while (m->l_motif_pos < (length - 1) || IS_MISMATCH(m->sequence[m->left_sequence]) )
		{
			shortenMatch(m, 1, direction);
		}
			
		break;
	case RIGHT:
		while (m->r_motif_pos > 0  || IS_MISMATCH(m->sequence[m->right_sequence]) )
		{
			shortenMatch(m, 1, direction);
		}
		break;
	}
	
}

void extendMatches(SeqStructPtr seqP, MatchPtr* matches)
{
	int length;
	char* motif;
	char status[1000];
	
	int j;

	int orig_size;
	int error_distance;
	
	MatchPtr match;
	
	memset(status, 0, 1000);
	updateStatus(status, 0, num_matches);
	
	for (j = 0; j < num_matches; j++)
	{
		if ((j%100) == 0) updateStatus(status, j, num_matches);
		match = matches[j];
		if (match->sequence_index != seqP->index) continue;
	
		motif = database_motifs[match->motif_index];
		length = strlen(motif);
		
		orig_size = match->size;
		
		//TODO: change orig_size/error distance if initially mismatches can be on the ends.
		
		// Go left
		error_distance = orig_size;
		recursiveMatch(seqP, match, LEFT, &error_distance);
			
		// Go right
		error_distance = orig_size;
		recursiveMatch(seqP, match, RIGHT, &error_distance);
		
		match->errors = match->insertions + match->deletions + match->substitutions;
		match->repeat_count = (double)(match->size - (match->insertions + match->substitutions)) / strlen(database_motifs[match->motif_index]);
		match->purity = 1.0 - ((double)(match->insertions + match->substitutions + match->deletions)
			/ (double) (match->size + match->deletions - match->insertions ));
		
		match->sequence[match->right_sequence+1]= '\0';
	}
	updateStatus(status, num_matches, num_matches);
		
}

void remove_match(int index)
{
	int i;
	deleteMatch(matches[index]);
	for (i=index; i < num_matches - 1; i++)
	{
		matches[i] = matches[i+1];
	}
	num_matches--;
	
}

void findMatches(SeqStructPtr seqP, char** repeats, int numRepeats) {
	
	int i, j;
	char status[1000];
	int *halt_until;
	int divisor = 100;

	memset(status, 0, 1000);
	halt_until = malloc(sizeof(int)*numRepeats);
	memset(halt_until, 0, sizeof(int)*numRepeats);
	
	updateStatus(status, 0, seqP->seqLen);
	
	if (numRepeats > 300) divisor = 10;
	
	for (i = 0; i < seqP->seqLen; i++) {
		if (i%divisor == 0)	updateStatus(status, i, seqP->seqLen);

		for (j=0; j < numRepeats; j++) {
			if (i < halt_until[j]) continue;
			matches[num_matches] = findMatch(seqP, i, j);
			
			// Deal with increasing the size of the matches array
			// if we pass the maximum.
			if (matches[num_matches]) {
				halt_until[j] = matches[num_matches]->end + 1;
				num_matches++;
				
				if (num_matches >= max_matches)
				{
					int step;
					step = (num_matches / MATCH_INCREMENT) + 1;
					matches = realloc(matches, sizeof(MatchPtr) * (MATCH_INCREMENT * step));
					if (!matches) exit(MEMORY_ERROR);
					max_matches = MATCH_INCREMENT * step;
					memset(matches+num_matches, 0, max_matches-num_matches);
				}
			}
				
			
		}
	}
	free(halt_until);
	updateStatus(status, seqP->seqLen, seqP->seqLen);
	
}

void find_compound_repeats(FileStructPtr fsp, Boolean compound_repeat,
Boolean degenerative_repeat, int compound_gap_size)
{
	int i, j;
	int old_num_matches;
	char status[1000];
//	SeqStructPtr seqP;
	MatchPtr cmp = NULL;

	memset(status, 0, 1000);
	old_num_matches = num_matches;
	
	updateStatus(status, 0, old_num_matches);
	
	for (i = 0; i < num_matches; i++) {
		
		MatchPtr m;
		m = matches[i];

		// backup in case the memory hasn't been zeroed, even
		// though it should be
		matches[num_matches] = 0;

		if (i%100 == 0)	updateStatus(status, i, old_num_matches);

		// If the current match is already a compound		
		if (m->num_sub_matches > 0) continue;
		
		if (!check_motif(m)) continue;
		
		if (!cmp)
		{
			cmp = malloc(sizeof(Match));
			memcpy(cmp, m, sizeof(Match));

			cmp->index = num_matches;
			cmp->truncated = 0;
			
			cmp->sub_matches = malloc(sizeof(struct Match**) * SUB_MATCH_INCREMENT);
			cmp->max_sub_matches = SUB_MATCH_INCREMENT;
			
			cmp->sub_matches[cmp->num_sub_matches] = (struct Match*)m;
			cmp->num_sub_matches++;
		}

		for (j=i+1; j < num_matches; j++) {
			// Think I'm making an assumption here that matches
			// will be in order.
			MatchPtr n;
			int diff;
			
			n = matches[j];
			
			// If the current match is already a compound.
			if (n->num_sub_matches > 0) continue;
			
			if (!check_motif(n)) continue;
			
			// What kinds of combined repeats do we want?
			if (compound_repeat && isDegenerate(cmp, n) && !degenerative_repeat) continue;
			if (!compound_repeat && !isDegenerate(cmp, n) && degenerative_repeat) continue;
			
			diff = n->start - cmp->end;
			if (diff < compound_gap_size && diff >= 0)
			{
				int k;
				
				i=j+1;
				m->in_compound = True;
				
				if (cmp->num_sub_matches >= cmp->max_sub_matches)
				{
					int step = 0;
					step = cmp->max_sub_matches / SUB_MATCH_INCREMENT;
					step++;
					cmp->max_sub_matches = SUB_MATCH_INCREMENT * step;
					cmp->sub_matches = realloc(cmp->sub_matches, sizeof(struct Match**) * cmp->max_sub_matches);
					if (!cmp->sub_matches)
						exit(MEMORY_ERROR);
				}
				cmp->sub_matches[cmp->num_sub_matches] = (struct Match*) n;
				cmp->total_gap_sizes += diff;
				cmp->num_sub_matches++;
				
				cmp->insertions += n->insertions;
				cmp->deletions += n->deletions;
				cmp->substitutions += n->substitutions;
				cmp->errors += n->errors;
				
				cmp->repeat_count += n->repeat_count;
				cmp->purity += n->purity;
				
				cmp->end = n->end;
				cmp->size += (n->size + diff - 1);
				
				add_secondary_motif(cmp, n->motif_index);
				
				// Add secondary_motifs
				for (k = 0; k < n->num_secondary_motifs; k++)
				{
					add_secondary_motif(cmp, n->secondary_motifs[k]);
				}
				
				// Join sequences
				for (k = 1; k < diff; k++)
				{
					cmp->sequence[cmp->right_sequence + 1] = '.';
					cmp->right_sequence++;
				}
				strcat(cmp->sequence + cmp->left_sequence, n->sequence+n->left_sequence);
				cmp->right_sequence = strlen(cmp->sequence+cmp->left_sequence) - 1;
				
				matches[num_matches] = cmp;
				
				n->in_compound = True;
				//cmp->is_degenerative = isDegenerate(cmp, n);
			} else if (diff > compound_gap_size)
			{
				j = num_matches;
			}
		}

		// Deal with increasing the size of the matches array
		// if we pass the maximum.
		if (matches[num_matches]) {
			num_matches++;
			if (num_matches >= max_matches)
			{
				int step;
				step = (num_matches / MATCH_INCREMENT) + 1;
				matches = realloc(matches, sizeof(MatchPtr) * (MATCH_INCREMENT * step));
				if (!matches) exit(MEMORY_ERROR);
				max_matches = MATCH_INCREMENT * step;
				memset(matches+num_matches, 0, max_matches-num_matches);
			}
		}

		if (cmp->num_sub_matches == 0) deleteMatch(cmp);
		else cmp->purity = cmp->purity / cmp->num_sub_matches;
		cmp = NULL;
	}
	updateStatus(status, old_num_matches, old_num_matches);
	
	
	
}

Boolean isDegenerate(MatchPtr m, MatchPtr n)
{
	int i;
	

	if (m->motif_index == n->motif_index) return True;
	else 
	{
		int m_length, n_length;
		
		// Andrew indicated he wanted degenerates to be the same length		
		m_length = strlen(database_motifs[m->motif_index]);
		n_length = strlen(database_motifs[n->motif_index]);
		if (m_length != n_length) return False;
		
		for (i = 0; i < m->num_secondary_motifs; i++)
		{
			if (m->secondary_motifs[i] == n->motif_index)
				return True;
		}
	}
	return False;
}

void add_secondary_motif(MatchPtr m, int m_index)
{
	if (m_index == -1)
	{
		printf("Error: repeat-finder tried to add a secondary motif with index -1\n");
	}
	
	
	if (m_index != m->motif_index)
	{
		if (m->num_secondary_motifs == 0) m->max_secondary_motifs = SECONDARY_MOTIF_INCREMENT;
		/*if (!m->secondary_motifs)
		{
			m->secondary_motifs = malloc(sizeof(int) * SECONDARY_MOTIF_INCREMENT);
			m->max_secondary_motifs = SECONDARY_MOTIF_INCREMENT;
			m->secondary_motifs[0] = m_index;
			m->num_secondary_motifs = 1;
			
		}
		else*/
		{
				
			int k = 0;
			while (k >= 0 && k < m->num_secondary_motifs)
			{
				if (m->secondary_motifs[k] == m_index) k = -1;
				else k++;
			}
			if (k >= 0)
			{
				if (m->num_secondary_motifs == m->max_secondary_motifs)
				{
					fprintf(stderr,"Match has too many secondary repeats associated with it."
					"Inform author - joel.pitt@gmail.com\n");
					/*int step = 0;
					step = m->max_secondary_motifs / SECONDARY_MOTIF_INCREMENT;
					step++;
					m->max_secondary_motifs = step * SECONDARY_MOTIF_INCREMENT;
					m->secondary_motifs = realloc(m->secondary_motifs, sizeof(int*)
						* m->max_secondary_motifs);*/
				}
				m->secondary_motifs[m->num_secondary_motifs] =
					m_index;
				m->num_secondary_motifs++;
			}
		}
	}
}

int compare_matches(const void* a, const void* b)
{
	MatchPtr n,m;
	int diff;
	n = *((MatchPtr*)a);
	m = *((MatchPtr*)b);
	
	diff = (n->start - m->start); // lower  is better
	if (!diff) diff = m->size - n->size; // higher
	if (!diff) diff = m->purity - n->purity; //higher
	if (!diff) diff = strlen(database_motifs[n->motif_index]) //lower
		- strlen(database_motifs[m->motif_index]);
	
	return diff;
}

void sort_matches()
{
	qsort(matches, num_matches, sizeof(MatchPtr),&compare_matches);
}

void shortenMatch(MatchPtr m, int amount, int direction)
{
	int length;
	int trim = 0;
	
	length = strlen(database_motifs[m->motif_index]);

	if (amount == 0) trim = 1;
	
	switch (direction)
	{
	case LEFT:
		while (amount > 0 || trim)
		{
			switch (m->sequence[m->left_sequence]){
			case '+':
				m->insertions--;
				m->start++;
				amount--;
				break;
			case '-':
				m->deletions--;
				MOTIF_RIGHT(m->l_motif_pos, length)
				break;
			case '*':
				m->substitutions--;
				MOTIF_RIGHT(m->l_motif_pos, length)
				m->start++;
				amount--;
				break;
			default:
				MOTIF_RIGHT(m->l_motif_pos, length)
				m->start++;
				amount--;
			}
			m->sequence[m->left_sequence] = '\0';
			m->left_sequence++;
			trim = 0;
		}
		break;
	case RIGHT:
		while (amount > 0 || trim)
		{
			switch (m->sequence[m->right_sequence]){
			case '+':
				m->insertions--;
				m->end--;
				amount--;
				break;
			case '-':
				m->deletions--;
				MOTIF_LEFT(m->r_motif_pos,length)
				break;
			case '*':
				m->substitutions--;
				MOTIF_LEFT(m->r_motif_pos,length)
				m->end--;
				amount--;
				break;
			default:
				MOTIF_LEFT(m->r_motif_pos,length)
				m->end--;
				amount--;
			}
			m->sequence[m->right_sequence] = '\0';
			m->right_sequence--;
			
			trim = 0;
		}
		break;
	}
	m->size = (m->end - m->start) + 1;
	m->errors = m->insertions + m->deletions + m->substitutions;
	m->repeat_count = (double) (m->size - (m->insertions + m->substitutions))
		/ strlen(database_motifs[m->motif_index]);
	m->purity = 1.0 - ((double)(m->insertions + m->substitutions + m->deletions)
			/ (double) (m->size + m->deletions - m->insertions ));
			//1.0 - ((double)(m->insertions + m->substitutions) / m->size);
	
	
	
	
	if (m->repeat_count < 2) m->size = 0;
	
}

int check_error_distance(MatchPtr m, int direction)
{
	int i;
	int distance = 0;
	int diff = 0;
	
	int l_distance = 0;
	int l_begin, l_end, l_sbegin, l_send;
		
	int begin = -1, end = -1;
	int sbegin = -1, send = -1;
	int pos = 0;
	
	diff = m->right_sequence - m->left_sequence + 1;
	
	for (i=0; i < diff && distance < (min_error_distance - 1); i++)
	{
	
		switch (direction)
		{
		case LEFT:
			if (IS_MISMATCH(m->sequence[m->left_sequence+i]))
			{
				if (m->sequence[m->left_sequence+i] == '+'
				|| m->sequence[m->left_sequence+i] == '*')
					pos++;
//				return distance;

// if distance is greater than l_distance update saved start and end 
				if (distance > l_distance)
				{
					l_begin = begin;
					l_end = end;
					l_sbegin = sbegin;
					l_send = send;
					l_distance = distance;
				}

				distance = 0;
				begin = -1;
				end = -1;
				sbegin = -1;
				send = -1;
			} else
			{
				if (begin==-1)
				{
					begin=m->start + pos;
					end=m->start + pos;
					sbegin=m->left_sequence+i;
				} else {
					end++;
				}
				
				send=m->left_sequence+i;
				pos++;
				distance++;	
			}
			break;
		case RIGHT:
			if (IS_MISMATCH(m->sequence[m->right_sequence-i]))
			{
				if (m->sequence[m->right_sequence-i] == '+'
				|| m->sequence[m->right_sequence-i] == '*')
					pos++;
//				return distance;

// if distance is greater than l_distance update saved start and end 
				if (distance > l_distance)
				{
					l_begin = begin;
					l_end = end;
					l_sbegin = sbegin;
					l_send = send;
					l_distance = distance;
				}

				distance = 0;
				begin = -1;
				end = -1;
				sbegin = -1;
				send = -1;
			} else
			{
				if (begin==-1)
				{
					begin=m->end - pos;
					end=m->end - pos;
					send=m->right_sequence-i;
				} else {
					begin--;
				}
				
				sbegin=m->right_sequence-i;
				pos++;
				distance++;	
				
			}
			break;
		}
	}
	
	if (distance > l_distance)
	{
		l_begin = begin;
		l_end = end;
		l_sbegin = sbegin;
		l_send = send;
		l_distance = distance;
	}
	
	if (l_distance < (min_error_distance - 1))
	{
		//int temp=0;
		
		if (l_begin != m->start)
		{
     		//fprintf(stderr," ldist=%d\n",l_distance);
     		//temp=1;
			shortenMatch(m, l_begin - m->start, LEFT);
		}
		if (l_end != m->end)
		{
     		//if (temp==0) fprintf(stderr," ldist=%d\n",l_distance);
			shortenMatch(m, m->end - l_end, RIGHT);
		}
		
		
	}
	
	return distance;
	

}

void trim_match(MatchPtr m)
{
	int length;
	length = strlen(database_motifs[m->motif_index]);
			
	while (m->l_motif_pos < (length - 1)
		|| IS_MISMATCH(m->sequence[m->left_sequence])
		|| IS_MISMATCH(m->sequence[m->left_sequence + 1]) )
	{
		shortenMatch(m, 0, LEFT);
	}
	
	while (m->r_motif_pos > 0
		|| IS_MISMATCH(m->sequence[m->right_sequence]) 
		|| IS_MISMATCH(m->sequence[m->right_sequence - 1]) )
	{
		shortenMatch(m, 0, RIGHT);
	}

	m->errors = m->insertions + m->deletions + m->substitutions;
	m->repeat_count = (double) (m->size - (m->insertions + m->substitutions))
		/ strlen(database_motifs[m->motif_index]);
	m->purity = 1.0 - ((double)(m->insertions + m->substitutions + m->deletions)
			/ (double) (m->size + m->deletions - m->insertions ));
}

int  remove_small_matches(int size)
{
	int i=0;
	int count=0;
	i = num_matches - 1;
	while( i >= 0 )
	{
		if (matches[i]->size <= size) {
			count++;
			remove_match(i);
		}
		i--;
//		else i++;
	}
	return count;
}

void dump_matches()
{
	int i = 0;
	
	fprintf(stdout,"\n----------DUMP------------\n");
	
	for (i=0; i < num_matches; i++)
	{
		dump_match(i);
	}
	
	fprintf(stdout,"--------END DUMP----------\n");
}

void dump_match(int i)
{
	MatchPtr m;
	m = matches[i];
	if (m->size < 2) fprintf(stdout,"*");
			
	fprintf(stderr,"%d\t%d\t%d\t%d\t%s\t%s\n",i,m->start,m->end,m->size,database_motifs[m->motif_index], m->sequence + m->left_sequence);
		
}

Boolean are_end_bases_pure(MatchPtr m, int direction, int length)
{
	int pos = 0;

	switch (direction) {
	case LEFT:
		while(pos < length)
		{
			if (IS_MISMATCH(m->sequence[m->left_sequence + pos]))
			{
				return 0;
			}
			pos++;
		}
		break;
	case RIGHT:
		while(pos < length)
		{
			if (IS_MISMATCH(m->sequence[m->right_sequence - pos]))
			{
				return 0;
			}
			pos++;
		}
	
		break;
	}
	
	return 1;
	
}
