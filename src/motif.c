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

#include "fasta.h"
#include "motif.h"

char* motifs[MAX_MOTIFS];
char* database_motifs[MAX_MOTIFS];
int database_motif_sizes[MAX_NUM_SIZES];
int num_sizes = 0;
int database_content[MAX_MOTIFS][4];
int num_database_motifs = 0;
int halt_until[MAX_MOTIFS];
int content[MAX_MOTIFS][4];
int num_motifs = 0;
// Combinatorial variables
int  comb_n;

char bases[] = { 'A', 'T', 'C', 'G' };

void freeMotifs(int collection)
{
	char ** coll;
	int * count;
	if (collection == DATABASE_MOTIFS) {
		coll = database_motifs;
		count = &num_database_motifs;
	} else {
		coll = motifs;
		count = &num_motifs;
	}
		
	int i;
	for (i=0; i < *count; i++)
	{
		free(coll[i]);
	}
	memset(coll, 0, MAX_MOTIFS);
	*count=0;
}

void makeMotifs(int n, int c, char* motif, int collection) {

	char ** coll;
	int * count;
	
	if (n-1 > MAX_NUM_SIZES)
	{
		fprintf(stderr, "Motifs >= 50 are not supported, because that is very large"
		" and I'm using a fixed array to represent what motifs are being searched for.\n");
		return;
	}

	if (collection == DATABASE_MOTIFS) {
		coll = database_motifs;
		count = &num_database_motifs;
		if (n==c){
			database_motif_sizes[num_sizes] = n;
			num_sizes++;
		}
	} else {
		coll = motifs;
		count = &num_motifs;
	}

	if (c > 0) {
		motif[c - 1] = bases[0]; makeMotifs(n, c - 1, motif, collection);
		motif[c - 1] = bases[1]; makeMotifs(n, c - 1, motif, collection);
		motif[c - 1] = bases[2]; makeMotifs(n, c - 1, motif, collection);
		motif[c - 1] = bases[3]; makeMotifs(n, c - 1, motif, collection);
	} else {
		int j;
		Boolean ok = True;

		// Check for repeats within the motif.
		for (j = 1; j <= n / 2 && ok != False; j++) {
			int k;
			Boolean tempOk = True;

			for (k = j; k < n; k += j) {
				if (strncmp(motif, motif + k, j) == 0) {
					if (tempOk == False || k == j) tempOk = False;
				} else tempOk = True;
			}

			ok = tempOk;

		}

		if (ok) {
			coll[*count] = malloc(sizeof(char) * (strlen(motif) + 1));
			strcpy(coll[*count], motif);
			for (j=0; j < strlen(motif); j++) {
				switch (collection) {
				case DATABASE_MOTIFS:
					if (motif[j] == bases[0]) database_content[*count][0]++;
					else if(motif[j] == bases[1]) database_content[*count][1]++;
					else if(motif[j] == bases[2]) database_content[*count][2]++;
					else if(motif[j] == bases[3]) database_content[*count][3]++;
					break;
				case FILTER_MOTIFS:
					if (motif[j] == bases[0]) content[*count][0]++;
					else if(motif[j] == bases[1]) content[*count][1]++;
					else if(motif[j] == bases[2]) content[*count][2]++;
					else if(motif[j] == bases[3]) content[*count][3]++;
					break;
				}
			}
			(*count)++;
			//fprintf(stderr, "\t\tusing motif %s \n", motifs[num_motifs - 1]);

		} //else fprintf(stderr, "skipping motif %s \n", motif);
	}
}

/*void pushMotifCollection(){
	memcpy(stored_motifs, motifs, MAX_MOTIFS);
	num_stored_motifs = num_motifs;
	memset(motifs, 0, MAX_MOTIFS);
	num_motifs = 0;
}
void popMotifCollection(){
	freeMotifs();
	memcpy(motifs, stored_motifs, MAX_MOTIFS);
	num_motifs = num_stored_motifs;
}*/

Boolean check_motif(MatchPtr m)
{
	int j, m_motif_size = 0;
	Boolean motif_present = False;
	if (!m->size) return motif_present;
	
	m_motif_size = strlen(database_motifs[m->motif_index]);
	
	for (j = 0; j < num_motifs && !motif_present; j++) {
		int k;
		if (strcmp(motifs[j], database_motifs[m->motif_index]) == 0)
			motif_present = True;
		else if (strlen(motifs[j]) == m_motif_size){
			for (k = 0; k < m->num_secondary_motifs; k++) {
				if (strcmp(motifs[j], database_motifs[m->secondary_motifs[k]]) == 0)
					motif_present = True;
			}
		}
	}
	return motif_present;
}
