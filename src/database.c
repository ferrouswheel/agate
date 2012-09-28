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

#include "database.h"
#include "match.h"
#include "motif.h"
#include "common.h"

#define BIN_SIZE 50

extern MatchPtr *matches;
extern char* database_motifs[];
extern int database_motif_sizes[MAX_NUM_SIZES];
extern int num_sizes;

extern int num_database_motifs;
extern int num_matches;
extern int max_matches;

extern int partial_ends;

Boolean trim_ends = 1;

//Debugging
Boolean debug_on = 0;

int **bins = NULL;
int *bin_sizes;
int *bin_max;
int bin_width = 0;
int nbins = 0;

int *priority_index;

int filter_database(SeqStructPtr seqP);
void write_database(FileStructPtr fsp);
void trim_database(SeqStructPtr seqP);
int cmp_index_filter(const void* a, const void* b);
int cmp_end_pos(const void* a, const void* b);
void dump_bin(int i);
void dump_bins();

void check_secondary_motifs()
{
	int i,j;
	
	for (i=0; i < num_matches; i++)
	{
		for (j=0; j < matches[i]->num_secondary_motifs; j++)
		{
			if (matches[i]->secondary_motifs[j] == -1)
			{
				printf("invalid secondary motif detected\n");
			}
		}
	}
}

int create_database(FileStructPtr fsp)
{
	SeqStructPtr seqP;
	int count;
	clock_t startTime, timer;

	fprintf(stdout,"Creating database");
	
/*	for (i = 1; i <= size; i++)
	{
		char* motif;
		motif = malloc(sizeof(char) * (i + 1));
		motif[i] = '\0';
					
		makeMotifs(i, i, motif, DATABASE_MOTIFS);
		free(motif);
	}*/
	fprintf(stdout," using %d motifs.\n", num_database_motifs);	
	
	initBuffer(fsp);
	matches = malloc(sizeof(MatchPtr)*MATCH_INCREMENT);
	startTime = clock();
	
	count = 1;
	while (alignSequence(fsp,count) == count)
	{
		if (count > 1) printf("\n Warning: more than one sequence in sequence file, "
		"this program will probably give incorrect results. The author apoligises for"
		" not having time to support this.\n");
		seqP = newSeqStruct(fsp);
		timer = clock();
		
		fprintf(stdout,"Processing sequence %d :\n", count);
		fprintf(stdout,"\tFinding pure matches\t- ");
		findMatches(seqP, database_motifs, num_database_motifs);
		//check_secondary_motifs();
		//dump_matches();
		fprintf(stdout,"\n\tExtending matches\t- ");
		extendMatches(seqP, matches);
		//check_secondary_motifs();
		//dump_matches();		
		if (trim_ends)
		{
			fprintf(stdout,"\n\tTrimming terminal mismatches\t- ");
			trim_database(seqP);
			//check_secondary_motifs();
			//dump_matches();
		}
		fprintf(stdout,"\n\tFiltering matches\t- ");
		filter_database(seqP);
		check_secondary_motifs();
		//dump_matches();
		//if (trim_ends)
		//{
			//fprintf(stdout,"\n\tRetrimming terminal mismatches\t- ");
			//trim_database(seqP);
			//check_secondary_motifs();
			//dump_matches();
		//}
		
		free((void *)seqP);
		count++;
		
		fprintf(stdout,"\n\tTook %4.2f seconds.\n", (float) (clock() - timer) / CLOCKS_PER_SEC);
	}
	fprintf(stdout, "Total time elapsed: %4.2f seconds.\n", (float) (clock() - startTime) / CLOCKS_PER_SEC);
	
	write_database(fsp);
		
	return num_matches;
}

void write_database(FileStructPtr fsp)
{

	FILE* fd;
	char data_fn[1024];
	int i, counter;
	
	strcpy(data_fn, fsp->fname);
	strcat(data_fn, ".database");
	
	fprintf(stdout, "Saving database to %s...", data_fn);
	
	fd = openOutFile(data_fn);
	
	counter=0;
	fprintf(fd, "%d ", num_sizes);
	for (i=0; i < num_sizes; i++)
	{
		fprintf(fd, "%d ", database_motif_sizes[i]);
	}
	fprintf(fd, "\n");
	
	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		int j;
		m = matches[i];
		m->index = i;
		
		fprintf(fd, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %d, %d, %d, %d, %d, %d, %d, %d, %s , %d, ",
		m->sequence_index,
		m->index,
		m->start,
		m->end,
		m->size,
		m->substitutions,
		m->deletions,
		m->insertions,
		m->truncated,
		m->repeat_count,
		m->purity,
		m->motif_index,
		m->content[0],
		m->content[1],
		m->content[2],
		m->content[3],
		m->content[4],
		m->left_sequence,
		m->right_sequence,
		m->sequence+m->left_sequence,
		m->num_secondary_motifs		
		);
		

		for (j=0; j < m->num_secondary_motifs; j++)
		{
			fprintf(fd,"%d, ", m->secondary_motifs[j]);
		}
		fprintf(fd,"\n");
	}
	
	
	fprintf(stdout, "saved!\n");
	
	fclose(fd);
}

int get_bin(int pos)
{
	int result;
	result = pos / bin_width;
	if (result >=nbins) result = nbins - 1;
	else if ( result < 0 ) result = 0;
	return result;
}

void remove_from_bin_pos(int bin, int pos)
{

	
//	if (bins[bin][pos] >= 8261 || bins[bin][pos] <= 8363)
//	{
//		i = 0;
//	}
	
	if (pos >= bin_sizes[bin])
	{
		printf("bin position out of range\n");
		exit(1);
	} else if (bin < 0 || bin > nbins)
	{
		printf("bin number out of range\n");
		exit(1);
	}
	
	bins[bin][pos] = (int) -1;
	
}

void remove_from_bin(int bin, int index)
{
	int i;
	for (i=0; i< bin_sizes[bin]; i++)
	{
		if (bins[bin][i] == index)
		{
			remove_from_bin_pos(bin, i);
			return;
		}
	}
	printf ("\nindex %d not present in bin %d\n", index, bin);
		
}


int reprioritise(MatchPtr m)
{
	MatchPtr n;
	int p_pos, p_i;
	int match_index = 0;
	
	// p_i is the index in priority_index for match m
	p_i = m->p_index;
	match_index = priority_index[p_i];

	if (p_i >= (num_matches - 1)) 
	{
		return p_i;
	}

	if (m->size == 0)
	{
		p_pos = p_i+1;
		while (p_pos < num_matches)
		{
			n = matches[priority_index[p_pos]];
			priority_index[p_pos-1] = priority_index[p_pos];
			n->p_index = p_pos - 1;
			
			p_pos++;
		}
		
	}
	else
	{
		p_pos = p_i + 1;
	
		
		while (cmp_index_filter(&match_index, &(priority_index[p_pos])) > 0)
		{
			n = matches[priority_index[p_pos]];
			priority_index[p_pos-1] = priority_index[p_pos];
			n->p_index = p_pos - 1;
			
			p_pos++;
			if (p_pos >= num_matches) break;
		}
	}
	p_pos--;
	priority_index[p_pos] = match_index;
	m->p_index = p_pos;
	return p_pos;
}

void process_overlaps(int m_i, SeqStructPtr seqP)//, int *end_index)
{
	int i,j,k;
	int m_start_bin = 0;
	int m_end_bin = 0;
	int direction = 0;
	MatchPtr m;
	m = matches[m_i];
	
	m_start_bin= get_bin(m->start);
	m_end_bin= get_bin(m->end);
	
	if (debug_on)
	{
		fprintf(stderr,"======================================"
		"\nProcessing overlaps for:\n\t");
		dump_match(m_i);
		for (i=m_start_bin; i <= m_end_bin; i++)
		{	
			dump_bin(i);
		}
	}
	
	for (i=m_start_bin; i <= m_end_bin; i++)
	{	
		for (j=0; j < bin_sizes[i]; j++)
		{
			int change = 0;
			MatchPtr n;
			int n_index;
			n_index = bins[i][j];
			if (n_index == -1) continue;
			n = matches[n_index];

			if (n->size < 2 || n == m) continue;
			if (matches[j]->sequence_index != seqP->index) continue;

			// Case 1: same location or n completely within m_l
			if (m->start <= n->start &&  n->end <= m->end)
			{
				int s_bin, e_bin;
				add_secondary_motif(m, n->motif_index);
				
				if (debug_on)
				{
					fprintf(stderr,"Overlap found\n\t");
					dump_match(n_index);
				}
				
				n->size = 0;
				
				s_bin = get_bin(n->start);
				e_bin = get_bin(n->end);
				
				for (k = s_bin; k <= e_bin; k++)
				{
					if (k == i) remove_from_bin_pos(k,j);
					else remove_from_bin(k, n_index);
					
				}
				
				change = 1;
			}
						
			// Case 2: No overlap - do nothing
//			if (n->end < m->start || n->start > m->end)
				
					
			// Case 3: the end of n overlaps the start of m
			else if (m->start <= n->end && n->start < m->start)
			{
				int diff;
				int bins_diff, end_bin;
				
				if (debug_on)
				{
					fprintf(stderr,"Overlap found\n\t");
					dump_match(n_index);
				}
				
				
				diff = n->end - m->start + 1;
				
				if (diff >= 2*strlen(database_motifs[n->motif_index])
				&& are_end_bases_pure(n,RIGHT,2*strlen(database_motifs[n->motif_index])))
				{
					add_secondary_motif(m, n->motif_index);
				}
				
				bins_diff = get_bin(n->end - diff);
				end_bin = get_bin(n->end);
				for (k = end_bin; k > bins_diff; k--)
				{
					if (k == i) remove_from_bin_pos(k,j);
					else remove_from_bin(k, n_index);
					
				}
				
				direction = RIGHT;				
				shortenMatch(n, diff, direction);
				change = 1;
				
			}
			
			// Case 4: the end of m overlaps the start of n
			else if (m->start < n->start && m->end >= n->start)
			{
				int diff;
				int bins_diff, start_bin;
				
				if (debug_on)
				{
					fprintf(stderr,"Overlap found\n\t");
					dump_match(n_index);
				}
				
				diff = m->end - n->start + 1;
				if (diff >= 2*strlen(database_motifs[n->motif_index])
				&& are_end_bases_pure(n,LEFT,2*strlen(database_motifs[n->motif_index])))
				{
					add_secondary_motif(m, n->motif_index);
				}
	
				bins_diff = get_bin(n->start+diff);
				start_bin = get_bin(n->start);
				for (k = start_bin; k < bins_diff; k++)
				{
					if (k == i) remove_from_bin_pos(k,j);
					else remove_from_bin(k, n_index);
				}
	
				direction = LEFT;
				shortenMatch(n, diff, direction);
				
				change = 1;
			}
			
			if (change)
			{
				int new_index = 0;
				
				if (n->size > 0)
				{
					// Make sure ends are trimmed.
					trim_match(n);
					// Check that a pure match segment of min.
					// length error distance is present
					check_error_distance(n, direction);
				}
				
				if (debug_on)
				{
					fprintf(stderr,"overlap now\n\t");
					dump_match(n_index);
					fprintf(stderr,"priority before = %d", n->p_index);
				}
				
				new_index = reprioritise(n);
				
				if (debug_on)
				{
					fprintf(stderr,", after = %d\n", new_index);
					dump_bin(i);
					
				}
								
				
			}
			
			
		}
	}
	
}

int* sort_indexes_by_priority(int* size_index, int from)
{
	int i;
		
	qsort(size_index + from, num_matches - from, sizeof(int),&cmp_index_filter);
	
	for (i = from; i < num_matches; i++)
	{
		matches[size_index[i]]->p_index = i;
	}
	
	return size_index;
	
}

int cmp_index_filter(const void* a, const void* b)
{
	// Only this comparison function needs to worry about order
	// and preference of different match variables.
	
	MatchPtr n,m;
	int diff;
	n = matches[*((int*)a)];
	m = matches[*((int*)b)];
	
	diff = m->size - n->size; // higher
	if (!diff) diff = m->purity - n->purity; //higher
	if (!diff) diff = strlen(database_motifs[n->motif_index]) //lower
		- strlen(database_motifs[m->motif_index]);
	if (!diff) diff = n->start - m->start; //lower
	
	return diff;
}

int* sort_indexes_by_end_pos(int* end_index, int from)
{
		
	qsort(end_index+from, num_matches-from, sizeof(int),&cmp_end_pos);
	
	return end_index;
	
}

int cmp_end_pos(const void* a, const void* b)
{

	MatchPtr n,m;
	int diff;
	n = matches[*((int*)a)];
	m = matches[*((int*)b)];
	
	diff = n->end - m->end;
	
	return diff;
}


// Called by put_matches_in_bins to put a match index in a bin
void add_match_to_bin(int m_index, int bin)
{
	bins[bin][bin_sizes[bin]] = m_index;
	bin_sizes[bin]++;
	
	if (bin_sizes[bin] >= (bin_max[bin] - 1))
	{
		int step = 0;
		step = bin_max[bin] / BIN_SIZE;
		step++;
		bins[bin] = (int*) realloc(bins[bin], sizeof(int) * step * BIN_SIZE);
		bin_max[bin] = step * BIN_SIZE;		
		
	}
	
}

void dump_bins()
{
	int i;
	fprintf(stderr,"\n");
	for (i=0; i< nbins; i++)
	{
		dump_bin(i);
	}
}

void dump_bin(int i)
{
	int j;
	fprintf(stderr,"bin %d - ",i);
	for (j=0; j < bin_sizes[i]; j++)
	{
		if (bins[i][j] == -1)
			fprintf(stderr,"* ");
		else
			fprintf(stderr,"%d ", bins[i][j]);
	}
	fprintf(stderr,"\n");
}

// Called to partition matches into regions
int** put_matches_in_bins(SeqStructPtr seqP)
{
	int i;
	
	if (bins){
		printf("Error: Bins already created...\n");
		exit(1);
	}
	
	nbins = num_matches / BIN_SIZE;
	if (nbins == 0) nbins = 1;
	bins = (int**) malloc(sizeof(int*) * nbins);
	bin_sizes = (int*) malloc(sizeof(int) * nbins);
	bin_max = (int*) malloc(sizeof(int) * nbins);
	bin_width = seqP->seqLen / nbins;
	
	if (!bin_width)
	{
		printf("Error: Bin width is zero\n");
		exit(1);
	}
	
	if (!bins || !bin_sizes || !bin_max)
	{
		printf("Error allocating bins\n");	
		exit(1);
	}
	
	for (i=0; i < nbins; i++)
	{
		bins[i] = (int*) malloc(sizeof(int) * BIN_SIZE);
		bin_sizes[i] = 0;
		bin_max[i] = BIN_SIZE;
	}
	
	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		int start_dest_bin = 0;
		int end_dest_bin = 0;
		int j;
		m = matches[i];
		
		start_dest_bin = get_bin(m->start);
		end_dest_bin = get_bin(m->end);;
		
		for (j=start_dest_bin; j <= end_dest_bin; j++)
		{
			add_match_to_bin(i, j);
		}
		
	}
	
	return bins;
}


int filter_database(SeqStructPtr seqP)
{
	char status[1000];
	int i;
	int removed_count = 0;
	//int *end_index;
	
	if (!priority_index) {
		priority_index=malloc(sizeof(int) * num_matches);
		if (!priority_index) fprintf(stderr, "\nFailed to malloc for priority index\n");
	}
	//end_index=malloc(sizeof(int) * num_matches);
	
	for (i=0; i<num_matches; i++)
	{
		priority_index[i]=i;
		//end_index[i]=i;
	}
	sort_indexes_by_priority(priority_index, 0);
	
	put_matches_in_bins(seqP);
	//dump_bins();
	
	memset(status, 0, 1000);
	updateStatus(status, 0, num_matches);
	
	//dump_matches();
	
	for (i=0; i < num_matches; i++)
	{
		int m_index;
		MatchPtr m;

		// Matches must be sorted and priority indexes recalculated each time
		// - optimise in future by only sorting remaining matches to apply filter to (>=i).
		//sort_matches();	
		//sort_indexes_by_end_pos(end_index, 0);
		
		m_index = priority_index[i];
		m = matches[m_index];
		
		if ((i%10) == 0) updateStatus(status, i, num_matches);
		if (m->sequence_index != seqP->index) continue;
		if (m->size <= 2) {
			i = num_matches;
			continue;
		}
		
		//if (m->start >= 182950 && m->start <= 183000)
		//	debug_on = 1;
		//else
		//	debug_on = 0;

		process_overlaps(m_index, seqP);//, end_index);
		
	
		//printf("--\nfiltering %d",i);
		//dump_matches();
		
		
		
		
	}
	
	removed_count += remove_small_matches(2);
	free(priority_index);
	//free(end_index);
	
	updateStatus(status, num_matches, num_matches);
	
	fprintf(stdout, "\n\t\tRemoved %d inferior matches\n",
	removed_count);
	
	return removed_count;
}

void trim_database(SeqStructPtr seqP)
{
	char status[1000];
	int i;
	
	// Specifics of trimming to be determined.
	// Currently removes mismatches in the terminal two bases if the motif size is 2.
	
	// Sort the matches
	sort_matches();	
	
	memset(status, 0, 1000);
	updateStatus(status, 0, num_matches);
	
	for (i=0; i < num_matches; i++)
	{

		MatchPtr m;
		//int length;
		
		m = matches[i];
		
		if ((i%100) == 0) updateStatus(status, i, num_matches);
		if (m->sequence_index != seqP->index) continue;
		
		if (m->size) {
			trim_match(m);
/*			length = strlen(database_motifs[m->motif_index]);
			
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
			m->purity = 1.0 - ((double)(m->insertions + m->substitutions) / m->size);*/
//			m->sequence[m->right_sequence+1] = '\0';
		}
	}
	
	i=0;
	while( i < num_matches )
	{
		if (matches[i]->size <= 2) {
//			removed_count++;
			remove_match(i);
		}
		else i++;
	}
	
	updateStatus(status, num_matches, num_matches);

	
}

int read_database(FileStructPtr fsp)
{
	FILE* fd;
	char data_fn[1024];
	int result = 18;
	int motif_sizes = 0;
	
	matches = malloc(sizeof(MatchPtr)*MATCH_INCREMENT);
	num_matches=0;
		
	strcpy(data_fn, fsp->fname);
	strcat(data_fn, ".database");
	
	fd = fopen(data_fn, "r");
	if (!fd)
	{
		printf("Failed to open existing database. Create one with -d\n");
		exit(FILE_ERROR);
	}

	fprintf(stdout,"Loading database %s...", data_fn);	
	if (!fscanf(fd, "%d", &motif_sizes)) {
		printf("\nFailed to open existing database. Create one with -d\n");
		exit(FILE_ERROR);
	}
	
	while (motif_sizes > 0) {
		int i;
	
		if (!fscanf(fd, "%d", &i)) {
			printf("\nFailed to open existing database. Create one with -d\n");
			exit(FILE_ERROR);
		} else {
			char* motif;
			motif = malloc(sizeof(char) * (i + 1));
			motif[i] = '\0';
						
			makeMotifs(i, i, motif, DATABASE_MOTIFS);
			free(motif);
			motif_sizes--;
		}
	}
	fprintf(stdout," with %d motifs.\n", num_database_motifs);	
	
	result=19;
	while (result == 19)
	{
		MatchPtr m;
		int j, num_secondary = 0;
		
		m = malloc(sizeof(Match));
		memset(m, 0, sizeof(Match));
		
		result = fscanf(fd, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %lf, %lf, %d, %d, %d, %d, %d, %d, %d, %d, ",
		&(m->sequence_index),
		&(m->index),
		&(m->start),
		&(m->end),
		&(m->size),
		&(m->substitutions),
		&(m->deletions),
		&(m->insertions),
		&(m->truncated),
		&(m->repeat_count),
		&(m->purity),
		&(m->motif_index),
		&(m->content[0]),
		&(m->content[1]),
		&(m->content[2]),
		&(m->content[3]),
		&(m->content[4]),
		&(m->left_sequence),
		&(m->right_sequence)
		);
		
		if (result != 19) {
			free (m);
			continue;
		}
		
		fscanf(fd,"%s , %d,",m->sequence,&num_secondary);
		m->left_sequence = 0;
		m->right_sequence = strlen(m->sequence);

		for (j=0; j < num_secondary; j++)
		{
			int m_index = 0;
			fscanf(fd,"%d,", &m_index);
			add_secondary_motif(m, m_index);
		}
		fscanf(fd,"\n");
		
		m->hotspot_index = -1;
//		m->left_sequence = 0;
//		m->right_sequence = strlen(m->sequence)-1;
		matches[num_matches] = m;
		// Deal with increasing the size of the matches array
		// if we pass the maximum.
		num_matches++;
		if (num_matches >= max_matches)
				{
					int step;
					step = (num_matches / MATCH_INCREMENT) + 1;
					matches = realloc(matches, sizeof(MatchPtr) * (MATCH_INCREMENT * step));
					max_matches = MATCH_INCREMENT * step;
				}
			
			
	}
	
	fclose(fd);
	
	//	matches = malloc(sizeof(MatchPtr)*MATCH_INCREMENT);
	return num_matches;
}
