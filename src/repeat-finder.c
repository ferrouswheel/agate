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

/* #define DEBUG_SPUTNIK 1 */

#include "common.h"
#include "fasta.h"
#include "hotspot.h"
#include "motif.h"
#include "match.h"
#include "database.h"

/* External variables from hotspot module */
extern Hotspot hotspots[MAX_HOTSPOTS];
extern int numHotspots;

extern FILE* hotspot_fd;
extern Boolean hotspot;

/* External variables from motif module */
extern char* motifs[MAX_MOTIFS];
extern char* database_motifs[MAX_MOTIFS];
extern Boolean database_motif_sizes[50];
extern int halt_until[MAX_MOTIFS];
extern int content[MAX_MOTIFS][4];
extern int num_motifs;
extern int num_database_motifs;
extern char bases[4];

/* External variables from match module */
extern int num_matches;
extern int max_matches;
extern MatchPtr *matches;
extern int min_error_distance;
extern Boolean trim_ends;

/* Parameters for constraining search */
double lower = 2.0, upper = 1000000000.0;
double lowerp = 0.0, upperp = 1.0;
double lower_gc = 0.0, upper_gc = 1.0;
int flank_length = 50;
int database_size = 0;

/* Parameters for exporting to fasta files */
#define FASTA_INDEX_STEP 1024
int* fasta_indexes = NULL;
int num_fasta_indexes = 0;
char* fasta_fn_out = 0;
int fasta_seq_index = 1;

/* Parameters for compound repeats */
int compound_gap_size = 5;
Boolean compound_repeat = False;
Boolean degenerative_repeat = False;

/* Parameter for the format to dump search results in */
Boolean fasta_results = False;
Boolean fasta_mask = False;
Boolean fasta_mask_all = False;

FILE* outf;
FILE* fastaf;

void tallyBaseContent(FileStructPtr fsp, MatchPtr* matches)
{
	int i,j;
	SeqStructPtr seqP;
	int seq_count = 1;

	alignSequence(fsp,seq_count);
	seqP = newSeqStruct(fsp);
	
	for (j = 0; j < num_matches; j++)
	{
		MatchPtr m;
		m = matches[j];
		
		if (m->sequence_index != seq_count) {
			free(seqP);
			seq_count=m->sequence_index;
			alignSequence(fsp,seq_count);
			seqP = newSeqStruct(fsp);
		}
		
		for (i = m->start ; i <= m->end; i++)
		{
			switch (getBaseAt(seqP,i))
			{
			case 'A':
				m->content[0]++;
				break;
			case 'T':
				m->content[1]++;
				break;
			case 'C':
				m->content[2]++;
				break;
			case 'G':
				m->content[3]++;
				break;
			default:
				m->content[4]++;
				break;
			}
		}
	}
	
	
}

void findHotspot(FileStructPtr fsp)
{
	SeqStructPtr seqP;
	int seq_count = 1;
		
	int j,i;
	
	alignSequence(fsp,seq_count);
	seqP = newSeqStruct(fsp);
	
	for (i=0; i < num_matches; i++) {
		MatchPtr m;
			
		m = matches[i];
		if (m->size <= 2) continue;
		if (compound_repeat || degenerative_repeat)
		{
			if (m->num_sub_matches == 0) continue;
		}
	
		for (j=0; j < numHotspots; j++) {
			int k;
			
			//matches[i]->hotspot_index = -1;
			
			for (k = hotspots[j].numIntervals-1; k >= 0; k--) {

				if (m->start >= hotspots[j].intervals[k][0] && m->end <= hotspots[j].intervals[k][1])
				{
					m->hotspot_index = j;
					m->hotspot_interval = k;
					
					hotspots[j].counts[k]++;
					
				}
			}
		}
	}
	
	for (j=0; j < numHotspots; j++)
	{
		for (i=0; i < hotspots[j].numIntervals; i++)
		{
			int k;
			for (k=hotspots[j].intervals[i][0];
				k <= hotspots[j].intervals[i][1]; k++)
			{
				switch (getBaseAtPos(seqP, k))
				{
					case 'A':
						hotspots[j].base_counts[i][0]++;
						break;
					case 'T':
						hotspots[j].base_counts[i][1]++;
						break;
					case 'C':
						hotspots[j].base_counts[i][2]++;
						break;
					case 'G':
						hotspots[j].base_counts[i][3]++;
						break;
					default:
						hotspots[j].base_counts[i][4]++;
				}
				
			}
				
		}
	}
	
}

void mark_matches_without_motifs()
{
	int i;
	
	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		m = matches[i];
		if (!degenerative_repeat && !compound_repeat)
		{
			int j;
			int motif_present = 0;
			for (j = 0; j < num_motifs && !motif_present; j++) {
				if (strcmp(database_motifs[m->motif_index], motifs[j]) == 0)
				{
					motif_present = 1;
				}
			}
			if (!motif_present) m->size = 0;
			
		} else {
			if (!check_motif(m)) m->size = 0;
		}
	}
}


void filterMatches(MatchPtr* matches)
{

	int i;
	
	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		m = matches[i];
		
		double match_gc = (double) (m->content[2] + m->content[3]) / (double) m->size;

		if (m->repeat_count < lower || m->repeat_count > upper)
		{
			m->size = 0;
		}
		else if (m->purity < lowerp || m->purity > upperp)
		{
			m->size = 0;
		}
		else if (match_gc < lower_gc || match_gc > upper_gc)
		{
			m->size = 0;

		/*else {
			int j;
			Boolean motif_present = False;
			for (j = 0; j < num_motifs && !motif_present; j++) {
				int k;
				MatchPtr n;
				if (strcmp(motifs[j], database_motifs[m->motif_index]) == 0)
					motif_present = True;
				else {
					for (k = 0; k < m->num_sub_matches; k++) {
						n = (MatchPtr) m->sub_matches[k];
						if (strcmp(motifs[j], database_motifs[n->motif_index]) == 0)
							motif_present = True;
					}
				}
				
				
					
			}
			if (!motif_present) m->size = 0;*/
				
		}
		
	}
		
}


void dumpResults(FileStructPtr fsp, MatchPtr* matches)
{
	char *hotspot_area[] = { "None", "ORF", "0-1kb", "1-2kb", "2-3kb", "3-4kb" };	
	int i;
	SeqStructPtr seqP;
	int seq_count = 1;
	
	char* flank_left;
	char* flank_right;
	char buffer[MAX_MATCH_LENGTH];	
	
	flank_left = malloc(sizeof(char) * (flank_length + 1));
	flank_right = malloc(sizeof(char) * (flank_length + 1));
	
	fprintf(outf, "Seq,ID,Motif, Motif Length, Repeat Count, Total Length, Start Pos,"
" End Pos, Purity, A, T, C, G, Insertions, Deletions, Substitutions,"
" Raw Sequence, Sequence, Left flank, Right flank, Truncated?, InCompound?, Compound?, Subsequences, Secondary Motifs, Hotspot, Hotspot Region\n");
	
	alignSequence(fsp,seq_count);
	seqP = newSeqStruct(fsp);

	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		int j;
		
		m = matches[i];
		if (m->size <= 2) continue;
		if (compound_repeat || degenerative_repeat)
		{
			if (m->num_sub_matches == 0) continue;
		}
		
		if (m->sequence_index != seq_count) {
			free(seqP);
			seq_count=m->sequence_index;
			alignSequence(fsp,seq_count);
			seqP = newSeqStruct(fsp);
		}
		
		baseNCopy(seqP, flank_left, m->start - flank_length, flank_length);
		flank_left[flank_length]='\0';
		baseNCopy(seqP, flank_right, m->end + 1, flank_length);
		flank_right[flank_length]='\0';
		baseNCopy(seqP, buffer, m->start, m->end - m->start + 1);
		buffer[m->end - m->start + 1]='\0';
		
		m->sequence[m->right_sequence + 1] = '\0';
				
		fprintf(outf, "%d,%d,%s, %d, %f, %d, %d, %d, %f, %f, %f, %f, %f, %d, %d, %d, %s, %s, %s, %s, %d, %d, %d, ",
		m->sequence_index,
		m->index,
		database_motifs[m->motif_index],
		strlen(database_motifs[m->motif_index]), m->repeat_count, m->size,
		m->start + 1, m->end + 1, m->purity,
		(float) m->content[0] / m->size,
		(float) m->content[1] / m->size,(float) m->content[2] / m->size,
		(float) m->content[3] / m->size,
		m->insertions,
		m->deletions,
		m->substitutions,
		buffer,
		m->sequence + m->left_sequence,
		flank_left,
		flank_right,
		m->truncated,
		m->in_compound,
		(m->num_sub_matches > 0));
		
		for (j=0; j < m->num_sub_matches; j++)
		{
			fprintf(outf, "%d", ((MatchPtr)m->sub_matches[j])->index);
			if (j < (m->num_sub_matches - 1)) fprintf(outf, ".");
		}
		fprintf(outf, ", ");
		
		for (j=0; j < m->num_secondary_motifs; j++)
		{
			if (strlen(database_motifs[m->secondary_motifs[j]]) == strlen(database_motifs[m->motif_index]))
				fprintf(outf, "%s ", database_motifs[m->secondary_motifs[j]]);
		}
		fprintf(outf, ", ");
		
		if (m->hotspot_index != -1)
		{
			fprintf(outf, "%s, %s,", hotspots[m->hotspot_index].name, hotspot_area[m->hotspot_interval + 1]);
		}
		
		fprintf(outf, "\n");
	}
}

void dumpFasta(FileStructPtr fsp, MatchPtr* matches)
{
	int i;
	SeqStructPtr seqP;
	int seq_count = 1;
	
	char* buffer;
	char* fname;
	
	fname = strdup(fsp->fname);
	for (i=0; i < strlen(fname); i++)
	{
		if (fname[i] == '.')
			fname[i] = '\0';
	}
	buffer = malloc(sizeof(char) * ((2 * flank_length) + MAX_MATCH_LENGTH + 1));

	//fprintf(outf, "Seq,ID,Motif, Motif Length, Repeat Count, Total Length, Start Pos,"
//" End Pos, Purity, A, T, C, G, Insertions, Deletions, Substitutions,"
//" Raw Sequence, Sequence, Left flank, Right flank, Truncated?, InCompound?, Compound?, Subsequences, Secondary Motifs, Hotspot, Hotspot Region\n");
	
	alignSequence(fsp,seq_count);
	seqP = newSeqStruct(fsp);

	for (i=0; i < num_matches; i++)
	{
		MatchPtr m;
		int j;
		int column = 0;
		int seq_length;
		
		m = matches[i];
		if (m->size <= 2) continue;
		if (compound_repeat || degenerative_repeat)
		{
			if (m->num_sub_matches == 0) continue;
		}
		
		if (m->sequence_index != seq_count) {
			free(seqP);
			seq_count=m->sequence_index;
			alignSequence(fsp,seq_count);
			seqP = newSeqStruct(fsp);
		}
		
		seq_length = m->end - m->start + 1 + (2 * flank_length);
		baseNCopy(seqP, buffer, m->start - flank_length, seq_length);
		buffer[seq_length]='\0';
		
		if (fasta_mask_all)
		{
			for (j=flank_length; j < (seq_length - flank_length); j++)
				buffer[j]='N';
		} else if (fasta_mask)
		{
			int counter = m->left_sequence;
			for (j=flank_length; counter <= m->right_sequence
				&& j < (seq_length - flank_length); counter++)
			{
				if (IS_MISMATCH(m->sequence[counter]))
				{
					switch (m->sequence[counter])
					{
					case '*':
						j++; 
						break;
					case '-':
						break;
					case '+':
						j++;
						break;
					}
				} else {
					buffer[j]='N';
					j++;
				}
				
			} 
		}

		fprintf(fastaf, ">%s|id %d|pos %d to %d|m %s|rc %.2f|purity %.2f|flank %d\n",
		fname, m->index, m->start+1, m->end+1, database_motifs[m->motif_index],
		m->repeat_count, m->purity, flank_length);

		for (j = 0; j < seq_length; j++)
		{
			fputc(buffer[j], fastaf);
			column++;
			if (column >= 80)
			{
				fputc('\n', fastaf);
				column = 0;
			}
		}
		
		fprintf(fastaf, "\n");
	}
}


int parse_args(int argc, char* argv[], FileStructPtr *fsp)
{
	int chromosome = 0;
	int i;
	
	if (argc < 2)
	{
		fprintf(stderr,"Usage: %s <options> <fasta format sequence file name>\n\n", argv[0]);
		fprintf(stderr,"Options:\n");
		fprintf(stderr,"    -----Database--------------------\n");
		fprintf(stderr,"\t-d <maxsize>\t- (re)create database, count motifs up to <maxsize>.\n");
		fprintf(stderr,"\t-e <dist>\t- (used with -d) Allow errors every <dist> nucleotides (default=5)\n");
		fprintf(stderr,"    -----FASTA export-----------------\n");
		fprintf(stderr,"\t-fx \t- Output filtered results as a fasta file (results.fsa)\n");
		fprintf(stderr,"\t-fm \t- (used with -fx) mask repeat, but not mismatches, with N\n");
		fprintf(stderr,"\t-fa \t- (used with -fx) mask entire repeat, including mismatches, with N\n");
		fprintf(stderr,"\t-i <indexes>\t- output the match with the index in FASTA form (-1 outputs *all* matches).\n");
		fprintf(stderr,"\t-is <sequence>\t- The index of the sequence in the original fasta file. \n");
		fprintf(stderr,"\t-io <file>\t- (used with -i) Filename of fasta file to create.\n");
		fprintf(stderr,"    -----Hotspots---------------------\n");		
		fprintf(stderr,"\t-h <file>\t- Sorts output into hotspots defined in <file>.\n");
		fprintf(stderr,"\t-c <chromosome>\t- specify chromosome number (for hotspots)\n");
		fprintf(stderr,"    -----Compound repeats-------------\n");
		fprintf(stderr,"\t-gap \t- Specify compound repeat gap size.\n");
		fprintf(stderr,"\t-cr \t- Enable compound repeats.\n");
		fprintf(stderr,"\t-dr \t- Enable degenerative repeats.\n");
		fprintf(stderr,"    -----Filtering--------------------\n");
		fprintf(stderr,"\t-f <size>\t- Output the <size> flanking bases to a match (Also can be used with fasta export)\n");
		fprintf(stderr,"\t-n <number>\t- search for motifs of size <number>.\n");
		fprintf(stderr,"\t-m <motif>\t- search for a particular motif.\n");
		fprintf(stderr,"\t-l <lower>\t- find at least <lower> repeats.\n");
		fprintf(stderr,"\t-u <upper>\t- ignore more than <upper> repeats.\n");
		fprintf(stderr,"\t-lg <lower>\t- lowest gc%% repeats to find\n");
		fprintf(stderr,"\t-ug <upper>\t- greatest gc%% repeats to find\n");
		fprintf(stderr,"\t-lp <lower>\t- lowest purity repeats to find\n");
		fprintf(stderr,"\t-up <upper>\t- greatest purity repeats to find\n");
		fprintf(stderr,"\t\t\t (Note: range of purities is constrained by error distance)\n");
						
		exit(1);
	}
	
	/* parse arguments */
	for(i = 1; i < argc; i++)
	{
		if(argv[i][0] == '-')
		{
			if((strcmp(argv[i], "-d") == 0) ||
			(strcmp(argv[i], "--database") == 0)) {
				int k;
				i++;
				if (sscanf(argv[i], "%d", &k) == 1)
				{
					if (k>database_size) database_size = k;
					char* motif;
					motif = malloc(sizeof(char) * (k + 1));
					motif[k] = '\0';
					
					makeMotifs(k, k, motif, DATABASE_MOTIFS);
					free(motif);
				}
					
			}
			else if((strcmp(argv[i], "-h") == 0) ||
			(strcmp(argv[i], "--hotspot") == 0)) {
				hotspot = True;
				i++;
				hotspot_fd = (FILE*) fopen(argv[i], "r");
				if (chromosome == 0) chromosome = 1;
				if (hotspot_fd == NULL) {
					printf("Couldn't open hotspots file: %s", argv[i]);
					exit(2);
				}
				
			}
			else if((strcmp(argv[i], "-n") == 0) ||
			(strcmp(argv[i], "--number") == 0)) {
				int k;
				i++;
				if (sscanf(argv[i], "%d", &k) == 1) {
					char* motif;
					motif = malloc(sizeof(char) * (k + 1));
					motif[k] = '\0';
					
					makeMotifs(k, k, motif, FILTER_MOTIFS);
					free(motif);
					
				}
			}
			else if((strcmp(argv[i], "-m") == 0) ||
			(strcmp(argv[i], "--motif") == 0)) {
				int j;
				i++;
				motifs[num_motifs] = malloc(sizeof(char) * (strlen(argv[i]) + 1));
				strcpy(motifs[num_motifs], argv[i]);
				for (j=0; j < strlen(motifs[num_motifs]); j++) {
					if (motifs[num_motifs][j] == bases[0]) content[num_motifs][0]++;
					else if (motifs[num_motifs][j] == bases[1]) content[num_motifs][1]++;
					else if (motifs[num_motifs][j] == bases[2]) content[num_motifs][2]++;
					else if (motifs[num_motifs][j] == bases[3]) content[num_motifs][3]++;
				}
				
				num_motifs++;
				
			}
			else if((strcmp(argv[i], "-l") == 0) ||
			(strcmp(argv[i], "--lower") == 0)) {
				i++;
				if (sscanf(argv[i], "%lf", &lower) != 1) exit(2);
			}
			else if((strcmp(argv[i], "-u") == 0) ||
			(strcmp(argv[i], "--upper") == 0)) {
				i++;
				if (sscanf(argv[i], "%lf", &upper) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-lp") == 0) {
				i++;
				if (sscanf(argv[i], "%lf", &lowerp) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-up") == 0) {
				i++;
				if (sscanf(argv[i], "%lf", &upperp) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-lg") == 0) {
				i++;
				if (sscanf(argv[i], "%lf", &lower_gc) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-ug") == 0) {
				i++;
				if (sscanf(argv[i], "%lf", &upper_gc) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-e") == 0) {
				i++;
				if (sscanf(argv[i], "%d", &min_error_distance) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-f") == 0) {
				i++;
				if (sscanf(argv[i], "%d", &flank_length) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-io") == 0) {
				i++;
				fasta_fn_out = malloc(sizeof(char) * (strlen(argv[i]) + 1));
				strcpy(fasta_fn_out, argv[i]);
			}
			else if(strcmp(argv[i], "-i") == 0) {
				int index;
				i++;
				if (sscanf(argv[i], "%d", &index) != 1) exit(2);
				if ((num_fasta_indexes % FASTA_INDEX_STEP) == 0) {
					int steps = 0;
					steps = num_fasta_indexes / FASTA_INDEX_STEP;
					steps++;
					if (steps > 1)
						fasta_indexes = realloc(fasta_indexes, sizeof(int) * FASTA_INDEX_STEP * steps);
					else
						fasta_indexes = malloc(sizeof(int) * FASTA_INDEX_STEP * steps);
				}
				fasta_indexes[num_fasta_indexes] = index;
				num_fasta_indexes++;
				
			}
			else if(strcmp(argv[i], "-is") == 0) {
				i++;
				if (sscanf(argv[i], "%d", &fasta_seq_index) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-gap") == 0) {
				i++;
				if (sscanf(argv[i], "%d", &compound_gap_size) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-cr") == 0) {
				compound_repeat = True;
			}
			else if(strcmp(argv[i], "-dr") == 0) {
				degenerative_repeat = True;
			}
			else if((strcmp(argv[i], "-c") == 0) ||
			(strcmp(argv[i], "--chromosome") == 0)) {
				i++;
				if (sscanf(argv[i], "%d", &chromosome) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-fx") == 0) {
				fasta_results = True;
			}
			else if(strcmp(argv[i], "-fa") == 0) {
				fasta_mask_all = True;
			}
			else if(strcmp(argv[i], "-fm") == 0) {
				fasta_mask = True;
			}
			
		}
		else { // Fasta filename:
			*fsp = openFile(argv[i]);
		}
	}
	
	return chromosome;
}

int export_to_fasta(FileStructPtr fsp)
{
		FILE* fasta_outf;
		int i;
		
		SeqStructPtr seqP;
		int seq_count = fasta_seq_index;
	
		char* rep_seq;
		rep_seq = malloc(sizeof(char) * (2 * flank_length + MAX_MATCH_LENGTH + 1));
	
		alignSequence(fsp,seq_count);
		seqP = newSeqStruct(fsp);
		
		if (!fasta_fn_out) fasta_fn_out = "repeats.fsa";
		fasta_outf = openOutFile(fasta_fn_out);
		
		for (i=0; i < num_fasta_indexes; i++)
		{
			int j;
			int index;
			
			index = fasta_indexes[i];
			if (index == -1) {
				for (j=0; j < num_matches; j++)
				{
					int length;
					index = matches[j]->index;
					length = matches[j]->end - matches[j]->start + 1 + (2 * flank_length);
					baseNCopy(seqP, rep_seq, matches[j]->start - flank_length, length);
					rep_seq[length] = '\0';
				
					fprintf(fasta_outf, ">Match %d, motif %s, flanking size %d\n", index, database_motifs[matches[j]->motif_index], flank_length);
					fprintf(fasta_outf, "%s\n", rep_seq);
				}
			}
			else {
				for (j=0; j < num_matches; j++)
				{
					if (matches[j]->index == index) {
						int length;
						length = matches[j]->end - matches[j]->start + 1 + (2 * flank_length);
						baseNCopy(seqP, rep_seq, matches[j]->start - flank_length, length);
						rep_seq[length] = '\0';
				
						fprintf(fasta_outf, "> Match %d, motif %s, flanking size %d\n", index, database_motifs[matches[j]->motif_index], flank_length);
						fprintf(fasta_outf, "%s\n", rep_seq);
					}
				}
			}
		}
		
		fclose(fasta_outf);
		return 0;
		
	}

int main(int argc, char* argv[])
{
	FileStructPtr fsp;
	
	int count;
	int i;
	int chromosome;
	clock_t startTime;
	
	// Clear all motifs
	for (i=0; i < MAX_MOTIFS; i++) {
		memset(content[i], 0, 4);
	}
	
	// clear motif size array
	memset(database_motif_sizes, 0, MAX_NUM_SIZES * sizeof(int));
	
	chromosome=parse_args(argc, argv, &fsp);
	
	if (database_size > 0) create_database(fsp);
	else read_database(fsp);
		
	if (!fsp) {
		printf ("No fasta sequence loaded.\n");
		exit(1);
	}
	
	// Give size = 0 to matches that do not have with specified motifs
	mark_matches_without_motifs();
	
	if (degenerative_repeat || compound_repeat)
	{
		printf("\tFinding compound repeats\t- ");
		find_compound_repeats(fsp, compound_repeat, degenerative_repeat, compound_gap_size);		
		printf("\n");
	}
	
	if (num_fasta_indexes) {
		export_to_fasta(fsp);
		exit(0);
	}
		
	if (num_motifs > 0)
	{
		fprintf(stdout,"Number of motifs to check = %d\n", num_motifs);
	}
	else
	{
		fprintf(stderr,"No motifs specified for searching (use -n or -m options)\n");
		exit(4);
	}
	
	if (hotspot) parseHotspots(hotspot_fd, chromosome);

	if (fasta_results)
		fastaf = openOutFile("results.fsa");
	
	outf = openOutFile("results.csv");
	
	startTime = clock();
	
	count = 1;
	
	tallyBaseContent(fsp, matches);
	filterMatches(matches);
	findHotspot(fsp);
	
	if (fasta_results)
		dumpFasta(fsp, matches);
	
	dumpResults(fsp, matches);
	
	if (hotspot) {
		int j;
			fprintf(outf, "\nTotal Counts and %%gc in each hotspot region:\n");
			fprintf(outf, "Hotspot, ORF, 0-1kb, 1-2kb, 2-3kb, 3-4kb, "
			"ORF %%gc, 0-1kb %%gc, 1-2kb %%gc, 2-3kb %%gc, 3-4kb %%gc,\n");
			
			for (j=0; j < numHotspots; j++) {
				int k;
				
				fprintf(outf,"%s", hotspots[j].name);
				for (k = 0; k < 5; k++) {
					if (k < hotspots[j].numIntervals)
						fprintf(outf,",%d", hotspots[j].counts[k]);
					else fprintf(outf,",");
				}
				for (k = 0; k < 5; k++) {
					if (k < hotspots[j].numIntervals) {
						int total;
						int gc;
						total = hotspots[j].base_counts[k][0] +
							hotspots[j].base_counts[k][1] +
							hotspots[j].base_counts[k][2] +
							hotspots[j].base_counts[k][3] +
							hotspots[j].base_counts[k][4];
						gc = hotspots[j].base_counts[k][2] +
							hotspots[j].base_counts[k][3];
						fprintf(outf,",%f", (double) gc / (double) total );
					} else {
						fprintf(outf,",");
					}
				}
				fprintf(outf,"\n");
				
			}
		}
/*	while (alignSequence(fsp,count) == count)
	{
		int j;
		
		seqP = newSeqStruct(fsp);
		timer = clock();
		
		fprintf(stdout,"processing sequence %d\n", count);
		fprintf(outf,"Sequence %d\n", count++);
		
		findRepeats(seqP, motifs, num_motifs);
		
		free((void *)seqP);
		
		if (hotspot) {
			fprintf(outf, "\nTotal Counts in each hotspot region:\n");
			fprintf(outf, "Hotspot, ORF, 0-1kb, 1-2kb, 2-3kb, 3-4kb\n");
			
			for (j=0; j < numHotspots; j++) {
				int k;
				
				fprintf(outf,"%s", hotspots[j].name);
				for (k = 0; k < hotspots[j].numIntervals; k++) {
					fprintf(outf,",%d", hotspots[j].counts[k]);
				}
				fprintf(outf,"\n");
			}
		}
		
		fprintf(stdout,"...took %4.2f seconds.\n", (float) (clock() - timer) / CLOCKS_PER_SEC);
	}
	fprintf(stdout, "Total time elapsed: %4.2f seconds.\n", (float) (clock() - startTime) / CLOCKS_PER_SEC);
	*/
	close((int)outf);
	
	return 0;
}







