/* #define DEBUG_SPUTNIK 1 */

#include "common.h"
#include "fasta.h"
#include "hotspot.h"
#include "motif.h"

FILE* outf;
enum of {normal, tabular} outputFormat;

void findRepeats(SeqStructPtr seqP, char** repeats, int numRepeats) {
	
	int i, j, k;
	int *sizes;
	char status[1000];
	int stat_length;
	
	char *hotspot_area[] = { "None", "ORF", "0-1kb", "1-2kb", "2-3kb", "3-4kb" };
	
	sizes = malloc(sizeof(int) * numRepeats);
	memset(halt_until, 0, MAX_MOTIFS);
	memset(status, 0, 1000);
	
	updateStatus(status, 0, seqP->seqLen);
	
	for (i = 0; i < seqP->seqLen; i++) {
		
		k = 0;
		
        int repeatCount = 0;
		int start = i + 1, end = i+1;
		int length;
		int longest = 0;
		
		if (i%100 == 0)	updateStatus(status, i, seqP->seqLen);
		
		for (j=0; j < numRepeats; j++) {
			int tempi;
			sizes[j] = 0;
			
			length = strlen(repeats[j]);
			
			if (halt_until[j] <= i) {
				tempi = i;
				
				while (baseNCompare(seqP, tempi, repeats[j], length) == 0)
				{
					tempi += length;
					sizes[j] += length;
					//repeatCount++;
				}
				halt_until[j] = tempi;
			}
		}
		
		for (j=0; j < numRepeats; j++) {
			
			if (sizes[j] > sizes[longest]) {
				
				longest = j;
			} else if (sizes[j] == sizes[longest]) {
				
				if (strlen(repeats[j]) < strlen(repeats[longest])) longest = j;
			}
			
		}
		length = strlen(repeats[longest]);
		repeatCount = sizes[longest] / length;
		end = i + sizes[longest] + 1;
		//i += sizes[longest];
		
		if (repeatCount >= lower &&
		repeatCount <= upper) {
			char *hotspot_name = NULL;
			int intervalIndex = 0;
			
			fprintf(outf, "%s, %d, %d, %d, %d, %d, %f, %f, %f, %f, ", repeats[longest], length, repeatCount, length * repeatCount, start, end,
			(float) content[longest][0] / length, (float) content[longest][1] / length,
			(float) content[longest][2] / length, (float) content[longest][3] / length);
			
			for (j=0; j < numHotspots; j++) {
				int k;
				
				for (k = hotspots[j].numIntervals-1; k >= 0; k--) {
					
					if (start >= hotspots[j].intervals[k][0] && end <= hotspots[j].intervals[k][1])
						
					{
						if (hotspot_name) free(hotspot_name);
						hotspot_name = malloc(sizeof(char) * (strlen(hotspots[j].name) + 1));
						strcpy(hotspot_name, hotspots[j].name);
						intervalIndex = k;
						hotspots[j].counts[k]++;
						
					}
				}
			}
			
			if (hotspot_name) {
				fprintf(outf, "%s, %s,", hotspot_name, hotspot_area[intervalIndex+1]);
				free(hotspot_name);
				hotspot_name = NULL;
			}
			
			fprintf(outf, "\n");
			
		}
	}
	
	updateStatus(status, seqP->seqLen, seqP->seqLen);
	
}

/* Find the frequency of how often a repeat sequence occurs,
in a repeat or otherwise. */
int findFrequency(SeqStructPtr seqP, char* repeat, int length)
{
	int frequency = 0;
	unsigned int i,j;
	
	for (i=0; i < seqP->seqLen; i++) {
		if (baseNCompare(seqP, i, repeat, length) == 0)
		{
			fprintf(outf, "%c%c,%i,%i,%i\n", repeat[0], repeat[1], length/2, i, i + length - 1 );
			frequency++;
			i += length;
		}
		
		
	}
	return frequency;
}

int main(int argc, char* argv[])
{
	SeqStructPtr seqP;
	FileStructPtr fsp;
	
	int count;
	int i;
	int chromosome = 0;
	int k;
	
	clock_t startTime, timer;
	
	if (argc < 2)
	{
        fprintf(stderr,"Usage: %s <options> <fasta format sequence file names>\n\n", argv[0]);
		fprintf(stderr,"Options:\n\t-h <hotspot file>\t- Sorts output into hotspots.\n");
		fprintf(stderr,"\t-n <number>\t- search for motifs of size <number>.\n");
		fprintf(stderr,"\t-m <motif>\t- search for a particular motif.\n");
		fprintf(stderr,"\t-l <lower>\t- find at least <lower> repeats.\n");
		fprintf(stderr,"\t-u <upper>\t- ignore more than <upper> repeats.\n");
		fprintf(stderr,"\t-c <chromosome>\t- specify chromosome number (for hotspots)\n");
		
        exit(1);
	}
	
	for (i=0; i < MAX_MOTIFS; i++) {
		memset(content[i], 0, 4);
	}
	
	/* parse arguments */
	for(i = 1; i < argc; i++)
	{
		if(argv[i][0] == '-')
		{
			if((strcmp(argv[i], "-h") == 0) ||
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
					
					makeMotifs(k, k, motif);
					
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
				if (sscanf(argv[i], "%d", &lower) != 1) exit(2);
			}
			else if((strcmp(argv[i], "-u") == 0) ||
			(strcmp(argv[i], "--upper") == 0)) {
				i++;
				if (sscanf(argv[i], "%d", &upper) != 1) exit(2);
			}
			else if((strcmp(argv[i], "-c") == 0) ||
			(strcmp(argv[i], "--chromosome") == 0)) {
				i++;
				if (sscanf(argv[i], "%d", &chromosome) != 1) exit(2);
			}
		}
		else { // Fasta filename:
			fsp = openFile(argv[i]);
		}
	}
		
	fprintf(stdout,"Number of motifs to check = %d\n", num_motifs);
	
	if (!fsp) {
		printf ("No fasta sequence loaded.\n");
		exit(1);
	}
	
	if (hotspot) {
		parseHotspots(hotspot_fd, chromosome);
	}
	
	initBuffer(fsp);
	
	openOutFile("results.csv");
	
    startTime = clock();
	
	count = 0;
	
	fprintf(outf, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n", "MOTIF", "MOTIF LENGTH", "REPEAT COUNT",
	"TOTAL LENGTH", "START POS", "END POS", "%A", "%T", "%C", "%G", "HS", "HS REGION");
	
	while (alignSequence(fsp,count) == count)
	{
		int j;
		
        timer = clock();
		
        fprintf(stdout,"processing sequence %d\n", count++);
		fprintf(outf,"Sequence %d\n", count);
		
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
	
	close((int)outf);
	
	return 0;
}








