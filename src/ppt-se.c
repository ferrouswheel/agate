
#include "fasta.h"
#include "hotspot.h"
#include "common.h"

FileStructPtr fsp;

extern FILE* hotspot_fd;
extern Boolean hotspot;

FILE* outf;
extern Hotspot hotspots[10000];
extern int numHotspots;
extern char bases[4];

int lower = 1, upper = 1000000000;
float error_freq = 10;
float lowerg = 0.0, upperg = 100.0;

float lowerg_se = 0.0, upperg_se = 100.0;
int lower_se = 1, upper_se = 1000000000;

char* findSymmetry(char* sequence) {
	
	int i;
	int length = strlen(sequence);
	int longest_se = 0;
	int start = 0, end = 0;
	
	char *se;
	
	for (i=0; i < length; i++) {
		int j = 1;
		int se_size = 1;
		
		// Odd sized SE:
		while ((i - j) >= 0 && (i+j) < length) {
			if (sequence[i-j] == sequence[i+j]) {
				j++;
				se_size += 2;
			} else break;
			
		}
		if (se_size > longest_se) {
			longest_se = se_size;
			j--;
			start = i - j;
			end = i + j;
		}
		
		// Even sized SE:
		j = 1;
		se_size = 0;
		while ((i - (j-1)) >= 0 && (i+j) < length) {
			if (sequence[i-(j-1)] == sequence[i+j]) {
				j++;
				se_size += 2;
			} else break;
			
		}
		if (se_size > longest_se) {
			longest_se = se_size;
			j--;
			start = i - (j-1);
			end = i + j;
		}
	}
	
	se = malloc(sizeof(char) * (longest_se + 1));
	strncpy(se, sequence + (start * sizeof(char)), longest_se);
	se[longest_se] = '\0';
	
	return se;
	
	
}

void findPPTs(SeqStructPtr seqP) {
	
    int i, j, k;
    
	char status[1000];
    int stat_length;
	
	char *hotspot_area[] = { "None", "ORF", "0-1kb", "1-2kb", "2-3kb", "3-4kb" };
	
    memset(status, 0, 1000);
	
    sprintf(status, "0 / %d", seqP->seqLen);
    printf("%s", status);
	
    for (i = 0; i < seqP->seqLen; i++) {
		
        k = 0;
		
        int purine = 0, pyramidine = 0;
        char base_string[20];
        int repeatCount = 0;
        int start = i + 1, end = i+1;
        int length;
        int countCG=0;
		
        if (i%100 == 0) {
            stat_length = strlen(status);
            while (k < stat_length) {
				putchar('\b');
				k++;
            }
            sprintf(status, "%d / %d", i, seqP->seqLen);
            printf("%s", status);
        }
		
	   	while (getBaseAt(seqP,i) == 'A' || getBaseAt(seqP,i) == 'G') {
			purine++;
			if (getBaseAt(seqP,i) == 'G') countCG++;
			i++;
      	}
		
        if (purine >= lower) {
            repeatCount = purine;
            strcpy(base_string, "AG");
            end = i;
        } else {
			start = i+1;
			
            while (getBaseAt(seqP,i) == 'C' || getBaseAt(seqP,i) == 'T') {
				// In case last nucleotide was a mismatch.
				pyramidine++;
				if (getBaseAt(seqP,i) == 'C') countCG++;
				i++;
				
			}
			
            end = i;
            repeatCount = pyramidine;
            strcpy(base_string, "CT");
		}
		i--;
		length = end - start + 1;
		
		if (repeatCount >= lower && repeatCount <= upper) {
			float ratioGC = 0.0;
			
			ratioGC = (float) countCG / (float) length * 100.0;
			if (ratioGC >= lowerg && ratioGC <= upperg) {
				char *hotspot_name = NULL;
				int intervalIndex = 0;
				
				char *sequence = NULL;
				char *se = NULL;
				sequence = malloc(sizeof(char) * (length + 1));
				//strncpy(sequence, seqP->seqStr + ((start - 1) * sizeof(char)), length);
				baseNCopy(seqP, sequence, start - 1, length);
				sequence[length] = '\0';
				
				se = findSymmetry(sequence);
				if (strlen(se) > lower_se && strlen(se) < upper_se) {
					int p = 0;
					int countGC_se = 0;
					float ratioGC_se = 0.0;
					for (p = 0; p < strlen(se); p++) {
						if (se[p] == 'G' || se[p] == 'C')
							countGC_se++;
					}
					ratioGC_se = (float) countGC_se / (float) strlen(se) * 100.0;
					
					if (ratioGC_se > lowerg_se && ratioGC_se < upperg_se) {
						
						fprintf(outf, "%s, %s, %d, %d, %d, %f, %s, %d, %f,", base_string, sequence, length, start, end, ratioGC, se, strlen(se), ratioGC_se);
						
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
			}
			
        }
    }
	
	k = 0;
	
	
	
	stat_length = strlen(status);
	while (k < stat_length) {
		putchar('\b');
		k++;
	}
	sprintf(status, "%d / %d", seqP->seqLen, seqP->seqLen);
	printf("%s", status);
	
	
}

int main(int argc, char* argv[])
{
	SeqStructPtr seqP;
	int count;
	int i;
	int chromosome = 0;
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
        fprintf(stderr,"\t-lg <lower GC%>\t- specify lowest GC content to report\n");
        fprintf(stderr,"\t-ug <upper GC%>\t- specify highest GC content to report\n");
		fprintf(stderr,"\t-selg <lower SE GC%>\t- specify lowest SE GC content to report\n");
		fprintf(stderr,"\t-seug <upper SE GC%>\t- specify highest SE GC content to report\n");
		fprintf(stderr,"\t-sel <lower>\t- specify shortest SE length to report\n");
		fprintf(stderr,"\t-seu <upper>\t- specify longest SE length to report\n");
		
        exit(1);
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
			else if(strcmp(argv[i], "-lg") == 0)
			{
				i++;
				if (sscanf(argv[i], "%f", &lowerg) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-ug") == 0)
			{
				i++;
				if (sscanf(argv[i], "%f", &upperg) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-e") == 0)
			{
				i++;
				if (sscanf(argv[i], "%f", &error_freq) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-seu") == 0)
			{
				i++;
				if (sscanf(argv[i], "%d", &upper_se) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-sel") == 0)
			{
				i++;
				if (sscanf(argv[i], "%d", &lower_se) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-seug") == 0)
			{
				i++;
				if (sscanf(argv[i], "%f", &upperg_se) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-selg") == 0)
			{
				i++;
				if (sscanf(argv[i], "%f", &lowerg_se) != 1) exit(2);
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
	
	if (fsp->fd < 0) {
		printf ("No fasta sequence loaded.\n");
		exit(1);
    }
	
	if (hotspot) {
		
		parseHotspots(hotspot_fd, chromosome);
		
	}
	
	
	
	initBuffer(fsp);
	
    openOutFile("ppt-se.csv");
	
    startTime = clock();
	
	count = 1;
	
	
	fprintf(outf, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n", "TYPE", "SEQUENCE", "LENGTH", "START POS", "END POS", "GC (%)", "SE", "SE LENGTH", "SE GC (%)", "HS", "HS REGION");
	
	while (alignSequence(fsp,count) == count)
    {
        int j;
		
		seqP = newSeqStruct(fsp);
        timer = clock();
		
        fprintf(stdout,"processing sequence %d\n", count);
        fprintf(outf,"Sequence %d\n", count);
		
        findPPTs(seqP);
		
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
		count++;
	}
	fprintf(stdout, "Total time elapsed: %4.2f seconds.\n", (float) (clock() - startTime) / CLOCKS_PER_SEC);
	
	
	close((int)outf);
	
	return 0;
}




