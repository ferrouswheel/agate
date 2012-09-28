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
#include "hotspot.h"
#include "common.h"

extern FILE* hotspot_fd;
extern Boolean hotspot;
extern Hotspot hotspots[MAX_HOTSPOTS];
extern int numHotspots;

extern char bases[4];

FILE* outf;

float error_gap = 10;
float lowerg = 0.0, upperg = 100.0;
int lower = 1, upper = 1000000000;

Boolean allow_mismatches = False;
Boolean find_symmetry  = False;

float lowerg_se = 0.0, upperg_se = 100.0;
int lower_se = 1, upper_se = 1000000000;
int flank_length = 50;

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
    	
	char *hotspot_area[] = { "None", "ORF", "0-1kb", "1-2kb", "2-3kb", "3-4kb" };
	
    memset(status, 0, 1000);
	
	updateStatus(status, 0, seqP->seqLen);
    
    i=0;
    while (i < seqP->seqLen) {
		
        k = 0;
		
		// Counts for the length of purine and pyramidine ppts
        int purine = 0, pyramidine = 0;
		// Storage for type of ppt
        char base_string[20];
        // repeatCount is set to count of whatever was found (purine/pyramidine)
        int repeatCount = 0;
        // Start and end of sequence in position values (1-end)
        int start = i + 1, end = i+1;
        // Length of repeat
        int length;
    
    	// Count number of GC bases in ppt
        int countCG=0;
		
		// the minimum size between mismatches in the ppt
		int gap = error_gap;
		// error flag - indicates whether to keep searching
		int error = 0;
		// the next index to start searching from
		int nextI = i;
		
		// tally of mismatch types
		int mmA = 0, mmT = 0, mmC = 0, mmG = 0, mmO = 0;
		
		// tally of total mismatches in ppt
		int error_count = 0;
		// percentage value of errors - will always be <= 1/gap;
        float error_rate = 0.0;
		
		// Update the progress meter
        if (i%100 == 0) updateStatus(status, i, seqP->seqLen);
		
		// If the current base is unknown skip past it
		if (getBaseAt(seqP,i) != 'A' && getBaseAt(seqP,i) != 'G' &&
			getBaseAt(seqP,i) != 'C' && getBaseAt(seqP,i) != 'T')
		{
			i=i+1;
			continue;
		}
		
		// While error flag still indicates we should continue searching
		while (error == 0) {
			// Gap set to zero because no Purines found since last error(or start)
			gap=0;
			// while the base at position i is a purine
        	while (getBaseAt(seqP,i) == 'A' || getBaseAt(seqP,i) == 'G') {
				purine++; // Increase match length
				gap++;    // Increase gap size since start/error
				// Update CG counter if needed
				if (getBaseAt(seqP,i) == 'G') countCG++;
				// Move the position i forward one.
				i++;
        	}
        	
        	// If mismatches are allowed and we've found at least 2 purines then
			// the mismatch may be allowed...
			if (allow_mismatches && purine >= 2) {
				if (error_count  > 0) nextI = i - gap - 1;
				else nextI = i;
				
				if ((gap+1) < error_gap && error_count > 0) {
					error = 1;
					i--;
				} else {
					error = 0;
					error_count++;
					// ...We to include it as part of the match
					if (getBaseAt(seqP,i) == 'C') {
						countCG++; //update cg counter
						mmC++; // and tally the mismatch
					} else if (getBaseAt(seqP,i) == 'T') mmT++;
					else mmO++;
					// Also increase the length of the match.
					purine++;
					
					i++;
				}

			} else if (!allow_mismatches && purine) {
				nextI = i;
				i--; // Without mismatching the final increment of i needs to be undone.
				error = 1;
			} else {
				nextI = i;
				error = 1;
			}
		
        }
		
		// Get rid of mismatch in terminal two nucleotides
		if (purine >= 2 ) {
			// Is last nucleotide a mismatch?
			while (getBaseAt(seqP,i) != 'A' && getBaseAt(seqP,i) != 'G') {
				if (getBaseAt(seqP,i) == 'C') {
					countCG--;
					mmC--;
				} else if (getBaseAt(seqP,i) == 'T') mmT--;
				else mmO--;
				
				i--;
				error_count--;
				purine--;
			}
			
			// What about second to last nucleotide?
			if (getBaseAt(seqP,i-1) != 'A' && getBaseAt(seqP,i-1) != 'G') {
				error_count--;
				
				if (getBaseAt(seqP,i) == 'G') countCG--;
				if (getBaseAt(seqP,i-1) == 'C') {
					countCG--;
					mmC--;
				} else if (getBaseAt(seqP,i-1) == 'T') mmT--;
				else mmO--;
					
				i = i - 2;
				purine = purine - 2;
				//nextI = nextI - 2;
			}
			
		}

		// If enough purines were found then work out info		
        if (purine >= lower) {
			
			repeatCount = purine;
            end = i+1;
			length = end - start + 1;
			
            if (repeatCount >= lower && repeatCount <= upper){
				strcpy(base_string, "AG");
			}
		// otherwise try matching pyramidines
        } else {
			i = nextI;
			
	    	start = i+1;
			
            error_rate = 0.0; error_count = 0;
            countCG = 0;
			gap = error_gap;
			error = 0;
			
			mmA = 0, mmT = 0, mmC = 0, mmG = 0, mmO = 0;
			
			while (error == 0) {
				gap = 0;
            	while (getBaseAt(seqP,i) == 'C' || getBaseAt(seqP,i) == 'T') {
					pyramidine++;
					gap++;
                	if (getBaseAt(seqP,i) == 'C') countCG++;
                	i++;
					
            	}

				if (allow_mismatches && pyramidine >= 2) {
					if (error_count  > 0) nextI = i - gap - 1;
					else nextI = i;
	
					if ((gap+1) < error_gap && error_count > 0) {
						error = 1;
						i--;
					} else {
						error = 0;
						error_count++;
						// ...We to include it as part of the match
						if (getBaseAt(seqP,i) == 'G') {
							countCG++; //update cg counter
							mmG++; // and tally the mismatch
						} else if (getBaseAt(seqP,i) == 'A') mmA++;
						else mmO++;
						// Also increase the length of the match.
						pyramidine++;
						
						i++;
					}
	
				} else if (!allow_mismatches && pyramidine) {
					nextI = i;
					i--; // Without mismatching the final increment of i needs to be undone.
					error=1;
				} else {
					nextI = i;
					error = 1;
				}
                

            }
			
			if (pyramidine >= 2 ) {
				while (getBaseAt(seqP,i) != 'C' && getBaseAt(seqP,i) != 'T') {
					if (getBaseAt(seqP,i) == 'G') {
						countCG--;
						mmG--;
					} else if (getBaseAt(seqP,i) == 'A') mmA--;
					else mmO--;
					
					i--; // In case last nucleotide is a mismatch
					error_count--;
					pyramidine--;
					
				}
				
				
				if (getBaseAt(seqP,i-1) != 'C' && getBaseAt(seqP,i-1) != 'T') {
					error_count--;
					if (getBaseAt(seqP,i) == 'C') countCG--;
					if (getBaseAt(seqP,i-1) == 'G') {
						countCG--;
						mmG--;
					} else if (getBaseAt(seqP,i-1) == 'A') mmA--;
					else mmO--;
					
					i = i - 2;
					pyramidine = pyramidine - 2;
					//nextI = nextI - 2;
					
				}
			}
			
			end = i +1;
            repeatCount = pyramidine;
			length = end - start + 1;
			
            if (repeatCount >= lower && repeatCount <= upper){// && (((double)error_count-1.0)/(double)length) <= (1.0/(double)error_freq)) {
				strcpy(base_string, "CT");
			}
			
		}

		
		if (repeatCount >= lower && repeatCount <= upper) { //&& (((double)error_count-1.0)/(double)length) <= (1.0/(double)error_freq)) {
				float ratioGC = 0.0;
				
				char *se = NULL;
				int countGC_se = 0;
				float ratioGC_se = 0.0;
				
				char *sequence = NULL;
				sequence = malloc(sizeof(char) * (length + 1));
				baseNCopy(seqP, sequence, (start - 1), length);
				sequence[length] = '\0';
				
				if (find_symmetry)
				{
					int p;
					se = findSymmetry(sequence);
					
					for (p = 0; p < strlen(se); p++) {
						if (se[p] == 'G' || se[p] == 'C')
							countGC_se++;
					}
					ratioGC_se = (float) countGC_se / (float) strlen(se) * 100.0;
				}
				
				if (!find_symmetry || // if not looking for symmetry, or
					(find_symmetry &&
					strlen(se) > lower_se && strlen(se) < upper_se &&
					ratioGC_se > lowerg_se && ratioGC_se < upperg_se))
					// or we are looking at an SE that statisfies conditions
				{
					ratioGC = (float) countCG / (float) length * 100.0;
					if (ratioGC >= lowerg && ratioGC <= upperg) {
						char *hotspot_name = NULL;
						int intervalIndex = 0;
						
						fprintf(outf, "%s, %s, %d, %d, %d, %f, ", base_string, sequence, length, start, end, ratioGC);
						if (allow_mismatches)
							fprintf(outf, "%d, %d, %d, %d, %d, %d, ", error_count, mmG, mmA, mmC, mmT, mmO);
						else fprintf(outf, ",,,,,,");
						
						if (find_symmetry)
						{
							fprintf(outf, "%s, %d, %f, ", se, strlen(se), ratioGC_se);
						}
						else fprintf(outf, ",,,");
						
						// Output flanks
						{
							char* flank_left;
							char* flank_right;	
						
							flank_left = malloc(sizeof(char) * (flank_length + 1));
							flank_right = malloc(sizeof(char) * (flank_length + 1));
							
							baseNCopy(seqP, flank_left, start - flank_length - 1, flank_length);
							flank_left[flank_length]='\0';
							baseNCopy(seqP, flank_right, end, flank_length);
							flank_right[flank_length]='\0';
							
							fprintf(outf, "%s, %s,", flank_left, flank_right);
						}
											
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
						} else {
							fprintf(outf, ",,");
						}
						
						fprintf(outf, "\n");
					}
				}
			//if (!allow_mismatches)
			i++;
        } else {
        	i = nextI;
		}
    }
	
    updateStatus(status, 0, seqP->seqLen);
}
int main(int argc, char* argv[])
{
	SeqStructPtr seqP;
	FileStructPtr fsp;
	
	int count;
	int i;
	int chromosome = 0;
	
	clock_t startTime, timer;
	
	if (argc < 2)
    {
        fprintf(stderr,"Usage: %s <options> <fasta format sequence file names>\n\n", argv[0]);
        fprintf(stderr,"Options:\n\t-h <hotspot file>    - Sorts output into hotspots.\n");
        fprintf(stderr,"\t-l <lower>           - find at least <lower> repeats.\n");
		fprintf(stderr,"\t-u <upper>           - ignore more than <upper> repeats.\n");
        fprintf(stderr,"\t-c <chromosome>      - specify chromosome number (for hotspots)\n");
        fprintf(stderr,"\t-lg <lower GC\%%>      - specify lowest GC content to report\n");
        fprintf(stderr,"\t-ug <upper GC\%%>      - specify highest GC content to report\n");
		fprintf(stderr,"\t-e <error rate>      - the block size for mismatches (default = don't allow mismatches)\n");
		fprintf(stderr,"\t-f <size>            - number of flanking nucleotides to output.\n");
		fprintf(stderr,"\t-se                  - report symmetrical elements only (uses exact matching)\n");
		fprintf(stderr,"\t-selg <lower SE GC%%> - specify lowest SE GC content to report\n");
		fprintf(stderr,"\t-seug <upper SE GC%%> - specify highest SE GC content to report\n");
		fprintf(stderr,"\t-sel <lower>         - specify shortest SE length to report\n");
		fprintf(stderr,"\t-seu <upper>         - specify longest SE length to report\n");
		
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
				if (sscanf(argv[i], "%f", &error_gap) != 1) exit(2);
				allow_mismatches = True;
			}
			else if(strcmp(argv[i], "-f") == 0)
			{
				i++;
				if (sscanf(argv[i], "%d", &flank_length) != 1) exit(2);
			}
			else if(strcmp(argv[i], "-se") == 0)
			{
				find_symmetry = True;
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
	
	if (!fsp) {
		printf ("No fasta sequence loaded.\n");
		exit(1);
    }
	
	if (hotspot) {
		
		parseHotspots(hotspot_fd, chromosome);
		
	}
	
	initBuffer(fsp);
	
    outf = openOutFile("ppt-mismatch.csv");
	
    startTime = clock();
	
	count = 1;
	
	fprintf(outf, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n", "TYPE", "SEQUENCE", "LENGTH", "START POS", "END POS", "GC (%)", "MISMATCHES", "mmG", "mmA", "mmC", "mmT", "mmOther","SE", "SE LENGTH", "SE GC (%)", "LEFT FLANK", "RIGHT FLANK" , "HS", "HS REGION");
	
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

