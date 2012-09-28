/* #define DEBUG_SPUTNIK 1 */


/*
  polyAT

  find repeats in fasta format seq file
  allows for indels, returns score.

  beta version.  caveat emptor.

  chrisa  29-Jul-94

  chris abajian
  University of Washington
  Dept. of Molecular Biotechnology  FJ-20
  Fluke Hall, Mason Road
  Seattle WA 98195
*/

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>

/* trivial defs */
#ifndef True
#define True 1
#endif
#ifndef False
#define False 0
#endif

typedef int Boolean;

/* size of buffer for reads. */
#define BUF_SIZE 1024*10   /* 10K */
/* max size of description line (begins with ">") */
#define MAX_DESCRIPTION_LEN 1024
/* max sequence length */
#define MAX_SEQUENCE_LEN 1024*2000 /* 800K */
/* max number of sequence chars dumped to line */
#define MAX_OUT_LINE_CHARS 60

/* for debugging only */
#define MAX_ERRCODES 1024

/* search params and definitions */
#define MIN_UNIT_LENGTH 1 /* start search with dinucleotide repeats */
/* will search for di, tri, tetra ... <n>nucleotide repeats up to
   this value for n */
#define MAX_UNIT_LENGTH 15  /* up to and including pentanucleotides */
/* this is the point score for each exact match */
#define EXACT_MATCH_POINTS 1
/* this is the point score for a mismatch, insertion or deletion */
#define ERROR_MATCH_POINTS -6
/* this is the minimum score required to be considered a match */
#define MATCH_MIN_SCORE 8
/* this is the low score at which we stop trying */
#define MATCH_FAIL_SCORE -1
/* this is the max recursion depth we try to recover errors */
#define MAX_RECURSION 5


char *repeatName[MAX_UNIT_LENGTH+1] =
{
   "***ERROR***",    /* bad programmer!  no latte! */
   "single nucleotide",
   "dinucleotide",
   "trinucleotide",
   "tetranucleotide",
   "pentanucleotide",
   "6",
   "7",
   "8",
   "9",
   "10",
   "11",
   "12",
   "13",
   "14",
   "15"
};


char readBuf[BUF_SIZE];
Boolean endOfFile;
int curBufLen;
int curBufPos;

int fd = -1;
FILE* hotspot_fd;
Boolean hotspot;

FILE* outf;
Boolean havePutBack;
char putBack;
//void (*dumpMatch); //(SeqStructPtr, MatchStructPtr, boolean);
enum of {normal, tabular} outputFormat;

#define MAX_HOTSPOTS 10000
typedef struct HS
{
	char name[MAX_DESCRIPTION_LEN];
	int intervals[10][2];
	int numIntervals;
	int counts[10][4]; //counts of each nucleotide in sequence

} Hotspot;
Hotspot hotspots[10000];
int numHotspots = 0;

//void dumpMatchNormal(SeqStructPtr, MatchStructPtr, Boolean);
//dumpMatch = dumpMatchNormal;
/* struct for indiv sequence in a file */
typedef struct ss
{
   char descStr[MAX_DESCRIPTION_LEN];
   char seqStr[MAX_SEQUENCE_LEN];
   unsigned int seqLen;
} SeqStruct, *SeqStructPtr;

typedef struct po
{
	int start;
	int end;
	int freq;
} MatchStruct, *MatchStructPtr;

#define MAX_MOTIFS 10000
char* motifs[MAX_MOTIFS];
int halt_until[MAX_MOTIFS];
int counts[MAX_MOTIFS][4]; // counts for atcg in motif

int total_count[4]; // counts for sequence nucleotides

int num_motifs = 0;
// Combinatorial variables
int  comb_n;

char bases[] = { 'A', 'T', 'C', 'G' };

int lower = 1, upper = 1000000000;

void makeMotifs(int n, int c, char* motif) {

	if (c > 0) {
		motif[c - 1] = bases[0]; makeMotifs(n, c - 1, motif);
		motif[c - 1] = bases[1]; makeMotifs(n, c - 1, motif);
		motif[c - 1] = bases[2]; makeMotifs(n, c - 1, motif);
		motif[c - 1] = bases[3]; makeMotifs(n, c - 1, motif);
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
			motifs[num_motifs] = malloc(sizeof(char) * strlen(motif));
			strcpy(motifs[num_motifs], motif);

			memset(counts[num_motifs], 0, 4);

			for (j=0; j < strlen(motif); j++) {
				if (motif[j] == bases[0]) counts[num_motifs][0]++;
				else if (motif[j] == bases[1]) counts[num_motifs][1]++;
				else if (motif[j] == bases[2]) counts[num_motifs][2]++;
				else if (motif[j] == bases[3]) counts[num_motifs][3]++;
			}
			num_motifs++;
			//fprintf(stderr, "\t\tusing motif %s \n", motifs[num_motifs - 1]);

		} else fprintf(stderr, "skipping motif %s \n", motif);
	}
}


/*
 ************************************************************
 * these routines are used to read and parse the fasta format
 * sequence file
 ************************************************************
 */

void fillBuf(void)
{
   size_t result;

   result = read(fd, (void *)readBuf, BUF_SIZE);
   if (result == -1)
     {
        fprintf(stderr,"error reading file! errno = %d\n",errno);
        exit(1);
     }
   else if (result == 0)
     endOfFile = True;
   else
     {
        curBufLen = result;
        curBufPos = 0;
     }
}  /* readBuf */


/* returns True on success */
Boolean getChar(char *achar)
{
   if (havePutBack)
     {
        *achar = putBack;
        havePutBack = False;
        return(True);
     }

   if (curBufPos == curBufLen)
     fillBuf();

   if (endOfFile)
     return (False);

   *achar = readBuf[curBufPos++];
   return (True);
}


void putCharBack(char c)
{
   havePutBack = True;
   putBack = c;
}


int openFile(char *fn)
{
   int FILE;
   /* open the specified file */
   FILE = open(fn, O_RDONLY);
   if (FILE == -1)
     {
        fprintf(stderr,"unable to open file %s\n", fn);
        exit(1);
     }
	return FILE;
}

void openOutFile(char *fn)
{
   /* open the specified file */
   outf = (FILE*) fopen(fn, "w");
   if (outf == NULL)
     {
        fprintf(stderr,"unable to open file %s\n", fn);
        exit(1);
     }
}

/* should call this once for each file read */
void initBuffer(void)
{
   /* initialize length and pointer */
   curBufPos = 0;
   curBufLen = 0;
   havePutBack = False;
   endOfFile = False;
}

void addCharToLine(char c, char *line, int *lineLen)
{
   if (*lineLen < MAX_DESCRIPTION_LEN)
     line[(*lineLen)++] = c;
   else
     fprintf(stderr,"warning: description line truncated\n");
}


/*
 *********************************************************************
 * these routines are (more) specific to reading the fasta file format
 *********************************************************************
 */


/*
 * pick up a non-blank line from the file, presumably description.
 * truncates all leading blanks and/or blank lines
 */
Boolean getNonBlankLine(char *line)
{
   Boolean stop, nonBlank;
   char c;
   int lineLen;

   lineLen = 0;
   stop = False;
   nonBlank = False;  /* will be set by any non whitespace char */
   while ((! endOfFile) && (! stop))
     if (getChar(&c))
       if (c == '\n')
         stop = nonBlank; /* stop if have anything. don't save eol char. */
       else
         if (nonBlank)
           /* add it to line no matter what */
           addCharToLine(c,line,&lineLen);
         else if ((c != ' ') && (c != '\t'))
           {
              /* only non whitespace will start the line */
              nonBlank = True;
              addCharToLine(c,line,&lineLen);
           }
}


/* load the sequence struct with comment line and bases */
SeqStructPtr getSeq(char *fname)
{
   SeqStructPtr newSeqP;
   Boolean endOfSeq;
   char c;

   if (endOfFile) return ((SeqStructPtr )0);   /* bombproofing */

   /* malloc a new seq */
   if (! (newSeqP = (SeqStructPtr )malloc(sizeof(SeqStruct)) ) )
     {
        fprintf(stderr,"unable to malloc() memory for sequence.\n");
        exit(1);
     }
   /* clear mem */
   memset( (void *)newSeqP, '\0', sizeof(SeqStruct));

   /* pick up description line */
   if (! getNonBlankLine(newSeqP->descStr) )
     {
        free(newSeqP);
        return ((SeqStructPtr )0);
     }

   /* did it start correctly ? */
   if (newSeqP->descStr[0] != '>')
     {
        fprintf(stderr,"format error in input file:  missing '>'\n");
        exit(1);
     }

   endOfSeq = False;
   while ((!endOfFile) && (!endOfSeq))
     {
        if (getChar(&c))
          {
             if (c == '>')
               {
                  /* hit new sequence */
                  endOfSeq = True;
                  putCharBack(c);
               }
             else if (((c >= 'A') && (c <= 'Z')) ||
                      ((c >= 'a') && (c <= 'z')))/* bogus test, chris */
               /* have nucleotide */
               newSeqP->seqStr[newSeqP->seqLen++] = toupper(c);
             else if ((c != '\n') && (c != ' ') && (c != '\t'))
               {
                  /* wierd shit in file.  bail. */
                  fprintf(stderr,"bad char in sequence, file %s: %c\n",fname,c);
                  exit(1);
               }
          }
     }

   if (! newSeqP->seqLen)
     {
        fprintf(stderr,"? Null sequence encountered in file %s (ignored)\n",fname);
        fprintf(stderr,"  %s\n", newSeqP->descStr);
        free(newSeqP);
        return ((SeqStructPtr )0);
     }

   return(newSeqP);
}  /* getSeq */


/* for debugging.  dump entire seq to stdout. */
#ifdef DEBUG_SPUTNIK
void dumpSeq(SeqStructPtr seqP)
{
   int i, charsOnLine;

   fprintf(stdout,"%s\n", seqP->descStr);
   fprintf(stdout,"Sequence (length = %d):\n", seqP->seqLen);
   i = 0;
   charsOnLine = 0;
   while (i < seqP->seqLen)
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        fprintf(stdout,"%c", seqP->seqStr[i++]);
     }
   fprintf(stdout,"\n");
} /* dumpSeq */
#endif /* DEBUG_SPUTNIK */

void findRepeats(SeqStructPtr seqP, char** repeats, int numRepeats) {

	int i, j, k;
	int *sizes;
        char status[1000];
	int stat_length;

	char *hotspot_area[] = { "None", "ORF", "0-1kb", "1-2kb", "2-3kb", "3-4kb" };

	sizes = malloc(sizeof(int) * numRepeats);

	//memset(halt_until, 0, MAX_MOTIFS);
	memset(status, 0, 1000);
	memset(total_count, 0, 4);
	for (i=0; i<numHotspots; i++) {
        	int k;
                for (k=0; k < hotspots[i].numIntervals; k++)
			memset(hotspots[i].counts[k], 0, 4);
	}

	sprintf(status, "0 / %d", seqP->seqLen);
        printf("%s", status);

	for (i = 0; i < seqP->seqLen; i++) {

		k = 0;

        if (i%100 == 0) {
                      stat_length = strlen(status);
                	while (k < stat_length) {
                		putchar('\b');
                        	k++;
                	}
                        sprintf(status, "%d / %d", i, seqP->seqLen);
                        printf("%s", status);
                }




		if (seqP->seqStr[i] == bases[0]) total_count[0]++;
		else if (seqP->seqStr[i] == bases[1]) total_count[1]++;
		else if (seqP->seqStr[i] == bases[2]) total_count[2]++;
		else if (seqP->seqStr[i] == bases[3]) total_count[3]++;

		for (j=0; j < numHotspots; j++) {
				int k;

				for (k = hotspots[j].numIntervals-1; k >= 0; k--) {

					if (((i+1) > hotspots[j].intervals[k][0] && (i+1) < hotspots[j].intervals[k][1]))
					{

						if (seqP->seqStr[i] == bases[0]) hotspots[j].counts[k][0]++;
						else if (seqP->seqStr[i] == bases[1]) hotspots[j].counts[k][1]++;
						else if (seqP->seqStr[i] == bases[2]) hotspots[j].counts[k][2]++;
						else if (seqP->seqStr[i] == bases[3]) hotspots[j].counts[k][3]++;

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

/* Find the frequency of how often a repeat sequence occurs,
   in a repeat or otherwise. */
int findFrequency(SeqStructPtr seqP, char* repeat, int length)
{
	int frequency = 0;
	unsigned int i,j;

	for (i=0; i < seqP->seqLen; i++) {
		if (strncmp(seqP->seqStr + i, repeat, length) == 0)
		{
			fprintf(outf, "%c%c,%i,%i,%i\n", repeat[0], repeat[1], length/2, i, i + length - 1 );
			frequency++;
			i += length;
		}
		/*for (j=0; j < length; j++) {
		   if (repeat[j] != seqP->seqStr[i])
			   break;
		}
		if (j == length) frequency++; */

	}
	return frequency;
}

void parseHotspots(FILE* fd, int chromosome) {

#define MAX_HS_LINE 2000
#define MAX_TOKEN 2000


	char line[MAX_HS_LINE]; // Max line size 2000.
	Boolean CS = False;
	char chr[10];
	int i;

	char c;
	int locater = 0;

    sprintf(chr, "chr%d", chromosome);
    printf("chr=%s\n", chr);

	while (!feof(fd) && locater != -1) {
		char token[MAX_TOKEN];
		Boolean readValue = False;

		i = 0;
		memset(line, '\0', MAX_HS_LINE);
		memset(token, '\0', MAX_TOKEN);

		while (!feof(fd) && (c = fgetc(fd)) != '\n') {
			line[i] = c;
			i++;
		}

		i=0;
		while ((c = line[i]) != '\t' && c != '\n' && c != '\r') {
			token[i] = line[i];
			i++;
		}
		token[i] = '\0';

        //printf("token=%s\n", token);
        switch (locater) {
		case 0:
			if (strcmp(chr, token) == 0) {
				locater = 1;
				printf("found chromosome %d in hotspot file\n", chromosome);
			break;
		case 1:
		//	locater = 2;
		//	break; // Don't need chromosome length
		case 2: // Header HS/CS
			if (strcmp("HS", token) == 0) {
				printf("Found HS\n");
				CS = False;
				locater = 3;
			}
			else if (strcmp("CS", token) == 0) {
				printf("Found CS\n");
				CS = True;
				locater = 4;
			}
			else if (token[0] == '\0') {
				locater = -1;
			}
			break;
		case 3:
		case 4:
			if (token[0] == '\0') {
				locater = 2;
			}
			else {

				int j, k;
				char numberToken[MAX_TOKEN];

				strcpy(hotspots[numHotspots].name, token);

				if (CS)
					strcat(hotspots[numHotspots].name, " CS");

				hotspots[numHotspots].numIntervals = 0;

				j = strlen(token);
				while (j < strlen(line)) {

					memset(numberToken, '\0', MAX_TOKEN);

					for (k=0; (c = line[j]) != '\t' && c != '\0'; k++) {
						numberToken[k] = c;
						//putchar(c);

						j++;
					}
					j++;

					if (numberToken[0] != '\0'){
						int value;

						if (sscanf(numberToken, "%d", &value)) {
							//printf("numberToken=%s\n", numberToken);
							if (readValue == False) {
								hotspots[numHotspots].intervals[hotspots[numHotspots].numIntervals][0] = value;
								//hotspots[numHotspots].counts[hotspots[numHotspots].numIntervals] = 0;

								readValue = True;
							} else {
								hotspots[numHotspots].intervals[hotspots[numHotspots].numIntervals][1] = value;
								//hotspots[numHotspots].counts[hotspots[numHotspots].numIntervals] = 0;

								readValue = False;
							}
							if (readValue == False) hotspots[numHotspots].numIntervals++;
						}
					}

				}

				numHotspots++;



			}
			break;





        }

        }

	}

	for (i=0; i < numHotspots; i++) {

		printf("%s:",hotspots[i].name);

		int j;
		for (j=0; j < hotspots[i].numIntervals; j++) {

			printf("%d-%d,",hotspots[i].intervals[j][0],hotspots[i].intervals[j][1]);

		}
		printf("\n");


	}

}

double calcExpected(int motif[4], int p[4], int N, int n, int t) {

	double expected = 0;

		double mult;
		double n_dash;

		mult = pow (p[0]/(double)N, motif[0]) * pow (p[1]/(double)N, motif[1]) *
			pow (p[2]/(double)N, motif[2]) * pow (p[3]/(double)N, motif[3]);

		n_dash = N - (t * n) - (2 * n) + 1;

		expected = pow(mult, t) * (1 - mult) * (n_dash * (1 - mult) + (2*n));



	return expected;

}

int main(int argc, char* argv[])
{
   SeqStructPtr seqP;
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

	 upper =10;

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

			memset(counts[num_motifs], 0, 4);

			for (j=0; j < strlen(argv[i]); j++) {
				if (argv[i][j] == bases[0]) counts[num_motifs][0]++;
				else if (argv[i][j] == bases[1]) counts[num_motifs][1]++;
				else if (argv[i][j] == bases[2]) counts[num_motifs][2]++;
				else if (argv[i][j] == bases[3]) counts[num_motifs][3]++;
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
		fd = openFile(argv[i]);
	}
   }


	fprintf(stdout,"Number of motifs to check = %d\n", num_motifs);

   if (fd < 0) {
   	printf ("No fasta sequence loaded.\n");
	exit(1);
	}

   if (hotspot) {

	parseHotspots(hotspot_fd, chromosome);

   }



   initBuffer();

	openOutFile("expected.csv");

    startTime = clock();

   count = 0;

   while (! endOfFile)
     if ((seqP = getSeq(argv[1])))
       {
		int j;

        timer = clock();

        fprintf(stdout,"processing sequence %d\n", count++);
		fprintf(outf,"Sequence %d\n", count);

		strcpy(hotspots[numHotspots].name, "Whole chromosome");
		hotspots[numHotspots].numIntervals = 1;
		hotspots[numHotspots].intervals[0][0] = 0;
		hotspots[numHotspots].intervals[0][1] = seqP->seqLen;
		memset(hotspots[numHotspots].counts[0], '\0', 4);
		numHotspots++;

		findRepeats(seqP, motifs, num_motifs);



		if (hotspot) {
		if (num_motifs < 100) {
                fprintf(outf, "Motifs:, ");

                for (i=0; i < num_motifs; i++) {
			fprintf(outf, "%s, ", motifs[i]);
		}
                }

		fprintf(outf, "\nTotal Expected (copynumbers %d to %d) Counts in each hotspot region:,, Expected frequency in region\n", lower, upper);
		fprintf(outf, "Hotspot, ORF, 0-1kb, 1-2kb, 2-3kb, 3-4kb\n");



		for (j=0; j < numHotspots; j++) {


			int k;

			//fprintf(outf,"%d, ", l);
			fprintf(outf,"%s, ", hotspots[j].name);

			for (k = 0; k < hotspots[j].numIntervals; k++) {
				double total;
                                total = 0;
				for (i=0; i < num_motifs; i++) {
					int l;
					for (l=lower; l <= upper; l++) {
						total += calcExpected(counts[i], hotspots[j].counts[k],
							hotspots[j].intervals[k][1] - hotspots[j].intervals[k][0], strlen(motifs[i]), l);
					}
				}
				fprintf(outf,"%.25f, ", total);
			}
			fprintf(outf,"\n");
		}



		fprintf(outf, "\nBase counts in each hotspot region:\n");
		fprintf(outf, "Base, Hotspot, ORF, 0-1kb, 1-2kb, 2-3kb, 3-4kb\n");
		for (j=0; j < numHotspots; j++) {
			int l;
			for (l=0; l < 4; l++) {
			int k;

			fprintf(outf,"%c, %s, ", bases[l], hotspots[j].name);
			for (k = 0; k < hotspots[j].numIntervals; k++) {
					fprintf(outf, "%d, ", hotspots[j].counts[k][l]);

				}
			fprintf(outf, "\n");


			}

		}

		}

		        free((void *)seqP);

		fprintf(stdout,"...took %4.2f seconds.\n", (float) (clock() - timer) / CLOCKS_PER_SEC);
       }
   fprintf(stdout, "Total time elapsed: %4.2f seconds.\n", (float) (clock() - startTime) / CLOCKS_PER_SEC);

	close((int)outf);

	return 0;
}








