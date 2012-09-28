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

#include <string.h>

#include "common.h"

void updateStatus(char* status, int value, int max)
{
	int stat_length, k = 0;
	
	stat_length = strlen(status);
	while (k < stat_length) {
		putchar('\b');
		k++;
	}
	sprintf(status, "%d / %d", value, max);

	printf("%s", status);
	fflush(stdout);
}

/* Find the frequency of how often a repeat sequence occurs,
   in a repeat or otherwise. */
int findFrequency(SeqStructPtr seqP, FILE* outf, char* repeat, int length)
{
    int frequency = 0;
    unsigned int i;

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

FILE* openOutFile(char *fn)
{
	FILE* outf;
	/* open the specified file */
	outf = (FILE*) fopen(fn, "w");
	if (outf == NULL)
	{
        fprintf(stderr,"unable to open file %s\n", fn);
        return NULL;
	}
	return outf;
}
