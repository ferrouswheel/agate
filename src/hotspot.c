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

#include "hotspot.h"

Hotspot hotspots[MAX_HOTSPOTS];
int numHotspots = 0;

FILE* hotspot_fd;
Boolean hotspot = 0;

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
    //printf("chr=%s\n", chr);

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
								hotspots[numHotspots].counts[hotspots[numHotspots].numIntervals] = 0;

								readValue = True;
							} else {
								hotspots[numHotspots].intervals[hotspots[numHotspots].numIntervals][1] = value;
								hotspots[numHotspots].counts[hotspots[numHotspots].numIntervals] = 0;

								readValue = False;
							}
							if (readValue == False) {
								int l;
								for (l=0; l < 5; l++)
									hotspots[numHotspots].base_counts[hotspots[numHotspots].numIntervals][l] = 0;
								hotspots[numHotspots].numIntervals++;
								
							}
						}
					}

				}

				numHotspots++;



			}
			break;





        }

        }

	}
/*
	for (i=0; i < numHotspots; i++) {

		printf("%s:",hotspots[i].name);

		int j;
		for (j=0; j < hotspots[i].numIntervals; j++) {

			printf("%d-%d,",hotspots[i].intervals[j][0],hotspots[i].intervals[j][1]);

		}
		printf("\n");


	}
*/
}

