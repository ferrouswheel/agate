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

#ifndef _REPEAT_MOTIF_H_
#define _REPEAT_MOTIF_H_

#include "common.h"
#include "match.h"

#define DATABASE_MOTIFS 0 
#define FILTER_MOTIFS 1
#define MAX_NUM_SIZES 50

/*typedef struct po
{
	int start;
	int end;
	int freq;
} MatchStruct, *MatchStructPtr;*/

#define MAX_MOTIFS 100000

void makeMotifs(int n, int c, char* motif, int collection);
void freeMotifs(int collection);
Boolean check_motif(MatchPtr m);

#endif // _REPEAT_MOTIF_H_
