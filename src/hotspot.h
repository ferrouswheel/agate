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

#ifndef _REPEAT_HOTSPOT_H_
#define _REPEAT_HOTSPOT_H_

#include "fasta.h"

#define MAX_HOTSPOTS 10000
typedef struct HS
{
	char name[MAX_DESCRIPTION_LEN];
	int intervals[10][2];
	int numIntervals;
	int counts[10];
	int base_counts[10][5];

} Hotspot;

void parseHotspots(FILE* fd, int chromosome);

#endif // _REPEAT_HOTSPOT_H_
