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

#ifndef _REPEAT_COMMON_H_
#define _REPEAT_COMMON_H_

#include <stdio.h>

#include "fasta.h"

#define MEMORY_ERROR -7
#define FILE_ERROR -6

void updateStatus(char* status, int value, int max);
int findFrequency(SeqStructPtr seqP, FILE* outf,  char* repeat, int length);

FILE* openOutFile(char *fn);

#endif // _REPEAT_COMMON_H_
