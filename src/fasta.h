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

/* Common file for loading (virtually) infinite sized FASTA files.
 * (c) 2005 Joel Pitt, Lincoln University, joel.pitt@gmail.com
 */

#ifndef _REPEAT_FASTA_H_
#define _REPEAT_FASTA_H_
 
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>

/* trivial defs */
#ifndef True
#define True 1
#endif
#ifndef False
#define False 0
#endif

/* size of buffer for reads. */
#define BUF_SIZE 1024*10   /* 10K */
/* max size of description line (begins with ">") */
#define MAX_DESCRIPTION_LEN 1024
/* sequence block length */
#define SEQUENCE_BLOCK_LEN (1024 * 2000) 
/* max number of sequence chars dumped to line */
#define MAX_OUT_LINE_CHARS 60

typedef char Boolean;

/* Structure referring to a fasta file which may contain more than one sequence*/
typedef struct fs
{
	int fd;
	char buf[BUF_SIZE];
	Boolean endOfFile;
	int curBufLen;
	int curBufPos;
	int curSeekEnd;
	int curBufSeekStart;
	int newlines; // number of newlines encountered in buffer up to curBufPos
	Boolean usesCRLF;
	Boolean havePutBack;
	char putBack;
	char *fname; 
	int seq_index;
} FileStruct, *FileStructPtr;

/* struct for an individual sequence in a file */
typedef struct ss
{
	int index;
	char descStr[MAX_DESCRIPTION_LEN];
	char seqStr[SEQUENCE_BLOCK_LEN];
	int blockNumber;
	long blockFileEndPosition;
	long seqStart;
	unsigned int seqLen;
	FileStructPtr fp;
} SeqStruct, *SeqStructPtr;

/* Openfile initialises a file structure and opens the file */
FileStructPtr openFile(char* filename);
/* If the file structure is aligned to the start of a sequence,
 * I.e. the decription line, then returns a sequence structure */
SeqStructPtr newSeqStruct(FileStructPtr fps);
/* Align to the Nth sequence in the file, returns the sequence
 * number aligned to. If there are not N sequences in the file
 * it aligns and returns the last sequence */
int alignSequence(FileStructPtr fp, int n);
/* Once a sequence structure is retrieved you can access the
 * position of any nucleotide in it, getBaseAtPos starts at 1
 * while getBaseAt starts at 0 */
char getBaseAtPos(SeqStructPtr seqP, int pos);
char getBaseAt(SeqStructPtr seqP, int pos);

int baseNCompare(SeqStructPtr seqP, int pos, char* motif, int motifLength);
int baseNCopy(SeqStructPtr seqP, char* dest, int start, int length);

/* --- Internal functions ---- */
void fillBuf(FileStructPtr fp);
Boolean getChar(char *achar, FileStructPtr fsp);
void initBuffer(FileStructPtr fsp);
/* Adds character to line, line is sizeof MAX_DESCRIPTION_LENGTH */
void addCharToLine(char c, char *line, int *lineLen);
void putCharBack(char c, FileStructPtr fsp);
void getNonBlankLine(char *line, FileStructPtr fs);

int getSeq(SeqStructPtr seqPtr);

void dumpSeq(SeqStructPtr seqP);

#endif // _REPEAT_FASTA_H_
