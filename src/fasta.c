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

//#define FASTA_TEST
#ifdef FASTA_TEST
int main(int argc, char* argv[])
{
	SeqStructPtr seqP;
	FileStructPtr fps;
	char base;
	int i = 1;

	int pos;

	fps = openFile(argv[1]);
	while (alignSequence(fps, i) == i) {
		//printf("%d : \nAfter alignment pos = %d, seekstart=%d\n",i,lseek(fps->fd, 0, SEEK_CUR),fps->curBufSeekStart);
		seqP = newSeqStruct(fps);
		//printf("After new seq. pos = %d, seekstart=%d\n",lseek(fps->fd, 0, SEEK_CUR),fps->curBufSeekStart);
	
		fprintf(stderr,"%d\n",seqP->seqLen);
	
		for (pos = 1; pos <= seqP->seqLen; pos++)
		{
			base = getBaseAtPos(seqP, pos);
			printf("%c", base);
			if ((pos % 60) == 0)
			{
				putchar('\n');
				fflush(stdout);
			}
	
		}
		printf("\n");
	
//		pos = 117;
//		base = getBaseAtPos(seqP, pos);
		free(seqP);
		i++;
	}
	/*printf("Base at pos %d = %c\n", pos, base);
	dumpSeq(seqP);

	pos = 20;
	base = getBaseAtPos(seqP, pos);
	printf("Base at pos %d = %c\n", pos, base);
	dumpSeq(seqP);*/


	return 0;
}
#endif

FileStructPtr openFile(char *fn)
{
	FileStructPtr fps;
	int FILE;

	if (! (fps = (FileStructPtr )malloc(sizeof(FileStruct)) ) )
	{
		fprintf(stderr,"unable to malloc() memory for file structure.\n");
		exit(1);
	}
	/* clear mem */
	memset( (void *)fps, '\0', sizeof(FileStruct));

	/* open the specified file */
	FILE = open(fn, O_RDONLY);
	if (FILE == -1)
    {
        fprintf(stderr,"unable to open file %s\n", fn);
        exit(1);
    }

	fps->fd = FILE;
	fps->fname = strdup(fn);
	initBuffer(fps);
	
	return fps;
}

void fillBuf(FileStructPtr fp)
{
   size_t result;
   
   /* Save the file seek position for where the buffer starts. */
   fp->curBufSeekStart = lseek(fp->fd, 0, SEEK_CUR);
   /* Fill the buffer */
   result = read(fp->fd, (void *) fp->buf, BUF_SIZE);
   if (result == -1)
   {
        fprintf(stderr,"error reading file! errno = %d\n",errno);
        exit(1);
   }
   else if (result == 0) fp->endOfFile = True;
   else
    {
		// This is a hack - should really be using binary
		// reads for random access.
		// Checks whether lseek result equals character position
		// if they don't then files uses CRLF for newline.
		// Stupid windows.
		if (lseek(fp->fd, 0, SEEK_CUR) != (result+fp->curSeekEnd))
			fp->usesCRLF = True;
		
		/* Save buffer file position end */
		fp->curSeekEnd = lseek(fp->fd, 0, SEEK_CUR);
		/* Init buffer properties. */
		fp->curBufLen = result;
		fp->curBufPos = 0;
    }
}  /* readBuf */

/* returns True on success */
Boolean getChar(char *achar, FileStructPtr fsp)
{
   if (fsp->havePutBack)
    {
        *achar = fsp->putBack;
        fsp->havePutBack = False;
        if (*achar == '\n') fsp->newlines++;
        return(True);
    }

   if (fsp->curBufPos == fsp->curBufLen)
    fillBuf(fsp);

   if (fsp->endOfFile)
    return (False);

   *achar = fsp->buf[fsp->curBufPos++];
   if (*achar == '\n') fsp->newlines++;
   return (True);
}

/* should call this once for each file read */
void initBuffer(FileStructPtr fsp)
{
   /* initialize length and pointer */
   fsp->curBufPos = 0;
   fsp->curBufLen = 0;
   fsp->endOfFile = False;
   fsp->havePutBack = False;
   fsp->endOfFile = False;
   
   fsp->curBufSeekStart = 0;
   fsp->curSeekEnd = 0;
   fsp->newlines = 0;
   
}

void addCharToLine(char c, char *line, int *lineLen)
{
   if (*lineLen < MAX_DESCRIPTION_LEN)
    line[(*lineLen)++] = c;
   else
    fprintf(stderr,"warning: description line truncated\n");
}

/*
 * pick up a non-blank line from the file, presumably description.
 * truncates all leading blanks and/or blank lines
 */
void getNonBlankLine(char* line, FileStructPtr fs)
{
   Boolean stop, nonBlank;
   char c;
   int lineLen;
   
   lineLen = 0;
   stop = False;
   nonBlank = False;  /* will be set by any non whitespace char */
   while ((! fs->endOfFile) && (! stop))
   {
   if (getChar(&c, fs)) {
    if (c == '\n') {// || c == '\r' || c == '\f') 
        stop = nonBlank; /* stop if have anything. don't save eol char. */
   } else {
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
   }
   }
}

int alignSequence(FileStructPtr fp, int n)
{
	char c;
	Boolean lineStart = True;
	int count = 0;
	long filePos = 0;
	
	lseek(fp->fd, 0, SEEK_SET);
	initBuffer(fp);

	while (count < n && !fp->endOfFile)
	{
		getChar(&c, fp);
		
		if (lineStart && (c == '>')  )
		{
			count++;
			if (count == n) break;
		}
		filePos++;
		
		if (c == '\n') {
			lineStart = True;
			fp->newlines++;
			if (fp->usesCRLF) filePos++;
		}
		else lineStart = False;

	}
	lseek(fp->fd, filePos, SEEK_SET);
	initBuffer(fp);
	fp->seq_index = count;
	return count;

}

SeqStructPtr newSeqStruct(FileStructPtr fsp)
{
	SeqStructPtr newSeqP;

	/* malloc a new seq */
	if (! (newSeqP = (SeqStructPtr )malloc(sizeof(SeqStruct)) ) )
	{
		fprintf(stderr,"unable to malloc() memory for sequence.\n");
		exit(1);
	}
	/* clear mem */
	memset( (void *)newSeqP, '\0', sizeof(SeqStruct));

	/* pick up description line */
	getNonBlankLine(newSeqP->descStr, fsp);
	/*if (! getNonBlankLine(newSeqP->descStr, fsp) )
	{
		free(newSeqP);
		return ((SeqStructPtr )0);
	}*/

	/* did it start correctly ? */
	if (newSeqP->descStr[0] != '>')
	{
		fprintf(stderr,"format error in input file:  missing '>'\n");
		exit(1);
	}
	newSeqP->seqStart = lseek(fsp->fd,fsp->curBufSeekStart + fsp->curBufPos + fsp->newlines,SEEK_SET);
	newSeqP->blockFileEndPosition = newSeqP->seqStart;
	newSeqP->blockNumber = 0;
	newSeqP->seqLen = -1;
	newSeqP->fp = fsp;
	newSeqP->index = fsp->seq_index;
	getSeq(newSeqP);

	return newSeqP;

}

/* load the sequence struct with comment line and bases */
int getSeq(SeqStructPtr seqPtr)
{
   SeqStructPtr newSeqP;
   FileStructPtr fsp;
   Boolean endOfSeq;
   long blockStartPos = 0;
   long counter = 0;
   char c;

   newSeqP = seqPtr;
   fsp = seqPtr->fp;
   initBuffer(fsp);

   lseek(fsp->fd, newSeqP->seqStart, SEEK_SET);
   counter = newSeqP->seqStart;
   blockStartPos = newSeqP->seqStart + (SEQUENCE_BLOCK_LEN * newSeqP->blockNumber);
   if (newSeqP->seqLen != -1) {
	   if ((SEQUENCE_BLOCK_LEN * newSeqP->blockNumber) >= newSeqP->seqLen) {
		   fprintf(stderr, "Block beginning (%d) greater than sequence length (%d)\n", (SEQUENCE_BLOCK_LEN * newSeqP->blockNumber), newSeqP->seqLen);
		   return -1;
	   }
   }


   endOfSeq = False;
   while ((! fsp->endOfFile) && (!endOfSeq))
    {
		if (getChar(&c, fsp))
        {
            if (c == '>')
            {
                /* hit new sequence */
                endOfSeq = True;
                putCharBack(c,fsp);
            }
            else if (((c >= 'A') && (c <= 'Z')) ||
                    ((c >= 'a') && (c <= 'z')))/* bogus test, chris */
			{
				unsigned int position = 0;
				position = counter - blockStartPos;
				/* have nucleotide */
				if (position >= 0 && position < SEQUENCE_BLOCK_LEN)
					newSeqP->seqStr[position] = toupper(c);
				counter++;
			}
            else if ((c != '\n') && (c != '\r') && (c != ' ') && (c != '\t'))
            {
                /* wierd shit in file.  bail. */
                fprintf(stderr,"bad char in sequence, file %s: %d\n",fsp->fname,c);
                exit(1);
            }
        }

    }
	newSeqP->seqLen = counter - newSeqP->seqStart;
	newSeqP->blockFileEndPosition = blockStartPos + SEQUENCE_BLOCK_LEN;

   if (! newSeqP->seqLen)
    {
        fprintf(stderr,"? Null sequence encountered in file %s (ignored)\n",fsp->fname);
        fprintf(stderr,"  %s\n", newSeqP->descStr);
        return -1;

    }

   return 0;
}  /* getSeq */

char getBaseAtPos(SeqStructPtr seqP, int pos)
{
	int block;
	int blockPos;
	pos = pos - 1;
	block = pos / SEQUENCE_BLOCK_LEN;
	blockPos = pos % SEQUENCE_BLOCK_LEN;

	if (pos >= seqP->seqLen) {
		//fprintf(stderr, "Position (%d) greater than sequence length (%d)\n", pos + 1, seqP->seqLen);
		return -1;
	}

	if (block != seqP->blockNumber)
	{
		seqP->blockNumber = block;
		getSeq(seqP);
	}


	return seqP->seqStr[blockPos];
}

char getBaseAt(SeqStructPtr seqP, int pos)
{
	int block;
	int blockPos;
	pos = pos;
	block = pos / SEQUENCE_BLOCK_LEN;
	blockPos = pos % SEQUENCE_BLOCK_LEN;

	if (pos >= seqP->seqLen || pos < 0) {
		//fprintf(stderr, "Position (%d) greater than sequence length (%d)\n", pos + 1, seqP->seqLen);
		return 0;
	}

	if (block != seqP->blockNumber)
	{
		seqP->blockNumber = block;
		getSeq(seqP);
	}


	return seqP->seqStr[blockPos];
}

int baseNCompare(SeqStructPtr seqP, int pos, char* motif, int motifLength)
{
	int i;
	int mismatches = 0;
	for (i=0; i < motifLength; i++)
	{
		if (getBaseAt(seqP, pos + i) != motif[i] ) mismatches += 1;
	}
	return mismatches;
}

int baseNCopy(SeqStructPtr seqP, char* dest, int start, int length)
{
	int i;
	char c;
	for (i=0; i < length; i++)
	{
		
		c = getBaseAt(seqP,start+i);
		if (c == 0)
			dest[i] = '-';
		else dest[i] = c;
		
	}
	return 0;
}

void putCharBack(char c, FileStructPtr fsp)
{
   fsp->havePutBack = True;
   fsp->putBack = c;
}

void dumpSeq(SeqStructPtr seqP)
{
   int i, charsOnLine;

   fprintf(stdout,"%s\n", seqP->descStr);
   fprintf(stdout,"Sequence (length = %d):\n", seqP->seqLen);
   fprintf(stdout,"Sequence file pos: %d\n", seqP->seqStart);
   fprintf(stdout,"Block number = %d\n", seqP->blockNumber);
   fprintf(stdout,"Block file pos end = %d\n", seqP->blockFileEndPosition);
   i = 0;
   charsOnLine = 0;
   while (i < SEQUENCE_BLOCK_LEN)
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
