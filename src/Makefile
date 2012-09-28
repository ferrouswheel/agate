CFLAGS=-Wall
LDFLAGS=
CC=gcc

PRGMS: repeat-finder ppt-mismatch expected
# ppt-se not needed because ppt-mismatch does this.

default: all

all: ${PRGMS}

#ppt-se: fasta.o ppt-se.o common.o hotspot.o

repeat-finder: fasta.o common.o hotspot.o motif.o match.o database.o repeat-finder.o
#	gcc -o repeat-finder.exe c:/mingw/lib/gcrt1.o repeat-finder.o fasta.o match.o database.o motif.o hotspot.o common.o -Lc:/mingw/lib/ -lgmon

ppt-mismatch: fasta.o common.o ppt-mismatch.o hotspot.o

expected: expected.c
	${CC} ${CFLAGS} expected.c -lm -o expected

match.o: match.c match.h fasta.h

database.o: database.h database.c fasta.h

clean:
	rm -f *.o

#fasta.o: fasta.c fasta.h


#ppt-se.o: ppt-se.c
