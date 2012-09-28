agate
=====

Optimal and non overlapping DNA repeat searcher

How to compile: See INSTALL file.

This repo contains source for four executables. They are probably not master
pieces of software engineering as I wrote them in my free time while doing
a PhD in an unrelated area, and the requirements changes substantially while
developing each tool:

- **repeat-finder**: Main program for finding optimal and exclusive repeating
  motifs using mismatches. Further description found below.
- **ppt-se**: Looks for poly-purine/poly-pyrimidine tracts (PPTs) with an
  internal axis of symmetry.  This includes, for example, AAAGAAA, GAGGAG. It
  has the additional commands "sel" and "seu" which specify the size limits of
  the symmetrical element.
- **ppt-mismatch**: Has the additional capability of finding PPT arrays with
  mismatches, the "e" command specified the block size within which only one
  mismatch is allowed. Mismatches are also disallowed within the terminal two
  nucleotides.
- **ppt**: Looks for PPTs eg AAAAAAAA, AGGGGA, GGGGG, CCCCCC, CCTCTTTCC,
  TTTTTT. "l" and "u" are size limiters for the overall tract size.
- **expected**: Calculates the expected number of repeats based on the
  equations in the paper: _de Wachter, Rupert (1981). The Number of Repeats
  Expected in Random Nucleic Acid Sequences and Found in Genes.  J. theor.
  Biol. 91, 71-98_

# Instructions for use

## Create a database of matches.

This is done with the "-d" option followed by a number indicating the maximum
size of motifs to look for. e.g.

    $ repeat-finder.exe -d 3 chr01.fsa

This will look for motifs of size 3 and below in the file "chr01.fsa".

This is the most time consuming part and once you have run it you shouldn't
have to run it again unless you want to change the maximum size of motifs or
you get a new version of the program.  The database filename is the same as the
sequence with ".database" appended. So for the above example, a database file
called "chr01.fsa.database" is created.

The database file *must* be the same as the sequence filename with ".database"
appended. So if you rename the sequence file, you must rename the database file
or regenerate it.

### Changing the error distance.

By default databases are created with a minimum error distance of 5.  If you
wish to alter this you can also specify an error distance when you create
a database using the "-e" option. e.g.

    $ repeat-finder.exe -d 3 -e 7 chr01.fsa

Sets the minimum error distance to be 7 nucleotides.

## Normal matching

Specify what matches you are interested in using the following:

```
-f <size>       - Output the <size> flanking bases
 to a match (default 50)
 
-n <number>     - search for motifs of size <number>.
You can use this option more than once, for example to
specify motifs of sizes 1 and 2 you can use "-n 2 -n 1"

-m <motif>      - search for a particular motif.
e.g "-m ATT". You can also use this more than once.

-l <lower>      - find matches with at least <lower> repeats.
-u <upper>      - ignore matches with more than <upper> repeats.

-lp <lower>     - lowest purity repeats to find
-up <upper>     - greatest purity repeats to find
(Note: range of purities is constrained by error distance)
```

Results are output to the file "results.csv". Be warned that this is
overwritten if it already exists so rename results you wish to keep.

## Compound Repeats

Use "-cr" to look for compound repeats, and "-gap" to specify the maximum gap
between repeats. If you only are interested in degenerate repeats then use
"-dr" as "-cr" looks for both compound and degenerate repeats.

You also must use the options for normal matches to specify which ones to
output. e.g.

    $ repeat-finder -cr -gap 4 -n 2 chr01.fsa

Looks for compound repeats with a gap of at most 4.  It then filters these for
only those that have at least one sub-repeat with a motif of size 2.

Results are output to the file "results.csv". The compound matches are reported
at the end, and also has a column to indicate the the indexes of the sub
repeats.

## Extracting to FASTA

To output to a fasta file use "-i" option along with the index of the match. If
the index is of a compound repeat you must also specify either -cr or -dr
(depending on what you used to get the index) and the same -gap number on the
command line.

You can use the index "-1" to indicate that all repeats should be output to
fasta.

The file created is by default called "repeats.fsa"

# The name

NOTE: This project used to be called repeatfinder and be [hosted on
sourceforge](http://repeatfinder.sf.net), but then they broke the drupal website
and I'm more fond of github. The reason for the name change is that there are
many other pieces of software called "repeatfinder" that do the same thing.
"Agate" is a [satellite/missle name](http://planet4589.org/space/misc/names.html)
and was the shortest name mostly composed of nucleotide letters ;-)

