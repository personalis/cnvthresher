Summary Information:
--------------------------------------------------------------------------------
CNVthresher is (C) 2015 by Personalis, Inc.
Author: Jason Harris <jason.harris@personalis.com>

Github:

License: CNVthresher is released as open source software under a
modified BSD license. See the file LICENSE in the same folder as this
README file.


Brief description:
--------------------------------------------------------------------------------
CNVthresher annotates CNV calls defined in an input GFF file, by
examining the corresponding BAM file(s) for evidence in support of
each call.  


Usage:
--------------------------------------------------------------------------------
Usage:  cnvthresher <options> BAMpattern GFFfile

Note: BAMpattern must contain a '*' in place of the chromosome name 
Options: 
    -o FILE   output annotated GFF file [default: adapted from input GFF; ends with .out.gff]
    -r FILE   FASTA file for the alignment reference
    -b INT    TigerCNV sampling size, in bp [default: 1000]
    -q        Quiet mode; suppress diagnostic console output
    -m        Calculate morphological annotations
    -k        Keep original breakpoints  [default: false; update when soft-clip/split-read edges found]
    -v        Print version string and exit
    -h        Print this help message and exit

Usage notes:

cnvthresher currently assumes that the BAM files are split by
chromosome, that all BAM files are in the same folder, and that the
chromosome name is encoded in the BAM filename.  The BAMpattern
argument should be the full path to the BAM files, including the
filename pattern with the chromosome name replaced with the wildcard
'*' character.  Note that you will likely need to escape the wildcard
('\*'), or enclose the entire BAMpattern in qoutes, to avoid argument
expansion by the shell interpreter.

The alignment reference FASTA file is optional; if it is omitted, then
the "fNonRef" stat will not be computed.


Installation notes:
--------------------------------------------------------------------------------
CNVthresher requires bamtools.

To get a copy of bamtools:
git clone git://github.com/pezmaster31/bamtools.git

Then follow the build instructions here:
https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

Once you have installed bamtools, you need to modify the BAMTOOLS
variable in the Makefile (in the same folder as this README file) to
point to your bamtools installation.

Then simply type "make" to compile the code.


Definition of cnvthresher's annotation tags:
--------------------------------------------------------------------------------
scbp1/scbp2		refined breakpoints based on soft-clip edges
scf1/scf2		fraction of reads with the same soft-clip edge
srbp1/srbp2		refined breakpoints based on split-read edges
srf			fraction of reads with the same split-read edge
fbiall			fraction of interior loci with biallelic allele distribution
funbal			fraction of biallelic loci that are unbalanced
ncov			mean normalized read depth in the feature
sdcov			std. dev. of normalized read depth in the feature
dcov			delta normalized depth across the breakpoints
sharp			edge sharpness parameter, ranges from 0 to 1
fgood			fraction of interior loci with normalized depth consistent with CNV call
mscore			Morphology score (aggregate of previous 5 morph. metrics)
ninsize			Number of bracketing read pairs with anomalous insert size (consistent with the CNV call)
fQ0			fraction of interior reads with Qmap=0
fQ0flank		fraction of flanking reads with Qmap=0
fNonRef			fraction of flanking mapped bases which mismatch the reference

