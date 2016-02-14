/*! \brief cnvthresher: verify CNV detections using BAM evidence
  ! \author Jason Harris <jason.harris@personalis.com>
  ! \copyright 2015 Personalis, Inc.
  !
  ! Given a GFF file of CNV detections, cnvthresher examines the
  ! aligned reads in the vicinity of the CNV (from a BAM file)
  ! to look for corroborating evidence in three categories:
  !
  ! + Edge evidence:
  !   look for aligned split-reads (for smaller features)
  !   or soft-clip edges (for larger features)
  !
  ! + Insert-size evidence
  !   look for read pairs with anomalous large insert sizes, 
  !   consistent with the size and position of the CNV
  !
  ! + Morphological evidence
  !   measure several morph. attributes: edge sharpness, 
  !   std. dev. of coverage, frac. of loci with normalized 
  !   coverage in expected range
  !
  ! + other:
  !   Fraction of reads in the feature with mapping Q = 0
  !   Fraction of loci in the feature with multiallelic base
  !   distribution
*/
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>           //for sort()
#include <unistd.h>
#include <ctime>               //for timestamp
#include <cmath>               //for sqrt()
#include <iomanip>             //for std::setprecision()

#include "api/BamReader.h"     //BamTools
#include "FileUtil.h"          //File IO convenience
#include "GFF.h"               //encapsulates an annotated GFF
#include "GRegion.h"           //one region in a GFF
#include "EdgeEvidence.h"      //handler for edge-evidence
#include "InsSizeEvidence.h"   //handler for insert-size evidence
#include "MorphEvidence.h"     //handler for morphological evidence
#include "Fasta.h"

using namespace std;
using namespace BamTools;

static const char *versionString = "1.0.1";

//! Print the version string
static int version();

//! Print version and usage information
static int usage();

//! Load BAM for the named chromosome into the BamTools reader
void loadBAM(BamReader &reader, string chrName, string bampattern, bool QuietMode);

//! Get the median insert size from an aligned BAM
size_t getInsertSize(BamReader &reader);

//! Variables related to the reference sequence
Fasta *Ref = NULL;
RefVector Vref;
uint16_t bamchr = 65535; //index into the list of chromosomes, of the currently-loaded BAM

bool QuietMode = false;

//! The main program loop
int main(int argc, char *argv[])
{
    string bampattern = "";
    string gffpath = "";
    string outgffpath = "";
    string fastafile = "";
    size_t binsize = 1000;
    size_t insertsize = 0;
    bool computeMorphMetrics = false;
    bool keepOriginalBreakpoints = false;
    //Parse command-line arguments
    int c;
    while ((c = getopt(argc, argv, "hvqmkb:o:r:")) >= 0) {
        switch(c) {
        case 'b': binsize = atoi(optarg); break;
        case 'o': outgffpath = optarg; break;
        case 'r': fastafile = optarg; break;
        case 'q': QuietMode = true; break;
        case 'm': computeMorphMetrics = true; break;
        case 'k': keepOriginalBreakpoints = true; break;
        case 'v': return version(); break;
        case 'h': return usage(); break;
        }
    }

    //Too few arguments?  Show usage and exit.
    if (argc <= optind+1) return usage();

    //Print version string
    version();

    //Set the input BAM filename pattern, and the input GFF filename
    bampattern = argv[optind];

    gffpath = argv[optind+1];

    //Make sure specified GFF exists
    if (! fexists(gffpath.c_str()))
    {
        cerr << "Error: specified GFF file not found: " << gffpath << endl;
        exit(1);
    }

    //Set the output GFF filename, if it wasn't set with -o
    if (outgffpath.length() == 0)
    {
        outgffpath = gffpath.substr(0, gffpath.length()-4) + ".out.gff";
    }

    //Load the Reference FASTA file
    if (fexists(fastafile.c_str()))
    {
        if (! QuietMode) cerr << GRegion::timestamp() << " Read reference FASTA: " << fastafile << endl;
        Ref = new Fasta(fastafile);
    }
    else
    {
        if (fastafile.length() > 0)
        {
            cerr << "Warning: reference FASTA file not found: " << fastafile << endl;
        }
        else
        {
            cerr << "Warning: no reference FASTA file provided." << endl;
        }
        cerr << "    will not compute fNonRef metric." << endl;
    }

    //Load the chromosome 1 BAM file now, so we can extract Reference data
    BamReader reader;
    loadBAM(reader, "1", bampattern, QuietMode);
    if (! reader.IsOpen())
    {
        loadBAM(reader, "chr1", bampattern, QuietMode);
        if (! reader.IsOpen())
        {
            cerr << "Error: could not open BAM for chromosome 1." << endl;
            exit(1);
        }
    }
    bamchr = 0;
    Vref = reader.GetReferenceData();

    //the GRegion class has two static string-vector members (TypeName and ChrName)
    //that need to be initialized here, using Vref
    GRegion::initTypeName();
    for (RefVector::iterator it=Vref.begin(); it != Vref.end(); ++it)
    {
        GRegion::ChrName.push_back((*it).RefName);
    }

    //Traverse some read pairs in the BAM to get the mean insert size
    insertsize = getInsertSize(reader);
    if (! QuietMode) cerr << GRegion::timestamp() << " median insert size in BAM: " << insertsize << endl;

    //Parse the input GFF
    if (! QuietMode) cerr << GRegion::timestamp() << " Parsing GFF: " << gffpath << endl;
    GFF G(gffpath);
    if (! QuietMode) cerr << GRegion::timestamp() << "  read " << G.nRegions() << " CNVs." << endl;

    //softLength is how much slop we allow in the detected CNV breakpoints
    //set to twice the detection algorithm's binning size, but not less than 1 kbp
    size_t softLength = max((size_t)1000, (size_t)(3*binsize));

    //Initialize objects for evidence gathering
    EdgeEvidence edge(softLength, insertsize, QuietMode, Ref);
    InsSizeEvidence inssize(softLength, insertsize, binsize, QuietMode);
    MorphEvidence morph(softLength, insertsize, QuietMode);

    //Loop over each region in the input GFF
    for (vector<GRegion*>::iterator it = G.regionList().begin(); it != G.regionList().end(); ++it)
    {
        GRegion *reg = *it;

        //Skip non-CNV features
        if (reg->type() != Deletion && reg->type() != Duplication)
        {
            if (! QuietMode)
            {
                cerr << "Skipping non-CNV feature @ "
                     << reg->chrName() << ":"
                     << reg->pos1() << "-" << reg->pos2() << endl;
            }
            continue;
        }

        if (! QuietMode)
        {
            cerr << endl << GRegion::timestamp() << " CNV @ "
                 << reg->chrName() << ":"
                 << reg->pos1() << "-" << reg->pos2() << endl;
        }
        
        if (reg->size() > 1e7)
        {
            cerr << "cnvthresher cannot currently handle features larger than 10 Mbp; skipping" << endl;
            continue;
        }

        //Load a new BAM if this region is on a different chromosome
        if (bamchr != reg->chrID())
        {
            loadBAM(reader, reg->chrName(), bampattern, QuietMode);
            bamchr = reg->chrID();
        }

        // measure edge-based evidence
        edge.traverseReads(reader, reg);
        edge.computeMetrics();
        edge.updateBreakpoints();
        //recompute metrics with new BPs
        // edge.computeMetrics();
        edge.addAttributes();

        //measure insert-size-based evidence
        inssize.traverseReads(reader, reg);
        inssize.computeMetrics();
        inssize.addAttributes();

        //measure morphological evidence
        if (computeMorphMetrics)
        {
            morph.traverseReads(reader, reg);
            morph.computeMetrics();
            morph.addAttributes();
        }
    } //end loop over regions

    G.add_evidence_attributes();

    //Write out annotated GFF
    G.writeGFF(outgffpath, ! keepOriginalBreakpoints);

    delete Ref;
}

size_t getInsertSize(BamReader &reader)
{
    size_t MaxPairs = 10000;
    vector<size_t> ISList;
    BamAlignment r;

    while (ISList.size() < MaxPairs && reader.GetNextAlignmentCore(r))
    {
        if (r.RefID == r.MateRefID //mate mapped to same chr
            && r.Position < r.MatePosition //this is the upstream read of a properly-mapped pair
            )
        {
            ISList.push_back(r.MatePosition - r.Position);
        }
    }

    //Return the median insert size
    sort(ISList.begin(), ISList.end());
    return ISList[ (size_t)(ISList.size()/2) ];
}

void loadBAM(BamReader &reader, string chrName, string bampattern, bool QuietMode)
{
    //Determine whether the BAM pattern contains a wildcard ("*" or "%CHROM%")
    //If not, then we just have a single BAM file
    size_t WildcardIndex = bampattern.find("%CHROM%");
    if (WildcardIndex != string::npos)
    {
        bampattern.replace(WildcardIndex, 7, "*");
    }

    WildcardIndex = bampattern.find("*");

    //Close the previously opened BAM file, if necessary
    if (reader.IsOpen())
    {
        reader.Close();
    }

    //Open the BAM file 
    string bamfile;
    if (WildcardIndex != string::npos)
    {
        bamfile = bampattern.substr(0, WildcardIndex)
            + chrName
            + bampattern.substr(WildcardIndex+1, string::npos);
    }
    else
    {
        //We were only given one BAM file to work with
        bamfile = bampattern;
    }

    if (! QuietMode)
        cerr << GRegion::timestamp() << " Opening BAM file: " << bamfile << endl;

    if (! reader.Open(bamfile))
    {
        cerr << "Could not open BAM file: " << bamfile << endl;
        cerr << "Message was: \n" << reader.GetErrorString() << endl;
        return;
    }

    string baifile = bamfile + ".bai";
    if (! fexists(baifile.c_str()))
    {
        baifile = bamfile;
        baifile.replace(bamfile.length()-4, string::npos, ".bai");
        if (! fexists(baifile.c_str()))
        {
            cerr << "Error: could not find index for BAM file: " << bamfile << endl;
            if (reader.IsOpen())
            {
                reader.Close();
            }
            return;
        }
    }

    bool result = reader.OpenIndex(baifile);
    if (! result)
    {
        cerr << "Error loading BAM index.  Message was:\n" << reader.GetErrorString() << endl;
            if (reader.IsOpen())
            {
                reader.Close();
            }
        return;
    }
}

int version() {
    fprintf(stderr, "cnvthresher version: %s\n", versionString);
    return 1;
}

int usage() {
    version();
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  cnvthresher <options> BAMpattern GFFfile\n\n");
    fprintf(stderr, "Note: BAMpattern must contain a '*' in place of the chromosome name\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "    -o FILE   output annotated GFF file [default: adapted from input GFF; ends with .out.gff]\n");
    fprintf(stderr, "    -r FILE   FASTA file for the alignment reference\n");
    fprintf(stderr, "    -b INT    TigerCNV sampling size, in bp [default: 1000]\n");
    fprintf(stderr, "    -q        Quiet mode; suppress diagnostic console output\n");
    fprintf(stderr, "    -m        Calculate morphological annotations\n");
    fprintf(stderr, "    -k        Keep original breakpoints  [default: false; update when soft-clip/split-read edges found]\n");
    fprintf(stderr, "    -v        Print version string and exit\n");
    fprintf(stderr, "    -h        Print this help message and exit\n");
    return 1;
}

