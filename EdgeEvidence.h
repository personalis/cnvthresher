/*
 * EdgeEvidence.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_EDGEEVIDENCE_H_
#define BAMCNV_EDGEEVIDENCE_H_

#include "BamEvidence.h"

using namespace BamTools;
using namespace std;

class Fasta;

//! \class EdgeEvidence
//! \brief handler for edge-based evidence in support of a CNV detection
//! True CNVs will have either aligned soft-clip events, or aligned 
//! split-read events, at their breakpoints (or maybe both).  
//! This class traverses the reads in the vicinity of the CNV detection,
//! looking for these aligned soft-clip/split-read edges, at or near 
//! the called CNV breakpoints.
//! NOTE: this class also measures the fraction of multallelic loci 
//! in the CNV detection, and the fraction of reads with Qmap=0
class EdgeEvidence: public BamEvidence 
{
public:
    EdgeEvidence(size_t softLength, size_t insertsize, bool qmode, Fasta *f);
    ~EdgeEvidence();
    
    //! Traverse the aligned reads in a BAM file, gather data 
    //! to search for edges and other evidence. 
    //! \param reader reference to the BamReader object
    //! \param reg pointer to the GRegion object (the called CNV to be analyzed)
    //! \note this overrides a virtual fucntion in the base class
    void traverseReads(BamReader &reader, GRegion *reg);

    //! Examine the data gathered by traverseReads(), identify 
    //! aligned soft-clip/split-read edges
    //! \note this overrides a virtual fucntion in the base class
    void computeMetrics();
    //! Add findings to the GRegion object's attributes list
    //! \note this overrides a virtual fucntion in the base class
    void addAttributes();
    //! If soft-clip of split-read edges were found, update the pos1/pos2 
    //! data in the GRegion object.
    void updateBreakpoints();

protected:
    //! Reset internal data for the new CNV region
    //! \note this overrides a virtual fucntion in the base class
    void resetForRegion(GRegion *reg);
    //! delete allocated memeory (called by resetForRegion() an the destructor)
    //! \note this overrides a virtual fucntion in the base class
    void clearData();
    
private:
    //! Identify split-read edges
    void find_split_read_edges();
    //! Identify soft-clip egdes
    void find_soft_clip_edges();
    //! Count the number of non-consensus basecalls among aligned clipped portions of reads
    //! soft-clip edge identification requires that the clipped sequences be coherent.
    size_t count_AlignedBases_mismatches(vector<string> &AlignedClippedBases);
    //! Count the number of loci in the detected CNV that have a multiallelic arrangement of bases
    void count_biallelic_loci();
    //! Count the number of loci in the flanking regions that are >90% non-Reference allele
    void count_nonref_loci();

    //! static member variables which I am using like globals or #defines
    static int SoftClipBufferSize;
    static float FMIN;

    uint16_t *cov; // stores coverage across the target region
    uint16_t *nA;  // stores the number of A calls across the target region
    uint16_t *nC;  // stores the number of C calls across the target region
    uint16_t *nG;  // stores the number of G calls across the target region
    uint16_t *nT;  // stores the number of T calls across the target region
    uint16_t *nSoftClip;  // stores the number of soft-clip events across the target region
    uint16_t *nNonRef; // stores the number of non-Reference calls across the flanking regions

    // SoftClipRead stores the aligned position and clipped sequence 
    // of all soft-clipped reads encountered in the target region
    vector<pair<uint32_t, string> > SoftClipRead;
    // SplitRead stores the aligned positions before/after the split, 
    // of all split reads encountered in the target regions
    vector<pair<uint32_t, uint32_t> > SplitRead;

    double fracQ0, fracQ0flank, fracNonRef_flank;
    size_t nread, nbiall, nunbalanced;
    size_t srbp1, srbp2, scbp1, scbp2;
    double srf, scf1, scf2;

    Fasta *Ref;
};

#endif //BAMCNV_EDGEEVIDENCE_H_
