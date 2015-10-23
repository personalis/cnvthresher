/*
 * BamEvidence.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_BAMEVIDENCE_H_
#define BAMCNV_BAMEVIDENCE_H_

#include "api/BamReader.h"

using namespace BamTools;

class GRegion;

//! class BamEvidence
//! This is the base class for objects responsible for traversing reads in a BAM
//! in the vicinity of a CNV detection, and gathering different kinds of evidence 
//! to support the detection.  
class BamEvidence {
public:
    //! ctor
    //! Set member variables according to the arguments
    BamEvidence(size_t softLength, size_t insertsize, bool qmode);
    //! dtor (empty)
    ~BamEvidence();

    virtual void traverseReads(BamReader &reader, GRegion *reg) = 0;
    virtual void computeMetrics() = 0;
    virtual void addAttributes() = 0;

protected:
    virtual void resetForRegion(GRegion *reg) = 0;
    virtual void clearData() = 0;

    inline bool quietMode() const { return m_QuietMode; }

    bool   m_QuietMode;
    size_t m_StartPos;    //! starting position of the target region (over which reads will be traversed)
    size_t m_EndPos;      //! ending position of the target region (over which reads will be traversed)
    size_t m_TargetSize;  //! size of the target region
    size_t m_SoftLength;  //! softening length ~ expected uncertainty in breakpoints
    size_t m_SoftLength2; //! softening length truncated at half the region size
    size_t m_InsertSize;  //! characteristic insert size of the read pairs
    GRegion *m_Region;    //! pointer to the CNV detection
};

#endif //BAMCNV_BAMEVIDENCE_H_
