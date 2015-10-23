/*
 * InsSizeEvidence.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_INSSIZEEVIDENCE_H_
#define BAMCNV_INSSIZEEVIDENCE_H_

#include "BamEvidence.h"

using namespace std;

/*! \class InsSizeEvidence
 *! \brief handler for insert-size evidence for a CNV call.
 *! CNVs imprint a signature of anomalous insert sizes on read pairs
 *! which happen to bracket the CNV breakpoints.  This class traverses
 *! reads in the upstream flank of the CNV, looking for reads whose 
 *! mates are mapped to the downstream flank.
 */

class InsSizeEvidence: public BamEvidence
{
public:
    InsSizeEvidence(size_t softLength, size_t insertsize, size_t binsize, bool qmode);
    ~InsSizeEvidence();

    void traverseReads(BamReader &reader, GRegion *reg);
    void computeMetrics();
    void addAttributes();

protected:
    void resetForRegion(GRegion *reg);
    void clearData();

private:
    vector<pair<uint32_t, uint32_t> > Del_ReadPair;
    vector<pair<uint32_t, uint32_t> > Dup_ReadPair;

    size_t m_MinInsertSize, m_MinInsertSize_Dup;
    size_t m_BinSize;
    size_t m_NumInsSize;
};

#endif //BAMCNV_INSSIZEEVIDENCE_H_
