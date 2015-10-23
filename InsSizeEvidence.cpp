/*
 * InsSizeEvidence.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <sstream>
#include "GRegion.h"
#include "InsSizeEvidence.h"

InsSizeEvidence::InsSizeEvidence(size_t softLength, size_t insertsize, size_t binsize, bool qmode) :
    BamEvidence::BamEvidence(softLength, insertsize, qmode)
{
    m_BinSize = binsize;
}

InsSizeEvidence::~InsSizeEvidence()
{
    clearData();
}

void InsSizeEvidence::resetForRegion(GRegion *reg)
{
    m_Region = reg;
    
    //Potentially modify m_softLength2 for small regions
    m_SoftLength2 = min(m_SoftLength, (size_t)(reg->size(true)/2));

    //For insert-size analysis, we only traverse reads in a +-softLength neighborhood around pos1
    //"true" argument to pos1()/pos2(): use positions updated by EdgeEvidence, if available
    //m_EndPos cannot exceed the midpoint of the region
    m_StartPos   = reg->pos1(true) - m_SoftLength2 - 2*m_InsertSize;
    m_EndPos     = min(reg->pos1(true) + m_SoftLength2 + 2*m_InsertSize, reg->pos1(true) + (size_t)(reg->size(true)/2));
    m_TargetSize = m_EndPos - m_StartPos + 1;

    //Set the minimum insert size to be considered consistent with the detected feature,
    //and significantly larger than the insert-size
    m_MinInsertSize = 0;
    if (m_Region->size(true) > m_InsertSize + 2 * m_BinSize)
    {
        m_MinInsertSize = m_Region->size(true) + m_InsertSize - 2 * m_BinSize;
    }

    m_MinInsertSize_Dup = 0;
    if (m_Region->size(true) > 2 * m_InsertSize + 2 * m_BinSize)
    {
        m_MinInsertSize_Dup = m_Region->size(true) - m_InsertSize - 2 * m_BinSize;
    }

    m_NumInsSize = 0;

    clearData();
}

void InsSizeEvidence::clearData()
{
    //Clear lists
    Del_ReadPair.clear();
    Dup_ReadPair.clear();
}

void InsSizeEvidence::traverseReads(BamReader &reader, GRegion *reg)
{
    BamAlignment r;

    resetForRegion(reg);

    //DEBUG
    if (! quietMode()) 
    {
        cerr << "Minimum insert size: " << m_MinInsertSize << "   " << m_MinInsertSize_Dup << endl;
    }

    if (m_MinInsertSize == 0 && m_MinInsertSize_Dup == 0)
    {
        if (! quietMode())
        {
            cerr << GRegion::timestamp() << " skipping insert-size analysis (feature is too small)" << endl;
        }
        return;
    }

    if (! quietMode())
    {
        cerr << GRegion::timestamp() << " Traversing reads for anomalous insert sizes: " << endl;
        cerr << "  Set target region: " << reg->chrName() << ":"
             << m_StartPos << "-" << m_EndPos << endl;
    }
    
    bool result = reader.SetRegion((int)m_Region->chrID(), (int)m_StartPos, 
                     (int)m_Region->chrID(), (int)m_EndPos);
    if (! result)
    {
        cerr << "Failed to set region: " << m_Region->chrID() << "  " << m_StartPos << "  " << m_EndPos << endl;
        cerr << "Error was: " << reader.GetErrorString() << endl; 
        exit(1);
    }

    while (reader.GetNextAlignmentCore(r))
    {
        //Deletions: look for large insert sizes among read pairs with:
        //  upstream read just outside feature's upstream bp
        //  downstream read just outside feature's downstream bp
        //  downstream read is within softLength of the downstream bp (not *too* far downstream)
        //  pair's insert size exceeds MinInsertSize
        if (m_Region->type(true) == Deletion 
            && m_MinInsertSize > 0 //Feature is big enough that we can detect anomalous insert sizes
            && r.RefID == r.MateRefID //mapped to same chr
            && r.Position  < r.MatePosition  //This is the upstream read
            && r.Position  < (int)m_Region->pos1(true)  //This read is upstream of feature
            && r.MatePosition > (int)m_Region->pos2(true) //Mate is downstream of feature
            && r.MatePosition <= (int)(m_Region->pos2(true) + m_SoftLength2 + m_InsertSize) //Mate is not TOO far downstream
            && r.MatePosition - r.Position >= (int)m_MinInsertSize //large insert size
            )
        {
            Del_ReadPair.push_back(make_pair(r.Position, r.MatePosition));
            
            if (! quietMode())
            {
                cerr << "  Large ins-size (Del):  " << r.Position << "   " << r.MatePosition 
                     << "  " << r.MatePosition - r.Position 
                     << endl;
            }
        }

        //Duplications: look for large insert sizes among read pairs with:
        //  upstream read just inside feature's upstream bp
        //  downstream read just inside feature's downstream bp
        //  pair's insert size exceeds MinInsertSize
        if (m_Region->type(true) == Duplication 
            && m_MinInsertSize_Dup > 0 //Feature is big enough that we can detect anomalous insert sizes
            && r.RefID == r.MateRefID //mapped to same chr
            && r.Position  < r.MatePosition  //This is the upstream read
            && r.Position  >= (int)m_Region->pos1(true) - 10 //This read is inside of feature
            && r.MatePosition <= (int)m_Region->pos2(true) + 10 //Mate is inside of feature
            && r.MatePosition - r.Position >= (int)m_MinInsertSize_Dup //large insert size
            )
        {
            Dup_ReadPair.push_back(make_pair(r.Position, r.MatePosition));
            
            if (! quietMode())
            {
                cerr << "  Large ins-size (Dup):  " << r.Position << "   " << r.MatePosition
                     << "  " << r.MatePosition - r.Position 
                     << endl;
            }
        }
    }
}

void InsSizeEvidence::computeMetrics()
{
    m_NumInsSize = max( Del_ReadPair.size(), Dup_ReadPair.size() );
}

void InsSizeEvidence::addAttributes()
{
    m_Region->set_n_insize(m_NumInsSize);
}
