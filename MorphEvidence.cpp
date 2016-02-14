/*
 * MorphEvidence.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <sstream>
#include <iomanip>
#include <cmath> //sqrt()
#include "GRegion.h"
#include "MorphEvidence.h"

MorphEvidence::MorphEvidence(size_t softLength, size_t insertsize, bool qmode) :
    BamEvidence(softLength, insertsize, qmode)
{
    cov = NULL;
}

MorphEvidence::~MorphEvidence()
{
    clearData();
}

void MorphEvidence::resetForRegion(GRegion *reg)
{
    m_Region = reg;

    //Potentially modify m_softLength2 for small regions
    m_SoftLength2 = min(m_SoftLength, (size_t)(reg->size(true)/2));

    //For morphological anaysis, we traverse reads throughout the region,
    //including flanks that are each SoftLength in size
    if (reg->pos1(true) < m_SoftLength2)
    {
        m_StartPos = 0;
    }
    else
    {
        m_StartPos = reg->pos1(true) - m_SoftLength2;
    }
    m_EndPos     = reg->pos2(true) + m_SoftLength2;
    m_TargetSize = m_EndPos - m_StartPos + 1;

    clearData();

    cov = new uint16_t[m_TargetSize];
    memset(cov, 0, m_TargetSize * sizeof(uint16_t));
}

void MorphEvidence::traverseReads(BamReader &reader, GRegion *reg)
{
    BamAlignment r;

    resetForRegion(reg);

    if (! quietMode())
    {
        cerr << GRegion::timestamp() << " Traversing reads for morphological evidence: " << endl;
        cerr << "  Set target region: " << m_Region->chrName() << ":"
             << m_StartPos << "-" << m_EndPos << endl;
    }
    bool result = reader.SetRegion((int)m_Region->chrID(), (int)m_StartPos,
                     (int)m_Region->chrID(), (int)m_EndPos);
    if (! result)
    {
        cerr << "Failed to set region: " << (int)m_Region->chrID() << "  " << m_StartPos << "  " << m_EndPos << endl;
        cerr << "Error was: " << reader.GetErrorString() << endl;
        exit(1);
    }
    //Loop over all reads in the target region
    //Accumulate coverage over the target region,
    while (reader.GetNextAlignmentCore(r))
    {
        //uncomment to exclude reads with mapping quality zero
        //if (r.MapQuality == 0)
        //{
        //    continue;
        //}

        int kRead = 0; //position along read from start of read
        int kRef = r.Position - m_StartPos; //position along target region

        r.BuildCharData();
        for (size_t i=0; i < r.CigarData.size(); ++i)
        {
            uint32_t oplen = r.CigarData[i].Length;
            switch(r.CigarData[i].Type) {
            case 'M':
                for (size_t ii = 0; ii < oplen; ++ii)
                {
                    if (kRef>=0 && kRef<(int)m_TargetSize)
                    {
                        cov[kRef]++;
                    }
                    kRead++;
                    kRef++;

                    //Stop processing this read if it overuns the size of the target region
                    if (kRef >= (int)m_TargetSize) break;
                }
                break;
            case 'I':
                kRead += oplen;
                break;
            case 'S':
                kRead += oplen;
                break;
            case 'H':
                kRead += oplen;
                break;
            case 'D':
                kRef += oplen;
                break;
            case 'N':
                kRef += oplen;
                break;
            }
        }
    }
}

void MorphEvidence::refine_breakpoints()
{
    //If we already have refined breakpoints from EdgeEvidence, just return
    if (m_Region->pos1_refined() == true && m_Region->pos2_refined() == true)
    {
        m_RefinedStart = m_Region->pos1(true);
        m_RefinedStop = m_Region->pos2(true);
        return;
    }

    //We must use the shape of the coverage profile to identify the true breakpoints.
    //1: measure mean coverage inside the rough breakpoints
    size_t iStart = m_Region->pos1(false) - m_StartPos;
    size_t iStop  = m_Region->pos2(false) - m_StartPos;
    size_t nskip  = max((size_t)1, size_t(m_Region->size(true)/1000));
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_InFeature, m_SDCov_InFeature);
    
    //2: measure mean coverage in flanking regions, softLength away from the rough breakpoints
    size_t dOffset = m_SoftLength2/2;
    size_t FlankSize = min((size_t)500, (size_t)(m_Region->size(true)/4));

    if (m_Region->pos1(false) - m_StartPos < FlankSize + dOffset)
    {
        iStart = 0;
        iStop  = FlankSize;
    }
    else
    {
        iStop = m_Region->pos1(false) - m_StartPos - dOffset;
        iStart = iStop - FlankSize + 1;
    }
    nskip  = 1;
    m_MeanCov_Flank1 = 0.0;
    double SDCov_Flank1   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_Flank1, SDCov_Flank1);

    if (m_EndPos - m_Region->pos2(false) < FlankSize + dOffset)
    {
        iStop = m_TargetSize - 1;
        iStart = iStop - FlankSize + 1;
    }
    else
    {
        iStart = m_Region->pos2(false) - m_StartPos + dOffset;
        iStop  = iStart + FlankSize - 1;
    }
    m_MeanCov_Flank2 = 0.0;
    double SDCov_Flank2   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_Flank2, SDCov_Flank2);

    //3: Identify midpoint coverage level between interior mean and flank mean on each side
    double MidCov1 = 0.5*(m_MeanCov_InFeature + m_MeanCov_Flank1);
    double MidCov2 = 0.5*(m_MeanCov_InFeature + m_MeanCov_Flank2);

    //Simplify the math: if the feature is a duplication (interior mean cov > flank cov)
    //we will multiply the cov by -1 so we can always assume interior cov < flank cov
    double sgn1 = 1.0;
    double sgn2 = 1.0;
    if (m_MeanCov_InFeature > m_MeanCov_Flank1) sgn1 = -1.0;
    if (m_MeanCov_InFeature > m_MeanCov_Flank2) sgn2 = -1.0;
 
    //4: Binary search to identify position where coverage crosses the midpoint.
    //   a: in most cases, the search area starts at the rough breakpoint 
    //      unless the coverage at the rough breakpoint is already beyond the midpoint.
    //      In that case, we redefine the search area (start halfway between CNV midpoint and rough breakpoint,
    //      end at rough breakpoint)
    //Refine start breakpoint
    if (m_Region->pos1(false) - m_StartPos > dOffset)
    {
        iStart = m_Region->pos1(false) - m_StartPos - dOffset;
        iStop  = m_Region->pos1(false) - m_StartPos;    
    }
    else
    {
        iStart = 0;
        iStop = dOffset;
    }
    if (sgn1*cov[iStart] > sgn1*MidCov1)
    {
        //Reset start position to be halfway between rough breakpoint and midpoint of the feature
        iStart = 0.5*(m_TargetSize/2 + m_Region->pos1(false)-m_StartPos);
        //Reset end position to be the rough breakpoint
        iStop = m_Region->pos1(false) - m_StartPos;

        //If the coverage at this new start position is still beyond the midpoint, 
        //set start to the midpoint of the feature and stop to the current start
        if (sgn1*cov[iStart] > sgn1*MidCov1)
        {
            iStop = iStart;
            iStart = m_TargetSize/2;
        }
    }
    m_RefinedStart = binary_search(iStart, iStop, MidCov1) + m_StartPos;

    //Refine stop breakpoint
    if (m_Region->pos2(false) + dOffset >= m_EndPos)
    {
        iStop = m_EndPos - m_StartPos - 1;
        iStart = iStop - dOffset;
    }
    else
    {
        iStart = m_Region->pos2(false) - m_StartPos;
        iStop  = m_Region->pos2(false) - m_StartPos + dOffset;
    }
    if (sgn2*cov[iStart] > sgn2*MidCov2)
    {
        //Reset start position to be halfway between rough breakpoint and midpoint of the feature
        iStart = 0.5*(m_TargetSize/2 + m_Region->pos2(false)-m_StartPos);
        //Reset end position to be the rough breakpoint
        iStop = m_Region->pos2(false) - m_StartPos;

        //If the coverage at this new start position is still beyond the midpoint, 
        //set start to the midpoint of the feature and stop to the current start
        if (sgn2*cov[iStart] > sgn2*MidCov2)
        {
            iStop = iStart;
            iStart = m_TargetSize/2;
        }
    }
    m_RefinedStop = binary_search(iStart, iStop, MidCov2) + m_StartPos;
}

//Return the index of the cov vector (between iStart and iStop) where 
//the value crosses cov0
size_t MorphEvidence::binary_search(size_t iStart, size_t iStop, double cov0)
{
    if (cov[iStart] == cov0) return iStart;
    if (cov[iStop]  == cov0) return iStop;

    //swap iStart and iStop if iStart is larger
    if (iStart > iStop)
    {
        size_t tmp = iStart;
        iStart = iStop;
        iStop = tmp;
    }

    //set sgn to make sure that cov increases over the search region
    double sgn = 1.0;
    if (cov[iStart] > cov[iStop]) sgn = -1.0;

    size_t itry = (size_t)(0.5*(iStart + iStop));
    size_t dpos = itry;

    while (dpos > 0)
    {
        if (itry > iStop) itry = iStop;
        if (itry < iStart) itry = iStart;

        if (cov[itry] == cov0) break;

        if (sgn*cov[itry] < sgn*cov0)
        {
            itry += dpos;
        }
        else
        {
            itry -= dpos;
        }

        if (dpos > 1)
        {
            dpos = (size_t)((float)dpos/2 + 0.5);
        }
        else
        {
            dpos = 0;
        }
    }

    return itry;
}

void MorphEvidence::computeMetrics()
{
    //Before measuring morphological metrics, it is essential to have the 
    //best possible breakpoints.  Especially if the detection algorithm 
    //doesn't inherently detect exact breakpoints (e.g., TigerCNV).
    refine_breakpoints();

    //1: Measure mean coverage in the feature.  
    //   If the feature is large, we can save time by subsampling
    size_t iStart = m_RefinedStart - m_StartPos;
    size_t iStop  = m_RefinedStop - m_StartPos;
    size_t nskip  = max((size_t)1, size_t(m_Region->size(true)/1000));
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_InFeature, m_SDCov_InFeature);

    //2: Measure mean coverage in a large flanking region on each side
    //   Flanks are 1/4 the size of the feature, up to a max of 500 bp
    //   Flanks are offset from the region bps by m_softLength2
    size_t dOffset = m_SoftLength2 / 2;
    size_t FlankSize = min((size_t)500, (size_t)(m_Region->size(true)/4));
    if (m_RefinedStart - m_StartPos < FlankSize + dOffset)
    {
        iStart = 0;
        iStop  = FlankSize;
    }
    else
    {
        iStop = m_RefinedStart - m_StartPos - dOffset;
        iStart = iStop - FlankSize + 1;
    }
    nskip  = 1;
    m_MeanCov_Flank1 = 0.0;
    double SDCov_Flank1   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_Flank1, SDCov_Flank1);

    if (m_EndPos - m_RefinedStop < FlankSize + dOffset)
    {
        iStop = m_TargetSize - 1;
        iStart = iStop - FlankSize + 1;
    }
    else
    {
        iStart = m_RefinedStop - m_StartPos + dOffset;
        iStop  = iStart + FlankSize - 1;
    }
    m_MeanCov_Flank2 = 0.0;
    double SDCov_Flank2   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, m_MeanCov_Flank2, SDCov_Flank2);

    //3: Measure mean coverage in a small flanking region on each side
    //   Edge flanks are 1/20 the size of the region, up to 100 bp
    //   Flanks are adjacent to the region bps
    size_t EdgeSize = min((size_t)100, (size_t)(m_Region->size(true)/20));
    iStop  = m_RefinedStart - m_StartPos - 1;
    iStart = iStop - EdgeSize + 1;
    nskip  = 1;
    double MeanCov_OuterEdge1 = 0.0;
    double SDCov_OuterEdge1   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_OuterEdge1, SDCov_OuterEdge1);

    iStart = m_RefinedStop - m_StartPos + 1;
    iStop  = iStart + EdgeSize - 1;
    double MeanCov_OuterEdge2 = 0.0;
    double SDCov_OuterEdge2   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_OuterEdge2, SDCov_OuterEdge2);

    //4: Measure mean coverage in large regions just inside each breakpoint
    //   Regions are 1/4 region size, up to 500 bp
    //   Regions are offset from region bps by 100 bp
    iStart = m_RefinedStart - m_StartPos + dOffset;
    iStop  = iStart + FlankSize - 1;
    double MeanCov_InnerRegion1 = 0.0;
    double SDCov_InnerRegion1   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_InnerRegion1, SDCov_InnerRegion1);

    iStop  = m_RefinedStop - m_StartPos - dOffset;
    iStart = iStop - FlankSize + 1;
    double MeanCov_InnerRegion2 = 0.0;
    double SDCov_InnerRegion2   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_InnerRegion2, SDCov_InnerRegion2);

    //5: Measure mean coverage in small regions just inside each breakpoint
    //   Regions are 1/10 the size of the region, up to 100 bp
    //   Regions are adjacent to the region bps
    iStart = m_RefinedStart - m_StartPos + 1;
    iStop  = iStart + EdgeSize - 1;
    double MeanCov_InnerEdge1 = 0.0;
    double SDCov_InnerEdge1   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_InnerEdge1, SDCov_InnerEdge1);

    iStop  = m_RefinedStop - m_StartPos - 1;
    iStart = iStop - EdgeSize + 1;
    double MeanCov_InnerEdge2 = 0.0;
    double SDCov_InnerEdge2   = 0.0;
    mean_and_sd_coverage(iStart, iStop, nskip, MeanCov_InnerEdge2, SDCov_InnerEdge2);


    //Compute DeltaCov: the mean difference between InnerRegion and Flank coverage 
    //(average taken between two breakpoints, unless one of the flanks is zero (Gap))
    m_DeltaCov = 0.0;
    double MeanCov_Flank = 0.0;
    if (m_MeanCov_Flank1 > 0 && m_MeanCov_Flank2 > 0)
    {
        m_DeltaCov = 0.5*(MeanCov_InnerRegion1-m_MeanCov_Flank1 + MeanCov_InnerRegion2-m_MeanCov_Flank2);
        MeanCov_Flank = 0.5*(m_MeanCov_Flank1 + m_MeanCov_Flank2);
    }
    else if (m_MeanCov_Flank1 > 0)
    {
        m_DeltaCov = MeanCov_InnerRegion1-m_MeanCov_Flank1;
        MeanCov_Flank = m_MeanCov_Flank1;
    }
    else if (m_MeanCov_Flank2 > 0)
    {
        m_DeltaCov = MeanCov_InnerRegion2-m_MeanCov_Flank2;
        MeanCov_Flank = m_MeanCov_Flank2;
    }

    //Compute Edge_DeltaCov: the mean difference between InnerEdge and OuterEdge coverage 
    //(average taken between two breakpoints, unless one of the flanks is zero (Gap))
    double Edge_DeltaCov = 0.0;
    if (MeanCov_OuterEdge1 > 0 && MeanCov_OuterEdge2 > 0)
    {
        Edge_DeltaCov = 0.5*(MeanCov_InnerEdge1-MeanCov_OuterEdge1 + MeanCov_InnerEdge2-MeanCov_OuterEdge2);
    }
    else if (MeanCov_OuterEdge1 > 0)
    {
        Edge_DeltaCov = MeanCov_InnerEdge1-MeanCov_OuterEdge1;
    }
    else if (MeanCov_OuterEdge2 > 0)
    {
        Edge_DeltaCov = MeanCov_InnerRegion2-MeanCov_OuterEdge2;
    }
    
    //Compute EdgeSharpness: the fraction Edge_DeltaCov/m_DeltaCov
    //(basically: how much of the coverage delta happens very close to the bps)
    m_EdgeSharpness = max(0.0, min(1.0, Edge_DeltaCov/m_DeltaCov));

    //Compute FracInRange: 
    //Compute the fraction of coverage values in the feature that lie within +-0.25
    //of the expected level, given the zygosity
    //(for dups and no-zyg dels, select an expected coverage level based on the mean coverage in the feature)
    //nCov0 is the coverage in the feature as a fraction of the flank coverage (to the nearest 0.5)
    double nCov0 = 0.5 * (double)((int)(2.0 * m_MeanCov_InFeature / MeanCov_Flank));
    if (m_Region->type(true) == Deletion)
    {
        if (m_Region->zygosity() == Homozygous || m_Region->zygosity() == Hemizygous)
        {
            nCov0 = 0.0;
        }
        else if (m_Region->zygosity() == Heterozygous)
        {
            nCov0 = 0.5;
        }
    }
    else if (m_Region->type(true) == Duplication && m_Region->zygosity() == Heterozygous)
    {
        nCov0 = 1.5;
    }

    double LowerLimit = max(0.0, MeanCov_Flank * (nCov0 - 0.25));
    double UpperLimit = MeanCov_Flank * (nCov0 + 0.25);
    iStart = max((size_t)0, m_RefinedStart - m_StartPos);
    iStop  = min(m_TargetSize - 1, m_RefinedStop - m_StartPos);
    size_t nInRange = 0;
    size_t ntotal = 0;
    for (size_t i=iStart; i < iStop; ++i)
    {
        ++ntotal;
        if (cov[i] >= LowerLimit && cov[i] <= UpperLimit)
        {
            ++nInRange;
        }
    }
    m_FracInRange = (double)nInRange/(double)ntotal;
}

//void MorphEvidence::compute_mscore()
//{
//    m_MSCORE = 0.0;
//    double cov0 = (double)((size_t)(m_MeanCov_InFeature + 0.5)); // all this just to round()
//    if (m_Region->type(true) == Deletion)
//    {
//        if (m_Region->zygosity() == Homozygous || m_Region->zygosity() == Hemizygous)
//        {
//            cov0 = 0.0;
//        }
//        else if (m_Region->zygosity() == Heterozygous)
//        {
//            cov0 = 1.0;
//        }
//    }
//
//    double dcov0 = cov0 - 2.0;
//
//    // normalized coverage in the feature (want it to be near cov0)
//    double f = 1.0 - 2.0*abs(m_MeanCov_InFeature - cov0);
//    m_MSCORE += min(max(f, 0.0), 1.0);
//
//    // std dev of coverage (smaller is better)
//    f = 1.2 - m_SDCov_InFeature;
//    m_MSCORE += min(max(f, 0.0), 1.0);
//
//    // delta coverage (want it to be near dcov0)
//    f = 1.1 - abs(dcov0 - m_DeltaCov);
//    m_MSCORE += min(max(f, 0.0), 1.0);
//
//    // edge sharpness
//    f = m_EdgeSharpness;
//    m_MSCORE += min(max(f, 0.0), 1.0);
//
//    // fraction of coverage points in expected range
//    f = m_FracInRange;
//    m_MSCORE += min(max(f, 0.0), 1.0);
//}

void MorphEvidence::addAttributes()
{
    m_Region->set_mean_cov_in_feature(m_MeanCov_InFeature);
    m_Region->set_sdev_cov_in_feature(m_SDCov_InFeature);
    m_Region->set_mean_cov_in_flank1(m_MeanCov_Flank1);
    m_Region->set_mean_cov_in_flank2(m_MeanCov_Flank2);
    m_Region->set_edge_sharpness(m_EdgeSharpness);
    m_Region->set_edge_delta(m_DeltaCov);
    m_Region->set_frac_in_range(m_FracInRange);
    //    m_Region->set_mscore(m_MSCORE);
}

void MorphEvidence::mean_and_sd_coverage(size_t iStart, size_t iStop, size_t nskip, double &meanCov, double &sdCov)
{
    if (iStart < 0)
    {
        cerr << "Warning: cannot access coverage data prior to Start position." << endl;
        cerr << "    requested iStart: " << iStart << endl;
        cerr << "    resetting iStart to 0." << endl;
        iStart = 0;
    }

    if (iStop > m_TargetSize)
    {
        cerr << "Warning: cannot compute coverage statistics beyond size of the target region." << endl;
        cerr << "    requested iStop: " << iStop << "    target region size: " << m_TargetSize << endl;
        cerr << "    resetting iStop to be equal to m_TargetSize - 1." << endl;
        iStop = m_TargetSize - 1;
    }

    //To compute the mean and SD, we need to accumulate two sums:
    //  sum of coverage values
    //  sum of coverage values, squared
    double sumCov = 0.0;
    double sumCov2 = 0.0;
    size_t nn = 0;
    for (size_t i=iStart; i < iStop; i += nskip)
    {
        sumCov += (double)cov[i];
        sumCov2 += (double)cov[i] * (double)cov[i];
        ++nn;
    }

    meanCov = sumCov / (double)nn;
    sdCov = sqrt( (sumCov2/(double)nn - meanCov * meanCov) * (double)nn / (double)(nn-1) );
}

void MorphEvidence::clearData()
{
    if (cov != NULL)
    {
        delete[] cov;
        cov = NULL;
    }
}
