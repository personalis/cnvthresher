/*
 * EdgeEvidence.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "EdgeEvidence.h"
#include "GRegion.h"
#include "Fasta.h"

int   EdgeEvidence::SoftClipBufferSize = 20;
float EdgeEvidence::FMIN = 0.80;

EdgeEvidence::EdgeEvidence(size_t softLength, size_t insertsize, bool qmode, Fasta *f) : 
    BamEvidence(softLength, insertsize, qmode)
{
    cov = NULL;
    nA = NULL;
    nC = NULL;
    nG = NULL;
    nT = NULL;
    nSoftClip = NULL;
    nNonRef = NULL;
    Ref = f;
}

EdgeEvidence::~EdgeEvidence()
{
    clearData();
}

void EdgeEvidence::resetForRegion(GRegion *reg)
{
    m_Region = reg;

    //For Edge-evidence analysis, we traverse reads throughout the region, 
    //including flanks that are each softLength in size
    //"true" argument to pos1()/pos2() means use bp-improved positions, if available
    m_StartPos   = reg->pos1(true) - m_SoftLength - 2*m_InsertSize;
    m_EndPos     = reg->pos2(true) + m_SoftLength + 2*m_InsertSize;
    m_TargetSize = m_EndPos - m_StartPos + 1;

    //Potentially modify m_softLength2 for small regions
    m_SoftLength2 = min(m_SoftLength, (size_t)(reg->size(true)/2));

    scbp1 = 0;
    scbp2 = 0;
    srbp1 = 0;
    srbp2 = 0;
    scf1 = 0.0;
    scf2 = 0.0;
    srf = 0.0;

    nread = 0; 
    nbiall = 0; 
    nunbalanced = 0;

    clearData();

    if (! quietMode())
    {
        cerr << m_StartPos << "  " << m_EndPos << " : " << m_TargetSize << endl;
    }

    cov = new uint16_t[m_TargetSize];
    nA  = new uint16_t[m_TargetSize];
    nC  = new uint16_t[m_TargetSize];
    nG  = new uint16_t[m_TargetSize];
    nT  = new uint16_t[m_TargetSize];
    nSoftClip = new uint16_t[m_TargetSize];
    nNonRef = new uint16_t[m_TargetSize];

    size_t nbytes = m_TargetSize * sizeof(uint16_t);
    memset(cov, 0, nbytes);
    memset(nA, 0, nbytes);
    memset(nC, 0, nbytes);
    memset(nG, 0, nbytes);
    memset(nT, 0, nbytes);
    memset(nSoftClip, 0, nbytes);
    memset(nNonRef, 0, nbytes);
}

void EdgeEvidence::traverseReads(BamReader &reader, GRegion *reg)
{
    BamAlignment r;

    resetForRegion(reg);

    if (! quietMode())
    {
        cerr << GRegion::timestamp() << " Traversing reads for BAM evidence: " << endl;
        cerr << "    Target region set to: " << m_Region->chrName() << ":"
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

    //Count number of reads in entire target area (region plus flanks)
    //Also count number of reads inside the feature region, and
    //the subset with mapQ==0
    nread = 0;
    fracQ0 = 0.0;
    fracQ0flank = 0.0;
    size_t nread_in_feature = 0;
    size_t nread_Q0 = 0;
    size_t nread_in_flanks = 0;
    size_t nread_Q0_flanks = 0;

    //loop over all reads in the target region
    while (reader.GetNextAlignmentCore(r))
    {
        //Count number of reads in the region, number of reads in the feature, and number of reads 
        //in the feature with mapping quality = 0
        nread++;
        if (r.Position >= (int)m_Region->pos1(true) && r.GetEndPosition() <= (int)m_Region->pos2(true))
        {
            nread_in_feature++;
            if (r.MapQuality == 0)
            {
                nread_Q0++;
            }
        }
        else
        {
            nread_in_flanks++;
            if (r.MapQuality == 0)
            {
                nread_Q0_flanks++;
            }
        }

        //Parse the CIGAR string; accumulate coverage for mapped positions
        int kRead = 0; //position along read from start of read
        int kRef = r.Position - m_StartPos; //position along target region
        //(note: kRef may be negative, if the read only partially overlaps)

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
                        if (r.QueryBases[kRead] == 'A') ++nA[kRef];
                        else if (r.QueryBases[kRead] == 'C') ++nC[kRef];
                        else if (r.QueryBases[kRead] == 'G') ++nG[kRef];
                        else if (r.QueryBases[kRead] == 'T') ++nT[kRef];

                        if (Ref != NULL)
                        {
                            if (r.QueryBases[kRead] != Ref->baseAt(m_Region->chrID(), kRef + m_StartPos)) ++nNonRef[kRef];
                        }
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
                if (kRef>=0 && kRef<(int)m_TargetSize)
                {
                    nSoftClip[kRef]++;
                    size_t nbases = min((int)oplen, EdgeEvidence::SoftClipBufferSize);
                    //If the CIGAR string started with S, then
                    //we want to reverse the clipped sequence
                    //string before storing it
                    string seq = r.QueryBases.substr(kRead, nbases);
                    if (i==0) //CIGAR started with S
                    {
                        string s = r.QueryBases.substr(kRead, oplen);
                        string srev = string(s.rbegin(), s.rend());
                        seq = srev.substr(0, nbases);
                    }
                    SoftClipRead.push_back(make_pair(kRef, seq));
                }
                kRead += oplen;
                break;
            case 'H':
                kRead += oplen;
                break;
            case 'D':
                //TODO: make these lists, and insert in Pos1 sorted order
                if (kRef + m_StartPos <= m_Region->pos1(true) + m_SoftLength2 &&
                    kRef + oplen + m_StartPos >= m_Region->pos2(true) - m_SoftLength2)
                {
                    SplitRead.push_back(make_pair(kRef + m_StartPos, kRef + oplen + m_StartPos));
                }
                kRef += oplen;
                break;
            case 'N':
                kRef += oplen;
                break;
            }
        }
    } //end of loop over all reads in target region

    if (nread_in_feature > 0)
    {
        fracQ0 = (double)nread_Q0/(double)nread_in_feature;
    }

    if (nread_in_flanks > 0)
    {
        fracQ0flank = (double)nread_Q0_flanks/(double)nread_in_flanks;
    }
}

void EdgeEvidence::computeMetrics()
{
    if (nread == 0) return;

    find_split_read_edges();

    find_soft_clip_edges();

    count_biallelic_loci();

    if (Ref != NULL)
    {
        count_nonref_loci();
    }
}

void EdgeEvidence::addAttributes()
{
    if (srbp1 != 0)
    {
        m_Region->set_srbp1(srbp1);
        m_Region->set_srbp2(srbp2);
        m_Region->set_srf(srf);
    }

    if (scbp1 != 0)
    {
        m_Region->set_scbp1(scbp1);
        m_Region->set_scf1(scf1);
    }

    if (scbp2 != 0)
    {
        m_Region->set_scbp2(scbp2);
        m_Region->set_scf2(scf2);
    }

    m_Region->set_n_biallelic(nbiall);
    m_Region->set_n_unbalanced(nunbalanced);
    m_Region->set_fQ0(fracQ0);
    m_Region->set_fQ0flank(fracQ0flank);
    
    if (Ref != NULL)
    {
        m_Region->set_fNonRef_flank(fracNonRef_flank);
    }
    else
    {
        m_Region->set_fNonRef_flank(0.0);
    }
}

void EdgeEvidence::updateBreakpoints()
{
    if (scbp1 != 0 && scbp2 != 0)
    {
        if (! quietMode())
        {
            fprintf(stderr, "  updating breakpoints:  (%zu -- %zu) ---> (%zu -- %zu)   [soft-clip edges]\n",
                    m_Region->pos1(), m_Region->pos2(), scbp1, scbp2);
        }
        m_Region->setPos1_bp(scbp1);
        m_Region->setPos2_bp(scbp2);
    }
    //else if (srbp1 != 0 && srbp1 != 0)
    //{
    //    if (! quietMode())
    //    {
    //        fprintf(stderr, "  updating breakpoints:  (%zu -- %zu) ---> (%zu -- %zu)   [split-read edges]\n",
    //                m_Region->pos1(), m_Region->pos2(), srbp1, srbp2);
    //    }
    //    m_Region->setPos1_bp(srbp1);
    //    m_Region->setPos2_bp(srbp2);
    //}
}

void EdgeEvidence::find_split_read_edges()
{
    //Make sure SplitRead list is sorted in order of increasing position
    sort(SplitRead.begin(), SplitRead.end());

    size_t lastPos = 0;
    size_t nSR = 0;
    size_t max_nSR = 0;
    vector<uint32_t>srpos2;
    for (size_t i=0; i < SplitRead.size(); ++i)
    {
        if (SplitRead[i].first == lastPos)
        {
            ++nSR;
            srpos2.push_back(SplitRead[i].second);
        }
        else
        {
            if (nSR > max_nSR)
            {
                max_nSR = nSR;
                size_t ii = SplitRead[i].first - m_StartPos;
                if (nSR >4 && 5*nSR > cov[ii])
                {
                    //Determine the consensus (plurality) Pos2 value
                    sort(srpos2.begin(), srpos2.end());
                    //Save time: most times, they will all be the same value
                    uint32_t srpos2_consensus = srpos2[0];
                    if (srpos2[0] != srpos2[srpos2.size()-1])
                    {
                        size_t nmax = 0;
                        vector<uint32_t>::iterator it = srpos2.begin();
                        while (it != srpos2.end())
                        {
                            size_t n = 0;
                            uint32_t currentval = *it;
                            while (it != srpos2.end() && *it == currentval)
                            {
                                n++;
                                it++;
                            }

                            if (n > nmax)
                            {
                                nmax = n;
                                srpos2_consensus = currentval;
                            }
                        }
                    }

                    srbp1 = lastPos;
                    srbp2 = srpos2_consensus;
                    srf   = (double)nSR/(double)cov[ii];
                }
            }

            nSR = 0;
            srpos2.clear();
            lastPos = SplitRead[i].first;
        }
    }
}

void EdgeEvidence::find_soft_clip_edges()
{
    //Identify aligned soft-clip edges in the target region that have coherent soft-clipped sequences
    //In case more than one significant soft-clip edge is identified, we will keep the one 
    //that is nearest the called breakpoint.

    //Sort soft-clip reads in order of increasing position
    sort(SoftClipRead.begin(), SoftClipRead.end());
    size_t iiSC = 0;
    size_t minOffset1 = 1e9;
    size_t minOffset2 = 1e9;
    for (size_t i=0; i < m_TargetSize; ++i)
    {
        if (nSoftClip[i] > 4 && 8*nSoftClip[i] > cov[i])
        {
            //Position i contains a pileup of soft-clip events.
            //collect all events at this position and examine the aligned soft-clipped sequences
            vector<string> AlignedClippedBases;
            for (size_t ii=0; ii<(size_t)EdgeEvidence::SoftClipBufferSize; ++ii)
            {
                AlignedClippedBases.push_back("");
            }

            size_t nClippedBases = 0;
            while (iiSC < SoftClipRead.size() && SoftClipRead[iiSC].first <= i)
            {
                uint32_t scPos = SoftClipRead[iiSC].first;
                string scSeq = SoftClipRead[iiSC].second;
                if (scPos == i)
                {
                    for (size_t kk=0; kk < min(scSeq.length(), (size_t)EdgeEvidence::SoftClipBufferSize); ++kk)
                    {
                        AlignedClippedBases[kk] = AlignedClippedBases[kk] + scSeq[kk];
                        ++nClippedBases;
                    }
                }
                ++iiSC;
            }

            //Count the number of times AligndClippedBases has a non-majority call
            size_t ndiff = count_AlignedBases_mismatches(AlignedClippedBases);

            if (! quietMode())
            {
                cerr << "  sc edge: " << i+m_StartPos << "  " << nSoftClip[i] << "  " << cov[i] 
                     << " : " << ndiff << "  " << nClippedBases << endl;
            }
            
            //Pass the soft clip if fewer than 5% of basecalls
            //in the aligned clipped reads were not the majority call
            if (20*ndiff < nClippedBases)
            {
                //Count it as a breakpoint if it is in
                //the upstream/downstream zone, AND
                //it is closer to the called breakpoint than the current best edge
                //each edge will be considered only for one breakpoint (whichever it is nearer)
                size_t offset1 = abs((int)(i+m_StartPos) - (int)(m_Region->pos1(true)));
                size_t offset2 = abs((int)(i+m_StartPos) - (int)(m_Region->pos2(true)));
                if (i + m_StartPos < m_Region->pos1(true) + m_SoftLength2
                    && offset1 < offset2 && offset1 < minOffset1)
                {
                    scbp1 = i+m_StartPos;
                    scf1  = (double)nSoftClip[i]/(double)cov[i-1];
                    minOffset1 = offset1;
                }
                if (i + m_StartPos > m_Region->pos2(true) - m_SoftLength2
                    && offset2 < offset1 && offset2 < minOffset2)
                {
                    scbp2 = i + m_StartPos;
                    scf2  = (double)nSoftClip[i]/(double)cov[i];
                    minOffset2 = offset2;
                }
            }
        }
    }

    if (! quietMode())
    {
        cerr << "  scbp1: " << scbp1 << "   scbp2: " << scbp2 << endl;
    }
}

size_t EdgeEvidence::count_AlignedBases_mismatches(vector<string> &AlignedClippedBases)
{
    size_t ndiff = 0;
    for (size_t ii=0; ii<(size_t)EdgeEvidence::SoftClipBufferSize; ++ii)
    {
        string s = AlignedClippedBases[ii];
        sort(s.begin(), s.end());

        if (s[s.length()-1] != s[0])
        {
            //More than one base present; count them
            vector<size_t> nb;
            size_t nA = s.find_first_of("C");
            if (nA == string::npos) nA = 0;
            nb.push_back(nA);
            size_t nC = s.find_first_of("G");
            if (nC == string::npos)
                nC = 0;
            else
                nC = nC - nA;
            nb.push_back(nC);
            size_t nG = s.find_first_of("T");
            if (nG == string::npos)
                nG = 0;
            else
                nG = nG - nA - nC;
            nb.push_back(nG);

            size_t nT = s.length() - nA - nC - nG;
            nb.push_back(nT);

            //accumulate ndiff by the total non-majority calls
            sort(nb.begin(), nb.end());
            ndiff = ndiff + nb[0] + nb[1] + nb[2];
        }
    }

    return ndiff;
}

void EdgeEvidence::count_biallelic_loci()
{
    nbiall = 0;
    nunbalanced = 0;
    for (size_t i=0; i<m_TargetSize; ++i)
    {
        //Multiallelic positions.  Require at least 12 reads,
        //and must be between called breakpoints
        if (cov[i] >= 12 && i + m_StartPos > m_Region->pos1(true)
            && i + m_StartPos < m_Region->pos2(true))
        {
            float fA = (float)nA[i]/(float)cov[i];
            float fC = (float)nC[i]/(float)cov[i];
            float fG = (float)nG[i]/(float)cov[i];
            float fT = (float)nT[i]/(float)cov[i];
            if (fA < EdgeEvidence::FMIN && fC < EdgeEvidence::FMIN && fG < EdgeEvidence::FMIN && fT < EdgeEvidence::FMIN)
            {
                //if (! quietMode())
                //{
                //    cerr << "  multiallelic:   "
                //         << i + m_StartPos << "    "
                //         << nA[i] << "  " << nC[i] << "  "
                //         << nG[i] << "  " << nT[i] << "  "
                //         << cov[i]
                //         << endl;
                //    
                //}
                
                nbiall++;

                //Identify the largest allele depth
                size_t nMax = nA[i];
                size_t imax = 0;
                if (nC[i] > nMax)
                {
                    nMax = nC[i];
                    imax = 1;
                }
                if (nG[i] > nMax)
                {
                    nMax = nG[i];
                    imax = 2;
                }
                if (nT[i] > nMax)
                {
                    nMax = nT[i];
                    imax = 3;
                }
    
                //Identify the second-largest allele depth
                size_t nMax2 = 0;
                if (imax != 0 && nA[i] > nMax2)
                {
                    nMax2 = nA[i];
                }
                if (imax != 1 && nC[i] > nMax2)
                {
                    nMax2 = nC[i];
                }
                if (imax != 2 && nG[i] > nMax2)
                {
                    nMax2 = nG[i];
                }
                if (imax != 3 && nT[i] > nMax2)
                {
                    nMax2 = nT[i];
                }
    
                //if nMax2/(nMax+nMax2) is "unbalanced", accumulate the nunbalanced counter
                if (nMax+nMax2 > 0)
                {
                    double f = (double)nMax2/(double)(nMax + nMax2);
                    if (f > 0.20 && f < 0.30)
                    {
                        nunbalanced++;
                    }
                }
            }
        }
    }
}

void EdgeEvidence::count_nonref_loci()
{
    size_t n_non_ref = 0;
    size_t n_in_flank = 0;
    fracNonRef_flank = 0.0;

    for (size_t i=0; i<m_TargetSize; ++i)
    {
        //Skip loci inside the region itself (so we only count nonref loci in the flanks)
        if (cov[i] < 10 || 
            (i + m_StartPos > m_Region->pos1(true) &&
             i + m_StartPos < m_Region->pos2(true))
            )
        {
            continue;
        }

        ++n_in_flank;
        double f = (double)nNonRef[i]/(double)cov[i];
        if (f >= 0.005)
        {
            ++n_non_ref;
        }
    }

    fracNonRef_flank = (double)n_non_ref/(double)n_in_flank;
}

void EdgeEvidence::clearData()
{
    SoftClipRead.clear();
    SplitRead.clear();

    if (cov != NULL)
    {
        delete[] cov;
        cov = NULL;
    }
    if (nA != NULL)
    {
        delete[] nA;
        nA = NULL;
    }
    if (nC != NULL)
    {
        delete[] nC;
        nC = NULL;
    }
    if (nG != NULL)
    {
        delete[] nG;
        nG = NULL;
    }
    if (nT != NULL)
    {
        delete[] nT;
        nT = NULL;
    }
    if (nSoftClip != NULL)
    {
        delete[] nSoftClip;
        nSoftClip = NULL;
    }
}

