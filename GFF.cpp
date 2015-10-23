/*
 * GFF.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "GRegion.h"
#include "GFF.h"

GFF::GFF()
{
}

GFF::GFF(const string &gff_file)
{
    fstream fgff(gff_file.c_str(), ios::in);
    bool inHeader = true;
    while (fgff.good())
    {
        string gff_line;
        getline(fgff, gff_line);
        //trim trailing whitespace
        size_t endpos = gff_line.find_last_not_of(" \t\n");
        if( endpos != string::npos )
        {
            gff_line = gff_line.substr( 0, endpos+1 );
        }

        //Ignore empty lines
        if (gff_line.length() == 0)
        {
            continue;
        }

        if (gff_line[0] != '#')
        {
            inHeader = false;
        }

        if (inHeader == true)
        {
            m_Header.push_back(gff_line);
        }
        else
        {
            if (gff_line[0] != '#')
            {
                GRegion *newRegion = new GRegion(gff_line);
                insertRegion(newRegion);
            }
        }
    }
}

GFF::~GFF()
{
    while(! m_RegionList.empty()) 
    {
        delete m_RegionList.back(); 
        m_RegionList.pop_back();
    }
}
        
void GFF::insertRegion(GRegion *newRegion)
{
    if (m_RegionList.empty())
    {
        m_RegionList.push_back(newRegion);
        return;
    }

    if (! areRegionsSorted())
    {
        cerr << "Error: GFF regions not sorted."  << endl;
        exit(1);
    }

    //Let's make a guess that the input GFF file is already sorted
    //(if so we can save some time)
    GRegion *lastRegion = m_RegionList[m_RegionList.size() - 1];
    if (newRegion->chrID() > lastRegion->chrID() || (newRegion->chrID() == lastRegion->chrID() && newRegion->pos1() >= lastRegion->pos1()))
    {
        m_RegionList.push_back(newRegion);
        return;
    }

    //Binary search to find sorted insertion index
    size_t ii = (size_t)(m_RegionList.size()/2);
    size_t delta = ii;

    while (delta > 0)
    {
        if (m_RegionList[ii]->chrID() == newRegion->chrID() && m_RegionList[ii]->pos1() == newRegion->pos1())
        {
            break;
        }

        if (m_RegionList[ii]->chrID() < newRegion->chrID())
        {
            size_t itry = ii + delta;
            if (itry >= m_RegionList.size())
            {
                itry = m_RegionList.size() - 1;
            }
            if (m_RegionList[itry]->chrID() <= newRegion->chrID())
            {
                ii = itry;
            }
        }
        else if (m_RegionList[ii]->chrID() > newRegion->chrID())
        {
            if (ii >= delta)
            {
                ii -= delta;
            }
            else
            {
                ii = 0;
            }
        }
        else //on the same chromosome
        {
            if (m_RegionList[ii]->pos1() < newRegion->pos1())
            {
                size_t itry = ii + delta;
                if (itry >= m_RegionList.size())
                {
                    itry = m_RegionList.size() - 1;
                }
                if (m_RegionList[itry]->pos1() < newRegion->pos1())
                {
                    ii = itry;
                }
            }
            else
            {
                if (ii >= delta)
                {
                    ii -= delta;
                }
                else
                {
                    ii = 0;
                }
            }
        }

        if (delta > 1)
        {
            delta = (size_t)((float)(delta + 1)/2);
        }
        else
        {
            delta = 0;
        }
    }

    while (ii < m_RegionList.size() && (m_RegionList[ii]->chrID() < newRegion->chrID() || (m_RegionList[ii]->chrID() == newRegion->chrID() && m_RegionList[ii]->pos1() < newRegion->pos1())))
    {
        ii++;
    }

    //DEBUG
    //cerr << "  " << ii << "  " << m_RegionList.size() << "  " << newRegion->chrName() << ":" << newRegion->pos1(true) << endl;

    if (ii >= m_RegionList.size())
    {
        m_RegionList.push_back(newRegion);
    }
    else
    {
        vector<GRegion*>::iterator regit = m_RegionList.begin() + ii;
        m_RegionList.insert(regit, newRegion);
    }
}

vector<uint8_t> GFF::chrIDs()
{
    vector<uint8_t> ids;
    for (vector<GRegion*>::iterator it = m_RegionList.begin(); it != m_RegionList.end(); it++)
    {
        ids.push_back((*it)->chrID());
    }
    return ids;
}

vector<size_t> GFF::startPositions()
{
    vector<size_t> startPos;
    for (vector<GRegion*>::iterator it = m_RegionList.begin(); it != m_RegionList.end(); it++)
    {
        startPos.push_back((*it)->pos1());
    }
    return startPos;
}

vector<size_t> GFF::stopPositions()
{
    vector<size_t> stopPos;
    for (vector<GRegion*>::iterator it = m_RegionList.begin(); it != m_RegionList.end(); it++)
    {
        stopPos.push_back((*it)->pos2());
    }
    return stopPos;
}

bool GFF::areRegionsSorted()
{
    uint8_t lastChr = 0;
    size_t lastPos = 0;
    for (vector<GRegion*>::iterator it = m_RegionList.begin(); it != m_RegionList.end(); it++)
    {
        if ((*it)->chrID() < lastChr) 
        {
            return false;
        }
        else if ((*it)->chrID() == lastChr)
        {
            if ((*it)->pos1() >= lastPos)
            {
                lastPos = (*it)->pos1();
            }
            else //pos1() < lastPos
            {
                return false;
            }
        }
        else //chrID() > lastChr
        {
            lastChr = (*it)->chrID();
            lastPos = (*it)->pos1();
        }
    }

    return true;
}

void GFF::add_evidence_attributes()
{
    //First determine robust coverage baseline by aggregating over 
    //flanking regions of all features on each chromosome
    compute_coverage_baseline();

    for (vector<GRegion*>::iterator itreg=m_RegionList.begin(); itreg != m_RegionList.end(); itreg++)
    {
        GRegion *reg = *itreg;
        //Divide by 2 so that normalized coverage represents local copy number
        //(i.e., we are assuming that the CovBaseline value represents 2 copies)
        double normFactor = m_CovBaseline[reg->chrID()] / 2;
        //Now that we have a proper normalization factor, we can compute MSCORE
        reg->compute_mscore(normFactor);

        ostringstream ss;
        ss << setiosflags(ios::fixed) << setprecision(3);

        if (reg->scbp1() > 0)
        {
            ss << "scbp1=" << reg->scbp1() << ";";
        }
        if (reg->scbp2() > 0)
        {
            ss << "scbp2=" << reg->scbp2() << ";";
        }
        if (reg->scbp1() > 0)
        {
            ss << "scf1=" << reg->scf1() << ";";
        }
        if (reg->scbp2() > 0)
        {
            ss << "scf2=" << reg->scf2() << ";";
        }
        if (reg->srbp1() > 0)
        {
            ss << "srbp1=" << reg->srbp1() << ";";
            ss << "srbp2=" << reg->srbp2() << ";";
            ss << "srf=" << reg->srf() << ";";
        }

        ss << setiosflags(ios::fixed) << setprecision(6);
        ss << "fbiall=" << reg->fbiall() << ";";
        ss << setiosflags(ios::fixed) << setprecision(3);
        ss << "funbal=" << reg->funbal() << ";";

        ss << "mscore=" << reg->mscore() << ";";
        ss << "ncov=" << reg->mean_cov_in_feature() / normFactor << ";";
        ss << "dcov=" << reg->edge_delta() / normFactor << ";";
        ss << "sdcov=" << reg->sdev_cov_in_feature() / normFactor << ";";
        ss << "sharp=" << reg->edge_sharpness() << ";";
        ss << "fgood=" << reg->frac_in_range() << ";";
        ss << "ninsize=" << reg->n_insize() << ";";
        ss << "fQ0=" << reg->fQ0() << ";";
        ss << "fQ0flank=" << reg->fQ0flank() << ";";
        ss << "fNonRef=" << reg->fNonRef() << ";";

        reg->addAttribute(ss.str());
    }
}

void GFF::compute_coverage_baseline()
{
    //Measure the median flank coverage value per chromosome
    vector<double> FlankCov;
    int lastchr = -1;
    for (vector<GRegion*>::iterator itreg=m_RegionList.begin(); itreg != m_RegionList.end(); itreg++)
    {
        GRegion *reg = *itreg;
        if (reg->chrID() != lastchr)
        {
            if (lastchr >= 0)
            {
                //Finished with previous chromosome
                //Add median flank coverage value to m_CovBaseline
                if (FlankCov.size() > 0)
                {
                    sort(FlankCov.begin(), FlankCov.end());
                    size_t ii = (size_t)(FlankCov.size() / 2); 
                    if (FlankCov.size() % 2 == 0)
                    {
                        m_CovBaseline[lastchr] = 0.5*(FlankCov[ii] + FlankCov[ii+1]);
                    }
                    else
                    {
                        m_CovBaseline[lastchr] = FlankCov[ii];
                    }
                }
                else
                {
                    m_CovBaseline[lastchr] = 0;
                }
            }

            FlankCov.clear();
            lastchr = reg->chrID();
        }

        //Aggregate flank coverage values
        if (reg->mean_cov_in_flank1() > 0)
        {
            FlankCov.push_back( reg->mean_cov_in_flank1() );
        }
        if (reg->mean_cov_in_flank2() > 0)
        {
            FlankCov.push_back( reg->mean_cov_in_flank2() );
        }
    }

    //Take care of the last chromosome
    if (FlankCov.size() > 0)
    {
        sort(FlankCov.begin(), FlankCov.end());
        size_t ii = (size_t)(FlankCov.size() / 2); 
        if (FlankCov.size() % 2 == 0)
        {
            m_CovBaseline[lastchr] = 0.5*(FlankCov[ii] + FlankCov[ii-1]);
        }
        else
        {
            m_CovBaseline[lastchr] = FlankCov[ii];
        }
    }
    else
    {
        m_CovBaseline[lastchr] = 0;
    }
}

void GFF::writeGFF(string gff_file, bool useBP=false)
{
    fstream fout(gff_file.c_str(), ios::out);

    for (vector<string>::iterator it=m_Header.begin(); it != m_Header.end(); it++)
    {
        fout << *it << endl;
    }

    for (vector<GRegion*>::iterator it=m_RegionList.begin(); it != m_RegionList.end(); it++)
    {
        GRegion *reg = *it;
        double score = reg->score();
        string sscore = ".";
        if (score != 0)
        {
            ostringstream ss;
            ss << setiosflags(ios::fixed) << setprecision(3) << score;
            sscore = ss.str();
        }

        fout << reg->chrName() << "\t" 
             << reg->source() << "\t" 
             << reg->typeString(useBP) << "\t"
             << reg->pos1(useBP) << "\t" 
             << reg->pos2(useBP) << "\t" 
             << sscore << "\t" 
             << reg->strandString() << "\t" 
             << reg->frameString() << "\t" 
             << reg->attributeString() << endl;
    }
}
