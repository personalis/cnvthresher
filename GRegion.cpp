/*
 * GRegion.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath> //abs(double)
#include <cstdlib>
#include <algorithm>

#include "GRegion.h"

vector<string> GRegion::ChrName;
vector<string> GRegion::TypeName;

GRegion::GRegion(const uint8_t &ichr, const size_t &pos1, const size_t &pos2, const FeatureType &type,
               const string &source, const string &attribute,
               const double score, const Strand strand, const Frame frame)
{
    m_ichr = ichr;
    m_start = pos1;
    m_stop  = pos2;
    m_type  = type;
    m_start_bp = pos1;
    m_stop_bp  = pos2;
    m_start_refined = false;
    m_stop_refined = false;
    m_type_bp  = type;
    m_source = source;
    m_score = score;
    m_strand = strand;
    m_frame = frame;
    m_Zyg = Zyg_Unknown;

    init_evidence_attributes();

    setAttributesFromString(attribute);
}

GRegion::GRegion(const string &gff_line)
{
    stringstream ss(gff_line);
    string sschr, stype, sscore, sstrand, sframe, sattrib;
    try {
        ss >> sschr >> m_source >> stype >> m_start >> m_stop;
    }
    catch (int e)
    {
        cerr << "Error: invalid GFF file?  Line was: " << gff_line << endl;
        exit(1);
    }
    m_start_bp = m_start;
    m_stop_bp = m_stop;

    //Assign chrID as index into ChrName
    m_ichr = 255;
    for (size_t i=0; i<ChrName.size(); ++i)
    {
        if (sschr.compare(ChrName[i]) == 0)
        {
            m_ichr = i;
            break;
        }
    }
    if (m_ichr == 255)
    {
        cerr << "Error: did not find '" << sschr << "' among known chromosome names." << endl;
        cerr << "       known names: ";
        for (size_t i=0; i<ChrName.size(); ++i)
        {
            cerr << ChrName[i] << "  ";
        }
        cerr << endl;
        exit(1);
    }

    setType(stype, false);
    m_type_bp = m_type;

    if (ss >> sscore)
    {
        //Note: m_score set to zero when the string is not a number (in particular "." is used when there is no score)
        if (sscore.compare(".") == 0)
        {
            m_score = 0;
        }
        else
        {
            m_score = atof(sscore.c_str());
        }

        if (ss >> sstrand)
        {
            if (sstrand.compare(".") == 0)
            {
                m_strand = NOSTRAND;
            }
            else if (sstrand.compare("+") == 0)
            {
                m_strand = FWD;
            }
            else if (sstrand.compare("-") == 0)
            {
                m_strand = REV;
            }
            else
            {
                cerr << "Warning: STRAND field must be one of '.', '+', or '-'.  Line was: " << gff_line << endl;
                m_strand = NOSTRAND;
            }

            if (ss >> sframe)
            {
                if (sframe.compare(".") == 0)
                {
                    m_frame = NOFRAME;
                }
                else if (sframe.compare("0") == 0)
                {
                    m_frame = INFRAME;
                }
                else if (sframe.compare("1") == 0)
                {
                    m_frame = OFF1;
                }
                else if (sframe.compare("2") == 0)
                {
                    m_frame = OFF2;
                }
                else
                {
                    cerr << "Warning: FRAME field must be one of '.', '0', '1', or '2'.  Line was: " << gff_line << endl;
                    m_frame = NOFRAME;
                }

                size_t len = 255;
                char buff[len];
                ss.getline(buff, len);
                sattrib = buff;
                //Remove leading whitespace
                size_t startpos = sattrib.find_first_not_of(" \t\n");
                if( startpos != string::npos )
                {
                    sattrib = sattrib.substr(startpos);
                    setAttributesFromString(sattrib);
                }
            }
        }
    }

    init_evidence_attributes();
}

void GRegion::init_evidence_attributes()
{
    m_scbp1 = 0;
    m_scbp2 = 0;
    m_srbp1 = 0;
    m_srbp2 = 0;
    m_nbiall = 0;
    m_nunbalanced = 0;
    m_ninsize = 0;

    m_scf1 = 0.0;
    m_scf2 = 0.0;
    m_srf = 0.0;
    m_fQ0 = 0.0;
    m_fQ0flank = 0.0;
    m_fNonRef_flank = 0.0;

    m_MeanCov_InFeature = 0.0;
    m_SDCov_InFeature = 0.0;
    m_MeanCov_Flank1 = 0.0;
    m_MeanCov_Flank2 = 0.0;
    m_EdgeSharpness = 0.0;
    m_DeltaCov = 0.0;
    m_FracInRange = 0.0;
    m_mscore = 0.0;
}

void GRegion::initTypeName()
{
    TypeName.push_back("Unknown");
    TypeName.push_back("Deletion");
    TypeName.push_back("Duplication");
    TypeName.push_back("Insertion");
    TypeName.push_back("Inversion");
    TypeName.push_back("Translocation");
}

void GRegion::compute_mscore(double normFactor)
{
    //If normFactor is zero (default), use average of the flank coverage values
    if (normFactor == 0)
    {
        //Divide flank coverage by 2 so normalized values represent local copy number
        if (m_MeanCov_Flank1 > 0 && m_MeanCov_Flank2 > 0)
        {
            normFactor = 0.5*(m_MeanCov_Flank1 + m_MeanCov_Flank2) / 2;
        }
        else if (m_MeanCov_Flank1 > 0)
        {
            normFactor = m_MeanCov_Flank1 / 2;
        }
        else if (m_MeanCov_Flank2 > 0)
        {
            normFactor = m_MeanCov_Flank2 / 2;
        }
        else
        {
            //Give up! Just use raw coverage numbers (bad)
            normFactor = 1.0;
        }
    }

    m_mscore = 0.0;

    //set normalized values of cov_in_feature, dcov, and sdcov
    double ncov  = mean_cov_in_feature() / normFactor;
    double sdcov = sdev_cov_in_feature() / normFactor;
    double dcov  = edge_delta() / normFactor;

    //set cov0, the "expected" normalized coverage given the CNV type
    //by default, just set to the rounded value of the normalized in-feature coverage
    //(only appropriate for duplications and deletions with no zygosity)
    double cov0 = (double)((size_t)(ncov + 0.5));
    if (type(true) == Deletion)
    {
        if (zygosity() == Homozygous || zygosity() == Hemizygous)
        {
            cov0 = 0.0;
        }
        else if (zygosity() == Heterozygous)
        {
            cov0 = 1.0;
        }
    }
    else if (type(true) == Duplication && zygosity() == Heterozygous)
    {
        cov0 = 3.0;
    }

    //dcov0 is the expected coverage delta, given cov0
    double dcov0 = cov0 - 2.0;

    // normalized coverage in the feature (want it to be near cov0)
    double f = 1.0 - 2.0*abs(ncov - cov0);
    m_mscore += min(max(f, 0.0), 1.0);

    // std dev of coverage (smaller is better)
    f = 1.2 - sdcov;
    m_mscore += min(max(f, 0.0), 1.0);

    // delta coverage (want it to be near dcov0)
    f = 1.1 - abs(dcov0 - dcov);
    m_mscore += min(max(f, 0.0), 1.0);

    // edge sharpness
    f = edge_sharpness();
    m_mscore += min(max(f, 0.0), 1.0);

    // fraction of coverage points in expected range
    f = frac_in_range();
    m_mscore += min(max(f, 0.0), 1.0);
}

string GRegion::attributeString() const
{
    string s = "";
    if (m_AttributeList.size() > 0)
    {
        for (vector<string>::const_iterator it = m_AttributeList.begin(); it != m_AttributeList.end(); ++it)
        {
            if (it == m_AttributeList.begin())
            {
                s = *it;
            }
            else
            {
                s = s + ";" + *it;
            }
        }
    }

    return s;
}

void GRegion::setType(string &stype, bool useBP)
{
    FeatureType new_type;
    vector<string>::iterator it = find(TypeName.begin(), TypeName.end(), stype);
    if (it != TypeName.end())
    {
        new_type = (FeatureType)(it - TypeName.begin());
    }
    else
    {
        //string not found in TypeName; default to Unknown
        cerr << "Warning: unrecognized SV type '" << stype << "'; using 'Unknown'." << endl;
        new_type = Unknown;
    }

    if (useBP)
    {
        m_type_bp = new_type;
    }
    else
    {
        m_type = new_type;
    }
}

string GRegion::typeString(bool useBP) const
{
    return TypeName[type(useBP)];
}

string GRegion::strandString() const
{
    string s;
    switch(strand())
    {
    case NOSTRAND: s = "."; break;
    case FWD: s = "+"; break;
    case REV: s = "-"; break;
    }

    return s;
}

string GRegion::frameString() const
{
    string s;
    switch(frame())
    {
    case NOFRAME: s = "."; break;
    case INFRAME: s = "0"; break;
    case OFF1: s = "1"; break;
    case OFF2: s = "2"; break;
    }

    return s;
}

void GRegion::setAttributesFromString(const string &attributeList)
{
    m_AttributeList.clear();
    addAttribute(attributeList);

    setZygosity(Zyg_Unknown);

    //Check for zygosity using Zyg tag
    size_t izyg = attributeList.find("Zyg=");
    if (izyg != string::npos)
    {
        izyg += 4;
        size_t iend = attributeList.find(";", izyg);
        string ssZyg = attributeList.substr(izyg, iend-izyg);
        size_t ii = ssZyg.find("omozyg");
        if (ii != string::npos)
        {
            setZygosity(Homozygous);
            return;
        }

        ii = ssZyg.find("eterozyg");
        if (ii != string::npos)
        {
            setZygosity(Heterozygous);
            return;
        }

        ii = ssZyg.find("emizyg");
        if (ii != string::npos)
        {
            setZygosity(Hemizygous);
            return;
        }
    }

    //Check for zygosity using GT tag
    size_t igt = attributeList.find("GT");
    if (igt != string::npos)
    {
        igt += 3;
        size_t iend = attributeList.find(";", igt);
        string ssGT = attributeList.substr(igt, iend-igt);
        if (ssGT.compare("0/0") == 0)
        {
            setZygosity(Homozygous);
        }
        else if (ssGT.compare("0/1") == 0)
        {
            setZygosity(Heterozygous);
        }
        else if (ssGT.compare("1") == 0)
        {
            setZygosity(Hemizygous);
        }
    }

    return;
}

void GRegion::addAttribute(const string &newAttribute)
{
    stringstream ss(newAttribute);
    while( ss.good() )
    {
        string w;
        getline( ss, w, ';' );
        if (w.length() > 0)
        {
            m_AttributeList.push_back( w );
        }
    }
}

const string GRegion::timestamp()
{
    time_t rawtime;
    struct tm * timeinfo;
    char *buffer = new char[80];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime(buffer,80,"%T",timeinfo);
    string st = buffer;
    st = "["+st+"]";
    delete[] buffer;
    return st;
}
