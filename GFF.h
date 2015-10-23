/*
 * GFF.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_GFF_H_
#define BAMCNV_GFF_H_

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

using namespace std;

class GRegion;

/*! \class GFF
 *! \brief Encapsulate the data in a GFF file
 *! \author Jason Harris
 *! \copyright 2014 Personalis, Inc.
 */
class GFF {
 public:
    //! default constructor
    //! Create an empty GFF object 
    GFF();

    //! constructor
    //! Initialize the GFF object from a file
    GFF(const string &gff_file);

    //! destructor
    ~GFF();

    //Accessors
    inline vector<string> header() const { return m_Header; }
    vector<GRegion*>& regionList() { return m_RegionList; }
    vector<size_t> startPositions();
    vector<size_t> stopPositions();
    vector<uint8_t> chrIDs();

    void setHeader(vector<string> h) { m_Header = h; }
    void addToHeader(const string &newLine) { m_Header.push_back(newLine); }
    //JH: Add if needed (note: should add in sorted order...optionally?)
    //void insertRegion(const string &gff_line);
    void insertRegion(GRegion *newRegion);

    inline size_t nRegions() { return m_RegionList.size(); }
    bool areRegionsSorted();

    void add_evidence_attributes();
    void writeGFF(string gff_file, bool useBP);

 private:
    void compute_coverage_baseline();

    map<int,double> m_CovBaseline;
    vector<GRegion*> m_RegionList;
    vector<string> m_Header;
};


#endif //BAMCNV_GFF_H_
