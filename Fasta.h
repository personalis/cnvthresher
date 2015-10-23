/*
 * Fasta.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_FASTA_H_
#define BAMCNV_FASTA_H_

#include <vector>
#include <stdlib.h>

using namespace std;

//! class Fasta
//! \brief Encapsulate a FASTA file, typically describing a reference sequence
//! Segment sequences are stored as strings. The collection of sequences
//! is stored in a string vector.
class Fasta
{
 public :
    //! Default ctor.  Create an empty Fasta object
    Fasta();
    //! Constructor.  Initialize a Fatsa object from a FATSA file
    Fasta(const string &fasta_file);
    //! dtor
    ~Fasta();

    //! \return the reference base at the given location
    //! \param iseg the index of the segment
    //! \param pos the zero-based index position to query
    inline char baseAt(const size_t &iseg, const size_t &pos) { return Ref[iseg][pos]; };
    inline char baseAt(const string &segName, const size_t &pos)
    {
        size_t iseg = segmentIndex(segName);
        return baseAt(iseg, pos);
    }

    //! \return a substring representing the reference sequence in a given region
    //! \param iseg the index of the segment
    //! \param pstart the zero-based index position at the start of the query
    //! \param pstop the zero-based index position at the end of the query
    inline string subsequence(const size_t &iseg, const size_t &pstart, const size_t &pstop)
    {
        return Ref[iseg].substr(pstart, pstop-pstart+1);
    }
    inline string subsequence(const string &segName, const size_t &pstart, const size_t &pstop)
    {
        size_t iseg = segmentIndex(segName);
        return subsequence(iseg, pstart, pstop);
    }

    inline size_t numberOfSegments() { return SegNames.size(); }

    size_t segmentIndex(const string &segName);
    inline string segmentName(const size_t &iseg) { return SegNames[iseg]; }
    inline size_t segmentLength(const size_t &iseg) { return Ref[iseg].length(); }
    inline size_t segmentLength(const string &segName)
    {
        size_t iseg = segmentIndex(segName);
        return segmentLength(iseg);
    }

 private:
    void initFromFile(const string &fasta_file);
    void addSegment(const string &line);
    inline void appendSequence(const size_t iseg, const string &line) { Ref[iseg] += line; }

    vector<string> Ref;
    vector<string> SegNames;
};

#endif //BAMCNV_FASTA_H_
