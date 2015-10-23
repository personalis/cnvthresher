/*
 * Fasta.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

#include "FileUtil.h"
#include "Fasta.h"

Fasta::Fasta()
{
    // (empty for now)
}

Fasta::Fasta(const string &fasta_file)
{
    //First call default ctor for any initialization of member data
    Fasta();

    initFromFile(fasta_file);
}

Fasta::~Fasta()
{
    // (empty for now)
}

void Fasta::initFromFile(const string &fasta_file)
{
    if (! fexists(fasta_file.c_str()))
    {
        cerr << "Specified FASTA file not found: " << fasta_file << endl;
        exit(1);
    }

    fstream fa;
    safeOpenStream(fa, fasta_file.c_str());

    size_t iseg = 0;
    while (fa.good())
    {
        string fa_line;
        getline(fa, fa_line);
        //trim trailing whitespace
        size_t endpos = fa_line.find_last_not_of(" \t\r\n");
        if (endpos != string::npos)
        {
            fa_line = fa_line.substr(0, endpos+1);
        }

        //Skip empty lines
        if (fa_line.length() == 0)
        {
            continue;
        }

        //Skip comment lines
        if (fa_line[0] == '#')
        {
            continue;
        }

        //If the line starts with ">", parse the segment name
        if (fa_line[0] == '>')
        {
            addSegment(fa_line);
            iseg = SegNames.size() - 1;
            continue;
        }

        //If we get here, we need to append the line 
        //to the current segment's sequence string.
        appendSequence(iseg, fa_line);
    }
}

void Fasta::addSegment(const string &line)
{
    size_t istart = 1;
    size_t istop = line.find(" ");
    string segName = line.substr(istart, istop-istart);
    //cerr << "  adding segment " << segName << endl;
    SegNames.push_back(segName);
    Ref.push_back("");
}

size_t Fasta::segmentIndex(const string &segName)
{
    vector<string>::iterator it = find(SegNames.begin(), SegNames.end(), segName);
    if (it == SegNames.end())
    {
        cerr << "Error: segment named '" << segName << "' not present in Reference." << endl;
        exit(1);
    }
    
    return (size_t)(it - SegNames.begin());
}
