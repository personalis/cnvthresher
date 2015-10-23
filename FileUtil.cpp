/*
 * FileUtil.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <string>
#include <iostream>
#include <fstream>
#include <cstring> //For strerror
#include <errno.h> //For errno
#include <stdlib.h>
#include "FileUtil.h"

using namespace std;

/**
  * \function getValueFromConfig
  * \parameter bampath path where th config file will be found
  * \parameter key the key name for which the value will be returned
  * \return the value associated with the specified key
  */
string getValueFromConfig(string &bampath, const char *key, bool useJobCfg)
{
    //The config file is named either .pbr.cfg (recent builds)
    //or .pbr.cfg (older builds)
    string cfgfile; 
    if (useJobCfg == true)
    {
        cfgfile = bampath+"/.job.cfg";
        if (!fexists(cfgfile.c_str()))
        {
            cerr << "Config file not found: " << cfgfile << endl;
            throw(21);
        }
    }
    else
    {
        cfgfile = bampath+"/.pbr.cfg";
        string cfgfile2 = bampath+"/.job.cfg";
        if (!fexists(cfgfile.c_str()))
        {
            cfgfile = cfgfile2;
            if (!fexists(cfgfile.c_str()))
            {
                cerr << "Config file not found: " << cfgfile << endl;
                throw(21);
            }
        }
    }

    fstream fcfg(cfgfile.c_str(), ios::in);
    fcfg.exceptions ( ifstream::failbit | ifstream::badbit );
    if (!fcfg.is_open())
    {
        cerr << "could not open config file for reading: " << cfgfile << endl;
        throw(21);
    }
    
    string sline;
    string svalue = "";
    try {
        while (getline(fcfg, sline))
        {
            //ignore comment lines
            if (sline[0]=='#') continue;
            
            //Find the line containing the specified key
            size_t ii = sline.find(key);
            if (ii != string::npos)
            {
                size_t istart = sline.find("=")+1;
                svalue = sline.substr(istart);
                break;
            }
        }
    } catch (ifstream::failure e) {
        if (!fcfg.eof() && !fcfg.bad()) //eof and bad are not a problem
        {
            cerr << "Error: " << e.what() << endl;
            throw(21);
        }
    }
    fcfg.close();

    if (svalue.length() == 0)
    {
        //If we tried .pbr.cfg, and .job.cfg also exists, 
        //try again with .job.cfg
        size_t ii = cfgfile.find(".pbr.cfg");
        if (ii != string::npos)
        {
            svalue = getValueFromConfig(bampath, key, true); //force .job.cfg
            if (svalue.length() > 0)
            {
                return svalue;
            }
        }
        
        cerr << "Error: specified key (" << key << ") was not found in config file: " << cfgfile << endl;
        throw(21);
    }

    return svalue;
}

void safeOpenStream(fstream &ff, string &filename, ios_base::openmode mode)
{
    safeOpenStream(ff, filename.c_str(), mode);
}

void safeOpenStream(fstream &ff, const char *filename, ios_base::openmode mode)
{
    ff.open(filename, mode);
    if (!ff.good())
    {
        cerr << "Could not open file " << filename
             << ": " << strerror(errno) << endl;
        exit(1);
    }
}

bool fexists(const char *filename)
{
    bool result = false;
    ifstream ifile(filename);
    result = ifile;
    if (result) ifile.close();
    return result;
}

