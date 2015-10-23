/*
 * FileUtil.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include <fstream>

using namespace std;

//Just some utility functions dealing with text file processing
string getValueFromConfig(string &bampath, const char *key, bool useJobCfg=false);

void safeOpenStream(fstream &ff, string &filename, ios_base::openmode mode=ios::in);
void safeOpenStream(fstream &ff, const char *filename, ios_base::openmode mode=ios::in);

bool fexists(const char *filename);

