/*
 * BamEvidence.cpp
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#include "BamEvidence.h"

BamEvidence::BamEvidence(size_t softLength, size_t insertsize, bool qmode)
{
    m_QuietMode   = qmode;
    m_SoftLength  = softLength;
    m_SoftLength2 = softLength;
    m_InsertSize  = insertsize;
    m_Region      = NULL;
}

BamEvidence::~BamEvidence()
{
}

