/*
 * MorphEvidence.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_MORPHEVIDENCE_H_
#define BAMCNV_MORPHEVIDENCE_H_

#include "BamEvidence.h"

/*! \class MorphEvidence
 *! \brief handler for morphological evidence in support of a CNV call
 *! CNVs have a particular coverage-profile morphology, which we 
 *! parameterize using a number of metrics:
 *! + Interior normalized coverage: should be near an integer value
 *! + std. dev. of interior normalized coverage
 *! + delta-coverage: change in normalized coverage across the breakpoints
 *! + sharpness: measure of how abruptly coverage changes across the breakpoints
 *! + fgood: fraction of interior loci with coverage consistent with CNV call
 */
class MorphEvidence: public BamEvidence
{
public:
    MorphEvidence(size_t softLength, size_t insertsize, bool qmode);
    ~MorphEvidence();

    void traverseReads(BamReader &reader, GRegion *reg);
    void computeMetrics();
    void addAttributes();

protected:
    void resetForRegion(GRegion *reg);
    void clearData();

private:
    void mean_and_sd_coverage(size_t iStart, size_t iStop, size_t nskip, double &meanCov, double &sdCov);
    void refine_breakpoints();
    size_t binary_search(size_t iStart, size_t iStop, double cov0);

    uint16_t *cov;

    size_t m_RefinedStart;
    size_t m_RefinedStop;
    double m_MeanCov_InFeature;
    double m_MeanCov_Flank1;
    double m_MeanCov_Flank2;
    double m_SDCov_InFeature;
    double m_EdgeSharpness;
    double m_DeltaCov;
    double m_FracInRange;
};

#endif //BAMCNV_MORPHEVIDENCE_H_
