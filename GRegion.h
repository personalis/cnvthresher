/*
 * GRegion.h
 *
 * This is part of CNVthresher.  Copyright (C) 2015 by Personalis, Inc.
 * Original author: Jason Harris <jason.harris@personalis.com>
 */
#ifndef BAMCNV_REGION_H_
#define BAMCNV_REGION_H_

#include <string>
#include <vector>
#include <stdint.h>

enum Zygosity {Zyg_Unknown=0, Homozygous=1, Heterozygous=2, Hemizygous=3};
enum FeatureType {Unknown=0, Deletion=1, Duplication=2, Insertion=3, Inversion=4, Translocation=5};
enum Strand {NOSTRAND, FWD, REV};
enum Frame {NOFRAME, INFRAME, OFF1, OFF2};

using namespace std;

/*! \class GRegion
 *! \brief Encapsulate a genomic region as defined in a GFF file
 */

class GRegion {
 public:
    //! static member for list of chromosome names
    static vector<string> ChrName;
    static vector<string> TypeName;

    //! constructor
    //! Each member initialized by an argument
    GRegion(const uint8_t &ichr, const size_t &pos1, const size_t &pos2, const FeatureType &type, 
           const string &source="", const string &attribute="", 
           const double score=0, const Strand strand=NOSTRAND, const Frame frame=NOFRAME);
    //! constructor
    //! Initialized by a GFF line (typical use case)
    GRegion(const string &gff_line);

    static const string timestamp();
    static void initTypeName();

    void compute_mscore(double normFactor=0);

    //Accessors
    inline uint8_t chrID() const { return m_ichr; }
    inline string chrName() const { return ChrName[m_ichr]; }
    inline size_t pos1(bool useBP=false) const { return useBP ? m_start_bp : m_start; }
    inline size_t pos2(bool useBP=false) const { return useBP ? m_stop_bp  : m_stop; }
    inline size_t size(bool useBP=false) const { return useBP ? m_stop_bp - m_start_bp + 1 : m_stop - m_start + 1; }
    inline bool pos1_refined() const { return m_start_refined; }
    inline bool pos2_refined() const { return m_stop_refined; }
    inline FeatureType type(bool useBP=false) const { return useBP ? m_type_bp : m_type; }
    inline Zygosity zygosity() const { return m_Zyg; }
    inline string source() const { return m_source; }
    inline vector<string> attributeList() const { return m_AttributeList; }
    inline double score() const { return m_score; }
    inline Strand strand() const { return m_strand; }
    inline Frame frame() const { return m_frame; }
    string typeString(bool useBP=false) const;
    string strandString() const;
    string frameString() const;
    string attributeString() const;

    //Getters for EdgeEvidence attributes
    inline size_t scbp1() const { return m_scbp1; }
    inline size_t scbp2() const { return m_scbp2; }
    inline double scf1() const { return m_scf1; }
    inline double scf2() const { return m_scf2; }
    inline size_t srbp1() const { return m_srbp1; }
    inline size_t srbp2() const { return m_srbp2; }
    inline double srf() const { return m_srf; }
    inline size_t n_biallelic() const { return m_nbiall; }
    inline size_t n_unbalanced() const { return m_nunbalanced; }
    inline double fbiall() const { return (double)n_biallelic()/(double)size(true); }
    inline double funbal() const { return n_biallelic() ? (double)n_unbalanced()/(double)n_biallelic() : 0.0; }
    inline double fQ0() const { return m_fQ0; }
    inline double fQ0flank() const { return m_fQ0flank; }
    inline double fNonRef() const { return m_fNonRef_flank; }

    //Getter for InsSizeEvidence attribute
    inline size_t n_insize() const { return m_ninsize; }

    //Getters for MorphEvidence attributes
    inline double mean_cov_in_feature() const { return m_MeanCov_InFeature; }
    inline double sdev_cov_in_feature() const { return m_SDCov_InFeature; }
    inline double mean_cov_in_flank1() const { return m_MeanCov_Flank1; }
    inline double mean_cov_in_flank2() const { return m_MeanCov_Flank2; }
    inline double edge_sharpness() const { return m_EdgeSharpness; }
    inline double edge_delta() const { return m_DeltaCov; }
    inline double frac_in_range() const { return m_FracInRange; }
    inline double mscore() const { return m_mscore; }
    
    //Setters
    inline void setChrID(const uint8_t &ichr) { m_ichr = ichr; }
    inline void setPos1(const size_t &pos1) { m_start = pos1; }
    inline void setPos2(const size_t &pos2) { m_stop = pos2; }
    inline void setType(const FeatureType &type) { m_type = type; }
    inline void setZygosity(const Zygosity &zyg) { m_Zyg = zyg; }
    inline void setPos1_bp(const size_t &pos1) { m_start_bp = pos1; m_start_refined = true; }
    inline void setPos2_bp(const size_t &pos2) { m_stop_bp = pos2; m_stop_refined = true; }
    inline void setType_bp(const FeatureType &type) { m_type_bp = type; }
    void setType(string &stype, bool useBP=false);
    inline void setSource(const string &source) { m_source = source; }
    void setAttributesFromString(const string &attributeList);
    void addAttribute(const string &newAttribute);
    inline void setScore(const double &score) { m_score = score; }
    inline void setStrand(const Strand &strand) { m_strand = strand; }
    inline void setFrame(const Frame &frame) { m_frame = frame; }

    //Setters for EdgeEvidence attributes
    inline void set_scbp1(const size_t new_bp) { m_scbp1 = new_bp; }
    inline void set_scbp2(const size_t new_bp) { m_scbp2 = new_bp; }
    inline void set_srbp1(const size_t new_bp) { m_srbp1 = new_bp; }
    inline void set_srbp2(const size_t new_bp) { m_srbp2 = new_bp; }
    inline void set_scf1(const double new_frac) { m_scf1 = new_frac; }
    inline void set_scf2(const double new_frac) { m_scf2 = new_frac; }
    inline void set_srf(const double new_frac) { m_srf = new_frac; }

    inline void set_n_biallelic(const size_t nn) { m_nbiall = nn; }
    inline void set_n_unbalanced(const size_t nn) { m_nunbalanced = nn; }
    inline void set_fQ0(const double fQ0) { m_fQ0 = fQ0; }
    inline void set_fQ0flank(const double fQ0flank) { m_fQ0flank = fQ0flank; }
    inline void set_fNonRef_flank(const double fNonRef) { m_fNonRef_flank = fNonRef; }
    //Setter for InsSizeEvidence attribute
    inline void set_n_insize(const size_t nn) { m_ninsize = nn; }

    //Setters for MorphEvidence attributes
    inline void set_mean_cov_in_feature(const double v) { m_MeanCov_InFeature = v; }
    inline void set_sdev_cov_in_feature(const double v) { m_SDCov_InFeature = v; }
    inline void set_mean_cov_in_flank1(const double v) { m_MeanCov_Flank1 = v; }
    inline void set_mean_cov_in_flank2(const double v) { m_MeanCov_Flank2 = v; }
    inline void set_edge_sharpness(const double v) { m_EdgeSharpness = v; }
    inline void set_edge_delta(const double v) { m_DeltaCov = v; }
    inline void set_frac_in_range(const double v) { m_FracInRange = v; }
    inline void set_mscore(const double v) { m_mscore = v; }

 private:
    void init_evidence_attributes();

    uint8_t m_ichr;
    size_t m_start,     m_stop;
    size_t m_start_bp,  m_stop_bp;
    bool m_start_refined, m_stop_refined;
    FeatureType m_type, m_type_bp;
    Zygosity m_Zyg;
    string m_source;
    vector<string> m_AttributeList;
    double m_score;
    Strand m_strand;
    Frame m_frame;

    //EdgeEvidence attributes
    size_t m_scbp1, m_scbp2, m_srbp1, m_srbp2;
    size_t m_nbiall, m_nunbalanced;
    double m_scf1, m_scf2, m_srf;
    double m_fQ0, m_fQ0flank, m_fNonRef_flank;

    //InsSize attribute
    size_t m_ninsize;

    //MorphEvidence attributes
    double m_MeanCov_InFeature, m_SDCov_InFeature, m_MeanCov_Flank1, m_MeanCov_Flank2;
    double m_EdgeSharpness, m_DeltaCov, m_FracInRange, m_mscore;
};

#endif //BAMCNV_REGION_H_
