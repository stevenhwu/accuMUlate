#pragma once
#ifndef __VariantVisitor_H_
#define __VariantVisitor_H_

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"

using namespace BamTools;
class VariantVisitor : public PileupVisitor{
public:
    VariantVisitor(const RefVector& bam_references,
            const SamHeader& header,
            const Fasta& idx_ref,
            GenomeData& all_the_data,
//                      const ModelParams& p,
            const SampleMap& samples,
            BamAlignment& ali,
            int qual_cut,
            int mapping_cut,
            double prob_cut):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references),
            m_header(header), m_samples(samples),
            m_qual_cut(qual_cut), m_ali(ali),
            m_all_the_data(all_the_data), m_prob_cut(prob_cut),
            m_mapping_cut(mapping_cut)
    { }
    ~VariantVisitor(void) { }
public:
    void Visit(const PileupPosition& pileupData);
private:
    RefVector m_bam_ref;
    SamHeader m_header;
    Fasta m_idx_ref;
    GenomeData& m_all_the_data;
    SampleMap m_samples;
    BamAlignment& m_ali;
//        ModelParams m_params;
    int m_qual_cut;
    int m_mapping_cut;
    double m_prob_cut;
    char current_base;
    uint64_t chr_index;
};


#endif //__VariantVisitor_H_
