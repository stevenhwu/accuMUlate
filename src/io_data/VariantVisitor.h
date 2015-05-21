#pragma once
#ifndef __VariantVisitor_H_
#define __VariantVisitor_H_

#include "api/BamReader.h"
//#include "bamtools_pileup_engine.h"
//#include "bamtools_fasta.h"

#include <mutations/model.h>
#include "parsers.h"

using namespace BamTools;
class VariantVisitor : public LocalBamToolsUtils::PileupVisitor{
public:
    VariantVisitor(const RefVector& bam_references,
            const SamHeader& header,
            const LocalBamToolsUtils::Fasta& idx_ref,
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
    virtual ~VariantVisitor(void) { }
public:
    virtual void Visit(const LocalBamToolsUtils::PileupPosition& pileupData);
private:
    LocalBamToolsUtils::Fasta m_idx_ref;
    RefVector m_bam_ref;
    SamHeader m_header;
    SampleMap m_samples;

//        ModelParams m_params;
    int m_qual_cut;
    BamAlignment& m_ali;
    GenomeData& m_all_the_data;
    double m_prob_cut;
    int m_mapping_cut;

    char current_base;
    uint64_t chr_index;
};


#endif //__VariantVisitor_H_
