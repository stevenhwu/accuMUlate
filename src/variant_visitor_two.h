/*
 * variant_visitor_two.h
 *
 *  Created on: 12/9/14
 *      Author: Steven Wu
 */

#pragma once
#ifndef __VARIANT_VISITOR_TWO_H_
#define __VARIANT_VISITOR_TWO_H_

#include <stdint.h>
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"

using namespace BamTools;

extern int global_count[10];


class VariantVisitorTwo : public PileupVisitor{

    static const char ZERO_CHAR = ((char) 0);

public:
    VariantVisitorTwo(const RefVector& bam_references,
            const SamHeader& header,
            const Fasta& idx_ref,
            GenomeData& all_the_data,
//                      const ModelParams& p,
//            const SampleMap& samples,
//            BamAlignment& ali,
            int qual_cut,
            int mapping_cut,
            double prob_cut);



    virtual ~VariantVisitorTwo(void) {
    }

public:
    virtual void Visit(const PileupPosition& pileupData);
    static bool filter_data(ReadDataVector &read_vector);
private:

    Fasta m_idx_ref;
    RefVector m_bam_ref;
    SamHeader m_header;
//    SampleMap m_samples;

//        ModelParams m_params;
    int m_qual_cut;
//    BamAlignment& m_ali;
    GenomeData& m_all_the_data;
    double m_prob_cut;
    int m_mapping_cut;

    char current_base;
    uint64_t chr_index;

    char qual_cut_char;

    std::string rg_tag;



//    unordered_map<std::string, std::string> map_tag_sample_two_stage;
    std::unordered_map<std::string, int> map_tag_sample;

    int quit = 0;
    int total_samele_count;

    int GetSampleIndex(std::string const &tag_data);


};


#endif //__VARIANT_VISITOR_TWO_H_
