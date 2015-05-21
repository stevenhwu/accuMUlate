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
#include "parsers.h"


#include <mutations/model.h>

#include "genome_data_stream.h"

using namespace BamTools;

extern int global_count[10];


class VariantVisitorTwo : public LocalBamToolsUtils::PileupVisitor{

    static const char ZERO_CHAR = ((char) 0);

public:
    VariantVisitorTwo(const RefVector &bam_references, const SamHeader &header, const LocalBamToolsUtils::Fasta &idx_ref,
            GenomeData &all_the_data, std::string binary_outfile, int qual_cut, int mapping_cut, double prob_cut);



    virtual ~VariantVisitorTwo(void) {
        gd_stream.close();
    }

public:
    virtual void Visit(const LocalBamToolsUtils::PileupPosition& pileupData);
    static bool filter_data(ReadDataVector &read_vector);

private:


    RefVector m_bam_ref;
    SamHeader m_header;
    LocalBamToolsUtils::Fasta m_idx_ref;
    GenomeData &m_all_the_data;
    int m_qual_cut;
    int m_mapping_cut;

    char current_base;
    char qual_cut_char;
    std::string rg_tag;
    std::unordered_map<std::string, int> map_tag_sample;


    double m_prob_cut;
//    SampleMap m_samples;
//        ModelParams m_params;
//    uint64_t chr_index;
//    unordered_map<std::string, std::string> map_tag_sample_two_stage;

private:
    int total_sample_count;

    int GetSampleIndex(std::string const &tag_data);


public:
    int getTotal_sample_count() const ;

    GenomeDataStream gd_stream;
};


#endif //__VARIANT_VISITOR_TWO_H_
