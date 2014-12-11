/*
 * variant_visitor_two.cc
 *
 *  Created on: 12/9/14
 *      Author: Steven Wu
 */


#include "variant_visitor_two.h"


int global_count = 0 ;


void VariantVisitorTwo::Visit(const PileupPosition &pileupData)
{

//    string chr = m_bam_ref[pileupData.RefId].RefName;
    uint64_t pos = pileupData.Position;
    ReadDataVector bcalls(total_samele_count, ReadData{{0, 0, 0, 0}});

    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);//Check: Room for improvemnt here?
//    string tag_id;

    uint16_t ref_base_idx = base_index2[(int)current_base];
    if (ref_base_idx < 4) { //TODO Model for bases at which reference is 'N'

    for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {//auto it2 = *it;
//    for (auto it2 : pileupData.PileupAlignments) {
//        auto *it = &it2;
//            if (include_site(*it, m_mapping_cut, m_qual_cut)) {
        int32_t pos_in_alignment = it->PositionInAlignment;
//        if (include_site_3(*it, m_mapping_cut, qual_cut_char)) {
        if (include_site_4( it->Alignment, pos_in_alignment, m_mapping_cut, qual_cut_char)) {


            string tag_data = it->Alignment.TagData;
            size_t start_index = tag_data.find(rg_tag);
            if(start_index != string::npos){
                start_index+=4;
            }
            else{
                size_t x = tag_data.find("RG");
                cout << "ERRER: "<< tag_data << "\t" << x << endl;
            }
            size_t end_index = tag_data.find(  ZERO_CHAR, start_index);
            string tag_id_2 = tag_data.substr(start_index, (end_index - start_index));
//
//
////            string tag_id;
////            it->Alignment.GetTag("RG", tag_id);
////            string sm = m_header.ReadGroups[tag_id].Sample;
////            if(tag_id != tag_id_2){
////                cout << tag_id << "\t" << tag_id_2 <<  "\t" << start_index << "\t" << end_index << "\n" << tag_data << endl;
////                exit(-1);
////            }
////            if(sm != sm2){
////                cout << sm << "\t" << sm2 << endl;
////                exit(-2);
////            }
//
////            string sm2 = map_tag_sample_two_stage[tag_id_2]; //TODO:catch exception
////            uint32_t sindex = m_samples[sm2]; //TODO check samples existed!
//
            int sindex = map_tag_sample[tag_id_2];

//            cout << tag_id_2 << endl;
//            cout << tag_id_2 << "\t" << sm2 << "\t" << sindex << "\t" << sindex2 << endl;


            uint16_t bindex = base_index2[(int) it->Alignment.QueryBases[pos_in_alignment]];

            if (bindex < 4) {
                global_count++;
                bcalls[sindex].reads[bindex] += 1;
            }

        }
    }


//                ModelInput d = {ref_base_idx, bcalls};
        m_all_the_data.push_back(ModelInput{ref_base_idx, bcalls});

    }
    else{//Hardly ever happen
//        quit ++;
//        bool x = m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
//        cout << "ref_base > 4: " << ref_base_idx << "===" << pos <<"="<< pileupData.RefId <<"=" << current_base << "=" << x << endl;
//        string t;
////        m_idx_ref.GetSequence(pileupData.RefId, 0, pos, t);
//        cout << quit << endl;
    }


}

VariantVisitorTwo::VariantVisitorTwo(const RefVector &bam_references, const SamHeader &header, const Fasta &idx_ref,
        GenomeData &all_the_data,
//        const SampleMap &samples, BamAlignment &ali,
        int qual_cut, int mapping_cut, double prob_cut) :
        PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references),
        m_header(header), //m_samples(samples),
        m_qual_cut(qual_cut), // m_ali(ali),
        m_all_the_data(all_the_data), m_prob_cut(prob_cut),
        m_mapping_cut(mapping_cut) {

    qual_cut_char = (char) (qual_cut + 33);

    rg_tag.push_back(ZERO_CHAR);
    rg_tag += "RGZ";


    SampleMap sample_map_temp;

    total_samele_count = 0;

    SamReadGroupDictionary dictionary = m_header.ReadGroups;
    size_t tag_count = dictionary.Size();

    for (auto it = dictionary.Begin(); it != dictionary.End(); it++) {
        if (it->HasSample()) {
            auto s = sample_map_temp.find(it->Sample);
            if (s == sample_map_temp.end()) { // not in there yet
                sample_map_temp[it->Sample] = total_samele_count;
                total_samele_count++;
            }
        }
    }

//        m_samples = sample_map_temp;

    map_tag_sample = unordered_map<std::string, int>(tag_count);
    for (auto dict = dictionary.Begin(); dict != dictionary.End(); dict++) {
        map_tag_sample.emplace(dict->ID, sample_map_temp[dict->Sample]);
//        map_tag_sample_two_stage.emplace(dict->ID, dict->Sample);
    }

//        cout << "===========================" << endl;
//        cout << map_tag_sample_two_stage.load_factor() << "\t" << map_tag_sample_two_stage.max_load_factor() << "\t" << map_tag_sample_two_stage.max_size() << endl;
//        cout << map_tag_sample_two_stage.bucket_count() << "\t" << map_tag_sample_two_stage.bucket_size(1) << "\t" << map_tag_sample_two_stage.size() << endl;
////
//        exit(29);
//



};


