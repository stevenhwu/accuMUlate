/*
 * variant_visitor_two.cc
 *
 *  Created on: 12/9/14
 *      Author: Steven Wu
 */



#include "variant_visitor_two.h"


int global_count[10];

void VariantVisitorTwo::Visit(const LocalBamToolsUtils::PileupPosition &pileupData) {

//    string chr = m_bam_ref[pileupData.RefId].RefName;
    uint64_t pos = pileupData.Position;

    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);//Check: Room for improvement here?

    uint16_t ref_base_idx = base_index2[(int) current_base];
    if (ref_base_idx < 4) { //TODO Model for bases at which reference is 'N'
        global_count[1]++;
        ReadDataVector bcalls(total_sample_count, ReadData{0});
        
        for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {//auto it2 = *it;

            int32_t pos_in_alignment = it->PositionInAlignment;
            if (include_site_4(it->Alignment, pos_in_alignment, m_mapping_cut, qual_cut_char)) {
                int sindex = GetSampleIndex(it->Alignment.TagData);

                uint16_t bindex = base_index2[(int) it->Alignment.QueryBases[pos_in_alignment]];

                if (bindex < 4) {
                    global_count[0]++;
                    bcalls[sindex].reads[bindex] += 1;
                }
            }
//13309742:m35
//14713917:m25
//15621555:m15
        }
        if( filter_data(bcalls) ) {
            global_count[2]++;
            m_all_the_data.push_back(ModelInput{0, ref_base_idx, bcalls});
            gd_stream.WriteModelInput(m_all_the_data.back());
        }


    }
    else {
        std::cout << "Invalid ref_base_idx: " << ref_base_idx << "\t" << current_base << std::endl;
    }


}

int VariantVisitorTwo::GetSampleIndex(std::string const &tag_data) {
    size_t start_index = tag_data.find(rg_tag);
    if (start_index != std::string::npos) {
                    start_index += 4;
                }
                else {
                    size_t x = tag_data.find("RG");//TODO: Check for (char)0, RG, Z
                    std::cout << "ERRER: " << tag_data << "\t" << x << std::endl;
                }
    size_t end_index = tag_data.find(ZERO_CHAR, start_index);
    std::string tag_id_2 = tag_data.substr(start_index, (end_index - start_index));
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
    return sindex;
}

VariantVisitorTwo::VariantVisitorTwo(const RefVector &bam_references, const SamHeader &header,
        const LocalBamToolsUtils::Fasta &idx_ref, GenomeData &all_the_data, std::string binary_outfile, int qual_cut, int mapping_cut, double prob_cut)
        :
        PileupVisitor(),  m_bam_ref(bam_references),
        m_header(header), m_idx_ref(idx_ref),
        m_all_the_data(all_the_data),
        m_qual_cut(qual_cut),
        m_mapping_cut(mapping_cut),
        m_prob_cut(prob_cut) {


    std::cout << "Init VariantVisitorTwo: " << m_qual_cut << "\t" << m_mapping_cut << "\t" << m_prob_cut << std::endl;
    qual_cut_char = (char) (m_qual_cut + 33);

    rg_tag.push_back(ZERO_CHAR);
    rg_tag += "RGZ";

    SampleMap sample_map_temp;
    total_sample_count = 0;
    SamReadGroupDictionary dictionary = m_header.ReadGroups;
    size_t tag_count = dictionary.Size();

    for (auto it = dictionary.Begin(); it != dictionary.End(); it++) {
        if (it->HasSample()) {
            auto s = sample_map_temp.find(it->Sample);
            if (s == sample_map_temp.end()) { // not in there yet
                sample_map_temp[it->Sample] = total_sample_count;
                total_sample_count++;
            }
        }
    }

//        m_samples = sample_map_temp;

    map_tag_sample = std::unordered_map<std::string, int>(tag_count);
    for (auto dict = dictionary.Begin(); dict != dictionary.End(); dict++) {
        map_tag_sample.emplace(dict->ID, sample_map_temp[dict->Sample]);
//        map_tag_sample_two_stage.emplace(dict->ID, dict->Sample);
    }
    for (auto sample : map_tag_sample) {
        std::cout << sample.first << "\t" << sample.second << std::endl;
    }

//        cout << "===========================" << endl;
//        cout << map_tag_sample_two_stage.load_factor() << "\t" << map_tag_sample_two_stage.max_load_factor() << "\t" << map_tag_sample_two_stage.max_size() << endl;
//        cout << map_tag_sample_two_stage.bucket_count() << "\t" << map_tag_sample_two_stage.bucket_size(1) << "\t" << map_tag_sample_two_stage.size() << endl;
////
//        exit(29);
//

    gd_stream = GenomeDataStream(binary_outfile, true);
    gd_stream.WriteHeader(total_sample_count);



};


bool VariantVisitorTwo::filter_data(ReadDataVector &read_vector) {
    return true;
    for (auto item : read_vector) {
        if (item.key==0){
            return false;
        }
    }
    return true;
    int pass_count = read_vector.size();
    for (auto item : read_vector) {
        if (item.key!=0){
            for (auto read : item.reads) {
                if(read>1){ //Simple filter, at least 1 base > than N==1 CHECK:
                    pass_count--;
                    break;
                }
            }
        }
    }

    if(pass_count < 3 ){// ALL samples must exist?? OR we allow some gaps. Use 3 (2 samples missing) for now?? CHECK:
//        std::cout << "PASS " << pass_count << std::endl;
        return true;
    }
//    std::cout << "FAIL " << pass_count << std::endl;
    return false;

}

int VariantVisitorTwo::getTotal_sample_count() const {
    return total_sample_count;
}
