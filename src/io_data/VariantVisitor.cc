#include "VariantVisitor.h"
#include "variant_visitor_two.h"


void VariantVisitor::Visit(const LocalBamToolsUtils::PileupPosition &pileupData) {
    std::string chr = m_bam_ref[pileupData.RefId].RefName;
    uint64_t pos = pileupData.Position;
    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    ReadDataVector bcalls(m_samples.size(), ReadData{0});
    std::string tag_id;

//    std::cout << "P:" << pileupData.PileupAlignments.size() << std::endl;

    for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {

        if (include_site(*it, m_mapping_cut, m_qual_cut)) {
            it->Alignment.GetTag("RG", tag_id);
            std::string sm = m_header.ReadGroups[tag_id].Sample;
            uint32_t sindex = m_samples[sm]; //TODO check samples existed!
            uint16_t bindex = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
            if (bindex < 4) {
                global_count[0]++;
                bcalls[sindex].reads[bindex] += 1;
            }
        }

    }


    uint16_t ref_base_idx = base_index(current_base);
    if (ref_base_idx < 4) { //TODO Model for bases at which reference is 'N'
        global_count[2]++;
            m_all_the_data.push_back(ModelInput{0, ref_base_idx, bcalls});
    }
    else{
        std::cout << "ELSE" << std::endl;
    }
}
