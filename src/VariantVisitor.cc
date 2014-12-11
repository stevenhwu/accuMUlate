#include "VariantVisitor.h"
#include "variant_visitor_two.h"


void VariantVisitor::Visit(const PileupPosition &pileupData) {
    string chr = m_bam_ref[pileupData.RefId].RefName;
    uint64_t pos = pileupData.Position;
    m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
    ReadDataVector bcalls(m_samples.size(), ReadData{{0, 0, 0, 0}});
    string tag_id;
    global_count2++;
    cout << "P:" << pileupData.PileupAlignments.size() << endl;
    for (auto it = begin(pileupData.PileupAlignments); it != end(pileupData.PileupAlignments); ++it) {
        global_count3++;
        if (include_site(*it, m_mapping_cut, m_qual_cut)) {
            it->Alignment.GetTag("RG", tag_id);
            string sm = m_header.ReadGroups[tag_id].Sample;
            uint32_t sindex = m_samples[sm]; //TODO check samples existed!
            uint16_t bindex = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
            cout << "B: " << bindex << "\t" << sindex << endl;
            if (bindex < 4) {
                global_count++;
                bcalls[sindex].reads[bindex] += 1;
            }
        }
    }


    uint16_t ref_base_idx = base_index(current_base);
    if (ref_base_idx < 4) { //TODO Model for bases at which reference is 'N'
//                ModelInput d = {ref_base_idx, bcalls};

        m_all_the_data.push_back(ModelInput{ref_base_idx, bcalls});
//                ModelInput d = {ref_base_idx, bcalls};
//                double prob_one = TetMAProbOneMutation(m_params,d);
//                double prob = TetMAProbability(m_params, d);
//                if(prob >= m_prob_cut){
//                     *m_ostream << chr << '\t'
//                                << pos << '\t'
//                                << current_base << '\t'
//                                << prob << '\t'
//                                << prob_one << '\t'
//                               << endl;
//                }
    }
    else{
        cout << "ELSE" << endl;
    }
}