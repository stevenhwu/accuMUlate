

#include <iostream>
#include "em_data_mutation_v1.h"


EmDataMutationV1::EmDataMutationV1(SequenceProb &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}


EmDataMutationV1::EmDataMutationV1(SiteGenotypes &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}

EmDataMutationV1::~EmDataMutationV1() {}


void EmDataMutationV1::UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) {

    summaryStat->Print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;
    //relpace with
//    EmSummaryStatMutation local;

    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);

    summaryStat->SetStats(std::vector<double> {stat_same, stat_diff});

//    std::cout << "call EMDATA2: " << stat_same << "\t" << stat_diff << "\t" << prob <<std::endl;
//    site->UpdateSummaryStat(sum_prob, local);
//    site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);
}


void EmDataMutationV1::UpdateEmModel(EmModel *em_model) {

//    cout << "In EmDataMutationV1, UpdateModel: " << endl;
//    EmModelMutationV1 *em = dynamic_cast<EmModelMutationV1*> (em_model);
    EmModelMutationV1 *em_model_mutation = static_cast<EmModelMutationV1 *> (em_model);
    EvolutionModel *evo_model = em_model_mutation->GetEvoModel();
    site.UpdateModel(*evo_model);
    //TODO: redo static_cast with visitor??? ?typeid?

}
//
//void EmDataMutationV1::UpdateEmModel(std::unique_ptr<EmModelMutationV1> &em_model) {
//
////    default_delete<EmModel> &deleter = em_model.get_deleter();
////    std::unique_ptr<EmDataMutationV1> em_model_mutation = static_unique_ptr_cast<EmModelMutationV1, EmModel, default_delete<EmModel> >(em_model);
//    EvolutionModel *evo_model = em_model->GetEvoModel();
//    site.UpdateModel(*evo_model);
//}
//
//
//void EmDataMutationV1::UpdateEmModel(unique_ptr<EmModel> &em_model) {
//    std::cout << "deal with unqiue pointer at update\n";
//    EmModelMutationV1 *em_model_mutation = static_cast<EmModelMutationV1*> (em_model.get());
////    EmModelMutationV1 *em_model_mutation2 = em_model->GetModel();
//    EvolutionModel *evo_model = em_model_mutation->GetEvoModel();
//    site.UpdateModel(*evo_model);
//}

void EmDataMutationV1::UpdateSummaryStat(double &prob, std::vector<double> &temp_stat) {

//    summaryStat->Print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;

    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);
    temp_stat[0] = stat_same;
    temp_stat[1] = stat_diff;
//    cout << "update with vector " << temp_stat[0] << "\t" << temp_stat[1] << std::endl;


}
