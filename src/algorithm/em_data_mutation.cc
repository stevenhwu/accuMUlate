

#include <iostream>
#include "em_data_mutation.h"
#include "em_model_mutation.h"


EmDataMutation::EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}


EmDataMutation::~EmDataMutation() {}


void EmDataMutation::UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) {


    summaryStat->print();
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


void EmDataMutation::UpdateEmModel(EmModel *em_model) {

//    cout << "In EmDataMutation, UpdateModel: " << endl;
//    EmModelMutation *em = dynamic_cast<EmModelMutation*> (em_model);
    EmModelMutation *em_model_mutation = static_cast<EmModelMutation*> (em_model);
    EvolutionModel *evo_model = em_model_mutation->GetEvoModel();
    site.UpdateModel(*evo_model);
    //TODO: redo static_cast with visitor??? ?typeid?

}
//
//void EmDataMutation::UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model) {
//
////    default_delete<EmModel> &deleter = em_model.get_deleter();
////    std::unique_ptr<EmDataMutation> em_model_mutation = static_unique_ptr_cast<EmModelMutation, EmModel, default_delete<EmModel> >(em_model);
//    EvolutionModel *evo_model = em_model->GetEvoModel();
//    site.UpdateModel(*evo_model);
//}
//
//
//void EmDataMutation::UpdateEmModel(unique_ptr<EmModel> &em_model) {
//    std::cout << "deal with unqiue pointer at update\n";
//    EmModelMutation *em_model_mutation = static_cast<EmModelMutation*> (em_model.get());
////    EmModelMutation *em_model_mutation2 = em_model->GetModel();
//    EvolutionModel *evo_model = em_model_mutation->GetEvoModel();
//    site.UpdateModel(*evo_model);
//}

void EmDataMutation::UpdateSummaryStat(double &prob, std::vector<double> &temp_stat) {


//    summaryStat->print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;
    //relpace with
//    EmSummaryStatMutation local;


    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);
    temp_stat[0] = stat_same;
    temp_stat[1] = stat_diff;
//    cout << "update with vector " << temp_stat[0] << "\t" << temp_stat[1] << std::endl;

//    summaryStat->SetStats(vector<double> {stat_same, stat_diff});
}
