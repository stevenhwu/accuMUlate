

#include <iostream>
#include "em_data_mutation.h"
#include "em_model_mutation.h"


EmDataMutation::EmDataMutation(SequenceProb &sequence_prob, EvolutionModel &evo_model){
    site = SiteProb(sequence_prob, evo_model);
}


EmDataMutation::~EmDataMutation() {}


void EmDataMutation::UpdateSummaryStat(double &prob, EmSummaryStat &summaryStat) {


    summaryStat.print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;
    //relpace with
//    EmSummaryStatMutation local;


    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);

    summaryStat.SetStats(vector<double> {stat_same, stat_diff});

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
    //FIXME: redo static_cast with visitor??? ?typeid?

}

void EmDataMutation::Test(double num) {
    cout << "In EmDataMutation, testing function: " << num << endl;
}
