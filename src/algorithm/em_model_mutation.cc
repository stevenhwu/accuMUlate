#include <iostream>
#include "em_model_mutation.h"

EmModelMutation::EmModelMutation(MutationModel &model) : mutation_model(model) {
//    mutation_model
}


void EmModelMutation::UpdateParameter(double param) {

//    MutationRate mutation_rate = evo_model->GetMutationRate();
    std::cout << "IN EmModelMutation: updateing: " << param << std::endl;
    mutation_model.UpdateExpBeta(param);


}


void EmModelMutation::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat) {


//    summaryStat->print();
    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;
    //relpace with
//    EmSummaryStatMutationV1 local;

    mutation_model.CalculateAncestorToDescendant(site_index, prob, stat_diff);
//    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);
    temp_stat[0] = 1-stat_diff;
    temp_stat[1] = stat_diff;
    std::cout << "update with vector " << temp_stat[0] << "\t" << temp_stat[1] << std::endl;

//    summaryStat->SetStats(vector<double> {stat_same, stat_diff});
}


void EmModelMutation::GetParameterInfo(){


}



size_t EmModelMutation::GetDataCount() {
    return mutation_model.GetSiteCount();
}
