#include <iostream>
#include "em_model_mutation.h"

EmModelMutation::EmModelMutation(MutationModel &model) : mutation_model(model) {
}


void EmModelMutation::UpdateParameter(double param) {

//    std::cout << "IN EmModelMutation: updateing: " << param << std::endl;
    mutation_model.UpdateExpBeta(param);


}


void EmModelMutation::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat) {

    prob = 0;
    double stat_same = 0;
    double stat_diff = 0;

    mutation_model.CalculateAncestorToDescendant(site_index, prob, stat_diff);
    temp_stat[0] = 1-stat_diff;
    temp_stat[1] = stat_diff;

}


void EmModelMutation::GetParameterInfo(){
    std::cout << "Error!! Not yet implemented" << std::endl;
    exit(41);
}



size_t EmModelMutation::GetDataCount() {
    return mutation_model.GetSiteCount();
}
