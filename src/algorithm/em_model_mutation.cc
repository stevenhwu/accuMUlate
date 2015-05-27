#include <iostream>
#include "em_model_mutation.h"

EmModelMutation::EmModelMutation(MutationModel &model) : mutation_model(model) {
}


void EmModelMutation::UpdateParameter(double param) {
    mutation_model.UpdateOneMinusExpBeta(param);
}


void EmModelMutation::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler) {

    prob = 0;
    log_likelihood_scaler = 0;
    double stat_diff = 0;

    mutation_model.CalculateAncestorToDescendant(site_index, prob, stat_diff, log_likelihood_scaler);
    temp_stat[0] = 1-stat_diff;
    temp_stat[1] = stat_diff;

}


size_t EmModelMutation::GetDataCount() {
    return mutation_model.GetSiteCount();
}


void EmModelMutation::GetParameterInfo(){
    std::cout << "Error!! Not yet implemented" << std::endl;
    std::exit(44);
}
