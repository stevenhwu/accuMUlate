#include <stddef.h>

/*
 * em_model_binary.cc.h
 *
 *  Created on: 12/6/14
 *      Author: Steven Wu
 */



#include <iostream>
#include "em_model_binomial.h"


size_t EmModelBinomial::GetDataCount() {
    return 0;
}

EmModelBinomial::EmModelBinomial(int n0, double prob0) : n(n0), prob(prob0){

//    MutationRate rate = evo_model->GetMutationRate();
//    cout << rate.prob << "\t" << rate.one_minus_p << endl;

}


void EmModelBinomial::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat) {
    std::cout << "Error!! Not yet implemented" << std::endl;
    exit(41);
}

void EmModelBinomial::UpdateParameter(double param) {
//    evo_model->UpdateMu(param);
//    evo_model->UpdateExpBeta(param);
//    MutationRate mutation_rate = evo_model->GetMutationRate();
    prob = param;
//    std::cout << "IN EmModelBinomial: updateing: " << param << "\t" << prob<< std::endl;
}


double EmModelBinomial::GetParameter() {
    return prob;
}

void EmModelBinomial::GetParameterInfo() {

}
