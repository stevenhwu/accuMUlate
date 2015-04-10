/*
 * em_model_binary.cc.h
 *
 *  Created on: 12/6/14
 *      Author: Steven Wu
 */



#include <iostream>
#include "em_model_binomial.h"



EmModelBinomial::EmModelBinomial(int n0, double prob0) : n(n0), prob(prob0){
}


void EmModelBinomial::UpdateSummaryStat(int site_index, double &prob, std::vector<double> &temp_stat, double &log_likelihood_scaler) {
    std::cout << "Error!! Not yet implemented" << std::endl;
    exit(44);
}

void EmModelBinomial::UpdateParameter(double param) {
    prob = param;
}


size_t EmModelBinomial::GetDataCount() {
    return 0;
}

double EmModelBinomial::GetParameter() {
    return prob;
}

//void EmModelBinomial::GetParameterInfo() {
//
//}
