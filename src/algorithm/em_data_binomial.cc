/*
 * em_data_binary.h
 *
 *  Created on: 12/6/14
 *      Author: Steven Wu
 */


#include "em_data_binomial.h"
#include "em_model_mutation.h"


EmDataBinomial::EmDataBinomial(std::vector<int> input){
    total_count = input.size();
    data = input;
    int count  = 0;
    for (auto item : input) {
        if(item==true){
            count ++;
        }
    }
}


EmDataBinomial::EmDataBinomial(int total0, int count0){
    total_count = total0;
    count = count0;
    count_negative = total_count - count;
}


EmDataBinomial::~EmDataBinomial() {}


void EmDataBinomial::UpdateSummaryStat(double &prob, std::unique_ptr<EmSummaryStat> &summaryStat) {


//    summaryStat.print();
    prob = 0;

    //relpace with
//    EmSummaryStatMutation local;

    prob = pow(binomial_prob, count)*pow((1-binomial_prob), count_negative);
//    std:: cout << "COUNT:" << count << "\t" << binomial_prob << "\t" << count_negative << std::endl;
//    double p1 = pow( binomial_prob, count);
//    double p2 = pow( (1-binomial_prob), count_negative);;
//    std::cout << p1 << "\t" << p2 << "\t\n";
//    exit(-100);
//    site.CalculateAncestorToDescendant(prob, stat_same, stat_diff);
//    double stat_same = prob * count;
//    double stat_diff = prob * count_negative;
    double stat_same = count;
    double stat_diff = count_negative;
    summaryStat->SetStats(vector<double> {stat_same, stat_diff});

    std::cout << "call EMDATA_Binomial: " << count << "\t" << binomial_prob << "\t" <<
        stat_same << "\t" << stat_diff << "\t" << prob <<std::endl;
//    site->UpdateSummaryStat(sum_prob, local);
//    site.CalculateAncestorToDescendant(sum_prob, stat_same, stat_diff);
}


void EmDataBinomial::UpdateEmModel(EmModel *em_model) {


//    EmModelMutation *em = dynamic_cast<EmModelMutation*> (em_model);
    EmModelBinomial *em_model_binomial = static_cast<EmModelBinomial*> (em_model);
    binomial_prob = em_model_binomial->GetParameter();
    cout << "In EmDataBinomial, UpdateModel: " << binomial_prob << endl;
//    site.UpdateModel(*evo_model);
    //TODO: redo static_cast with visitor??? ?typeid?

}

void EmDataBinomial::UpdateEmModel(std::unique_ptr<EmModelBinomial> &em_model) {

}

void EmDataBinomial::UpdateEmModel(std::unique_ptr<EmModelMutation> &em_model) {

}
