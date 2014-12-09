/*
 * em_algorithm.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */



#include "em_algorithm_binomial.h"
#include "em_algorithm.h"
#include "em_summary_stat_mutation.h"

#include <iostream>


EmAlgorithm::EmAlgorithm(int num_category0, std::vector <std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0) :
        num_category(num_category0), em_data_ptr(&data_ptr), em_model0(&em_model0) {



}


EmAlgorithm::EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr):
        em_data_ptr(&data_ptr), em_model_ptr(&model_ptr) {
    num_category = em_model_ptr->size();
    std::cout << "New Construct: " << num_category << std::endl;

}



void EmAlgorithm::Init() {

    site_count = em_data_ptr->size();
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);
    InitialiseProportion();

    InitialiseParameters();

    InitialiseSummaryStat();

}

void EmAlgorithm::ExpectationStep() {


    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->Reset();
        em_model[r]->UpdateParameter(parameters[r]); //exp_beta
    }

    for (size_t s = 0; s < site_count; ++s) {

        double sum_prob = 0;


        for (size_t r = 0; r < num_category; ++r) {

//            std::cout << "em_count: " << " num_category: " << r << " site: " << s << std::endl;
//            std::cout << "em_count: " << " num_category: " << r << " site: " << s << std::endl;
//em_model_ptr->at(r)
            em_data_ptr->at(s)->UpdateEmModel( em_model[r].get() );

            em_data_ptr->at(s)->UpdateSummaryStat(sum_prob, temp_stats[r]);
            //FIXME: Should make model take the data, lot's of refactor required to do this for mutation model.
            //CHECK: Should be possible to avoid static_cast

//            all_em_stats[r]->UpdateSumWithProportion(proportion[r], em_stat_local_single);

            all_probs(r, s) = sum_prob;
        }

        double sum = all_probs.col(s).sum();


        for (size_t r = 0; r < num_category; ++r) {

            double prob = all_probs(r,s) / sum;
            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
        }

    }


}

void EmAlgorithm::ExpectationStep2() {
    printf("here\n");
    vector<vector<double>> temp_stat = std::vector<std::vector<double>>(num_category);
    for (size_t r = 0; r < num_category; ++r) {
        temp_stat[r] = vector<double> (2);
        all_em_stats[r]->Reset();
        em_model_ptr->at(r)->UpdateParameter(parameters[r]); //exp_beta
    }

    for (size_t s = 0; s < site_count; ++s) {

        double sum_prob = 0;

        for (size_t r = 0; r < num_category; ++r) {
            std::cout << "em_count: " << " num_category: " << r << " site: " << s << std::endl;

//            em_model0->UpdateParameter(parameters[r]); //exp_beta

            (*em_data_ptr)[s]->UpdateEmModel(em_model_ptr->at(r));

            (*em_data_ptr)[s]->UpdateSummaryStat(sum_prob, temp_stat[r]);
            //TODO: Should make model take the data, lot's of refactor required to do this for mutation model.

//            all_em_stats[r]->UpdateSumWithProportion(proportion[r], em_stat_local_single);

            all_probs(r, s) = sum_prob;
        }

        double sum = all_probs.col(s).sum();

//        auto sum_col = all_probs.colwise().sum();
//        auto prop_col = all_probs.row(0) / sum_col;
//        double proportion_sum = prop_col.sum();
//        proportion[0] = proportion_sum / site_count;
//        proportion[1] = 1 - proportion[0];


        for (size_t r = 0; r < num_category; ++r) {

            double prob = all_probs(r,s) / sum;
//            std::cout << "===============nem prob: "<< prob << std::endl;
            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stat[r]);
        }

    }
//    CalculateProportion();


}
void EmAlgorithm::ExpectationStep_Old() {

    for (size_t r = 0; r < num_category; ++r) {
        em_model0->UpdateParameter(parameters[r]); //exp_beta

//        all_em_stats[r]->Reset();

        for (size_t s = 0; s < site_count; ++s) {
            std::cout << "em_count: "  << " num_category: " << r << " site: " << s << std:: endl;
            double sum_prob = 0;

            (*em_data_ptr)[s]->UpdateEmModel(em_model0);

            (*em_data_ptr)[s]->UpdateSummaryStat(sum_prob, em_stat_local_single);
            //TODO: Should make model take the data, lot's of refactor required to do this for mutation model.

//            all_em_stats[r]->UpdateSumWithProportion(proportion[r], em_stat_local_single);

            all_probs(r, s) = sum_prob;
        }
    }
//    CalculateProportion();

    for (size_t r = 0; r < num_category; ++r) {
        em_model0->UpdateParameter(parameters[r]); //exp_beta

//        all_em_stats[r]->Reset();

        for (size_t s = 0; s < site_count; ++s) {
            std::cout << "em_count: "  << " num_category: " << r << " site: " << s << std:: endl;
            double sum_prob = 0;

            (*em_data_ptr)[s]->UpdateEmModel(em_model0);

            (*em_data_ptr)[s]->UpdateSummaryStat(sum_prob, em_stat_local_single);
            //TODO: Should make model take the data, lot's of refactor required to do this for mutation model.

            all_em_stats[r]->UpdateSumWithProportion(proportion[r], em_stat_local_single);

//            all_probs(r, s) = sum_prob;

        }



    }
}

void EmAlgorithm::MaximizationStep() {
    CalculateProportion();


    for (size_t r = 0; r < num_category; ++r) {

        double new_parameter = all_em_stats[r]->MaximiseStats();
//        double new_mu = mutation_prob.ConvertExpBetaToMu(new_exp_beta);

//        double new_one_minus_mu = mutation_prob.ConvertExpBetaToMu(-(new_one_minus_exp_beta - 1));


        parameters[r] = new_parameter;

    }


//    printf("======= NEM_MU_r: %zu \tMu: %.5e %.5f =expBeta=  %.5f %.5f \t =prop= %.5f %.5f "
//                    "\t =stat= %.5f %.5f \n", r,
//            new_mu, new_one_minus_mu,
//            new_exp_beta, new_one_minus_exp_beta,
//            proportion[0], proportion[1],
//            all_stats_same[0], all_stats_same[1]

}

void EmAlgorithm::CalculateProportion() {

    auto sum_col = all_probs.colwise().sum();
    auto prop_col = all_probs.row(0) / sum_col;
    double proportion_sum = prop_col.sum();
    proportion[0] = proportion_sum / site_count;
    proportion[1] = 1 - proportion[0];
    std::cout << "Update Proportion: " << "\t" << proportion[0] << "\t" << proportion[1] << std::endl;


    if (num_category != 2) {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
    }

}


void EmAlgorithm::InitialiseProportion() {
    double default_proportion = 1.0 / num_category;
    proportion = std::vector<double>(num_category, default_proportion);
}

std::vector<double> EmAlgorithm::GetParameters() {
    return parameters;
}

vector<double> EmAlgorithm::GetProportion() {
    return proportion;
}