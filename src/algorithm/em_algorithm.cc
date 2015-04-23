/*
 * em_algorithm.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */


#include "em_algorithm_mutation.h"
#include <iostream>
#include "em_algorithm_mutation_v1.h"
#include "em_algorithm_binomial.h"
#include "em_algorithm.h"
#include "em_summary_stat_mutation.h"

#include <stddef.h>
#include <glob.h>

const double EM_CONVERGE_THRESHOLD = 1e-10;
const double EM_CONVERGE_RATIO_THRESHOLD = 1e-10;
const int EM_MAX_ITE = 50;
int VERBOSE_ITE = 100;
int LOG_ITE = 10;

EmAlgorithm::EmAlgorithm(MutationModelMultiCategories &model_multi0) : model_multi(model_multi0){
    num_category = model_multi.GetCategoriesCount();
    std::cout << "ModelMultiCategories Construct: " << num_category << "\t One data, one model" << std::endl;

}


EmAlgorithm::EmAlgorithm(int num_category0, std::vector <std::unique_ptr<EmData>> &data_ptr, EmModel &em_model0) :
        num_category(num_category0), em_data_ptr(&data_ptr), em_model0(&em_model0) {

    std::cout << "Old Construct: " << num_category << "\t One data, one model" << std::endl;
}


EmAlgorithm::EmAlgorithm(std::vector<std::unique_ptr<EmData>> &data_ptr, std::vector<std::unique_ptr<EmModel>> &model_ptr):
        em_data_ptr(&data_ptr), em_model_ptr(&model_ptr) {
    num_category = em_model_ptr->size();

    std::cout << "New Construct: " << num_category << std::endl;
}


EmAlgorithm::EmAlgorithm(std::vector<std::unique_ptr<EmModel>> &model_ptr):
        em_model_ptr(&model_ptr) {
    num_category = em_model_ptr->size();

    std::cout << "New Construct model only: " << num_category << std::endl;
}

void EmAlgorithm::InitWithModel() {

    if (num_category != 2) {
        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
        exit(222);
    }

    std::vector<int> temp_data_count(num_category);
    for (size_t i = 0; i < num_category; ++i) {
        temp_data_count[i] = (*em_model_ptr)[i]->GetDataCount();
    }
    for (size_t i = 1; i < num_category; ++i) {
        if(temp_data_count[i] != temp_data_count[i-1]){
            std::cout << "ERROR: Data count not equal\n" << temp_data_count[i] << "\t" << temp_data_count[i-1] << std::endl;
            exit(210);
        }
    }
    site_count = temp_data_count[0];
    std::cout << "Site_count: " << site_count << std::endl;
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);
    parameters = std::vector<double>(num_category);
    cache_parameters = std::vector<double>(num_category, 0);

    InitialiseProportion();

    InitialiseParameters();

    InitialiseSummaryStat();

}


void EmAlgorithm::InitWithData() {

    if (num_category != 2) {
        std::cout << "Not yet implemented for more than 2 categories: input_cat: " << "\t" << num_category<< std::endl;
        exit(222);
    }

    site_count = em_data_ptr->size();
    std::cout << "Site_count: " << site_count << std::endl;
    all_probs = Eigen::ArrayXXd::Zero(num_category, site_count);
    parameters = std::vector<double>(num_category);
    cache_parameters = std::vector<double>(num_category, 0);

    InitialiseProportion();

    InitialiseParameters();

    InitialiseSummaryStat();

}

void EmAlgorithm::ExpectationStepModel() {


    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->Reset();
        em_model[r]->UpdateParameter(parameters[r]); //exp_beta
    }
    double total = 0;
    for (size_t s = 0; s < site_count; ++s) {
        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {

            ExpectationStepCustom(s, r, sum_prob, temp_stats[r]);

            //FIXME: Should make model take the data, lot's of refactor required to do this for mutation model.
            //CHECK: Should be possible to avoid static_cast
            all_probs(r, s) = proportion[r] * sum_prob;
        }

        double sum = all_probs.col(s).sum();
        double total_site = 0;
        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;
            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
            total_site += all_probs(r,s);
        }
        total += log(total_site);

    }
    fflush(stdout);
}



void EmAlgorithm::ExpectationStepModelPtr() {


    UpdateEmParameters();
    double log_likelihood_scaler = 0;
    log_likelihood = 0;
    for (size_t s = 0; s < site_count; ++s) {
        double sum_prob = 0;
        for (size_t r = 0; r < num_category; ++r) {
            (*em_model_ptr)[r]->UpdateSummaryStat(s, sum_prob, temp_stats[r], log_likelihood_scaler);
            all_probs(r, s) = proportion[r] * sum_prob;
        }
        log_likelihood = log_likelihood + log(all_probs(0, s)+all_probs(1, s))+log_likelihood_scaler;
        double sum = all_probs.col(s).sum();

        for (size_t r = 0; r < num_category; ++r) {
            double prob = all_probs(r,s) / sum;
            all_em_stats[r]->UpdateSumWithProportion(prob, temp_stats[r]);
        }
    }
//    for (size_t r = 0; r < num_category; ++r) {
//        all_em_stats[r]->Print();
//    }
//    std::cout << "==Ln: " << log_likelihood << std::endl;
//    printf("Ln: %.40f\n", log_likelihood);

}

void EmAlgorithm::UpdateEmParameters() {
    for (size_t r = 0; r < num_category; ++r) {
        all_em_stats[r]->Reset();
        em_model_ptr->at(r)->UpdateParameter(parameters[r]); //exp_beta
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

    if (num_category != 2) {
        std::cout << "Not yet implemented for more than 2 categories" << std::endl;
        exit(222);
    }

    auto sum_col = all_probs.colwise().sum();
    auto prop_col = all_probs.row(0) / sum_col;
    double proportion_sum = prop_col.sum();
    proportion[0] = proportion_sum / site_count;
    proportion[1] = 1 - proportion[0];

//    std::cout << sum_col << std::endl;
//    std::cout << prop_col << std::endl;
//    std::cout << proportion_sum << std::endl;
//    std::cout << "Update Proportion: " << "\t" << proportion[0] << "\t" << proportion[1] << std::endl;

}


void EmAlgorithm::InitialiseProportion() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> lower(0.8,1);
    double lower_proportion = lower(gen);

    double default_proportion = 1.0 / num_category;
    lower_proportion = default_proportion;

    proportion = {{1-lower_proportion, lower_proportion}};


}

std::vector<double> EmAlgorithm::GetParameters() {
    return parameters;
}

std::vector<double> EmAlgorithm::GetProportion() {
    return proportion;
}

bool EmAlgorithm::EmStoppingCriteria(int ite) {

    double sum_diff = 0;

    sum_ratio = 0;
    for (size_t r = 0; r < num_category; ++r) {
        double diff = fabs(parameters[r] - cache_parameters[r]);
        sum_diff += diff;
        sum_ratio += (diff / cache_parameters[r]);
        cache_parameters[r] = parameters[r];
    }

    if ( (ite % VERBOSE_ITE) == 0) {
        std::cout << "Ite: " << ite << " sum_diff: " << sum_diff << "\tsum_ratio: " << sum_ratio << std::endl;
        PrintSummary();
    }

    if ( (ite % LOG_ITE) == 0) {
//        std::cout << "Ite: " << ite << " sum_diff: " << sum_diff << "\tsum_ratio: " << sum_ratio << std::endl;
//        PrintSummary();
        LogEmSummary(ite);
    }

    if (ite == EM_MAX_ITE){
        std::cout <<"============ DONE. IN DEBUG mode, fix at " << EM_MAX_ITE << " ites ======= " << sum_diff << " Total ite:" << ite << "\n";
        return false;
    }
//    if (sum_diff < EM_CONVERGE_THRESHOLD) {
//        std::cout <<"============ DONE (diff ~~ 0) ======= Diff:" << sum_diff << "\tRatio:" << sum_ratio  <<" Total ite:" << ite << "\n";
//        return false;
//    }
    if( sum_ratio < EM_CONVERGE_RATIO_THRESHOLD){
        std::cout <<"============ DONE (ratio < THRESHOLD) ======= Diff: " << sum_diff << "\tRatio:" << sum_ratio << " Total ite:" << ite << "\n";
        return true;
        return false;
    }

    if( std::isnan(sum_ratio) ){
//        std::cout <<"============ FAIL (ratio == -nan) ======= Diff: " << sum_diff << "\tRatio:" << sum_ratio << " Total ite:" << ite << "\n";
        return true;
        return false;
    }
    return true;
}

void EmAlgorithm::PrintSummary(){
//    std::cout << "Ite: " << ite << " sum_diff: " << sum_diff << "\tsum_ratio: " << sum_ratio << std::endl;
    printf("EM Summary: Ln:%.10f\nParameters: ", log_likelihood.load()); //TODO: Add this back
     for (auto &item :parameters) {
        printf("%.15e\t", item);
    }

    printf("\nProportions: ");
    for (auto &item :proportion) {
        printf("%.15e\t", item);
    }

    printf("\n================= END SUMMARY ============\n");
    fflush(stdout);
}

std::string EmAlgorithm::GetEMSummary(){
    std::string out;// = "EM Summary\tParameters0 ";
    char temp[1000];
    for (auto &item :parameters) {
        sprintf(temp, "%.12e\t", item);
        out.append(temp);
    }
    for (auto &item :proportion) {
        sprintf(temp, "%.12e\t", item);
        out.append(temp);
    }
//    sprintf(temp, "%.12f\t%.12e\t", log_likelihood, sum_ratio);//TODO: add this back
    out.append(temp);
    return out;

}

void EmAlgorithm::SetOutfilePrefix(const std::string & infile) {
    em_logger.SetOutifle(infile);

    std::string header;
    char temp[1000];
    header.append("Ite\t");
    for (int i = 0; i < num_category; ++i) {
        sprintf(temp, "Parameter_%d\t", i);
        header.append(temp);
    }
    for (int i = 0; i < num_category; ++i) {
        sprintf(temp, "Proportion_%d\t", i);
        header.append(temp);
    }
    header.append("Likelihood\tRatio\t");
    em_logger.Header(header);

}

void EmAlgorithm::Run() {

    std::cout << "Init ";
    PrintSummary();
    RunEM();
    em_logger.Stop();

}

void EmAlgorithm::LogEmSummary(int ite) {

    std::string summary = GetEMSummary();
    em_logger.LogLine(ite, summary);
}

