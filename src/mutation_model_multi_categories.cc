/*
 * mutation_model_multi_categories.cc
 *
 *  Created on: 4/10/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/10/15.
//

#include "mutation_model_multi_categories.h"

#include <strings.h>



////std::vector<SiteGenotypesIndex> MutationModelMultiCategories::a = std::vector<SiteGenotypesIndex>();
//std::vector<SiteGenotypesIndex> MutationModelMultiCategories::all_sequence_prob_index;
//
//std::vector<HaploidProbs> MutationModelMultiCategories::convert_index_key_to_haploid;
//std::vector<DiploidProbsIndex10> MutationModelMultiCategories::convert_index_key_to_diploid_10;
//
//std::vector<double> MutationModelMultiCategories::convert_index_key_to_haploid_scaler;
//std::vector<double> MutationModelMultiCategories::convert_index_key_to_diploid_10_scaler;
//
//std::vector<HaploidProbs> MutationModelMultiCategories::convert_index_key_to_haploid_unnormalised = std::vector<HaploidProbs> ();
//std::vector<DiploidProbsIndex10> MutationModelMultiCategories::convert_index_key_to_diploid_10_unnormalised;
//
////    std::vector<DiploidProbs> MutationModelMultiCategories::convert_index_key_to_diploid;
//    std::array<DiploidProbs, 4> MutationModelMultiCategories::ref_diploid_probs;


MutationModelMultiCategories::MutationModelMultiCategories(int num_categories, EvolutionModel &evo_model0)
        :categories_count(num_categories){

    MutationProb mutation_prob = evo_model0.GetMutationProb();
    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();
    evo_model = &evo_model0;

    transition_matrix_a_to_d_single = evo_model->GetTranstionMatirxAToD();
    mutation_rate_single = evo_model->GetMutationRate();

}

MutationModelMultiCategories::MutationModelMultiCategories(SequencingFactory &factory) {

    convert_index_key_to_haploid = factory.RemoveConvertIndexKeyToHaploid();
    convert_index_key_to_diploid_10 = factory.RemoveConvertIndexKeyToDiploidIndex10();

    convert_index_key_to_haploid_scaler = factory.RemoveConvertIndexKeyToHaploidScaler();
    convert_index_key_to_diploid_10_scaler = factory.RemoveConvertIndexKeyToDiploidIndex10Scaler();

//    all_sequence_prob_index = std::move(all);
    all_sequence_prob_index = factory.MoveSiteGenotypeIndexVector();
    site_count = all_sequence_prob_index.size();
    descendant_count = all_sequence_prob_index[0].GetDescendantCount();
    std::cout << "site_count: " << "\t" << site_count << "\tDescendant_count: " << descendant_count  << std::endl;
    std::cout << "Assuming all data have the same number of descendants. If not, reimplement this!!." << std::endl;

//    MutationModelMultiCategories::convert_index_key_to_haploid_unnormalised = factory.RemoveConvertIndexKeyToHaploidUnnormalised();
//    MutationModelMultiCategories::convert_index_key_to_diploid_10_unnormalised = factory.RemoveConvertIndexKeyToDiploidIndex10Unnormalised();



//    MutationModelMultiCategories::convert_index_key_to_haploid = factory.RemoveConvertIndexKeyToHaploidUnnormalised();
//    MutationModelMultiCategories::convert_index_key_to_diploid_10 = factory.RemoveConvertIndexKeyToDiploidIndex10Unnormalised();

//    MutationModelMultiCategories::convert_index_key_to_haploid_unnormalised = factory.RemoveConvertIndexKeyToHaploid();
//    MutationModelMultiCategories::convert_index_key_to_diploid_10_unnormalised = factory.RemoveConvertIndexKeyToDiploidIndex10();

}


int MutationModelMultiCategories::GetCategoriesCount() const {
    return categories_count;
}

int MutationModelMultiCategories::GetSiteCount() const {
    return site_count;
}


void MutationModelMultiCategories::MoveSequenceProb(std::vector<SiteGenotypesIndex> &&all) {
//    std::vector<SiteGenotypesIndex> &local = all;
    all_sequence_prob_index = std::move(all);
    site_count = all_sequence_prob_index.size();
    descendant_count = all_sequence_prob_index[0].GetDescendantCount();
    std::cout << "site_count: " << "\t" << site_count << "\tDescendant_count: " << descendant_count  << std::endl;
    std::cout << "Assuming all data have the same number of descendants. If not, reimplement this!!." << std::endl;


//    InitCache();//TODO: call this somewhere elseS
}


void MutationModelMultiCategories::UpdateExpBeta(int category_index, double expBeta) {

    evo_model->UpdateExpBeta(expBeta);

    transition_matrix_a_to_d_single = evo_model->GetTranstionMatirxAToD();
    mutation_rate_single = evo_model->GetMutationRate();

    UpdateCache(category_index);

}
//void MutationModelMultiCategories::SummaryIndexToHaploid() {
//    std::cout << "sum index2Hap" << std::endl;
//    for (int i = 0; i < 5; ++i) {
//        std::cout << convert_index_key_to_haploid[i].format(nice_row) << std::endl;
//    }
//    std::cout << convert_index_key_to_haploid[convert_index_key_to_haploid.size()-1].format(nice_row) << std::endl;
//    std::cout << convert_index_key_to_haploid[convert_index_key_to_haploid.size()-2].format(nice_row) << std::endl;
//}

void MutationModelMultiCategories::InitCache(int category_index) {
    int index_size = convert_index_key_to_haploid.size();
//    std::cout << map_rd_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << "\t" << index_size << std::endl;
//    std::cout << "INDEX size: " << convert_index_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << "\t" << map_rd_key_to_haploid.size() << "\t" <<
//            index_size << std::endl;

    for (int k = 0; k < 10; ++k) {
        std::vector<std::pair<double, double>> temp;
        temp.assign(index_size, std::make_pair(0,0));
        cache_read_data_to_all_index_rev[k] = temp;
    }

    UpdateCache(category_index);

}


void MutationModelMultiCategories::UpdateCache(int category_index) {


    std::array<double, 4> temp_base_prob;
    for (int b = 0; b < BASE_COUNT; ++b) {
        temp_base_prob[b] = mutation_rate_single * frequency_prior[b];
    }



    for (size_t i = 0; i < convert_index_key_to_haploid.size(); ++i) {
        auto &genotype = convert_index_key_to_haploid[i];

        double summary_stat_diff= 0;
        for (int b = 0; b < BASE_COUNT; ++b) {
            summary_stat_diff += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
        }


        for (int k = 0; k < ANCESTOR_COUNT; ++k) {
            int index16 = LookupTable::index_converter_10_to_16[k];
            double prob_reads_d_given_a = 0;

            for (int b = 0; b < BASE_COUNT; ++b) {
                double prob = transition_matrix_a_to_d_single(index16, b) * genotype[b];
                prob_reads_d_given_a += prob;
//                //stat_same += genotype[b] * evo_model.GetMutationProb().mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[k][b];
            }

            cache_read_data_to_all_index_rev[k][i].first = prob_reads_d_given_a;
            cache_read_data_to_all_index_rev[k][i].second = summary_stat_diff / prob_reads_d_given_a;
        }
    }
}


void MutationModelMultiCategories::CalculateAncestorToDescendant(int site_index, double &prob_reads, double &all_stats_diff, double &log_likelihood_scaler) {

    prob_reads = 0;
    all_stats_diff = 0;

    const auto &descendant_genotypes_index = all_sequence_prob_index[site_index].GetDescendantIndex();
    uint32_t anc_genotype_index = all_sequence_prob_index[site_index].GetAncestorIndex();
    auto &ancestor_genotype_10 =  convert_index_key_to_diploid_10[anc_genotype_index];
//    auto &cache =  convert_index_key_to_diploid[anc_genotype_index];

    log_likelihood_scaler = convert_index_key_to_diploid_10_scaler[anc_genotype_index]; //TODO: constant for size, ++ memory vs ++ time??
    for (auto &genotypes_index : descendant_genotypes_index) {
        log_likelihood_scaler += convert_index_key_to_haploid_scaler[genotypes_index];
    }
    log_likelihood_scaler=0;
    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;
    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        CacheLoopDesAll2(index10, descendant_genotypes_index, prod_prob_ancestor, summary_stat_diff_ancestor);//Uses this one

//        NoCacheCalculateDes(site_index, index10, prod_prob_ancestor, t, summary_stat_diff_ancestor);
//        CacheLoopDesAll(site_index, index10, prod_prob_ancestor, summary_stat_diff_ancestor);
//        CacheLoopDesAll3(index10, aa, prod_prob_ancestor, summary_stat_diff_ancestor);

        double prob_reads_given_a = ancestor_genotype_10[index10] * prod_prob_ancestor;
        prob_reads += prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;

    }
    all_stats_diff /= prob_reads;
    all_stats_diff /= descendant_count;//TODO: need this to auto calculate stat_same. sum to 1

#ifdef DEBUG7
    if(site_index >0) {
        double likelihood = 0;
        auto &ancestor_genotype_10_unnormalised =  convert_index_key_to_diploid_10_unnormalised[anc_genotype_index];
        for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
            double prob = 0.0;
            double temp = 0.0;
            CalculateLikelihoodUnnormalised(site_index, index10, prob, temp);
//        CacheLoopDesAll2(index10, descendant_genotypes_index, prod_prob_ancestor, summary_stat_diff_ancestor);//Uses this one
            double prob_reads_given_a = ancestor_genotype_10_unnormalised[index10] * prob;
            likelihood += prob_reads_given_a;
//            std::cout << (ancestor_genotype_10_unnormalised[index10]/ancestor_genotype_10[index10]) << "\t";

        }
//        std::cout << "" << std::endl;
        double temp_s = 1;
        double temp_2 = 0;
        std::vector<uint32_t> const &descendant_index = all_sequence_prob_index[site_index].GetDescendantIndex();
        for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
            HaploidProbs pune = convert_index_key_to_haploid_unnormalised[descendant_index[d]];//UNTESTED:
            HaploidProbs pp = convert_index_key_to_haploid[descendant_index[d]];//UNTESTED:
//            for (int i = 0; i < pp.size(); ++i) {
//                std::cout << (pune[i]/pp[i]) << "\t" ;
//            }std::cout << "" << std::endl;
            temp_s *= (pune[0]/pp[0]);
            temp_2 += convert_index_key_to_haploid_scaler[descendant_index[d]];
        }

        temp_2 += convert_index_key_to_diploid_10_scaler[anc_genotype_index];

//        std::cout << temp_s << std::endl;
        temp_s *= (ancestor_genotype_10_unnormalised[0]/ancestor_genotype_10[0]);
//        std::cout << temp_s << std::endl;
        temp_s = 1/temp_s;
//        std::cout << temp_s << std::endl;
//            std::cout << site_index << "\t" << prob_reads << "\t" << likelihood << "\t" << (prob_reads/likelihood) << "\t" << std::endl;
        if( log(prob_reads)+ temp_2 - log(likelihood) > 1e-13 ) {
            std::cout << site_index << "\t" << log(prob_reads) << "\t" << (log(prob_reads)+ temp_2) << "\t" <<
                    log(likelihood) << "\t" << ( (log(prob_reads)+ temp_2) - log(likelihood)) << "\t" << log(temp_s) << "\t" << temp_2 << std::endl;
        }//        std::cout << temp_2 << std::endl;
//        std::exit(4);
        if(log_likelihood_scaler == temp_2){
            std::cout << log_likelihood_scaler << "\t" << temp_2 << std::endl;
        }
    }
#endif


#ifdef DEBUG5
    	double stat_same = 0;
        double all_stats_same = 0;

        prob_reads = 0;
        all_stats_diff = 0;
        for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10)
        {
            NoCacheCalculateDes(site_index, index10, prod_prob_ancestor, stat_same, summary_stat_diff_ancestor);
            double prob_reads_given_a = all_ancestor_genotype_prior[site_index][index10] * prod_prob_ancestor;
            prob_reads += prob_reads_given_a;
            all_stats_same += stat_same*prob_reads_given_a;
            all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;
		}
        all_stats_diff /= prob_reads * descendant_count;
        all_stats_same /= prob_reads * descendant_count;
#endif

}




void MutationModelMultiCategories::CacheLoopDesAll2(int anc_index, const std::vector<uint32_t> &aa, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;

    auto &vv = cache_read_data_to_all_index_rev[anc_index];
//    for (int d = 0; d < descendant_count; ++d) {
    for (auto &item : aa) {

        std::pair<double, double> &cp = vv[item];
        product_prob_given_ancestor *= cp.first;
        summary_stat_diff_ancestor += cp.second;
    }
}


void MutationModelMultiCategories::NoCacheCalculateDes(int site_index, int a, double &product_prob_given_ancestor, double &stat_same, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    stat_same = 0;
    std::vector<uint32_t> const &descendant_index = all_sequence_prob_index[site_index].GetDescendantIndex();
    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;

//        HaploidProbs prob_reads_given_descent = all_sequence_prob[site_index].GetDescendantGenotypes(d);
        HaploidProbs prob_reads_given_descent = convert_index_key_to_haploid[descendant_index [d]];//UNTESTED:
        CalculateOneDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        product_prob_given_ancestor *= sum_over_probs;

    }
}

void MutationModelMultiCategories::CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs &prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;


    for (int b = 0; b < BASE_COUNT; ++b) {

        double prob = transition_matrix_a_to_d_single(index16, b) * prob_reads_given_descent[b];
//        double prob = cache_data_transition[t][anc_index10][b];
        prob_reads_d_given_a += prob;

        summary_stat_same += prob_reads_given_descent[b] * (1.0- mutation_rate_single) * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate_single * frequency_prior[b];


    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

}

void MutationModelMultiCategories::CalculateLikelihoodUnnormalised(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    double stat_same = 0;
    std::vector<uint32_t> const &descendant_index = all_sequence_prob_index[site_index].GetDescendantIndex();
    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;

//        HaploidProbs prob_reads_given_descent = all_sequence_prob[site_index].GetDescendantGenotypes(d);
        HaploidProbs prob_reads_given_descent = convert_index_key_to_haploid_unnormalised[descendant_index [d]];//UNTESTED:
        CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        product_prob_given_ancestor *= sum_over_probs;

    }
}

