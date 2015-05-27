/*
 * mutation_model_multi_categories.cc
 *
 *  Created on: 4/10/15
 *      Author: Steven Wu
 */



#include "mutation_model_multi_categories.h"


MutationModelMultiCategories::MutationModelMultiCategories(int num_categories, EvolutionModel &evo_model0, SequencingFactory &factory)
    :categories_count(num_categories){

    evo_model = &evo_model0;

    MutationProb mutation_prob = evo_model->GetMutationProb();
    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();

    transition_matrix_a_to_d_single.reserve(categories_count);
    mutation_rate_single.reserve(categories_count);
    for (int i = 0; i < categories_count; ++i) {
        transition_matrix_a_to_d_single.push_back(evo_model->GetTranstionMatirxAToD());
        mutation_rate_single.push_back(evo_model->GetMutationRate());
    }

    convert_index_key_to_haploid = factory.RemoveConvertIndexKeyToHaploid();
    convert_index_key_to_diploid_10 = factory.RemoveConvertIndexKeyToDiploidIndex10();

    convert_index_key_to_haploid_scaler = factory.RemoveConvertIndexKeyToHaploidScaler();
    convert_index_key_to_diploid_10_scaler = factory.RemoveConvertIndexKeyToDiploidIndex10Scaler();

    all_sequence_prob_index = factory.RemoveSiteGenotypeIndexVector();

    site_count = all_sequence_prob_index.size();
//    descendant_count = all_sequence_prob_index[0].GetDescendantCount();

    for (int i = 0; i < categories_count; ++i) {
        InitCache(i);
        UpdateCache(i);
    }

    std::cout << "Site_count: " << "\t" << site_count << "\tDescendant_count varies! "<< std::endl;
//    std::cout << "Assuming all data have the same number of descendants. If not, reimplement this!!." << std::endl;

}

int MutationModelMultiCategories::GetCategoriesCount() const {
    return categories_count;
}

int MutationModelMultiCategories::GetSiteCount() const {
    return site_count;
}


void MutationModelMultiCategories::UpdateOneMinusExpBeta(int category_index, double oneMinusExpBeta) {

    evo_model->UpdateOneMinusExpBeta(oneMinusExpBeta);

    transition_matrix_a_to_d_single[category_index] = evo_model->GetTranstionMatirxAToD();
    mutation_rate_single[category_index] = evo_model->GetMutationRate();

    UpdateCache(category_index);

}

void MutationModelMultiCategories::InitCache(int category_index) {
    size_t index_size = convert_index_key_to_haploid.size();

    cache_read_data_to_all_index_rev.reserve(categories_count);
    for (int i = 0; i < categories_count; ++i) {
        std::array<std::vector<std::pair<double, double>>, 10> temp_array;
        for (int k = 0; k < 10; ++k) {
            std::vector<std::pair<double, double>> temp(index_size, std::make_pair(0, 0) );
//            temp.assign(index_size, std::make_pair(0, 0));
            temp_array[k] = temp;
        }
        cache_read_data_to_all_index_rev.push_back(temp_array);
    }

    UpdateCache(category_index);

}


void MutationModelMultiCategories::UpdateCache(int category_index) {


    std::array<double, 4> temp_base_prob;
    for (int b = 0; b < BASE_COUNT; ++b) {
        temp_base_prob[b] = mutation_rate_single[category_index] * frequency_prior[b];
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
                double prob = transition_matrix_a_to_d_single[category_index](index16, b) * genotype[b];
                prob_reads_d_given_a += prob;
//                //stat_same += genotype[b] * evo_model.GetMutationProb().mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[k][b];
            }

            cache_read_data_to_all_index_rev[category_index][k][i].first = prob_reads_d_given_a;
            cache_read_data_to_all_index_rev[category_index][k][i].second = summary_stat_diff / prob_reads_d_given_a;
        }
    }
}


void MutationModelMultiCategories::CalculateAncestorToDescendant(int category_index, int site_index, double &prob_reads,
                                                                 double &all_stats_diff, double &log_likelihood_scaler) {

    prob_reads = 0;
    all_stats_diff = 0;

    int site_descendant_count = all_sequence_prob_index[site_index].GetDescendantCount();

    const auto &descendant_genotypes_index = all_sequence_prob_index[site_index].GetDescendantIndex();
    uint32_t anc_genotype_index = all_sequence_prob_index[site_index].GetAncestorIndex();
    auto &ancestor_genotype_10 =  convert_index_key_to_diploid_10[anc_genotype_index];
//    auto &cache =  convert_index_key_to_diploid[anc_genotype_index];

    log_likelihood_scaler = convert_index_key_to_diploid_10_scaler[anc_genotype_index]; //TODO: constant for size, ++ memory vs ++ time??
    for (auto &genotypes_index : descendant_genotypes_index) {
        log_likelihood_scaler += convert_index_key_to_haploid_scaler[genotypes_index];
    }
//    log_likelihood_scaler=0;
    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;

//    std::cout << site_descendant_count << ":";
//#define NO_CACHE
#ifdef NO_CACHE
    double no_cache_prob = 1;
    double no_cache_stat_diff = 0;
    double no_cache_stat_same = 0;

    double no_cache_prob_reads = 0;
    double no_cache_all_stats_diff = 0;
    double no_cache_all_stats_same = 0;
#endif


    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        CacheLoopDesAll2(category_index, index10, descendant_genotypes_index, prod_prob_ancestor,
                         summary_stat_diff_ancestor);//Uses this one

        double prob_reads_given_a = ancestor_genotype_10[index10] * prod_prob_ancestor;
        prob_reads += prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor * prob_reads_given_a;
    }
#ifdef NO_CACHE
    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        NoCacheCalculateDes(category_index, site_index, index10, no_cache_prob, no_cache_stat_same, no_cache_stat_diff);
        double no_cache_prob_reads_given_a = ancestor_genotype_10[index10] * no_cache_prob;
        no_cache_prob_reads += no_cache_prob_reads_given_a;
        no_cache_all_stats_diff += no_cache_stat_diff * no_cache_prob_reads_given_a;
        no_cache_all_stats_same += no_cache_stat_same * no_cache_prob_reads_given_a;
    }
//    std::cout << (no_cache_all_stats_diff +no_cache_all_stats_same )/no_cache_prob_reads << "\t"
//    std::cout << site_descendant_count << "\t" << no_cache_all_stats_diff << "\t" << no_cache_all_stats_same  << std::endl;
//    std::cout << prob_reads << "==" << no_cache_prob_reads  << "\t" << all_stats_diff <<"=="<< no_cache_all_stats_diff<< std::endl;

    if ((summary_stat_diff_ancestor - no_cache_stat_diff) > 1e-20) {
        std::cout << "STAT:" << site_index << ":" << category_index << ":" << summary_stat_diff_ancestor << "!=" <<
        no_cache_stat_diff << "\t" << (summary_stat_diff_ancestor - no_cache_stat_diff) << std::endl;
    }
    if ((prob_reads - no_cache_prob_reads) > 1e-18) {
        std::cout << "P:" << site_index << ":" << category_index << ":" << prob_reads << "==" << no_cache_prob_reads <<
        std::endl;
    }
    if ((all_stats_diff - no_cache_all_stats_diff) > 1e-18) {
        std::cout << "S:" << site_index << ":" << category_index << ":" << all_stats_diff << "==" <<
        no_cache_all_stats_diff << std::endl;
    }
#endif


    all_stats_diff /= prob_reads;
//    all_stats_diff /= descendant_count;//NOTE: need this to auto calculate stat_same. sum to 1
//    std::cout << "\t" << all_stats_diff << std::endl;
//    all_stats_diff /= site_descendant_count;



}




void MutationModelMultiCategories::CacheLoopDesAll2(int category_index, int anc_index, const std::vector<uint32_t> &aa,
                                                    double &product_prob_given_ancestor,
                                                    double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;

    auto &vv = cache_read_data_to_all_index_rev[category_index][anc_index];
//    for (int d = 0; d < descendant_count; ++d) {
//    std::cout << aa.size() ;
    for (auto &item : aa) {

        std::pair<double, double> &cp = vv[item];
        product_prob_given_ancestor *= cp.first;
        summary_stat_diff_ancestor += cp.second;
    }
}


void MutationModelMultiCategories::NoCacheCalculateDes(int categories_index, int site_index, int anc_index10,
                                                       double &product_prob_given_ancestor, double &summary_stat_same_ancestor,
                                                       double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    summary_stat_same_ancestor = 0;
    std::vector<uint32_t> const &descendant_index = all_sequence_prob_index[site_index].GetDescendantIndex();
//    std::cout << " "<< descendant_index.size() ;
    for (int d = 0; d < descendant_index.size(); ++d) {
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;

//        HaploidProbs prob_reads_given_descent = all_sequence_prob[site_index].GetDescendantGenotypes(d);
        HaploidProbs prob_reads_given_descent = convert_index_key_to_haploid[descendant_index [d]];//UNTESTED:
        CalculateOneDescendantGivenAncestor(categories_index, anc_index10, prob_reads_given_descent, sum_over_probs, summary_stat_same,
                                            summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        summary_stat_same_ancestor += summary_stat_same;
        product_prob_given_ancestor *= sum_over_probs;

    }
}

void MutationModelMultiCategories::CalculateOneDescendantGivenAncestor(int category_index, int anc_index10,
                                                                       HaploidProbs &prob_reads_given_descent,
                                                                       double &prob_reads_d_given_a,
                                                                       double &summary_stat_same,
                                                                       double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;


    for (int b = 0; b < BASE_COUNT; ++b) {
//        transition_matrix_a_to_d_single[category_index]()
        double prob = transition_matrix_a_to_d_single[category_index](index16, b) * prob_reads_given_descent[b];
//        double prob = cache_data_transition[t][anc_index10][b];
        prob_reads_d_given_a += prob;

        summary_stat_same += prob_reads_given_descent[b] * (1.0- mutation_rate_single[category_index]) * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate_single[category_index] * frequency_prior[b];

    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

}

void MutationModelMultiCategories::CalculateLikelihoodUnnormalised(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    double stat_same = 0;
    std::vector<uint32_t> const &descendant_index = all_sequence_prob_index[site_index].GetDescendantIndex();
    for (int d = 0; d < descendant_index.size(); ++d) {
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;

//        HaploidProbs prob_reads_given_descent = all_sequence_prob[site_index].GetDescendantGenotypes(d);
        HaploidProbs prob_reads_given_descent = convert_index_key_to_haploid_unnormalised[descendant_index [d]];//UNTESTED:
//        CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        product_prob_given_ancestor *= sum_over_probs;

    }
}

int MutationModelMultiCategories::GetDescendantCount(int site_index) {
    return all_sequence_prob_index[site_index].GetDescendantCount();
}
