/*
 * site_prob.cc.cc
 *
 *  Created on: Dec 1, 2014
 *      Author: Steven Wu
 */

#include <iostream>

#include "site_prob.h"


const int DEBUG = 0;

SiteProb::SiteProb(SiteGenotypes &sequence_prob,
		 MutationProb &mutation_prob, EvolutionModel &evo_model) {

    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();
    UpdateMuProb(mutation_prob);
    UpdateTransitionMatrix(evo_model);

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
    descendant_count = all_descendant_genotypes.size();
    std::cout << "Deprecated method" << std::endl;
    std::exit(99);
}


SiteProb::SiteProb(SiteGenotypes &sequence_prob, EvolutionModel &evo_model) {
    MutationProb mutation_prob = evo_model.GetMutationProb() ;
    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
    descendant_count = all_descendant_genotypes.size();

//    all_descendant_diff_stats.resize(descendant_count);
//    all_descendant_diff_stats2.resize(descendant_count, 0);
    all_descendant_diff_stats2.assign(descendant_count, 0);
//    master_prob.resize(16, std::vector<std::array<double, 4>>(descendant_count) );
//    master_prob2.resize(16*descendant_count, std::array<double, 4> {{0,0,0,0,}} );

    UpdateModel(evo_model);


//    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
//        int index16 = LookupTable::index_converter_10_to_16[index10];
//        for (int d = 0; d < descendant_count; ++d) {
//            int offset = index16*d;
//            for (int b = 0; b < BASE_COUNT; ++b) {
////                master_prob[index16][d][b] = transition_matrix_a_to_d(index16, b) * all_descendant_genotypes[d][b];
//                master_prob2[offset+d][b] = transition_matrix_a_to_d(index16, b) * all_descendant_genotypes[d][b];
//            }}}


//    for (int i = 0; i < descendant_count; ++i) {
//        for (int b = 0; b < BASE_COUNT; ++b) {
//            all_descendant_diff_stats[i][b] = all_descendant_genotypes[i][b] * frequency_prior_mutation_rate[b];
//        }
//    }

//    double a,b,c;
//    CalculateAncestorToDescendant(a,b,c);
}


SiteProb::SiteProb(SequenceProb &sequence_prob, EvolutionModel &evo_model) {
    MutationProb mutation_prob = evo_model.GetMutationProb();
    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
    descendant_count = all_descendant_genotypes.size();

    all_descendant_diff_stats2.assign(descendant_count, 0);

    UpdateModel(evo_model);
}


    SiteProb::~SiteProb() {
}

//int g_count = 0;
//int g_count2 = 0;
void SiteProb::UpdateModel(EvolutionModel &evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
    mutation_rate = evo_model.GetMutationRate();

    for (int b = 0; b < BASE_COUNT; ++b) {
        frequency_prior_mutation_rate[b] = mutation_rate * frequency_prior[b];
    }
//TODO: time trial
//    all_descendant_diff_stats2.resize(descendant_count,0);
//    all_descendant_diff_stats2.assign(descendant_count,0);
    std::fill(all_descendant_diff_stats2.begin(), all_descendant_diff_stats2.end(), 0);
    for (int i = 0; i < descendant_count; ++i) {

        for (int b = 0; b < BASE_COUNT; ++b) {
//            all_descendant_diff_stats[i][b] = all_descendant_genotypes[i][b] * frequency_prior_mutation_rate[b];
            all_descendant_diff_stats2[i] += all_descendant_genotypes[i][b] * frequency_prior_mutation_rate[b];
        }
    }



//    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
//        int index16 = LookupTable::index_converter_10_to_16[index10];
//        int offset = index16*descendant_count;
//        for (int d = 0; d < descendant_count; ++d) {
//            int offset2 = offset+d;
//            for (int b = 0; b < BASE_COUNT; ++b) {
////                master_prob[index16][d][b] = transition_matrix_a_to_d(index16, b) * all_descendant_genotypes[d][b];
//                master_prob2[offset2][b] = transition_matrix_a_to_d(index16, b) * all_descendant_genotypes[d][b];
//            }
//        }
//    }
}


void SiteProb::UpdateMuProb(MutationProb mutation_prob){

    mutation_rate = mutation_prob.GetMutationRate();
    for (int b = 0; b < BASE_COUNT; ++b) {
        frequency_prior_mutation_rate[b] = mutation_rate * frequency_prior[b];
    }
}

void SiteProb::UpdateTransitionMatrix(EvolutionModel &evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
}


double total_sum2 = 0;
void SiteProb::CalculateAncestorToDescendant(double &prob_reads, double &all_stats_same, double &all_stats_diff) {

    prob_reads = 0;
    all_stats_same = 0;
    all_stats_diff = 0;


    double summary_stat_same_ancestor = 0;
    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;

    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {

        int index16 = LookupTable::index_converter_10_to_16[index10];

        CalculateAllDescendantGivenAncestor(index16, prod_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);

        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor;
            //TODO: IF need speed, cache this
        prob_reads += prob_reads_given_a;
//        std::cout << prob_reads_given_a << "\t" << ancestor_genotypes[index16] << "\t" <<  ancestor_prior[index10] << "\t" <<  prod_prob_ancestor<< "\t" <<  std::endl;
//        all_stats_same += summary_stat_same_ancestor*prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;


//        std::cout << "S: " << prob_reads_given_a << "\t" << (summary_stat_diff_ancestor) << "\t" << prob_reads << "\t" << all_stats_diff  << std::endl;

        if (DEBUG>1) {
            std::cout << "==A: " << index10 << " " << index16 << " " << LookupTable::genotype_lookup_10[index10] << " " <<
                    ancestor_genotypes[index16] << "\t";
            std::cout << "Same: " << all_stats_same << "\tDiff:" << all_stats_diff <<
            "\t" << summary_stat_same_ancestor << "\t" <<summary_stat_diff_ancestor <<
            "\tProb:" << prod_prob_ancestor  << "\t" << prob_reads_given_a << std::endl ;
        }

    }
//    std::exit(3);
//    all_stats_same /= prob_reads;
    all_stats_diff /= prob_reads;
    all_stats_diff /= descendant_count;
    all_stats_same = 1 - all_stats_diff;
//    all_stats_same = descendant_count- all_stats_diff;
    if(DEBUG>0){
        std::cout << "summaryALL\tSame:" << all_stats_same << "\tDiff:" << all_stats_diff << "\t" << (all_stats_diff + all_stats_same) << std::endl;
        std::cout << "total_sum2: "<< total_sum2 << "\tProb: " << prob_reads <<std::endl;

    }

}

void SiteProb::CalculateAllDescendantGivenAncestor(int index16, double &product_prob_given_ancestor,
        double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor) {


    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    summary_stat_same_ancestor = 0;

    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 0;
//        HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now
//        std::array<double, 4> all_diff_stats = all_descendant_diff_stats[d];
        CalculateOneDescendantGivenAncestor(index16, d, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
//        summary_stat_same_ancestor += summary_stat_same;

        product_prob_given_ancestor *= sum_over_probs;
        if (DEBUG>0) {
            HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now
            std::cout << "====D: " << d << "\t Sum:" <<
                    sum_over_probs << "\t" << product_prob_given_ancestor << "\t" <<
                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << std::endl;
        }
//        std::cout << sum_over_probs << "\t" << std::endl;
    }
//    std::cout << std::endl;

}
//void SiteProb::CalculateOneDescendantGivenAncestor(int index16, int des_index, HaploidProbs &prob_reads_given_descent, std::array<double, 4> &all_diff_stats, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {
void SiteProb::CalculateOneDescendantGivenAncestor(int index16, int des_index, double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;

//prob_reads_given_descent = {0.4,0.1,0.1,0.4};
//    prob_reads_given_descent = {0.03,0.03,0.04,0.9};
//    prob_reads_given_descent = {0.25,0.25,0.25,0.25};
//    prob_reads_given_descent = {0.3,0.2,0.2,0.3};

    for (int b = 0; b < BASE_COUNT; ++b) {
        double prob = transition_matrix_a_to_d(index16, b) * all_descendant_genotypes[des_index][b];

//        double prob = master_prob2[index16][b];
        prob_reads_d_given_a += prob;

//        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
//        summary_stat_diff += prob_reads_given_descent[b] * frequency_prior_mutation_rate[b];

//        summary_stat_diff += all_descendant_diff_stats[des_index][b];
        //TODO: Cache these as lookup as well

//        total_sum2 += summary_stat_diff + summary_stat_same;
//        total_sum2 += mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[index16][b] + mutation_rate.prob * frequency_prior[b];
        if (DEBUG > 3) {
//            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[index16][b];
//            double t2 = prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
//            std::cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << std::endl;//t1 << "\t" << t2 <<std::endl;
        }
    }
//    summary_stat_same /= prob_reads_d_given_a;
//    if (summary_stat_diff != all_descendant_diff_stats2[des_index]) {
//       std::cout << summary_stat_diff << "\t" <<  all_descendant_diff_stats2[des_index] << std::endl;
//    }

    summary_stat_diff = all_descendant_diff_stats2[des_index] / prob_reads_d_given_a;
//    summary_stat_diff /= prob_reads_d_given_a;

    if (DEBUG>2) {
        std::cout << index16 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << std::endl;
    }
//    std::cout << index16 << "\t" << prob_reads_d_given_a << "\t" << summary_stat_diff << std::endl;
}



