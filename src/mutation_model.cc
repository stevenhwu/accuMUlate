/*
 * mutation_model.cc
 *
 *  Created on: 12/15/14
 *      Author: Steven Wu
 */



#include "mutation_model.h"




double total_sum_5 = 0;
int DEBUG = 0;


MutationModel::MutationModel(
        MutationProb mutation_prob, EvolutionModel evo_model) {

    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();
    UpdateMuProb(mutation_prob);
    UpdateTransitionMatrix(evo_model);

//    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
//    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
//    descendant_count = all_descendant_genotypes.size();
descendant_count = 7;
}


void MutationModel::UpdateModel(EvolutionModel &evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
    mutation_rate = evo_model.GetMutationRate();
//    cout << "new mutation rate: " << mutation_rate.prob  << "\t" << mutation_rate.one_minus_p <<
//            "\t" << thisCount << endl;

}

void MutationModel::UpdateMuProb(MutationProb mutation_prob){

    mutation_rate = mutation_prob.GetMutationRate();

}

void MutationModel::UpdateTransitionMatrix(EvolutionModel evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
}



void MutationModel::CalculateAncestorToDescendant(double &prob_reads, double &all_stats_same, double &all_stats_diff) {

    prob_reads = 0;
    all_stats_same = 0;
    all_stats_diff = 0;


    double summary_stat_same_ancestor = 0;
    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;

    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        int index16 = LookupTable::index_converter_10_to_16[index10];
//        std::cout << a << "\t" ;

//        for (int b = 0; b < BASE_COUNT; ++b) {
//            std::cout  << transition_matrix_a_to_d(index16, b) << "\t";
//        }
//        std:: cout << std::endl;

        CalculateAllDescendantGivenAncestor(index10, prod_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);

        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor;
        //TODO: IF need speed, cache this
        prob_reads += prob_reads_given_a;


//        all_stats_same += summary_stat_same_ancestor*prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;

//        if (DEBUG>1) {
//            std::cout << "==A: " << index10 << " " << index16 << " " << LookupTable::genotype_lookup_10[index10] << " " <<
//                    ancestor_genotypes[index16] << "\t";
//            std::cout << "Same: " << all_stats_same << "\tDiff:" << all_stats_diff <<
//                    "\t" << summary_stat_same_ancestor << "\t" <<summary_stat_diff_ancestor <<
//                    "\tProb:" << prod_prob_ancestor  << "\t" << prob_reads_given_a << std::endl ;
//        }

    }


//    all_stats_same /= prob_reads;
    all_stats_diff /= prob_reads;

    if(DEBUG>0){
        std::cout << "summaryALL\tSame:" << all_stats_same << "\tDiff:" << all_stats_diff << "\t" << (all_stats_diff + all_stats_same) << std::endl;
        std::cout << "total_sum_5: "<< total_sum_5 << "\tProb: " << prob_reads <<std::endl;

    }

}


void MutationModel::CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor,
        double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor=1;
    double product_prob_given_ancestor1 = 1;
    double summary_stat_diff_ancestor1 = 0;
    double product_prob_given_ancestor2 = 1;
    double summary_stat_diff_ancestor2 = 0;

//    summary_stat_same_ancestor = 0;


//    NoCacheCalculateDes(a, product_prob_given_ancestor1, summary_stat_diff_ancestor1);//-10
//    product_prob_given_ancestor = product_prob_given_ancestor1;
//    summary_stat_diff_ancestor = summary_stat_diff_ancestor1;

//    CacheLoopMap(a, product_prob_given_ancestor2, summary_stat_diff_ancestor2);//7~
    CacheLoopDes(a, product_prob_given_ancestor2, summary_stat_diff_ancestor2);//7~
    product_prob_given_ancestor = product_prob_given_ancestor2;
    summary_stat_diff_ancestor = summary_stat_diff_ancestor2;

//    if( (product_prob_given_ancestor1 - product_prob_given_ancestor2) > (product_prob_given_ancestor1/1e10) ){
//        std::cout << "Diff prob" << "\t" << product_prob_given_ancestor1 << "\t" << product_prob_given_ancestor2 << "\t" << (product_prob_given_ancestor1/1e8) << std::endl;
//    }
//    if( (summary_stat_diff_ancestor1 - summary_stat_diff_ancestor2) > (summary_stat_diff_ancestor1/1e10) ){
//        std::cout << "Diff stat" << "\t" << summary_stat_diff_ancestor1 << "\t" << summary_stat_diff_ancestor2 << "\t" << (summary_stat_diff_ancestor1/1e8) << std::endl;
//    }


}

void MutationModel::NoCacheCalculateDes(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;
        HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now
        CalculateOneDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        product_prob_given_ancestor *= sum_over_probs;
    }
}

void MutationModel::CacheLoopMap(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
//    for (auto item : map) {
//    if(vec.size() == descendant_count){
//
//        for (int d = 0; d < descendant_count; ++d) {
//            auto key = all_descendant_data[d].key;
//
//            auto cache = cache_read_data_to_all_2[key][a];
//            double cache_prob = cache[0];
//            double cache_stat_diff = cache[1];
//            summary_stat_diff_ancestor += cache_stat_diff;
//            product_prob_given_ancestor *= cache_prob;
//        }
//
//    }
//    else {

        for (auto item : vec) {

            auto cache = cache_read_data_to_all_2[item.first][a];
            double cache_prob = cache[0];
            double cache_stat_diff = cache[1];

            if (item.second == 1) {

                summary_stat_diff_ancestor += cache_stat_diff;
                product_prob_given_ancestor *= cache_prob;
            }
            else {
                summary_stat_diff_ancestor += (cache_stat_diff * item.second);
                for (int i = 0; i < item.second; ++i) {
                    product_prob_given_ancestor *= cache_prob;;
                }
            }

        }
//    }
}

void MutationModel::CacheLoopDes(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    for (int d = 0; d < descendant_count; ++d) {
        auto key = all_descendant_data[d].key;

        auto cache = cache_read_data_to_all_2[key][a];
//        sum_over_probs = cache.first;
//        summary_stat_diff = cache.second;
        product_prob_given_ancestor *= cache[0];
        summary_stat_diff_ancestor += cache[1];

    }
}


void MutationModel::CacheLoopDesAll(int index, int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    auto read = all_sequence_prob[index].GetDescendantReadData();
    for (int d = 0; d < descendant_count; ++d) {
        auto key = read[d].key;

        auto cache = cache_read_data_to_all_2[key][a];
//        sum_over_probs = cache.first;
//        summary_stat_diff = cache.second;
        product_prob_given_ancestor *= cache[0];
        summary_stat_diff_ancestor += cache[1];

    }
}

void MutationModel::CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;


    for (int b = 0; b < BASE_COUNT; ++b) {
        double p = prob_reads_given_descent[b];
        double prob = transition_matrix_a_to_d(index16, b) * p;
//        double prob = cache_data_transition[t][anc_index10][b];
        prob_reads_d_given_a += prob;

//        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
//        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
        summary_stat_diff += p * mutation_rate.prob * frequency_prior[b];
//        summary_stat_diff += cache_data[p][b] * mutation_rate.prob;

//        auto o1 =  p * mutation_rate.prob * frequency_prior[b];
//        auto c1 =  cache_data[p][b] * mutation_rate.prob;
//
//        auto c2 = cache_data_transition[p][anc_index10][b];
//        if( (o1 - c1) > (o1/100000000.0)){
//            std::cout << "NOT equal stat: " << b << "\t" << o1 << "\t" << c1 << std::endl;
//        }
//        if(prob != c2){
//            std::cout << "NOT equal prob: " << b << "\t" << prob << "\t" << c2 << std::endl;
//        }
//        TODO: Cache these as lookup as well
//
        if (DEBUG>3) {
            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
            std::cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << std::endl;//t1 << "\t" << t2 <<std::endl;
        }
    }
//    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

    if (DEBUG>2) {
        std::cout << anc_index10 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << std::endl;
    }
}

void MutationModel::AddCache(std::unordered_map<double, std::array<double, 4>> &allDMaps4,
        std::unordered_map<double,  Std2DArray>  &c2) {
    std::cout << "ADD\n";
    cache_data = allDMaps4;
    cache_data_transition = c2;
}


void MutationModel::CalculateLikelihoodAll(std::vector<SequenceProb> &all) {
    all_sequence_prob = all;
//    descendant_count

}

void MutationModel::CalculateLikelihood(SequenceProb &sequence_prob) {

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
    all_descendant_data = sequence_prob.GetDescendantReadData();
//    all_descendant_data = *sequence_prob.GetDescendantReadData3();
//    all_descendant_data = sequence_prob.GetDescendantReadData2();
//    sequence_prob.GetDescendantReadDataCOPY(all_descendant_data);
    descendant_count = sequence_prob.GetDescendantCount();

    map = sequence_prob.temp_map;
    vec = sequence_prob.condense_genotype;

//    CalculateAncestorToDescendant(a,b,c);
}

void MutationModel::AddCache2(
//        std::unordered_map<uint64_t, std::array<double, 10>> prob,
//        std::unordered_map<uint64_t, std::array<double, 10>> stat,
//        std::unordered_map<uint64_t, std::pair<std::array<double, 10>, std::array<double, 10> > > all,
        std::unordered_map<uint64_t, std::array<std::array<double, 2>, 10> > &all_2){
//    cache_read_data_to_prob = prob;
//    cache_read_data_to_stat = stat;
//    cache_read_data_to_all = all;
    cache_read_data_to_all_2 = all_2;
}


void MutationModel::CalculateLikelihoodOriginal(SequenceProb &sequence_prob, double &prob, double &stat_same, double &stat_diff) {

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
//    all_descendant_data = sequence_prob.GetDescendantReadData();
    descendant_count = all_descendant_genotypes.size();

//    map = sequence_prob.temp_map;
//    vec = sequence_prob.condense_genotype;

    CalculateAncestorToDescendant(prob, stat_same, stat_diff);

}

void MutationModel::CalculateAncestorToDescendantOrig(double &prob_reads, double &all_stats_same, double &all_stats_diff) {

    prob_reads = 0;
    all_stats_same = 0;
    all_stats_diff = 0;


    double summary_stat_same_ancestor = 0;
    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;

    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        int index16 = LookupTable::index_converter_10_to_16[index10];

        CalculateAllDescendantGivenAncestorOrig(index10, prod_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);

        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor;
        prob_reads += prob_reads_given_a;

//        all_stats_same += summary_stat_same_ancestor*prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;

        if (DEBUG>1) {
            std::cout << "==A: " << index10 << " " << index16 << " " << LookupTable::genotype_lookup_10[index10] << " " <<
                    ancestor_genotypes[index16] << "\t";
            std::cout << "Same: " << all_stats_same << "\tDiff:" << all_stats_diff <<
                    "\t" << summary_stat_same_ancestor << "\t" <<summary_stat_diff_ancestor <<
                    "\tProb:" << prod_prob_ancestor  << "\t" << prob_reads_given_a << std::endl ;
        }

    }

//    all_stats_same /= prob_reads;
    all_stats_diff /= prob_reads;

    if(DEBUG>0){
        std::cout << "summaryALL\tSame:" << all_stats_same << "\tDiff:" << all_stats_diff << "\t" << (all_stats_diff + all_stats_same) << std::endl;
        std::cout << "total_sum_5: "<< total_sum_5 << "\tProb: " << prob_reads <<std::endl;

    }

}


void MutationModel::CalculateAllDescendantGivenAncestorOrig(int a, double &product_prob_given_ancestor,
        double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor) {


    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    summary_stat_same_ancestor = 0;

    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 0;
        HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now

        CalculateOneDescendantGivenAncestorOrig(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

//        summary_stat_diff_ancestor += summary_stat_diff;//_d[d];
        summary_stat_same_ancestor += summary_stat_same;//_d[d];
        product_prob_given_ancestor *= sum_over_probs;

        if (DEBUG>0) {
            std::cout << "====D: " << d << "\t Sum:" <<
                    sum_over_probs << "\t" << product_prob_given_ancestor << "\t" <<
                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << std::endl;
        }

    }


}



void MutationModel::CalculateOneDescendantGivenAncestorOrig(int anc_index10, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;


    for (int b = 0; b < BASE_COUNT; ++b) {
        double p = prob_reads_given_descent[b];
//        double prob = transition_matrix_a_to_d(index16, b) * prob_reads_given_descent[b];
        double prob = transition_matrix_a_to_d(index16, b) * p;
        prob_reads_d_given_a += prob;

//        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
//        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
        summary_stat_diff += p * mutation_rate.prob * frequency_prior[b];
//
        if (DEBUG>3) {
            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
            std::cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << std::endl;//t1 << "\t" << t2 <<std::endl;
        }
    }
//    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

    if (DEBUG>2) {
        std::cout << anc_index10 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << std::endl;
    }
}



void MutationModel::CalculateOneDescendantGivenAncestorCache(int anc_index10, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;

    for (int b = 0; b < BASE_COUNT; ++b) {
        double p = prob_reads_given_descent[b];
        double prob = transition_matrix_a_to_d(index16, b) * p;
//        double prob = cache_data_transition[t][anc_index10][b];
        prob_reads_d_given_a += prob;

//        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
//        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
        summary_stat_diff += p * mutation_rate.prob * frequency_prior[b];
//        summary_stat_diff += cache_data[p][b] * mutation_rate.prob;

        auto o1 =  p * mutation_rate.prob * frequency_prior[b];
        auto c1 =  cache_data[p][b] * mutation_rate.prob;

        auto c2 = cache_data_transition[p][anc_index10][b];
        if( (o1 - c1) > (o1/100000000.0)){
            std::cout << "NOT equal stat: " << b << "\t" << o1 << "\t" << c1 << std::endl;
        }
        if(prob != c2){
            std::cout << "NOT equal prob: " << b << "\t" << prob << "\t" << c2 << std::endl;
        }

    }
//    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

}