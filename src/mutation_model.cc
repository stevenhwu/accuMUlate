/*
 * mutation_model.cc
 *
 *  Created on: 12/15/14
 *      Author: Steven Wu
 */



#include <stdint.h>
#include <strings.h>

#include "mutation_model.h"


int DEBUG = 0;

MutationModel::MutationModel(EvolutionModel &evo_model0) {

    MutationProb mutation_prob = evo_model0.GetMutationProb();
    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();
    evo_model = &evo_model0;

    transition_matrix_a_to_d = evo_model->GetTranstionMatirxAToD();
    mutation_rate = evo_model->GetMutationRate();

}

int MutationModel::GetSiteCount() const {
    return site_count;
}


void MutationModel::AddSequenceProb(std::vector<SiteGenotypes> &all) {
    all_sequence_prob = all;
    site_count = all.size();
    descendant_count = all_sequence_prob[0].GetDescendantCount();
    std::cout << "Assuming all data have the same number of descendants. If not, reimplement this!!." << std::endl;
    all_ancestor_genotype_prior.resize(all_sequence_prob.size());

    for (size_t i = 0; i < all_sequence_prob.size(); ++i) {
        auto ancestor_genotypes = all_sequence_prob[i].GetAncestorGenotypes();

        std::array<double, 10> ancestor_genotype_prior;
        for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
            int index16 = LookupTable::index_converter_10_to_16[index10];
            ancestor_genotype_prior[index10] = ancestor_genotypes[index16] * ancestor_prior[index10];
        }
        all_ancestor_genotype_prior[i] = ancestor_genotype_prior;
    }

    InitCache();
}


void MutationModel::UpdateExpBeta(double expBeta) {

    evo_model->UpdateExpBeta(expBeta);

    transition_matrix_a_to_d = evo_model->GetTranstionMatirxAToD();
    mutation_rate = evo_model->GetMutationRate();

//    clock_t t_start = clock();
//    for (int i = 0; i < 1e5; ++i) {

        UpdateCache();
//    }
//    clock_t t_end = clock();
//    std::cout << "MicroB time cache: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << std::endl;

}



void MutationModel::InitCache() {

    std::unordered_map<uint64_t, int> count;

//    cache_read_data_to_all.reserve(10000000);
//    cache_read_data_to_all.rehash( 10000000);
//    cache_read_data_to_all_index.reserve(all_sequence_prob.size());
//    cache_read_data_to_all_index = std::vector<std::array<std::pair<double, double>, 10> >(all_sequence_prob.size());

    int index = 0;
    for (size_t i = 0; i < all_sequence_prob.size(); ++i) {
        SiteGenotypes &item = all_sequence_prob[i];
        for (int j = 0; j < item.GetDescendantCount(); ++j) {
            ReadData rd = item.GetDescendantReadData(j);
            auto rd_key = rd.key;
            count[rd_key]++;
            auto find_key = map_rd_to_index.find(rd_key);

            if(find_key == map_rd_to_index.end()){
//                std::cout << item.GetDescendantIndex(j) << "\t" << index << std::endl;
                map_rd_key_to_haploid[rd_key]= item.GetDescendantGenotypes(j);

//                std::array<std::array<double, 2>, 10> temp;
                std::array<std::pair<double, double>, 10> temp;


                //                cache_read_data_to_all_index.push_back( temp);
                cache_read_data_to_all_index.push_back(temp) ;
                map_index_key_to_haploid.push_back( item.GetDescendantGenotypes(j) );
                map_rd_to_index[rd_key] = index;
                index++;
//                std::cout << map_rd_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << std::endl;
//                std::cout << "\t" << rd.reads[0] << "\t" << rd.reads[1] << "\t" <<rd.reads[2] << "\t" <<rd.reads[3] << "\t" << rd_key << std::endl;

            }
            item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
        }

    }
    std::cout << map_rd_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << "\t" << index << std::endl;

    for (int k = 0; k < 10; ++k) {
        std::vector<std::pair<double, double>> temp;
        temp.assign(index, std::make_pair(0,0));
//        std::vector<std::array<double, 2>> temp;
//        for (int i = 0; i < index; ++i) {
//            std::array<double, 2> t2;
//            temp.push_back(t2);
//        }
        cache_read_data_to_all_index_rev[k] = temp;
    }

    for (size_t i = 0; i < all_sequence_prob.size(); ++i) {
        SiteGenotypes &item = all_sequence_prob[i];
//        for (int j = 0; j < item.GetDescendantCount(); ++j) {
////            std::cout << item.GetDescendantIndex(j) << "\t" ;
//            if(item.GetDescendantIndex(j) == -1){
//                std::cout << i << "\t" << j << "\t" <<  item.GetDescendantIndex(j) << "\t" << std::endl;
//            }
//        }
//        item.SortIndex();
//        std::cout << std::endl;

    }
    summary_stat_diff_vec = std::vector<double>(map_index_key_to_haploid.size());
        UpdateCache();
//    std::exit(6);
}
void MutationModel::UpdateCache() {


    std::array<double, 4> temp_base_prob;
    for (int b = 0; b < BASE_COUNT; ++b) {
        temp_base_prob[b] = mutation_rate * frequency_prior[b];
    }

//    cache_read_data_to_all_index[index];// = temp;
//    map_index_key_to_haploid[index];// = item.GetDescendantGenotypes(j);


//    for (int i = 0; i < map_index_key_to_haploid.size(); ++i) {
//
////        auto rd_key = item.first;
////        auto genotype = item.second;
//        auto &genotype = map_index_key_to_haploid[i];
//        double summary_stat_diff = 0;
//        for (int b = 0; b < BASE_COUNT; ++b) {
//            summary_stat_diff += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
//        }
//
//        auto &cache_all = cache_read_data_to_all_index[i];
//
//        for (int k = 0; k < ANCESTOR_COUNT; ++k) {
//            int index16 = LookupTable::index_converter_10_to_16[k];
//            double prob_reads_d_given_a = 0;
//
//            for (int b = 0; b < BASE_COUNT; ++b) {
//                double prob = transition_matrix_a_to_d(index16, b) * genotype[b];
//                prob_reads_d_given_a += prob;
////                //stat_same += genotype[b] * evo_model.GetMutationProb().mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[k][b];
//            }
//
////            cache_all[k][0] = prob_reads_d_given_a;
////            cache_all[k][1] = summary_stat_diff / prob_reads_d_given_a;
////            cache_read_data_to_all_index_rev[k][i][0]= prob_reads_d_given_a;
////            cache_read_data_to_all_index_rev[k][i][1]= summary_stat_diff / prob_reads_d_given_a;
//
////            cache_all[k].first = prob_reads_d_given_a;
////            cache_all[k].second = summary_stat_diff / prob_reads_d_given_a;
//            cache_read_data_to_all_index_rev[k][i].first = prob_reads_d_given_a;
//            cache_read_data_to_all_index_rev[k][i].second = summary_stat_diff / prob_reads_d_given_a;
//        }
//
//
//    }
//    for (int i = 0; i < map_index_key_to_haploid.size(); ++i) {
//        auto &genotype = map_index_key_to_haploid[i];
//        summary_stat_diff_vec[i] = 0;
//        for (int b = 0; b < BASE_COUNT; ++b) {
//            summary_stat_diff_vec[i] += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
//        }
//    }
//
//
//    for (int k = 0; k < ANCESTOR_COUNT; ++k) {
//        auto &temp = cache_read_data_to_all_index_rev[k];
//        for (int i = 0; i < map_index_key_to_haploid.size(); ++i) {
//
//            auto &genotype = map_index_key_to_haploid[i];
////            double summary_stat_diff = 0;
////            for (int b = 0; b < BASE_COUNT; ++b) {
////                summary_stat_diff += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
////            }
//
//
//            int index16 = LookupTable::index_converter_10_to_16[k];
//            double prob_reads_d_given_a = 0;
//
//            for (int b = 0; b < BASE_COUNT; ++b) {
//                double prob = transition_matrix_a_to_d(index16, b) * genotype[b];
//                prob_reads_d_given_a += prob;
//    //                //stat_same += genotype[b] * evo_model.GetMutationProb().mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[k][b];
//            }
//
//            temp[i].first = prob_reads_d_given_a;
//            temp[i].second = summary_stat_diff_vec[i] / prob_reads_d_given_a;
//            }
//    }

//    std::cout << cache_read_data_to_all.load_factor() << "\t" << cache_read_data_to_all.bucket_count() << "\t" << cache_read_data_to_all.max_bucket_count()<< std::endl;
//      for (unsigned i=0; i<cache_read_data_to_all.bucket_count(); ++i) {
//    std::cout << "bucket #" << i << " has " << cache_read_data_to_all.bucket_size(i) << " elements.\n";
//  }
//     for (auto& x: cache_read_data_to_all) {
//    std::cout << "Element [" << x.first << ":"<< "]";
//    std::cout << " is in bucket #" << cache_read_data_to_all.bucket (x.first) << std::endl;
//  }
//    std::exit(3);

//    std::vector<double> summary_stat_diff_vec(map_index_key_to_haploid.size());
//    std::fill(summary_stat_diff_vec.begin(), summary_stat_diff_vec.end(), 0);
//    for (int i = 0; i < map_index_key_to_haploid.size(); ++i) {
//        auto &genotype = map_index_key_to_haploid[i];
////        summary_stat_diff_vec[i] = 0;
//        for (int b = 0; b < BASE_COUNT; ++b) {
//            summary_stat_diff_vec[i] += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
//        }
//    }


    for (size_t i = 0; i < map_index_key_to_haploid.size(); ++i) {
        auto &genotype = map_index_key_to_haploid[i];

        double summary_stat_diff= 0;
        for (int b = 0; b < BASE_COUNT; ++b) {
            summary_stat_diff += genotype[b] * temp_base_prob[b];//p * evo_model.GetMutationProb().mutation_rate.prob * frequency_prior[b];
        }


        for (int k = 0; k < ANCESTOR_COUNT; ++k) {
            int index16 = LookupTable::index_converter_10_to_16[k];
            double prob_reads_d_given_a = 0;

            for (int b = 0; b < BASE_COUNT; ++b) {
                double prob = transition_matrix_a_to_d(index16, b) * genotype[b];
                prob_reads_d_given_a += prob;
//                //stat_same += genotype[b] * evo_model.GetMutationProb().mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[k][b];
            }

//            cache_all[k][0] = prob_reads_d_given_a;
//            cache_all[k][1] = summary_stat_diff / prob_reads_d_given_a;
//            cache_read_data_to_all_index_rev[k][i][0]= prob_reads_d_given_a;
//            cache_read_data_to_all_index_rev[k][i][1]= summary_stat_diff / prob_reads_d_given_a;

//            cache_all[k].first = prob_reads_d_given_a;
//            cache_all[k].second = summary_stat_diff / prob_reads_d_given_a;
            cache_read_data_to_all_index_rev[k][i].first = prob_reads_d_given_a;
            cache_read_data_to_all_index_rev[k][i].second = summary_stat_diff / prob_reads_d_given_a;
//            cache_read_data_to_all_index_rev[k][i].second = summary_stat_diff_vec[i] / prob_reads_d_given_a;
        }


    }


}


void MutationModel::CalculateAncestorToDescendant(int site_index, double &prob_reads, double &all_stats_diff) {

    prob_reads = 0;
    all_stats_diff = 0;

    double summary_stat_diff_ancestor = 0;
    double prod_prob_ancestor = 1;

    const std::vector<int> &aa = all_sequence_prob[site_index].GetDescendantIndex();
//    std::vector< std::array<std::array<double, 2>, 10>* > cd;

//    for (int i = 0; i < descendant_count; ++i) {
//        cd.push_back( &cache_read_data_to_all_index[aa[i]]);
//    }


    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
//        int index16 = LookupTable::index_converter_10_to_16[index10];
//        CalculateAllDescendantGivenAncestor(index10, prod_prob_ancestor, summary_stat_diff_ancestor);

//        clock_t t_start = clock();
//        for (int i = 0; i < 1e9; ++i) {

//            CacheLoopDesAll(site_index, index10, prod_prob_ancestor, summary_stat_diff_ancestor);
        CacheLoopDesAll2(index10, aa, prod_prob_ancestor, summary_stat_diff_ancestor);
//        CacheLoopDesAll3(index10, aa, prod_prob_ancestor, summary_stat_diff_ancestor);
//            01:13    2:7-8    3:10-11
//        }
//        clock_t t_end = clock();
//        std::cout << "MicroB time: " << (t_end - t_start) / CLOCKS_PER_SEC << "\t" << (t_end - t_start) << std::endl;

        double prob_reads_given_a = all_ancestor_genotype_prior[site_index][index10] * prod_prob_ancestor;
        prob_reads += prob_reads_given_a;
        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;

//        prob_reads += cache_read_data_to_all_index_rev[index10][aa[index10]].second;
//        all_stats_diff += cache_read_data_to_all_index_rev[index10][aa[index10]].first;
//        prob_reads += cache_read_data_to_all_index_rev[index10][index10].second;
//        all_stats_diff += cache_read_data_to_all_index_rev[index10][index10].first;

#ifdef DEBUG1
        std::cout << "==A: " << index10 << " " << all_ancestor_genotype_prior[site_index][index10] << "\t" << "Diff:" <<
        		all_stats_diff << "\t" << summary_stat_diff_ancestor << "\tProb:" << prod_prob_ancestor  << "\t" << prob_reads_given_a << std::endl ;
#endif
    }
    all_stats_diff /= prob_reads;
    all_stats_diff /= descendant_count;

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

//
//
//void MutationModel::CacheLoopDes(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
//    product_prob_given_ancestor = 1;
//    summary_stat_diff_ancestor = 0;
//    for (int d = 0; d < descendant_count; ++d) {
//        auto key = all_descendant_data[d].key;
//
//        auto cache = cache_read_data_to_all[key][a];
////        sum_over_probs = cache.first;
////        summary_stat_diff = cache.second;
//        product_prob_given_ancestor *= cache[0];
//        summary_stat_diff_ancestor += cache[1];
//
//    }
//}

void MutationModel::CacheLoopDesAll4(int anc_index, std::vector<std::array<std::array<double, 2>, 10>*> &cd, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;

//    auto &vv = cache_read_data_to_all_index_rev[anc_index];
    for (int d = 0; d < descendant_count; ++d) {

        std::array<double, 2> &cache = (*(cd[d]) )[anc_index];
        product_prob_given_ancestor *= cache[0];
        summary_stat_diff_ancestor += cache[1];

    }
}


void MutationModel::CacheLoopDesAll3(int anc_index, const std::vector<int> &aa, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;


//    for (int d = 0; d < descendant_count; ++d) {
    for (auto &item : aa) {

//        std::array<double, 2> &cache = cache_read_data_to_all_index[aa[d]][anc_index];
//        product_prob_given_ancestor *= cache[0];
//        summary_stat_diff_ancestor += cache[1];
        std::pair<double, double> cp = cache_read_data_to_all_index[item][anc_index];
        product_prob_given_ancestor *= cp.first;
        summary_stat_diff_ancestor += cp.second;


    }
}



void MutationModel::CacheLoopDesAll2(int anc_index, const std::vector<int> &aa, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {

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


void MutationModel::CacheLoopDesAll(int site_index, int anc_index, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {

    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;

//    auto &t = all_sequence_prob[site_index];
//    auto key = all_sequence_prob[site_index].GetDescendantReadDataKey(0);

//    auto &aa = cache_read_data_to_all2[anc_index];
    const std::vector<int> &aa = all_sequence_prob[site_index].GetDescendantIndex();
    auto &vv = cache_read_data_to_all_index_rev[anc_index];
    for (int d = 0; d < descendant_count; ++d) {

//        auto key = all_sequence_prob[site_index].GetDescendantReadDataKey(d);
//        auto &cache = cache_read_data_to_all[key][anc_index];
//
////        auto key = aa.GetDescendantIndex(d);
////        std::array<double, 2> cache;
//
        std::pair<double, double> &cp = cache_read_data_to_all_index[aa[d]][anc_index];
        product_prob_given_ancestor *= cp.first;
        summary_stat_diff_ancestor += cp.second;

    }
//    std::cout << std::endl;
}


//
//void MutationModel::CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
//
////    product_prob_given_ancestor=1;
////    double product_prob_given_ancestor1 = 1;
////    double summary_stat_diff_ancestor1 = 0;
////    double product_prob_given_ancestor2 = 1;
////    double summary_stat_diff_ancestor2 = 0;
//
////    summary_stat_same_ancestor = 0;
//
//
////    NoCacheCalculateDes(a, product_prob_given_ancestor1, summary_stat_diff_ancestor1);//-10
//
////    CacheLoopMap(a, product_prob_given_ancestor2, summary_stat_diff_ancestor2);//7~
//    CacheLoopDes(a, product_prob_given_ancestor, summary_stat_diff_ancestor);//7~
//
//    //    product_prob_given_ancestor = product_prob_given_ancestor1;
////    summary_stat_diff_ancestor = summary_stat_diff_ancestor1;
//
////    product_prob_given_ancestor = product_prob_given_ancestor2;
////    summary_stat_diff_ancestor = summary_stat_diff_ancestor2;
//
////    if( (product_prob_given_ancestor1 - product_prob_given_ancestor2) > (product_prob_given_ancestor1/1e10) ){
////        std::cout << "Diff prob" << "\t" << product_prob_given_ancestor1 << "\t" << product_prob_given_ancestor2 << "\t" << (product_prob_given_ancestor1/1e8) << std::endl;
////    }
////    if( (summary_stat_diff_ancestor1 - summary_stat_diff_ancestor2) > (summary_stat_diff_ancestor1/1e10) ){
////        std::cout << "Diff stat" << "\t" << summary_stat_diff_ancestor1 << "\t" << summary_stat_diff_ancestor2 << "\t" << (summary_stat_diff_ancestor1/1e8) << std::endl;
////    }
//
//
//}
//
void MutationModel::NoCacheCalculateDes(int site_index, int a, double &product_prob_given_ancestor, double &stat_same, double &summary_stat_diff_ancestor) {
    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    stat_same = 0;
    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 1;
        HaploidProbs prob_reads_given_descent = all_sequence_prob[site_index].GetDescendantGenotypes(d);
        CalculateOneDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;
        product_prob_given_ancestor *= sum_over_probs;

        stat_same += summary_stat_same;
    }
}

void MutationModel::CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs &prob_reads_given_descent,
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

        summary_stat_same += prob_reads_given_descent[b] * (1.0-mutation_rate) * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
        summary_stat_diff += p * mutation_rate * frequency_prior[b];

        if (DEBUG>3)
        {
            double t1 = prob_reads_given_descent[b] * (1.0-mutation_rate) * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate * frequency_prior[b];
            std::cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << std::endl;//t1 << "\t" << t2 <<std::endl;
        }
    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

    if (DEBUG>2) {
        std::cout << anc_index10 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << std::endl;
    }
}


//void MutationModel::CacheLoopMap(int a, double &product_prob_given_ancestor, double &summary_stat_diff_ancestor) {
//    product_prob_given_ancestor = 1;
//    summary_stat_diff_ancestor = 0;
////    for (auto item : map) {
////    if(vec.size() == descendant_count){
////
////        for (int d = 0; d < descendant_count; ++d) {
////            auto key = all_descendant_data[d].key;
////
////            auto cache = cache_read_data_to_all[key][a];
////            double cache_prob = cache[0];
////            double cache_stat_diff = cache[1];
////            summary_stat_diff_ancestor += cache_stat_diff;
////            product_prob_given_ancestor *= cache_prob;
////        }
////
////    }
////    else {
//
//    for (auto item : vec) {
//
//        auto cache = cache_read_data_to_all[item.first][a];
//        double cache_prob = cache[0];
//        double cache_stat_diff = cache[1];
//
//        if (item.second == 1) {
//
//            summary_stat_diff_ancestor += cache_stat_diff;
//            product_prob_given_ancestor *= cache_prob;
//        }
//        else {
//            summary_stat_diff_ancestor += (cache_stat_diff * item.second);
//            for (int i = 0; i < item.second; ++i) {
//                product_prob_given_ancestor *= cache_prob;;
//            }
//        }
//
//    }
////    }
//}

//
//void MutationModel::AddCache(std::unordered_map<double, std::array<double, 4>> &allDMaps4,
//        std::unordered_map<double,  Std2DArray>  &c2) {
//    std::cout << "ADD\n";
//    cache_data = allDMaps4;
//    cache_data_transition = c2;
//}
//
//
//void MutationModel::CalculateLikelihood(SequenceProbV1 &sequence_prob) {
//
//    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
//    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
//    all_descendant_data = sequence_prob.GetDescendantReadData();
////    all_descendant_data = *sequence_prob.GetDescendantReadData3();
////    all_descendant_data = sequence_prob.GetDescendantReadData2();
////    sequence_prob.GetDescendantReadDataCOPY(all_descendant_data);
//    descendant_count = sequence_prob.GetDescendantCount();
//
//    map = sequence_prob.temp_map;
//    vec = sequence_prob.condense_genotype;
//
////    CalculateAncestorToDescendant(a,b,c);
//}


//
//void MutationModel::CalculateLikelihoodOriginal(SequenceProbV1 &sequence_prob, double &prob, double &stat_same, double &stat_diff) {
//
//    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
//    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
////    all_descendant_data = sequence_prob.GetDescendantReadData();
//    descendant_count = all_descendant_genotypes.size();
//
////    map = sequence_prob.temp_map;
////    vec = sequence_prob.condense_genotype;
//
//    CalculateAncestorToDescendant(prob, stat_same, stat_diff);
//
//}
//
//void MutationModel::CalculateAncestorToDescendantOrig(double &prob_reads, double &all_stats_same, double &all_stats_diff) {
//
//    prob_reads = 0;
//    all_stats_same = 0;
//    all_stats_diff = 0;
//
//
//    double summary_stat_same_ancestor = 0;
//    double summary_stat_diff_ancestor = 0;
//    double prod_prob_ancestor = 1;
//
//    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
//        int index16 = LookupTable::index_converter_10_to_16[index10];
//
//        CalculateAllDescendantGivenAncestorOrig(index10, prod_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);
//
//        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor;
//        prob_reads += prob_reads_given_a;
//
////        all_stats_same += summary_stat_same_ancestor*prob_reads_given_a;
//        all_stats_diff += summary_stat_diff_ancestor*prob_reads_given_a;
//
//        if (DEBUG>1) {
//            std::cout << "==A: " << index10 << " " << index16 << " " << LookupTable::genotype_lookup_10[index10] << " " <<
//                    ancestor_genotypes[index16] << "\t";
//            std::cout << "Same: " << all_stats_same << "\tDiff:" << all_stats_diff <<
//                    "\t" << summary_stat_same_ancestor << "\t" <<summary_stat_diff_ancestor <<
//                    "\tProb:" << prod_prob_ancestor  << "\t" << prob_reads_given_a << std::endl ;
//        }
//
//    }
//
////    all_stats_same /= prob_reads;
//    all_stats_diff /= prob_reads;
//
//    if(DEBUG>0){
//        std::cout << "summaryALL\tSame:" << all_stats_same << "\tDiff:" << all_stats_diff << "\t" << (all_stats_diff + all_stats_same) << std::endl;
//        std::cout << "total_sum_5: "<< total_sum_5 << "\tProb: " << prob_reads <<std::endl;
//
//    }
//
//}
//
//
//void MutationModel::CalculateAllDescendantGivenAncestorOrig(int a, double &product_prob_given_ancestor,
//        double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor) {
//
//
//    product_prob_given_ancestor = 1;
//    summary_stat_diff_ancestor = 0;
//    summary_stat_same_ancestor = 0;
//
//    for (int d = 0; d < descendant_count; ++d) {//TODO: Check descendant info, merge some of them together
//        double summary_stat_same = 0;
//        double summary_stat_diff = 0;
//        double sum_over_probs = 0;
//        HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now
//
//        CalculateOneDescendantGivenAncestorOrig(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);
//
////        summary_stat_diff_ancestor += summary_stat_diff;//_d[d];
//        summary_stat_same_ancestor += summary_stat_same;//_d[d];
//        product_prob_given_ancestor *= sum_over_probs;
//
//        if (DEBUG>0) {
//            std::cout << "====D: " << d << "\t Sum:" <<
//                    sum_over_probs << "\t" << product_prob_given_ancestor << "\t" <<
//                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
//                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << std::endl;
//        }
//
//    }
//
//
//}
//
//
//
//void MutationModel::CalculateOneDescendantGivenAncestorOrig(int anc_index10, HaploidProbs prob_reads_given_descent,
//        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {
//
//    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
//    prob_reads_d_given_a = 0;
//    summary_stat_same = 0;
//    summary_stat_diff = 0;
//
//
//    for (int b = 0; b < BASE_COUNT; ++b) {
//        double p = prob_reads_given_descent[b];
////        double prob = transition_matrix_a_to_d(index16, b) * prob_reads_given_descent[b];
//        double prob = transition_matrix_a_to_d(index16, b) * p;
//        prob_reads_d_given_a += prob;
//
//        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
////        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
//        summary_stat_diff += p * mutation_rate.prob * frequency_prior[b];
////
//        if (DEBUG>3) {
//            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
//            double t2 = prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
//            std::cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << std::endl;//t1 << "\t" << t2 <<std::endl;
//        }
//    }
//    summary_stat_same /= prob_reads_d_given_a;
//    summary_stat_diff /= prob_reads_d_given_a;
//
//    if (DEBUG>2) {
//        std::cout << anc_index10 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << std::endl;
//    }
//}
//
//
//
//void MutationModel::CalculateOneDescendantGivenAncestorCache(int anc_index10, HaploidProbs prob_reads_given_descent,
//        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {
//
//    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
//    prob_reads_d_given_a = 0;
//    summary_stat_same = 0;
//    summary_stat_diff = 0;
//
//    for (int b = 0; b < BASE_COUNT; ++b) {
//        double p = prob_reads_given_descent[b];
//        double prob = transition_matrix_a_to_d(index16, b) * p;
////        double prob = cache_data_transition[t][anc_index10][b];
//        prob_reads_d_given_a += prob;
//
////        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
////        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
//        summary_stat_diff += p * mutation_rate.prob * frequency_prior[b];
////        summary_stat_diff += cache_data[p][b] * mutation_rate.prob;
//
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
//
//    }
////    summary_stat_same /= prob_reads_d_given_a;
//    summary_stat_diff /= prob_reads_d_given_a;
//
//}

//void MutationModel::AddCache2(std::unordered_map<uint64_t, std::array<std::array<double, 2>, 10> > &all_2) {
//    cache_read_data_to_all = all_2;
//
//}

