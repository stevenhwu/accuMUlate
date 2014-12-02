/*
 * SiteProb.cc
 *
 *  Created on: Dec 1, 2014
 *      Author: Steven Wu
 */
#include <iostream>
//#include <tkDecls.h>
#include "SiteProb.h"


//
//double TetMAProbability(const ModelParams &params, const ModelInput site_data) {
//	MutationMatrix m = MutationAccumulation(params, false);
//	MutationMatrix mt = MutationAccumulation(params, true);
//
//	MutationMatrix mn = m - mt;
//	DiploidProbs pop_genotypes = DiploidPopulation(params, site_data.reference);
//
//	auto data = site_data.all_reads.begin();
//	DiploidProbs ancestor_genotypes = DiploidSequencing(params, site_data.reference,
//			*data);
//	ancestor_genotypes *= pop_genotypes;
//	DiploidProbs num_genotypes = ancestor_genotypes;
//	for (++data; data != site_data.all_reads.end(); ++data) {
//		HaploidProbs p = HaploidSequencing(params, site_data.reference, *data);
//
//		ancestor_genotypes *= (m.matrix() * p.matrix()).array();
//		num_genotypes *= (mn.matrix() * p.matrix()).array();
//		exit(-1);
//	}
//
//	//cerr << "\n" << ancestor_genotypes/ancestor_genotypes.sum() << endl;
//
//	return 1.0 - num_genotypes.sum() / ancestor_genotypes.sum();
//}

const int DEBUG = 0;

SiteProb::SiteProb(SequenceProb sequence_prob,
		 MutationProb mutation_prob, EvolutionModel evo_model) {



    ancestor_prior = mutation_prob.GetAncestorPrior();
    frequency_prior = mutation_prob.GetFrequencyPrior();
    UpdateMuProb(mutation_prob);
    UpdateTransitionMatrix(evo_model);

    ancestor_genotypes = sequence_prob.GetAncestorGenotypes();
    all_descendant_genotypes = sequence_prob.GetDescendantGenotypes();
    descendant_count = all_descendant_genotypes.size();

}



SiteProb::~SiteProb() {
}


void SiteProb::UpdateMuProb(MutationProb mutation_prob){

    mutation_rate = mutation_prob.GetMutationRate();

}

void SiteProb::UpdateTransitionMatrix(EvolutionModel evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
}


double total_sum2 = 0;
void SiteProb::CalculateAncestorToDescendant(double &sum_prob, double &stat_same, double &stat_diff) {

    double prob_reads = 0;
    double sum_all_stats_same = 0;
    double sum_all_stats_diff = 0;
    double summary_stat_same_ancestor[ANCESTOR_COUNT];
    double summary_stat_diff_ancestor[ANCESTOR_COUNT];
    double prod_prob_ancestor[ANCESTOR_COUNT];

    for (int a = 0; a < ANCESTOR_COUNT; ++a) {
        int index10 = a;
        int index16 = LookupTable::index_converter_10_to_16[a];

//summary_stat_same_ancestor[a]=0;
        CalculateAllDescendantGivenAncestor(a, prod_prob_ancestor[a], summary_stat_same_ancestor[a], summary_stat_diff_ancestor[a]);

        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  prod_prob_ancestor[a];
//            prob_reads_given_a = log(ancestor_genotypes(index16)) + log( ancestor_prior[index10]) + prod_prob_ancestor[a];
        prob_reads += prob_reads_given_a;

//        summary_stat_same_ancestor[a] *= prob_reads_given_a;
//        summary_stat_diff_ancestor[a] *= prob_reads_given_a;

        sum_all_stats_same += summary_stat_same_ancestor[a]*prob_reads_given_a;
        sum_all_stats_diff += summary_stat_diff_ancestor[a]*prob_reads_given_a;
        //            sum_all_stats += prob_reads_given_a * summary_stat_ancestor[a];

        //        prod_prob_ancestor[a] = sum(sum_prob_d);
        //        summary_stat_ancestor[a] =
        if (DEBUG>1) {
            cout << "==A: " << a << " " << index16 << " " << LookupTable::genotype_lookup_10[a] << " " <<
                    ancestor_genotypes[index16] << "\t";
            cout << "Same: " << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff <<
            "\t" << summary_stat_same_ancestor[a] << "\t" <<summary_stat_diff_ancestor[a] <<
            "\tProb:" << prod_prob_ancestor[a]  << "\t" << prob_reads_given_a << endl ;

        }


    }

    sum_all_stats_same /= prob_reads;
    sum_all_stats_diff /= prob_reads;

//    cout << ancestor_prior<< end;
    if(DEBUG>0){
        cout << "summaryALL\tSame:" << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff << "\t" << (sum_all_stats_diff + sum_all_stats_same) << endl << endl;
        cout << "total_sum2: "<< total_sum2 << "\tProb: " << prob_reads <<endl;

    }
    stat_same = sum_all_stats_same;
    stat_diff = sum_all_stats_diff;
    sum_prob = prob_reads;
}


void SiteProb::CalculateAllDescendantGivenAncestor(int a, double &product_prob_given_ancestor,
        double &summary_stat_same_ancestor, double &summary_stat_diff_ancestor) {


    product_prob_given_ancestor = 1;
    summary_stat_diff_ancestor = 0;
    summary_stat_same_ancestor = 0;
//    descendant_count = 1; // FIXME
    for (int d = 0; d < descendant_count; ++d) {
        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double sum_over_probs = 0;
        HaploidProbs prob_reads_given_descent = all_descendant_genotypes[d]; //Fixed value for now

//        CalculateOneDescendantGivenAncestor(a, prob_reads_given_descent, sum_prob_d[d], summary_stat_same_d[d], summary_stat_diff_d[d]);

        CalculateOneDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor += summary_stat_diff;//_d[d];
        summary_stat_same_ancestor += summary_stat_same;//_d[d];

//        sum_prob_ancestor[a] += log(sum_over_probs);
        product_prob_given_ancestor *= sum_over_probs;
        if (DEBUG>2) {
            cout << "====D: " << d << "\t Sum:" <<
                    sum_over_probs << "\t" << product_prob_given_ancestor << "\t" <<
                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << endl;
        }

    }
//    sum_prob_ancestor[a] = exp(sum_prob_ancestor[a]);

}

void SiteProb::CalculateOneDescendantGivenAncestor(int anc_index10, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[anc_index10];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;

//prob_reads_given_descent = {0.4,0.1,0.1,0.4};
//    prob_reads_given_descent = {0.03,0.03,0.04,0.9};
//    prob_reads_given_descent = {0.25,0.25,0.25,0.25};
//    prob_reads_given_descent = {0.3,0.2,0.2,0.3};

    for (int b = 0; b < BASE_COUNT; ++b) {

        double prob = transition_matrix_a_to_d(index16, b) * prob_reads_given_descent[b];
        prob_reads_d_given_a += prob;

        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
//
//        summary_stat_same += prob * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[a][b];
//        summary_stat_diff += prob * mutation_rate.mu * frequency_prior[b];

//        summary_stat_same += prob_reads_given_descent[b]; //* mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[a][b];
//        summary_stat_diff += prob_reads_given_descent[b]; //* mutation_rate.mu * frequency_prior[b];

//        total_sum2 += summary_stat_diff + summary_stat_same;
//        total_sum2 += mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b] + mutation_rate.prob * frequency_prior[b];
        if (DEBUG>3) {
            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_p * LookupTable::summary_stat_same_lookup_table[anc_index10][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate.prob * frequency_prior[b];
            cout << "======Loop base: " << b << "\t" << "\tP:" << prob <<"\tReadGivenD:"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << endl;//t1 << "\t" << t2 <<endl;
        }
    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;

    if (DEBUG>2) {
        cout << anc_index10 << "\t" <<summary_stat_same << "\t" << summary_stat_diff << endl;
    }
}


void SiteProb::PrintReads(ReadData data) {
    printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}




