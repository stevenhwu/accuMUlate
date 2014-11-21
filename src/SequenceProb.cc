/*
 * SequencingProb.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
//#include <tkDecls.h>
#include "SequenceProb.h"


//
//double TetMAProbability(const ModelParams &params, const ModelInput site_data) {
//	MutationMatrix m = MutationAccumulation(params, false);
//	MutationMatrix mt = MutationAccumulation(params, true);
//
//	MutationMatrix mn = m - mt;
//	DiploidProbs pop_genotypes = DiploidPopulation(params, site_data.reference);
//
//	auto data = site_data.all_reads.begin();
//	DiploidProbs anc_genotypes = DiploidSequencing(params, site_data.reference,
//			*data);
//	anc_genotypes *= pop_genotypes;
//	DiploidProbs num_genotypes = anc_genotypes;
//	for (++data; data != site_data.all_reads.end(); ++data) {
//		HaploidProbs p = HaploidSequencing(params, site_data.reference, *data);
//
//		anc_genotypes *= (m.matrix() * p.matrix()).array();
//		num_genotypes *= (mn.matrix() * p.matrix()).array();
//		exit(-1);
//	}
//
//	//cerr << "\n" << anc_genotypes/anc_genotypes.sum() << endl;
//
//	return 1.0 - num_genotypes.sum() / anc_genotypes.sum();
//}


const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;

Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const bool DEBUG = true;

SequenceProb::SequenceProb(const ModelParams &model_params,
		const ModelInput site_data, MutationProb muProb) {

	params = model_params;

	UpdateMuProb(muProb);

	data = site_data;
	pop_genotypes = DiploidPopulation(site_data.reference);

    descendant_count = site_data.all_reads.size()-1;

    ancestor = site_data.all_reads[0];
	anc_genotypes = DiploidSequencing(ancestor);
	anc_genotypes *= pop_genotypes;
//	DiploidProbs num_genotypes = anc_genotypes;

    for (int i = 0; i < descendant_count; ++i) {
        ReadData data = site_data.all_reads[i+1];
        all_descendant.push_back(data);

        HaploidProbs p = HaploidSequencing(data);
        all_hap.push_back(p);

        HaploidProbs normalised = p;///p.sum();
        all_normalised_hap.push_back(normalised);
    }

	UpdateLikelihood();

}

SequenceProb::~SequenceProb() {
}


void SequenceProb::CountReadToGenotype(){

//    DiploidPopulation(site_data.reference);

    uint16_t ref_count[4];
//    for (auto base : acgt) {
    for (int i = 0; i < 4; ++i) {
        ref_count[0] = ancestor.reads[0];
    }

//    = it.reads[0];
//    ref_count
//            + it.reads[1]+ it.reads[2]+ it.reads[3];
    int match_count[4];
    int mismatch_count[4];
    for (auto descendant : all_descendant) {

        for (int i = 0; i < 4; ++i) {
            int delta = fabs( ref_count[i] = descendant.reads[i]);
            mismatch_count[i] = delta;
            match_count[i] = 0;//??
        }
    }
}



void SequenceProb::UpdateMuProb(MutationProb muProb){

	mutation = muProb.GetMutation();
	non_mutation = muProb.GetNonMutation();
    mutation_rate = muProb.GetMutationRate();

    ancestor_prior = muProb.GetAncestorPrior();
    frequency_prior = muProb.GetFrequencyPrior();
    beta = muProb.GetBeta();


    UpdateLikelihood();
}

DiploidProbs SequenceProb::DiploidPopulation(int ref_allele) {
	ReadData d;
	DiploidProbs result;

	double alphas[4];
	for(int i : {0,1,2,3}){
		alphas[i] = params.theta*params.nuc_freq[i];
//		printf("A:%f\n",alphas[i]);
	}
//	printf("ref:%d\n", ref_allele);
	for(int i : {0,1,2,3}) {
		for(int j=0;j<i;++j) {
			d.key = 0;
			d.reads[ref_allele] = 1;
			d.reads[i] += 1;
			d.reads[j] += 1;
//			printf("  %d %d: %u %u %u %u\n", i, j, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
			result[i+j*4] = DirichletMultinomialLogProbability(alphas, d);
			result[j+i*4] = result[i+j*4];
		}
		d.key = 0;
		d.reads[ref_allele] = 1;
		d.reads[i] += 2;
//		printf("D:%d %d: %u %u %u %u\n", i, i, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
		result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
	}
//	std::cout << result << std::endl;
	return result.exp();//TODO: find out why no -maxCoef here??
}

DiploidProbs SequenceProb::DiploidSequencing(ReadData data) {
	DiploidProbs result;
	double alphas_total = (1.0 - params.phi_diploid) / params.phi_diploid;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j = 0; j < i; ++j) {
			double alphas[4];
			for (int k : { 0, 1, 2, 3 }) {
				if (k == i || k == j)
					alphas[k] = (0.5 - params.error_prob / 3.0) * alphas_total;
				else
					alphas[k] = (params.error_prob / 3.0) * alphas_total;
			}
			result[i * 4 + j] = DirichletMultinomialLogProbability(alphas,
                    data);
			result[j * 4 + i] = result[i * 4 + j];
		}
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - params.error_prob) * alphas_total;
			else
				alphas[k] = params.error_prob / 3.0 * alphas_total;
		}
		result[i * 4 + i] = DirichletMultinomialLogProbability(alphas, data);
	}
	double scale = result.maxCoeff();
	return (result - scale).exp();
//    result -= scale;
//    result = NormaliseLogArray(result);
//    return result;
}

HaploidProbs SequenceProb::HaploidSequencing(ReadData data) {
	HaploidProbs result;
	double alphas_total = (1.0 - params.phi_haploid) / params.phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - params.error_prob) * alphas_total;
			else
				alphas[k] = params.error_prob / 3.0 * alphas_total;
		}
		result[i] = DirichletMultinomialLogProbability(alphas, data);
	}

	double scale = result.maxCoeff();
	return (result - scale).exp();
//    result -= scale;
//    result = NormaliseLogArray(result);
//    return result;
}

template <class T> //FIXME: make sure it only take Eigen::Array class
T SequenceProb::NormaliseLogArray(T result) {
    result = result.exp();
    double normalise = result.sum();
    result /=  normalise;
    return result;
//    double scale = result.maxCoeff();
//	return (result - scale).exp();
}


void SequenceProb::UpdateLikelihood() {

//	auto data = site_data.all_reads.begin();
//	DiploidProbs anc_genotypes = DiploidSequencing(params, site_data.reference, *data);
//	anc_genotypes *= pop_genotypes;
	DiploidProbs anc_genotypes2 = anc_genotypes;
	DiploidProbs num_genotypes = anc_genotypes;
	for (auto t : all_hap) {

		anc_genotypes2 *= (mutation.matrix() * t.matrix()).array();
		num_genotypes *= (non_mutation.matrix() * t.matrix()).array();

	}

	likelihood = 1- num_genotypes.sum() / anc_genotypes2.sum();
	//cerr << "\n" << anc_genotypes/anc_genotypes.sum() << endl;

	likelihood = anc_genotypes2.sum();

}





double SequenceProb::GetLikelihood() {
	return likelihood;
}

//ModelInput GetData();
ModelInput SequenceProb::GetData(){
	return data;
}

HaploidProbs SequenceProb::GetDescendantToReads(int descent_index) {

    return all_normalised_hap[descent_index];
}

DiploidProbs SequenceProb::GetAncGenotypesToReads() {
    return anc_genotypes;
}


//int index_vector2[10][14] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
//        {0,1,1,1,2,2,2,2,2,2},// A
//        {2,1,2,2,0,1,1,2,2,2},// C
//        {2,2,1,2,2,1,2,0,1,2},// G
//        {2,2,2,1,2,2,1,2,1,0} // T
//
//};
double total_sum = 0;
void SequenceProb::CalculateAncestorToDescendant() {
    cout << "==================In Calc AtoD====================\n";
    for (int i = 0; i < descendant_count; ++i) {

        PrintReads(data.all_reads[i]);
    }

//    mu = 0.9;
//    prob_reads_given_descent[0] = 0.2;
//    prob_reads_given_descent[1] = 0.2;
//    prob_reads_given_descent[2] = 0.4;
//    prob_reads_given_descent[3] = 0.2;

//    for (int r = 0; r < 2; ++r) {
cout << anc_genotypes << endl;

    double prob_reads;
    double sum_all_stats_same;
    double sum_all_stats_diff;
    double summary_stat_same_ancestor[ANCESTOR_COUNT];
    double summary_stat_diff_ancestor[ANCESTOR_COUNT];
    double sum_prob_ancestor[ANCESTOR_COUNT];

        for (int a = 0; a < ANCESTOR_COUNT; ++a) {
            int index10 = a;
            int index16 = LookupTable::index_converter_10_to_16[a];
            if (DEBUG) {
                cout << "==Loop A: " << a << "\t" << index16 << "\t" << LookupTable::genotype_lookup_10[a] << "\t" <<
                        anc_genotypes(index16) << endl;
            }
//summary_stat_same_ancestor[a]=0;
            CalculateAllDescendantGivenAncestor(a, summary_stat_same_ancestor, summary_stat_diff_ancestor, sum_prob_ancestor);

anc_genotypes(index16) = 1.0/16.0 ;//-> sum up to 1 FIXME!!
            double prob_reads_given_a = anc_genotypes(index16) * ancestor_prior[index10] *  sum_prob_ancestor[a];// TODO: check ? P(R_Y | X) and P(R_X | X)
//            prob_reads_given_a = log(anc_genotypes(index16)) + log( ancestor_prior[index10]) + sum_prob_ancestor[a];// TODO: check ? P(R_Y | X) and P(R_X | X)

            summary_stat_same_ancestor[a] *= prob_reads_given_a;
            summary_stat_diff_ancestor[a] *= prob_reads_given_a;
            prob_reads += prob_reads_given_a;

            sum_all_stats_same += summary_stat_same_ancestor[a];
            sum_all_stats_diff += summary_stat_diff_ancestor[a];
            //            sum_all_stats += prob_reads_given_a * summary_stat_ancestor[a];

            //        sum_prob_ancestor[a] = sum(sum_prob_d);
            //        summary_stat_ancestor[a] =
            if (DEBUG) {

                cout << "==summary: Same:\t" << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff <<
                "====\t" << sum_prob_ancestor[a] << "\t" << anc_genotypes[index16] << "\t" << ancestor_prior[index10] << "\t" << prob_reads_given_a << endl << endl;

            }


        }
    sum_all_stats_same /= prob_reads;
    sum_all_stats_diff /= prob_reads;

//    cout << ancestor_prior<< end;
cout << "summaryALL\tSame:" << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff << "\t" << (sum_all_stats_diff + sum_all_stats_same) << endl << endl;
cout << "total_sum: "<<total_sum<<endl;
    exit(99);
    summary_stat_same_ancestor;
    summary_stat_diff_ancestor;// sum(stat*prb/denom)
}



void SequenceProb::CalculateAllDescendantGivenAncestor(int a, double summary_stat_same_ancestor[], double summary_stat_diff_ancestor[], double sum_prob_ancestor[]) {

//    std::vector<double> sum_prob_d(descendant_count);
//    std::vector<double> summary_stat_diff_d(descendant_count);
//    std::vector<double> summary_stat_same_d(descendant_count);;


    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double sum_over_probs = 0;
sum_prob_ancestor[a] = 1;

    for (int d = 0; d < descendant_count; ++d) {
        HaploidProbs prob_reads_given_descent = GetDescendantToReads(d); //Fixed value for now

//        CalculateDescendantGivenAncestor(a, prob_reads_given_descent, sum_prob_d[d], summary_stat_same_d[d], summary_stat_diff_d[d]);

        CalculateDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor[a] += summary_stat_diff;//_d[d];
        summary_stat_same_ancestor[a] += summary_stat_same;//_d[d];

        sum_prob_ancestor[a] *= sum_over_probs;

        if (DEBUG) {
            cout << "====D: " << d << "\t Sum:" <<
                    sum_over_probs << "\t" << sum_prob_ancestor[a] << "\t" <<
                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << endl;
        }

    }

}

void SequenceProb::CalculateDescendantGivenAncestor(int a, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[a];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;

//prob_reads_given_descent = {0.1,0.2,0.3,0.4};

    for (int b = 0; b < BASE_COUNT; ++b) {
//prob_reads_given_descent[b] = 0.25;
        double prob = transition_matrix_a_to_d(index16, b) * prob_reads_given_descent[b];
        prob_reads_d_given_a += prob;

        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.mu * frequency_prior[b];
//        total_sum += summary_stat_diff + summary_stat_same;
        total_sum += mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b] + mutation_rate.mu * frequency_prior[b];
        if (DEBUG) {

            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate.mu * frequency_prior[b];

//            cout << "======Loop base: " << b << "\t" << "\t" << prob <<"\t"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << endl;//t1 << "\t" << t2 <<endl;
          }
    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;
}


void SequenceProb::PrintReads(ReadData data) {
    printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}


double SequenceProb::CalculateExpectedValueForMu(Array10D summary_stat_AtoD) {

//    anc_genotypes;
    double sum = 0;
    double nume_sum = 0;
    for (int i = 0; i < 4; ++i) {
        for(int j=i;j< 4;++j) {
            int index10 = LookupTable::index_converter_16_to_10[i][j];
            int index16 = i*4+j;
            double denominator = anc_genotypes(index16)* frequency_prior[index10];
            sum += denominator;

//            nume_sum += denominator * summary_stat_AtoD(index10); // TODO: add P(X) freq_1 X freq_2
//            cout << i<<":"<<j<< " " << index10 << ":\t" << anc_genotypes(index16) << "\t" << summary_stat_AtoD(index10) << endl;
        }
    }

    cout << nume_sum << "\t" << sum << "\t" << (nume_sum/sum) << endl;
    return nume_sum/sum;


}

double SequenceProb::Maximisation(double summery_stat) {

    double max_hat = -beta * log(1- summery_stat / beta);
    return max_hat;

}

void SequenceProb::UpdateTransitionMatrix(EvolutionModel evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
}
