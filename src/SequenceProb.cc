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


const int ANCESTOR_COUNT = 10;
const int BASE_COUNT = 4;

Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const bool DEBUG = true;

SequenceProb::SequenceProb(const ModelParams &model_params,
		const ModelInput site_data, MutationProb muProb) {

//	params = model_params;

    phi_haploid = model_params.phi_haploid;
    phi_diploid = model_params.phi_diploid;
    error_prob = model_params.error_prob;
    theta = model_params.theta;

    UpdateMuProb(muProb);

//	data = site_data;

    int i=0;
    DiploidProbs pop_genotypes = DiploidPopulation(site_data.reference);
    ancestor_data = site_data.all_reads[i];
	ancestor_genotypes = DiploidSequencing(ancestor_data);
	ancestor_genotypes *= pop_genotypes;

//	DiploidProbs num_genotypes = ancestor_genotypes;

    for (i=1; i < site_data.all_reads.size(); ++i) {
        ReadData data = site_data.all_reads[i];
        all_descendant_data.push_back(data);

        HaploidProbs p = HaploidSequencing(data);
        all_descendant_genotypes.push_back(p);

    }
    descendant_count = site_data.all_reads.size()-1;


}

SequenceProb::~SequenceProb() {
}


void SequenceProb::CountReadToGenotype(){

//    DiploidPopulation(site_data.reference);

    uint16_t ref_count[4];
//    for (auto base : acgt) {
    for (int i = 0; i < 4; ++i) {
        ref_count[0] = ancestor_data.reads[0];
    }

//    = it.reads[0];
//    ref_count
//            + it.reads[1]+ it.reads[2]+ it.reads[3];
    int match_count[4];
    int mismatch_count[4];
    for (auto descendant : all_descendant_data) {

        for (int i = 0; i < 4; ++i) {
            int delta = fabs( ref_count[i] = descendant.reads[i]);
            mismatch_count[i] = delta;
            match_count[i] = 0;//??
        }
    }
}



void SequenceProb::UpdateMuProb(MutationProb muProb){

//	mutation = muProb.GetMutation();
//	non_mutation = muProb.GetNonMutation();
    mutation_rate = muProb.GetMutationRate();

    ancestor_prior = muProb.GetAncestorPrior();
    frequency_prior = muProb.GetFrequencyPrior();
//    exp_beta = muProb.GetBeta();


    UpdateLikelihood();
}

DiploidProbs SequenceProb::DiploidPopulation(int ref_allele) {
	ReadData d;
	DiploidProbs result;

	double alphas[4];
	for(int i : {0,1,2,3}){

        alphas[i] = theta * frequency_prior[i];
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
    double alphas_total = (1.0 - phi_diploid) / phi_diploid;
	for (int i : { 0, 1, 2, 3 }) {
        for (int j = 0; j < i; ++j) {
			double alphas[4];
			for (int k : { 0, 1, 2, 3 }) {
				if (k == i || k == j)
					alphas[k] = (0.5 - error_prob / 3.0) * alphas_total;
				else
					alphas[k] = (error_prob / 3.0) * alphas_total;
			}
			result[i * 4 + j] = DirichletMultinomialLogProbability(alphas, data);
			result[j * 4 + i] = result[i * 4 + j];
		}
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - error_prob) * alphas_total;
			else
				alphas[k] = error_prob / 3.0 * alphas_total;
		}
		result[i * 4 + i] = DirichletMultinomialLogProbability(alphas, data);
	}


	double scale = result.maxCoeff();
//	return (result - scale).exp();
//    result -= scale;
//    result = NormaliseLogArray(result);
    return result.exp();
}

HaploidProbs SequenceProb::HaploidSequencing(ReadData data) {
	HaploidProbs result;
    
    double alphas_total = (1.0 - phi_haploid) / phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
		double alphas[4];
		for (int k : { 0, 1, 2, 3 }) {
			if (k == i)
				alphas[k] = (1.0 - error_prob) * alphas_total;
			else
				alphas[k] = error_prob / 3.0 * alphas_total;
		}
		result[i] = DirichletMultinomialLogProbability(alphas, data);
	}

	double scale = result.maxCoeff();
//	return (result - scale).exp();
//    result -= scale;
//    result = NormaliseLogArray(result);
    return result.exp();
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
//	DiploidProbs ancestor_genotypes = DiploidSequencing(params, site_data.reference, *data);
//	ancestor_genotypes *= pop_genotypes;
	DiploidProbs anc_genotypes2 = ancestor_genotypes;
	DiploidProbs num_genotypes = ancestor_genotypes;
	for (auto t : all_descendant_genotypes) {

//		anc_genotypes2 *= (mutation.matrix() * t.matrix()).array();
//		num_genotypes *= (non_mutation.matrix() * t.matrix()).array();

	}

	likelihood = 1- num_genotypes.sum() / anc_genotypes2.sum();
	//cerr << "\n" << ancestor_genotypes/ancestor_genotypes.sum() << endl;

	likelihood = anc_genotypes2.sum();

}





double SequenceProb::GetLikelihood() {
	return likelihood;
}

//ModelInput GetData();
//ModelInput SequenceProb::GetData(){
//	return data;
//}

HaploidProbs SequenceProb::GetDescendantGenotypes(int descent_index) {
    if(descent_index >= descendant_count || descent_index < 0){
        //TODO: throw error
    }
    return all_descendant_genotypes[descent_index];
}

DiploidProbs SequenceProb::GetAncestorGenotypes() {
    return ancestor_genotypes;
}


//int index_vector2[10][14] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
//        {0,1,1,1,2,2,2,2,2,2},// A
//        {2,1,2,2,0,1,1,2,2,2},// C
//        {2,2,1,2,2,1,2,0,1,2},// G
//        {2,2,2,1,2,2,1,2,1,0} // T
//
//};
double total_sum = 0;
void SequenceProb::CalculateAncestorToDescendant(double &stat_same, double &stat_diff) {
if(DEBUG){
    cout << "==================In Calc AtoD====================\n";
    for (int i = 0; i < descendant_count; ++i) {

        PrintReads(all_descendant_data[i]);
    }
    for (auto f : ancestor_prior) {
        cout << f << "\t";
    }
    cout << endl;
//    mu = 0.9;
//    prob_reads_given_descent[0] = 0.2;
//    prob_reads_given_descent[1] = 0.2;
//    prob_reads_given_descent[2] = 0.4;
//    prob_reads_given_descent[3] = 0.2;

//    for (int r = 0; r < 2; ++r) {

    cout << ancestor_genotypes << endl;
}
    double prob_reads = 0;
    double sum_all_stats_same = 0;
    double sum_all_stats_diff = 0;
    double summary_stat_same_ancestor[ANCESTOR_COUNT];
    double summary_stat_diff_ancestor[ANCESTOR_COUNT];
    double sum_prob_ancestor[ANCESTOR_COUNT];

        for (int a = 0; a < ANCESTOR_COUNT; ++a) {
            int index10 = a;
            int index16 = LookupTable::index_converter_10_to_16[a];
            if (DEBUG) {
                cout << "==Loop A: " << a << "\t" << index16 << "\t" << LookupTable::genotype_lookup_10[a] << "\t" <<
                        ancestor_genotypes[index16] << endl;
            }
//summary_stat_same_ancestor[a]=0;
            CalculateAllDescendantGivenAncestor(a, summary_stat_same_ancestor, summary_stat_diff_ancestor, sum_prob_ancestor);

            double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[index10] *  sum_prob_ancestor[a];
//            prob_reads_given_a = log(ancestor_genotypes(index16)) + log( ancestor_prior[index10]) + sum_prob_ancestor[a];
            prob_reads += prob_reads_given_a;

            summary_stat_same_ancestor[a] *= prob_reads_given_a;
            summary_stat_diff_ancestor[a] *= prob_reads_given_a;

            sum_all_stats_same += summary_stat_same_ancestor[a];
            sum_all_stats_diff += summary_stat_diff_ancestor[a];
            //            sum_all_stats += prob_reads_given_a * summary_stat_ancestor[a];

            //        sum_prob_ancestor[a] = sum(sum_prob_d);
            //        summary_stat_ancestor[a] =
            if (DEBUG) {

                cout << "==summary: Same: " << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff <<
                " ==\t" << sum_prob_ancestor[a] << "\t" << ancestor_genotypes[index16] << "\t" << ancestor_prior[index10] << "\t" << prob_reads_given_a << endl << endl;

            }


        }
    sum_all_stats_same /= prob_reads;
    sum_all_stats_diff /= prob_reads;

//    cout << ancestor_prior<< end;
cout << "summaryALL\tSame:" << sum_all_stats_same << "\tDiff:" << sum_all_stats_diff << "\t" << (sum_all_stats_diff + sum_all_stats_same) << endl << endl;
cout << "total_sum: "<<total_sum<<endl;
    stat_same = sum_all_stats_same;
    stat_diff = sum_all_stats_diff;
//    exit(99);
//    summary_stat_same_ancestor;
//    summary_stat_diff_ancestor;// sum(stat*prb/denom)
}



void SequenceProb::CalculateAllDescendantGivenAncestor(int a, double summary_stat_same_ancestor[], double summary_stat_diff_ancestor[], double sum_prob_ancestor[]) {

//    std::vector<double> sum_prob_d(descendant_count);
//    std::vector<double> summary_stat_diff_d(descendant_count);
//    std::vector<double> summary_stat_same_d(descendant_count);;


    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double sum_over_probs = 0;
sum_prob_ancestor[a] = 1;
    summary_stat_diff_ancestor[a] = 0;
    summary_stat_same_ancestor[a] = 0;
    descendant_count = 1; // FIXME
    for (int d = 0; d < descendant_count; ++d) {
        HaploidProbs prob_reads_given_descent = GetDescendantGenotypes(d); //Fixed value for now

//        CalculateDescendantGivenAncestor(a, prob_reads_given_descent, sum_prob_d[d], summary_stat_same_d[d], summary_stat_diff_d[d]);

        CalculateDescendantGivenAncestor(a, prob_reads_given_descent, sum_over_probs, summary_stat_same, summary_stat_diff);

        summary_stat_diff_ancestor[a] += summary_stat_diff;//_d[d];
        summary_stat_same_ancestor[a] += summary_stat_same;//_d[d];

//        sum_prob_ancestor[a] += log(sum_over_probs);
        sum_prob_ancestor[a] *= sum_over_probs;
        if (DEBUG) {
            cout << "====D: " << d << "\t Sum:" <<
                    sum_over_probs << "\t" << sum_prob_ancestor[a] << "\t" <<
                    "\tSame:" << summary_stat_same << "\tDiff:" << summary_stat_diff << "\t" <<
                    " BASE FREQ: " << prob_reads_given_descent.format(nice_row) << endl;
        }

    }
//    sum_prob_ancestor[a] = exp(sum_prob_ancestor[a]);

}

void SequenceProb::CalculateDescendantGivenAncestor(int a, HaploidProbs prob_reads_given_descent,
        double &prob_reads_d_given_a, double &summary_stat_same, double &summary_stat_diff) {

    int index16 = LookupTable::index_converter_10_to_16[a];
    prob_reads_d_given_a = 0;
    summary_stat_same = 0;
    summary_stat_diff = 0;

//prob_reads_given_descent = {0.4,0.1,0.1,0.4};
//    prob_reads_given_descent = {0.03,0.03,0.04,0.9};
//    prob_reads_given_descent = {0.3,0.2,0.2,0.3};
    for (int b = 0; b < BASE_COUNT; ++b) {
//prob_reads_given_descent[b] = 0.25;
        double prob = transition_matrix_a_to_d(index16, b) * prob_reads_given_descent[b];
        prob_reads_d_given_a += prob;

        summary_stat_same += prob_reads_given_descent[b] * mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
        summary_stat_diff += prob_reads_given_descent[b] * mutation_rate.mu * frequency_prior[b];
//
//        summary_stat_same += prob * mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
//        summary_stat_diff += prob * mutation_rate.mu * frequency_prior[b];

//        summary_stat_same += prob_reads_given_descent[b]; //* mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
//        summary_stat_diff += prob_reads_given_descent[b]; //* mutation_rate.mu * frequency_prior[b];

//        total_sum += summary_stat_diff + summary_stat_same;
        total_sum += mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b] + mutation_rate.mu * frequency_prior[b];
        if (DEBUG) {

            double t1 = prob_reads_given_descent[b] * mutation_rate.one_minus_mu * LookupTable::summary_stat_same_lookup_table[a][b];
            double t2 = prob_reads_given_descent[b] * mutation_rate.mu * frequency_prior[b];

            cout << "======Loop base: " << b << "\t" << "\t" << prob <<"\t"<< prob_reads_given_descent[b] << "\t T1:" << t1 << "\t T2:" << t2 <<"\t SAME:"<<summary_stat_same << "\t" << summary_stat_diff << endl;//t1 << "\t" << t2 <<endl;
          }
    }
    summary_stat_same /= prob_reads_d_given_a;
    summary_stat_diff /= prob_reads_d_given_a;
}


void SequenceProb::PrintReads(ReadData data) {
    printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}


double SequenceProb::CalculateExpectedValueForMu(Array10D summary_stat_AtoD) {

//    ancestor_genotypes;
    double sum = 0;
    double nume_sum = 0;
    for (int i = 0; i < 4; ++i) {
        for(int j=i;j< 4;++j) {
            int index10 = LookupTable::index_converter_16_to_10[i][j];
            int index16 = i*4+j;
            double denominator = ancestor_genotypes[index16]* frequency_prior[index10];
            sum += denominator;

//            nume_sum += denominator * summary_stat_AtoD(index10); //
//            cout << i<<":"<<j<< " " << index10 << ":\t" << ancestor_genotypes(index16) << "\t" << summary_stat_AtoD(index10) << endl;
        }
    }

    cout << nume_sum << "\t" << sum << "\t" << (nume_sum/sum) << endl;
    return nume_sum/sum;


}

double SequenceProb::Maximisation(double summery_stat) {
    double beta = 0;
    double max_hat = -beta * log(1- summery_stat / beta);
    return max_hat;

}

void SequenceProb::UpdateTransitionMatrix(EvolutionModel evo_model) {
    transition_matrix_a_to_d = evo_model.GetTranstionMatirxAToD();
}
