/*
 * SequencingProb.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
#include "SequenceProb.h"
//#include "models/JC69.h"



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


Eigen::IOFormat nice_row(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");

const bool DEBUG = true;

SequenceProb::SequenceProb(const ModelParams &model_params,
		const ModelInput site_data, MutationProb muProb) {

	params = model_params; //TODO check copy method/ref

    double beta0 = 0;
    for (int i = 0; i < 4; ++i) {
        beta0 += params.nuc_freq[i] * params.nuc_freq[i];
        for (int j = i; j < 4; ++j) {
            int index10 = index_convert_16_to_10[i][j];
            frequency_prior[index10] = params.nuc_freq[i] * params.nuc_freq[j];
            if(i != j){
                frequency_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }
    beta = 1-beta0;


	UpdateMuProb(muProb);
	data = site_data;
	pop_genotypes = DiploidPopulation(site_data.reference);

    int descendant_count = site_data.all_reads.size();

    ancestor = site_data.all_reads[0];
	anc_genotypes = DiploidSequencing(ancestor);
	anc_genotypes *= pop_genotypes;
//	DiploidProbs num_genotypes = anc_genotypes;

    for (int i = 1; i < descendant_count; ++i) {
        ReadData data = site_data.all_reads[i];
        all_descendant.push_back(data);

        HaploidProbs p = HaploidSequencing(data);
        all_hap.push_back(p);

        HaploidProbs normalised = p/p.sum();
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
    mu = muProb.GetMu();
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


int index_vector[4][10] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
        {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
        {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
        {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
        {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T
};
double index_summary[4][10] ={
        {0, 1, 1, 1, 2, 2, 2, 2, 2, 2},// A
        {2, 1, 2, 2, 0, 1, 1, 2, 2, 2},// C
        {2, 2, 1, 2, 2, 1, 2, 0, 1, 2},// G
        {2, 2, 2, 1, 2, 2, 1, 2, 1, 0} // T

//        {0, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1},// A
//        {1, 0.5, 1, 1, 0, 0.5, 0.5, 1, 1, 1},// C
//        {1, 1, 0.5, 1, 1, 0.5, 1, 0, 0.5, 1},// G
//        {1, 1, 1, 0.5, 1, 1, 0.5, 1, 0.5, 0} // T
};
/*
10 cat version
  AA AC AG AT CC CG CT GG GT TT

  descnt = A
      A  C  G  T
 *AA == !! !! !!
  AC =! =! !! !!
  AG =! !! =! !!
  AT =! !! !! =!
 *CC !! == !! !!
  CG !! =! =! !!
  CT !! =! !! =!
 *GG !! !! == !!
  GT !! !! =! =!
 *TT !! !! !! ==


*/
int index_16_to_10_single[16] = {
    0,1,2,3,
    1,4,5,6,
    2,5,7,8,
    3,6,8,9,
};
int index_10_to_16[10] = {//??
        0,1,2, 3,
          5,6, 7,
            10,11,
               15
};

string genotype_lookup[10]= {
        "AA", "AC", "AG", "AT",
        "CC", "CG", "CT",
        "GG", "GT", "TT"
};

//int index_vector2[10][14] = { //Same as summary stat now, number of mismatches between D (o) and A(m,f)
//        {0,1,1,1,2,2,2,2,2,2},// A
//        {2,1,2,2,0,1,1,2,2,2},// C
//        {2,2,1,2,2,1,2,0,1,2},// G
//        {2,2,2,1,2,2,1,2,1,0} // T
//
//};

Array10D SequenceProb::CalculateAncestorToDescendant(Array4D prob_reads_given_descent, EvolutionModel *error_model) {
    cout << "==================In Calc AtoD====================\n";

    double likelihood = 0;
//    mu = 0.9;
//    prob_reads_given_descent[0] = 0.2;
//    prob_reads_given_descent[1] = 0.2;
//    prob_reads_given_descent[2] = 0.4;
//    prob_reads_given_descent[3] = 0.2;
    if (DEBUG) {
        cout << "MU: " << mu <<" BASE FREQ: " << prob_reads_given_descent.format(nice_row) << endl;
    }

    //JC69
//    array<double, 3> conditional_prob_jc = (JC69) error_model->GetConditionalProbSpecial();
//    double three_prob[3] = {0,0,0};

    MutationMatrix conditional_prob = error_model->GetConditionalProb();
    double sum_over_d = 0;
    double sum_stat = 0;
    double all_comb[10][4];
//    S_t(R, X ) = sum_j [ sum_stat / sum  ]
    double overall_summary[10];
    Array10D summary_stat;
    
    for (int a = 0; a < 10; ++a) {//10 or 16??
        sum_over_d = 0;
        sum_stat = 0;
        if(DEBUG) {
            cout << "==Loop A: " << a << "\t" << genotype_lookup[a] << endl;
        }

        for (int d = 0; d < 4; ++d) {
            int index16 = index_10_to_16[a];
             // 16 to 10 to 16
            all_comb[a][d] = conditional_prob(index16, d) * prob_reads_given_descent[d];
            sum_over_d += all_comb[a][d];

            double summary_given_a_d = all_comb[a][d] * index_summary[d][a];
            sum_stat += summary_given_a_d;
            if(DEBUG){
                cout << "==== Loop D: " << d << "\t" << all_comb[a][d] << "\t" << summary_given_a_d <<endl;
            }
        }
        overall_summary[a] = (sum_stat/ sum_over_d);
        summary_stat[a] =  (sum_stat/ sum_over_d);
    }
    if(DEBUG) {
        cout << summary_stat.format(nice_row) << endl << "================END=============\n\n";
    }

    return summary_stat;

}



void SequenceProb::PrintReads(ReadData data) {
    printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}


double SequenceProb::CalculateExpectedValueForMu(Array10D summary_stat_AtoD) {

    anc_genotypes;
    double sum = 0;
    double nume_sum = 0;
    for (int i = 0; i < 4; ++i) {
        for(int j=i;j< 4;++j) {
            int index10 = index_convert_16_to_10[i][j];
            int index16 = i*4+j;
            double denominator = anc_genotypes(index16)* frequency_prior[index10];
            sum += denominator;

            nume_sum += denominator * summary_stat_AtoD(index10); // TODO: add P(X) freq_1 X freq_2
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
