/*
 * SequencingProb.cc
 *
 *  Created on: Nov 7, 2014
 *      Author: Steven Wu
 */
#include <iostream>
#include "SequenceProb.h"

#include <string>

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

        all_norm_hap.push_back(normalised);
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
//			result[i+j*4] = 5;
//			result[j+i*4] = 3;
		}
		d.key = 0;
		d.reads[ref_allele] = 1;
		d.reads[i] += 2;
//		printf("D:%d %d: %u %u %u %u\n", i, i, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
		result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
	}
//	std::cout << result << std::endl;
//	Eigen::Matrix4d m;
//	for (int i = 0; i < result.rows(); ++i) {
//		m(i)= result[i];
//	}
//	std::cout << m<< std::endl;
	return result.exp();
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

//
//
//    ReadData data = all_descendant[descent_index];
//    HaploidProbs prob = all_hap[descent_index];
//    double likelihood = 0;
////    PrintReads(data);
//    for (int b = 0; b < 4; ++b) {
//
//    }

    HaploidProbs cond_likelihood = all_norm_hap[descent_index];
    return cond_likelihood;
//    std:: cout << descent_index << "\t" << prob.size()<< "\n" <<prob <<"\n==================================\n" << endl;



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
//index_summary = index_vector;


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

Array10D SequenceProb::CalculateAncestorToDescendant(Array4D prob_reads_given_descent) {
cout << "In Calc AtoD\n";
//    ReadData data = all_descendant[descent_index];
//    HaploidProbs prob = all_hap[descent_index];
    double likelihood = 0;
    cout << "BASE FREQ\n";
    cout << prob_reads_given_descent << "\t" << prob_reads_given_descent[0]<<
                "\t" << prob_reads_given_descent[1]<<
                "\t" << prob_reads_given_descent[2]<<
                "\t" << prob_reads_given_descent[3]<< endl;
//mu = 0.9;
//
////    for (int i = 0; i < 4; ++i) {
////        prob_reads_given_descent[i] = 0;
////    }
//    prob_reads_given_descent[0] = 0.1;
//    prob_reads_given_descent[1] = 0.1;
//    prob_reads_given_descent[2] = 0.7;
//    prob_reads_given_descent[3] = 0.1;

    double emu = computeExpFourThirdMu(mu);//FIXME still JC, change to F81 soon
    double three_prob[3] = {probThreeEqual(emu), probOneEqual(emu), probNotEqual(emu)};
    double sum = 0;
    double sum_stat = 0;
    double all_comb[4][10];
//    S_t(R, X ) = sum_j [ sum_stat / sum  ]
    double overall_summary[10];
    Array10D summary_stat;
    for (int a = 0; a < 10; ++a) {
        sum = 0;
//        cout << "==Loop A: "<<a << "\t" << genotype_lookup[a] <<endl;
        for (int b = 0; b < 4; ++b) {
//        prob_reads_given_descent[b] = 0.25;
//            cout << "====stat: b=" << b << "\tMu: " << mu << endl;
//
//        int index[10] = index_vector[b];

            all_comb[b][a] = three_prob[index_vector[b][a]] * prob_reads_given_descent[b];
            sum += all_comb[b][a];

//            cout << all_comb[b][a] << "\t";
        }
//        cout << endl;
        sum_stat = 0;
        for (int b = 0; b < 4; ++b) {
            all_comb[b][a]  *= index_summary[b][a];
            sum_stat += all_comb[b][a];
//            cout << b << ":" << a << "(" << all_comb[b][a] << ")" << index_summary[b][a] << "\t";
        }
//        cout << endl << ;
//        cout << endl << sum << "\t" << sum_stat << "\t" << (sum_stat/sum) << endl;
        overall_summary[a] = (sum_stat/sum);
        summary_stat[a] =  (sum_stat/sum);
    }

    cout << "\n========\n" ;
    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");
    cout << summary_stat.format(CommaInitFmt) << endl;
//    for (int i = 0; i < 10; ++i) {
//        cout << overall_summary[i] << " ";
//    }
    cout << endl;

    return summary_stat;
//        for (int i = 0; i < 4; ++i) {//anc_1
//        for (int j = 0; j <=i ; ++j) {//anc_2
//            cout << i << "\t" << j << "\t" << sum << endl;
//            for (int k = 0; k < 4; ++k) {
//                index_vector[k];
//
//
//            }
//        }
//    }
//    cout << "Result: " << sum << endl;

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
}


//FIXME, JC atm, no pi yet, add that and turn into F81
double SequenceProb::computeExpFourThirdMu(double mu) {
    return exp(-four_third*mu);
}

double SequenceProb::probThreeEqual(double emu) {
    return quarter + three_quarter* emu;
}

double SequenceProb::probOneEqual(double emu) {
    return quarter + quarter* emu;
}


double SequenceProb::probNotEqual(double emu) {
    return quarter - quarter* emu;
}

void SequenceProb::PrintReads(ReadData data) {
    printf("%d %d %d %d\n", data.reads[0], data.reads[1], data.reads[2], data.reads[3]);

}

DiploidProbs SequenceProb::getAncGenotypes() {
    return anc_genotypes;
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
