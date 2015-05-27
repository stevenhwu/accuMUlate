#include <iostream>
#include "F81.h"

F81::~F81() {
}

F81::F81(double mu) : EvolutionModel(MutationProb(mu)) {
//	cout << "C_F81_mu\n";
	freqs = mu_prob.GetFrequencyPrior();
	UpdateTransitionMatrix();
}

F81::F81(MutationProb &mutation_prob) : EvolutionModel(mutation_prob) {

	freqs = mu_prob.GetFrequencyPrior();
	UpdateTransitionMatrix();
}

F81::F81(F81 const &self) : EvolutionModel(self.GetMutationProb()) {
//	std::cout << "Call F81 Copy cons\n";
	mu_prob = self.GetMutationProb();
	freqs = mu_prob.GetFrequencyPrior();
	UpdateTransitionMatrix();
}

//
//F81::F81(double mu, Array4D freq) : EvolutionModel(mu) {
//	cout << "C_F81_array\n";
//    this->freqs = freq;
//	beta0 = MutationProb::CalculateBeta0(freqs);
//    UpdateTransitionMatrix();
//}
//
//F81::F81(double mu, vector<double> freq_v) : EvolutionModel(mu) {
//	cout << "C_F81_vector\n";
//	for (size_t i = 0; i < freqs.size(); ++i) {
//		freqs[i] = freq_v[i];
//	}
//
//	beta0 = MutationProb::CalculateBeta0(freqs);
//	UpdateTransitionMatrix();
//}


void F81::UpdateOneMinusExpBeta(double exp_beta) {
//FIXME: clarify expBeta or 1-expBeta, shoulde be 1-expBeta here
	double mu = mu_prob.ConvertOneMinusExpBetaToMu(exp_beta);
	UpdateMu(mu);
}



void F81::UpdateTransitionMatrix() {

    exp_beta = mu_prob.GetExpBeta();

	double one_minus_exp_beta_half = (1.0 - exp_beta)*0.5;

	for (int i : { 0, 1, 2, 3 }) {
		double prob = freqs[i] * one_minus_exp_beta_half;
		m.col(i) = Eigen::Vector4d::Constant(prob);
	}
	m.diagonal() += Eigen::Vector4d::Constant (exp_beta*0.5);//moved  *0.5 here


	for (int l = 0; l < 16; ++l) {
		auto index4_4 = LookupTable::index_converter_16_to_4_4[l];
		for (int k : { 0, 1, 2, 3 }) {
			transition_matrix_a_to_d(l, k) = (m(index4_4[0], k)+ m(index4_4[1], k));
		}
	}

//	std::cout << "F81:\t" << mu << "\t" << exp_beta << std::endl;

//    Eigen::Matrix4d m;
//	for (int i : { 0, 1, 2, 3 }) {
//        double prob = freqs[i] * one_minus_exp_beta;
//        for (int j : { 0, 1, 2, 3 }) {
//            m(j, i) = prob;
//			//MatrixXd::Constant(3,3,1.2)) replace with
//        }
//        m(i, i) += exp_beta;
//    }

//
//	Eigen::Matrix4d m;
//	for (int i : { 0, 1, 2, 3 }) {
//		double prob = freqs[i] * one_minus_exp_beta;
//		m.col(i) = Eigen::Vector4d::Constant(prob);
//	}
//	m.diagonal() += Eigen::Vector4d::Constant (exp_beta);;
//
//
//	//TODO: use col row to assign them
//    transition_matrix_a_to_d = MutationMatrix::Zero();
//    for (int i : { 0, 1, 2, 3 }) {
//        for (int j : { 0, 1, 2, 3 }) {
//            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
//            for (int k : { 0, 1, 2, 3 }) {
//                transition_matrix_a_to_d(index16, k) += 0.5 * m(i, k);
//                transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);
//            }
//        }
//    }
    
//    std::cout << "Calculate F81:\n" << mu <<  "\t" << beta0 << "\t" << exp_beta << "\n" <<
//             freqs[0] << freqs[1] << "\n" << m << "\n========================\n" << transition_matrix_a_to_d <<
//            "\n============================\n" << endl;

}




/* Old code from other places
    mutation = MutationAccumulation(params, false);
	MutationMatrix mt = MutationAccumulation(params, true);
	non_mutation = mutation - mt;

MutationMatrix MutationProb::MutationAccumulation2(bool and_mut) {

	using namespace Eigen;
	double exp_beta = CalculateExpBeta();
//	printf("beta0:%.10f\n", beta0); //0.9999999852
	Eigen::Matrix4d m;
	Eigen::Matrix4d m2 = Eigen::Matrix4d::Zero();


	Matrix3d m3 = Matrix3d::Identity();
	m3.row(1) = Vector3d(4,5,6);
//	std::cout << m3 << std::endl;

	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = frequency_prior[i] * (1.0 - exp_beta);
		}
		m(i, i) += exp_beta;
	}
//	element wise operation
	for (int i : { 0, 1, 2, 3 }) {
		double c = frequency_prior[i] * (1.0 - exp_beta);
//		c = (i+1) * 1.1;
		m2.row(i) = Vector4d::Constant(c);
	}
	m2.diagonal() += Vector4d::Constant(exp_beta);
//	std::cout << m << std::endl;
//	std::cout << "\n" << m2 << std::endl;

	//cerr << m << endl;
	MutationMatrix result;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			for (int k : { 0, 1, 2, 3 }) {
				result(i * 4 + j, k) = 0.0;
				if (!and_mut || i != k)
					result(i * 4 + j, k) += 0.5 * m(i, k);
				if (!and_mut || j != k)
					result(i * 4 + j, k) += 0.5 * m(j, k);
			}
		}
	}

	MutationMatrix result2 = MutationMatrix::Zero();;
		for (int i : { 0, 1, 2, 3 }) {
			for (int j : { 0, 1, 2, 3 }) {
				for (int k : { 0, 1, 2, 3 }) {
//					result(i * 4 + j, k) = 0.0;
					if (!and_mut || i != k)
						result2(i * 4 + j, k) += 0.5 * m2(i, k);
					if (!and_mut || j != k)
						result2(i * 4 + j, k) += 0.5 * m2(j, k);
				}
			}
		}
//	printf("%d\n", and_mut);
//	std::cout << result << std::endl;
//	std::cout << "\n" << std::endl;

//	std::cout << (result2 - result).sum() << "<-Should be 0" << std::endl;



	return result;
}

MutationMatrix MutationProb::MutationAccumulation(const ModelParams &params, bool and_mut) {

	double exp_beta = CalculateExpBeta();
//	printf("beta0:%.10f\n", beta0); //0.9999999852
	Eigen::Matrix4d m;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			m(i, j) = frequency_prior[i] * (1.0 - exp_beta);
		}
		m(i, i) += exp_beta;
	}
//	std::cout << m << std::endl;

	//cerr << m << endl;
	MutationMatrix result;
	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			for (int k : { 0, 1, 2, 3 }) {
				result(i * 4 + j, k) = 0.0;
				if (!and_mut || i != k)
					result(i * 4 + j, k) += 0.5 * m(i, k);
				if (!and_mut || j != k)
					result(i * 4 + j, k) += 0.5 * m(j, k);
			}
		}
	}
//	printf("%d\n", and_mut);
//	std::cout << result << std::endl;

	return result;
}




 */



void F81::OtherTransitionMatrix() {

	exp_beta = mu_prob.GetExpBeta();

	double one_minus_exp_beta_half = (1.0 - exp_beta) * 0.5;

	for (int i : {0, 1, 2, 3}) {
		double prob = freqs[i] * one_minus_exp_beta_half;
		m.col(i) = Eigen::Vector4d::Constant(prob);
	}
	m.diagonal() += Eigen::Vector4d::Constant(exp_beta * 0.5);//moved  *0.5 here

	//TODO: Need this to calculate zero mutation and estimate P(M=1)
	MutationMatrix transition_matrix_mutation_a_to_d = MutationMatrix::Zero();
	for (int l = 0; l < 16; ++l) {
		auto index4_4 = LookupTable::index_converter_16_to_4_4[l];

		for (int k : {0, 1, 2, 3}) {
			for (auto item : index4_4) {
				if( item != k ){
					transition_matrix_mutation_a_to_d(l, k) += m(item, k);
				}
			}
//			if( index4_4[0] != k ){
//				transition_matrix_mutation_a_to_d(l, k) += m(index4_4[0], k)
//			}
//			transition_matrix_a_to_d(l, k) = (m(index4_4[0], k) + m(index4_4[1], k));
//			l = i*4+j
//			if(i != k)
//				result(,k) += 0.5*m(i,k);
//			if( j != k)
//				result(i*4+j,k) += 0.5*m(j,k);
		}
	}
}
std::unique_ptr<EvolutionModel> F81::Clone() const {
//	printf("call F81 clone function\n");
	return std::unique_ptr<EvolutionModel> (new F81(*this) );
}

EvolutionModel* F81::Clone2() const {
//	printf("call F81 clone function raw pointer\n");
	return new F81(*this) ;
}

double F81::m1() {

//	Eigen::Matrix4d m;
	double one_minus_exp_beta_half = (1.0 - exp_beta)*0.5;
//	for (int i : { 0, 1, 2, 3 }) {
//		double prob = freqs[i] * one_minus_exp_beta;
//		for (int j : { 0, 1, 2, 3 }) {
//			m(j, i) = prob;
//			//MatrixXd::Constant(3,3,1.2)) replace with
//		}
//		m(i, i) += exp_beta;
//	}

	for (int i : { 0, 1, 2, 3 }) {
		double prob = freqs[i] * one_minus_exp_beta_half;
		m.col(i) = Eigen::Vector4d::Constant(prob);
//		m.col(i) = std::move(Eigen::Vector4d::Constant(prob));
	}
	m.diagonal() += Eigen::Vector4d::Constant (exp_beta*0.5);//moved  *0.5 here


//	transition_matrix_a_to_d = MutationMatrix::Zero();
//	for (int i : { 0, 1, 2, 3 }) {
//		for (int j : { 0, 1, 2, 3 }) {
//			int index16 = LookupTable::index_converter_4_4_to_16[i][j];
//			for (int k : { 0, 1, 2, 3 }) {
////				transition_matrix_a_to_d(index16, k) = 0;
//				transition_matrix_a_to_d(index16, k) = (m(i, k)+ m(j, k))*0.5;
////				transition_matrix_a_to_d(index16, k) *= 0.5;
////				transition_matrix_a_to_d(index16, k) += 0.5 * m(i, k);
////				transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);
//			}
//		}
//	}

	for (int l = 0; l < 16; ++l) {
			auto index4_4 = LookupTable::index_converter_16_to_4_4[l];
			for (int k : { 0, 1, 2, 3 }) {
//				transition_matrix_a_to_d(index16, k) = 0;
				transition_matrix_a_to_d(l, k) = (m(index4_4[0], k)+ m(index4_4[1], k));
//				transition_matrix_a_to_d(index16, k) *= 0.5;
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(i, k);
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);
		}
	}



//	transition_matrix_a_to_d *= 0.5;

//	transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);


	return transition_matrix_a_to_d.sum();

}

double F81::m2() {


	double one_minus_exp_beta = 1.0 - exp_beta;

//	for (int i : { 0, 1, 2, 3 }) {
//		double prob = freqs[i] * one_minus_exp_beta;
//		for (int j : { 0, 1, 2, 3 }) {
//			m(j, i) = prob;
//			//MatrixXd::Constant(3,3,1.2)) replace with
//		}
//		m(i, i) += exp_beta;
//	}
	for (int i : { 0, 1, 2, 3 }) {
		double prob = freqs[i] * one_minus_exp_beta;
		m.col(i) = Eigen::Vector4d::Constant(prob);
	}
	m.diagonal() += Eigen::Vector4d::Constant (exp_beta);;


//	transition_matrix_a_to_d = MutationMatrix::Zero();
//	for (int i : { 0, 1, 2, 3 }) {
//		for (int j : { 0, 1, 2, 3 }) {
//			int index16 = LookupTable::index_converter_4_4_to_16[i][j];
//			for (int k : { 0, 1, 2, 3 }) {
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(i, k);
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);
//			}
//		}
//	}


	for (int i : { 0, 1, 2, 3 }) {
		for (int j : { 0, 1, 2, 3 }) {
			int index16 = LookupTable::index_converter_4_4_to_16[i][j];
			for (int k : { 0, 1, 2, 3 }) {
//				transition_matrix_a_to_d(index16, k) = 0;
				transition_matrix_a_to_d(index16, k) = (m(i, k)+ m(j, k))*0.5;
//				transition_matrix_a_to_d(index16, k) *= 0.5;
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(i, k);
//				transition_matrix_a_to_d(index16, k) += 0.5 * m(j, k);
			}
		}
	}

	return transition_matrix_a_to_d.sum();
}
