#include "gtest/gtest.h"
#include "evolution_models/F81.h"

class F81Test : public ::testing::Test {

public:
    const double ERROR_THRESHOLD = 1e-15;

protected:
    const std::vector<double> freq_equal {0.25,0.25,0.25,0.25};     //beta0= 1.333333
    const std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    const double freq_equal_beta0 = (1.0 + (1.0 / 3.0));
    const double freq_not_equal_beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );

    MutationProb mu_prob_equal;
    MutationProb mu_prob_not_equal;


    double mu = (-log(0.8)) / freq_equal_beta0; //p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    double mu_not_equal = (-log(0.8)) / freq_not_equal_beta0; //p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8

    MutationMatrix expected_freq_equal_matrix = MutationMatrix::Zero();

    virtual void SetUp() {
        ModelParams params_equal;
        ModelParams params_not_equal;

        params_equal = {0.01, freq_equal, mu, 0.01, 0.01, 0.01};
        params_not_equal = {0.01, freq_not_equal, mu_not_equal, 0.01, 0.01, 0.01};

        mu_prob_equal = MutationProb(params_equal);
        mu_prob_not_equal = MutationProb(params_not_equal);

        // A -> A     = 0.2*0.25 + 0.8 = 0.85   => 0.85 * 0.5 = 0.425
        // A -> C/G/T = 0.2*0.25       = 0.05   => 0.05 * 0.5 = 0.025

        double s2 = 0.85; // same_same = 0.425 + 0.425
        double sd = 0.45; // same_diff = 0.425 + 0.025
        double d2 = 0.05; // diff_diff = 0.025 + 0.025

        //                                     A   C   G   T
        expected_freq_equal_matrix.row(0)  << s2, d2, d2, d2; //AA
        expected_freq_equal_matrix.row(1)  << sd, sd, d2, d2; //AC
        expected_freq_equal_matrix.row(2)  << sd, d2, sd, d2; //AG
        expected_freq_equal_matrix.row(3)  << sd, d2, d2, sd; //AT
        expected_freq_equal_matrix.row(4)  << sd, sd, d2, d2; //CA
        expected_freq_equal_matrix.row(5)  << d2, s2, d2, d2; //CC
        expected_freq_equal_matrix.row(6)  << d2, sd, sd, d2; //CG
        expected_freq_equal_matrix.row(7)  << d2, sd, d2, sd; //CT
        expected_freq_equal_matrix.row(8)  << sd, d2, sd, d2; //GA
        expected_freq_equal_matrix.row(9)  << d2, sd, sd, d2; //GC
        expected_freq_equal_matrix.row(10) << d2, d2, s2, d2; //GG
        expected_freq_equal_matrix.row(11) << d2, d2, sd, sd; //GT
        expected_freq_equal_matrix.row(12) << sd, d2, d2, sd; //TA
        expected_freq_equal_matrix.row(13) << d2, sd, d2, sd; //TC
        expected_freq_equal_matrix.row(14) << d2, d2, sd, sd; //TG
        expected_freq_equal_matrix.row(15) << d2, d2, d2, s2; //TT


    };
};


TEST_F(F81Test, F81EqualDoubleConstructorTest) {

    F81 model (mu);
    MutationProb prob = model.GetMutationProb();
    double expected_exp_beta = 0.8;
    ASSERT_DOUBLE_EQ(mu, prob.GetMu());
    ASSERT_DOUBLE_EQ(expected_exp_beta, prob.GetExpBeta());

    double mu_rate = prob.GetMutationRate();
//    ASSERT_DOUBLE_EQ(expected_exp_beta, mu_rate.one_minus_p);
    ASSERT_DOUBLE_EQ(1-expected_exp_beta, mu_rate);

    MutationMatrix matrix = model.GetTranstionMatirxAToD();

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
            for (int k = 0; k < 4; ++k) {
                ASSERT_DOUBLE_EQ(expected_freq_equal_matrix(index16, k), matrix(index16, k));
            }
        }
    }
}

TEST_F(F81Test, F81EqualTest) {

    F81 model (mu_prob_equal);
    MutationMatrix matrix = model.GetTranstionMatirxAToD();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
            for (int k = 0; k < 4; ++k) {
                ASSERT_DOUBLE_EQ(expected_freq_equal_matrix(index16, k), matrix(index16, k));
            }
        }
    }
}


TEST_F(F81Test, F81NotEqualTest) {

    F81 model (mu_prob_not_equal);
    MutationMatrix matrix = model.GetTranstionMatirxAToD();
    MutationMatrix expected_matrix = MutationMatrix::Zero();
    // A -> A     = 0.2*0.1  + 0.8 = 0.82   => 0.82 * 0.5 = 0.41
    // C -> C     = 0.2*0.2  + 0.8 = 0.84   => 0.84 * 0.5 = 0.42
    // G -> G     = 0.2*0.3  + 0.8 = 0.86   => 0.86 * 0.5 = 0.43
    // T -> T     = 0.2*0.4  + 0.8 = 0.88   => 0.88 * 0.5 = 0.44

    // X -> A/C/G/T = 0.2*0.1/.2/.3/.4    = 0.02/.04/,06/,0.8   => 0.x * 0.5 = 0.01/.02/.03/.04

    double s[4] = {0.41, 0.42, 0.43, 0.44};
    double d[4] = {0.01, 0.02, 0.03, 0.04};

    //                             A          C          G          T
    expected_matrix.row(0)  << s[0]+s[0], d[1]+d[1], d[2]+d[2], d[3]+d[3]; //AA
    expected_matrix.row(1)  << s[0]+d[0], d[1]+s[1], d[2]+d[2], d[3]+d[3]; //AC
    expected_matrix.row(2)  << s[0]+d[0], d[1]+d[1], d[2]+s[2], d[3]+d[3]; //AG
    expected_matrix.row(3)  << s[0]+d[0], d[1]+d[1], d[2]+d[2], d[3]+s[3]; //AT

    expected_matrix.row(4)  << d[0]+s[0], s[1]+d[1], d[2]+d[2], d[3]+d[3]; //CA
    expected_matrix.row(5)  << d[0]+d[0], s[1]+s[1], d[2]+d[2], d[3]+d[3]; //CC
    expected_matrix.row(6)  << d[0]+d[0], s[1]+d[1], d[2]+s[2], d[3]+d[3]; //CG
    expected_matrix.row(7)  << d[0]+d[0], s[1]+d[1], d[2]+d[2], d[3]+s[3]; //CT

    expected_matrix.row(8)  << d[0]+s[0], d[1]+d[1], s[2]+d[2], d[3]+d[3]; //GA
    expected_matrix.row(9)  << d[0]+d[0], d[1]+s[1], s[2]+d[2], d[3]+d[3]; //GC
    expected_matrix.row(10) << d[0]+d[0], d[1]+d[1], s[2]+s[2], d[3]+d[3]; //GG
    expected_matrix.row(11) << d[0]+d[0], d[1]+d[1], s[2]+d[2], d[3]+s[3]; //GT

    expected_matrix.row(12) << d[0]+s[0], d[1]+d[1], d[2]+d[2], s[3]+d[3]; //TA
    expected_matrix.row(13) << d[0]+d[0], d[1]+s[1], d[2]+d[2], s[3]+d[3]; //TC
    expected_matrix.row(14) << d[0]+d[0], d[1]+d[1], d[2]+s[2], s[3]+d[3]; //TG
    expected_matrix.row(15) << d[0]+d[0], d[1]+d[1], d[2]+d[2], s[3]+s[3]; //TT

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
            for (int k = 0; k < 4; ++k) {
                ASSERT_DOUBLE_EQ(expected_matrix(index16, k), matrix(index16, k));
            }
        }
    }
}


TEST_F(F81Test, F81NotEqualUpdateTest) {

    F81 model (mu_prob_not_equal);
    double new_mu = (-log(0.4)) / freq_not_equal_beta0; //p = 1-exp(beta*mu) = 0.6, => (1-p) = 0.4
    model.UpdateMu(new_mu);

    MutationProb prob = model.GetMutationProb();
    double expected_exp_beta = 0.4;
    ASSERT_DOUBLE_EQ(new_mu, prob.GetMu());
    ASSERT_DOUBLE_EQ(expected_exp_beta, prob.GetExpBeta());

//    MutationRate mu_rate = prob.GetMutationRate();
//    ASSERT_DOUBLE_EQ(expected_exp_beta, mu_rate.one_minus_p);
//    ASSERT_DOUBLE_EQ(1-expected_exp_beta, mu_rate.prob);
    double mu_rate = prob.GetMutationRate();
//    ASSERT_DOUBLE_EQ(expected_exp_beta, mu_rate.one_minus_p);
    ASSERT_DOUBLE_EQ(1-expected_exp_beta, mu_rate);

    MutationMatrix matrix = model.GetTranstionMatirxAToD();
    MutationMatrix expected_matrix = MutationMatrix::Zero();
    // A -> A     = 0.6*0.1  + 0.4 = 0.46   => 0.46 * 0.5 = 0.23
    // C -> C     = 0.6*0.2  + 0.4 = 0.52   => 0.52 * 0.5 = 0.26
    // G -> G     = 0.6*0.3  + 0.4 = 0.58   => 0.58 * 0.5 = 0.29
    // T -> T     = 0.6*0.4  + 0.4 = 0.64   => 0.64 * 0.5 = 0.32

    // X -> A/C/G/T = 0.6*0.1/.2/.3/.4    = 0.06/.12/,18/,24   => 0.x * 0.5 = 0.03/.06/.09/.12

    double s[4] = {0.23, 0.26, 0.29, 0.32};
    double d[4] = {0.03, 0.06, 0.09, 0.12};

    //                             A          C          G          T
    expected_matrix.row(0)  << s[0]+s[0], d[1]+d[1], d[2]+d[2], d[3]+d[3]; //AA
    expected_matrix.row(1)  << s[0]+d[0], d[1]+s[1], d[2]+d[2], d[3]+d[3]; //AC
    expected_matrix.row(2)  << s[0]+d[0], d[1]+d[1], d[2]+s[2], d[3]+d[3]; //AG
    expected_matrix.row(3)  << s[0]+d[0], d[1]+d[1], d[2]+d[2], d[3]+s[3]; //AT

    expected_matrix.row(4)  << d[0]+s[0], s[1]+d[1], d[2]+d[2], d[3]+d[3]; //CA
    expected_matrix.row(5)  << d[0]+d[0], s[1]+s[1], d[2]+d[2], d[3]+d[3]; //CC
    expected_matrix.row(6)  << d[0]+d[0], s[1]+d[1], d[2]+s[2], d[3]+d[3]; //CG
    expected_matrix.row(7)  << d[0]+d[0], s[1]+d[1], d[2]+d[2], d[3]+s[3]; //CT

    expected_matrix.row(8)  << d[0]+s[0], d[1]+d[1], s[2]+d[2], d[3]+d[3]; //GA
    expected_matrix.row(9)  << d[0]+d[0], d[1]+s[1], s[2]+d[2], d[3]+d[3]; //GC
    expected_matrix.row(10) << d[0]+d[0], d[1]+d[1], s[2]+s[2], d[3]+d[3]; //GG
    expected_matrix.row(11) << d[0]+d[0], d[1]+d[1], s[2]+d[2], d[3]+s[3]; //GT

    expected_matrix.row(12) << d[0]+s[0], d[1]+d[1], d[2]+d[2], s[3]+d[3]; //TA
    expected_matrix.row(13) << d[0]+d[0], d[1]+s[1], d[2]+d[2], s[3]+d[3]; //TC
    expected_matrix.row(14) << d[0]+d[0], d[1]+d[1], d[2]+s[2], s[3]+d[3]; //TG
    expected_matrix.row(15) << d[0]+d[0], d[1]+d[1], d[2]+d[2], s[3]+s[3]; //TT


    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int index16 = LookupTable::index_converter_4_4_to_16[i][j];
            for (int k = 0; k < 4; ++k) {
                ASSERT_DOUBLE_EQ(expected_matrix(index16, k), matrix(index16, k));
            }
        }
    }
}