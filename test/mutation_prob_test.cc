#include <mutations/model.h>
#include <mutations/mutation_prob.h>
#include "gtest/gtest.h"


class MutationProbTest : public ::testing::Test {
protected:
    int foo;
    ModelParams params_equal;
    ModelParams params_not_equal;

    std::vector<double> freq_equal {0.25,0.25,0.25,0.25};
    std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    double mu = 0.1;
    double mu_4 = 0.0001;

    virtual void SetUp() {
        params_equal = {
                0,//vm["theta"].as<double>(),
                freq_equal,//vm["nfreqs"].as<vector< double> >(),
                mu,//vm["mu"].as<double>(),
                0,//vm["seq-error"].as<double>(),
                0,//vm["phi-haploid"].as<double>(),
                0,//vm["phi-diploid"].as<double>(),
        };

        params_not_equal = {0, freq_not_equal, mu_4, 0, 0, 0};
    }
};




TEST_F(MutationProbTest, EqualFreqsInit) {

    MutationProb muProb = MutationProb(params_equal);
    Array4D d = muProb.GetFrequencyPrior();
    ASSERT_EQ(0.25, d[0]);
    ASSERT_EQ(0.25, d[1]);
    ASSERT_EQ(0.25, d[2]);
    ASSERT_EQ(0.25, d[3]);

    Array10D a10D = muProb.GetAncestorPrior();
    ASSERT_EQ(0.25*0.25, a10D[0]);
    ASSERT_EQ(0.25*0.25*2, a10D[1]);
    ASSERT_EQ(0.25*0.25*2, a10D[2]);
    ASSERT_EQ(0.25*0.25*2, a10D[3]);
    ASSERT_EQ(0.25*0.25, a10D[4]);
    ASSERT_EQ(0.25*0.25*2, a10D[5]);
    ASSERT_EQ(0.25*0.25*2, a10D[6]);
    ASSERT_EQ(0.25*0.25, a10D[7]);
    ASSERT_EQ(0.25*0.25*2, a10D[8]);
    ASSERT_EQ(0.25*0.25, a10D[9]);

    ASSERT_EQ(0.1, muProb.GetMu());
    ASSERT_DOUBLE_EQ(1.333333333333333, muProb.GetBeta0());
    double beta = muProb.GetExpBeta();

    ASSERT_EQ(0.87517331904294748401, beta);
    ASSERT_EQ(0.87517331904294748, beta); // ***474 fail, needs 474(8)
    ASSERT_DOUBLE_EQ(0.875173319042947, beta); // **4294 fail, needs **4294(7)


    double rate = muProb.GetMutationRate(); // exp(-1/ (1- 0.25^2*4) * 0.1)
    double expected_beta = 0.87517331904294748401;
//    ASSERT_DOUBLE_EQ(expected_beta, rate.one_minus_p);
    ASSERT_DOUBLE_EQ(1-expected_beta, rate);
    
}



TEST_F(MutationProbTest, NotEqualFreqsUpdateMu) {


    double new_mu = 0.001;
    MutationProb muProb = MutationProb(params_not_equal);
    muProb.UpdateMu( new_mu);

    Array4D d = muProb.GetFrequencyPrior();
    ASSERT_EQ(0.1, d[0]);
    ASSERT_EQ(0.2, d[1]);
    ASSERT_EQ(0.3, d[2]);
    ASSERT_EQ(0.4, d[3]);

    Array10D a10D = muProb.GetAncestorPrior();
    ASSERT_EQ(0.1*0.1, a10D[0]);
    ASSERT_EQ(0.1*0.2*2, a10D[1]);
    ASSERT_EQ(0.1*0.3*2, a10D[2]);
    ASSERT_EQ(0.1*0.4*2, a10D[3]);
    ASSERT_EQ(0.2*0.2, a10D[4]);
    ASSERT_EQ(0.2*0.3*2, a10D[5]);
    ASSERT_EQ(0.2*0.4*2, a10D[6]);
    ASSERT_EQ(0.3*0.3, a10D[7]);
    ASSERT_EQ(0.3*0.4*2, a10D[8]);
    ASSERT_EQ(0.4*0.4, a10D[9]);

    ASSERT_EQ(0.001, muProb.GetMu());
    ASSERT_DOUBLE_EQ(1.4285714285714286031, muProb.GetBeta0());

    double beta = muProb.GetExpBeta();
    double expected_beta = 0.99857244849385662366; //exp(-1/ (1- sum((1:4/10)^2)) * 0.001)

    ASSERT_DOUBLE_EQ(expected_beta, beta);

    double rate = muProb.GetMutationRate();
//    ASSERT_DOUBLE_EQ(expected_beta, rate.one_minus_p);
    ASSERT_DOUBLE_EQ(1-expected_beta, rate);

}



