/*
 * em_algorithm_mutation_test.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */

#include "em_model_mutation_v1.h"

#include "em_data_mutation_v1.h"
#include "em_algorithm_mutation_old.h"


class EmAlgorithmMutationTest : public ::testing::Test {
protected:
//
//    EmModelMutationV1 em_model_mutation;
    int total_sample_count = 12;
    ModelParams params_equal;
    ModelParams params_not_equal;

    std::vector<double> freq_equal {0.25,0.25,0.25,0.25};     //beta0= 1.333333
    std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    double mu = 0.1;
    double mu_4 = 0.0001;

    std::random_device rd;
    std::mt19937 e2;
    std::uniform_int_distribution<int> uniform_dist;

    virtual void SetUp() {
        params_equal = {
                0.01,//vm["theta"].as<double>(),
                freq_equal,//vm["nfreqs"].as<vector< double> >(),
                mu,//vm["mu"].as<double>(),
                0.01,//vm["seq-error"].as<double>(),
                0.01,//vm["phi-haploid"].as<double>(),
                0.01,//vm["phi-diploid"].as<double>(),
        };

        params_not_equal = {0.01, freq_not_equal, mu_4, 0.01, 0.01, 0.01};

        e2 = std::mt19937 (rd());
        uniform_dist = std::uniform_int_distribution<int> (0, 5);


    }
};


TEST_F(EmAlgorithmMutationTest, EmAlgorithmMutationTest1) {


//    MutationProb mutation_prob = MutationProb(0.1);
//    MutationProb mutation_prob = MutationProb(0.1);
    F81 evo_model0(1);

//    std::vector<SequenceProb> sp;
    std::vector<std::unique_ptr<EmData>> em_site_data;
//    size_t site_count = 10
    size_t fake_count = 100;
    double fake_prop = 0.01;
    size_t fake_diff_count = fake_count * fake_prop;
    for (size_t s = 0; s<fake_count; ++s) {

        ReadDataVector bcalls(total_sample_count, ReadData {0});
        uint16_t ref_index = 0;
        for (int i = 0; i < total_sample_count; ++i) {
            bcalls[i].key = 0;
            for (int j = 0; j < 4; ++j) {
                bcalls[i].reads[j] = (uint16_t) uniform_dist(e2);
            }
            bcalls[i].reads[0] = (uint16_t) 100 + uniform_dist(e2);
        }
        int i = 0;
        if (s < fake_diff_count) {
//            ref_index=3;
            bcalls[0].reads[0] = (uint16_t) uniform_dist(e2);
            bcalls[0].reads[3] = (uint16_t) 100 + uniform_dist(e2);
//            for (int i = 1; i < total_sample_count; ++i) {
//                bcalls[i].reads[1] = (uint16_t) (100 + uniform_dist(e2));
//                bcalls[i].reads[2] = (uint16_t) (100 + uniform_dist(e2));
//                bcalls[i].reads[0] = (uint16_t) uniform_dist(e2);
//            }
        }

        SequenceProb sp1 = SequenceProb(ModelInput{ref_index, bcalls}, params_equal);

        for (auto item : sp1.GetDescendantGenotypes()) {
            std::cout << item.format(nice_row) << std::endl;
        }
//        exit(34);
        em_site_data.emplace_back(new EmDataMutationV1(sp1, evo_model0));
    }


    EmModelMutationV1 em_model0(evo_model0);



//    EmAlgorithmMutationV1 em_alg (2, site_prob, model, em_site_data, em_model0);
    EmAlgorithmMutation em_alg2(2, em_site_data, em_model0);
    em_alg2.Run();

    em_alg2.PrintSummary();

exit(35);

}

