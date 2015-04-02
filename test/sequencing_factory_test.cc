//
// Created by steven on 4/2/15.
//

#include "gtest/gtest.h"
#include "sequencing_factory.h"
#include <random>
class SequenceFactorTest : public ::testing::Test {

public:
    constexpr static double ERROR_THRESHOLD = 1e-10;

protected:
    const std::vector<double> freq_equal {0.25,0.25,0.25,0.25};     //beta0= 1.333333
    const std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    ModelParams params_equal;
    ModelParams params_not_equal;


    double mu = 0.1;
    double mu_4 = 0.0001;
    ModelInput base_custom3;//1 anc 3 des
    ModelInput base_custom1;//1 anc 1 des
    ModelInput base_custom10;//1 anc 10 des

    virtual void SetUp() {
        params_equal = {0.01, freq_equal, mu, 0.01, 0.01, 0.01};
        params_not_equal = {0.01, freq_not_equal, mu_4, 0.01, 0.01, 0.01};

        base_custom1.reference = 0;
        ReadData r;
        for (int j = 0; j < 4; ++j) {
            r.reads[j] = (uint16_t) (j + 1);
        }
        base_custom1.all_reads.push_back(r);//anc
        base_custom1.all_reads.push_back(r);//des


        base_custom3.reference = 0;
        for (int i = 0; i < (3+1); ++i) {
            r.key=0;
            for (int j = 0; j < 4; ++j) {
                r.reads[j] = (uint16_t) (j + i + 1);
            }
            base_custom3.all_reads.push_back(r);
        }


        base_custom10.reference = 0;
        for (int i = 0; i < (1000+1); ++i) {
            r.key=0;
            for (int j = 0; j < 4; ++j) {
                r.reads[j] = (uint16_t) (j + i + 1);
            }
            base_custom10.all_reads.push_back(r);
        }
    }
};



void SimulateGenomeData(GenomeData &genome_data, int descendant_count, size_t fake_sample_count, double fake_prop) {

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_int_distribution<uint16_t> uniform_dist(0, 5);
    std::uniform_int_distribution<uint16_t> uniform3(0, 3);

    size_t fake_diff_count = fake_sample_count * fake_prop;
    descendant_count++;// for ancestor

    for (size_t s = 0; s < fake_sample_count; ++s) {

        ReadDataVector bcalls(descendant_count, ReadData{0});
        uint16_t ref_index = uniform3(e2);
        for (int i = 0; i < descendant_count; ++i) {
            bcalls[i].key=0;
            for (int j = 0; j < 4; ++j) {
                bcalls[i].reads[j] = uniform_dist(e2);
            }
            bcalls[i].reads[ref_index] = 20 + uniform_dist(e2);
        }


        if(s < fake_diff_count) { //diff ref
            uint16_t old_ref = ref_index;
            do{
                ref_index = uniform3(e2);
            }while (ref_index == old_ref);
            for (int j = 0; j < 4; ++j) {
                bcalls[0].reads[j] = uniform_dist(e2);
            }
            bcalls[0].reads[ref_index] = 20 + uniform_dist(e2);
        }

        genome_data.emplace_back( ModelInput{ref_index, bcalls} ) ;
    }

}

void ASSERT_GENOTYPES(HaploidProbs expected, HaploidProbs data){
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], SequenceFactorTest::ERROR_THRESHOLD );
    }
}
void ASSERT_GENOTYPES(DiploidProbs expected, DiploidProbs data){
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], data[i], SequenceFactorTest::ERROR_THRESHOLD );
    }
}

TEST_F(SequenceFactorTest, TestInit){



    GenomeData genome_data;
    int descendant_count = 5;
    int site_count = 1000;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0.5);

//    for (auto data: genome_data) {
//        std::cout << data.reference << std::endl;
//        for (auto readdata : data.all_reads) {
//            auto r = readdata.reads;
//            std::cout << r[0] << "\t" << r[1] << "\t" <<r[2] << "\t" <<r[3] << "\t" <<std::endl;
//        }
//        std::cout << std::endl;
//    }

//    sequencing_factory_v1.CreateSequenceProbV1(sp1, base_counts);
//    cout << "Time init seq v1: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;
//    cout <<  std::numeric_limits<double>::epsilon() << endl;
    std::vector<SequenceProb> sp;
    SequencingFactory sequencing_factory (params_not_equal);
    sequencing_factory.CreateSequenceProbsVector(sp, genome_data);
//    cout << "Time init seq latest: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;

    std::vector<SequenceProb> sp_expected;
    sp_expected.reserve(genome_data.size());
    for (auto data : genome_data) {
        sp_expected.emplace_back(data, params_not_equal);
    }

    for (int i = 0; i < sp_expected.size(); ++i) {

        ModelInput &site_data = genome_data[i];

        SequenceProb &factory = sp[i];
        SequenceProb &expected = sp_expected[i];
        auto foctory_genotype = factory.GetAncestorGenotypes();
        auto expected_genotype = expected.GetAncestorGenotypes();
        ASSERT_GENOTYPES(expected_genotype, foctory_genotype);
        for (int j = 0; j < expected.GetDescendantCount(); ++j) {
            auto foctory_genotype = factory.GetDescendantGenotypes(j);
            auto expected_genotype = expected.GetDescendantGenotypes(j);

            ASSERT_GENOTYPES(expected_genotype, foctory_genotype);
        }
    }


}
