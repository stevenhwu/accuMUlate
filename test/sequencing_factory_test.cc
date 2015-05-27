//
// Created by steven on 4/2/15.
//


#include <mutations/sequencing_factory.h>
#include "gtest/gtest.h"
#include "unittest_utils.h"

class SequenceFactorTest : public ::testing::Test {

public:

protected:
    SequenceFactorTest() : base_custom1(0, 0), base_custom10(0, 0) {}

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





TEST_F(SequenceFactorTest, TestInit){


    GenomeData genome_data;
    int descendant_count = 5;
    int site_count = 1000;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0.5);
    GenomeData genome_data2 = genome_data;
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
    std::vector<SiteGenotypesIndex> sp;
    SequencingFactory sequencing_factory (params_not_equal);
    sequencing_factory.CreateSequenceProbsVector(sp, genome_data);
//    cout << "Time init seq latest: " << ((clock() - t1) / CLOCKS_PER_SEC) << "\t" << (clock() - t1) << endl;

    std::vector<SequenceProb> sp_expected;
    sp_expected.reserve(genome_data2.size());
    for (auto data : genome_data2) {
        sp_expected.emplace_back(data, params_not_equal);
    }
    MutationProb mutation_prob (params_not_equal);
    Array10D expected_ancestor_prior = mutation_prob.GetAncestorPrior();


    std::vector<DiploidProbsIndex10> convertIndexToDiploidIndex10 = sequencing_factory.RemoveConvertIndexKeyToDiploidIndex10();
    std::vector<HaploidProbs> convertIndexToHaploid = sequencing_factory.RemoveConvertIndexKeyToHaploid();
    
    for (int i = 0; i < sp_expected.size(); ++i) {

        ModelInput &site_data = genome_data2[i];

        SiteGenotypesIndex &factory = sp[i];
        auto &expected = sp_expected[i];
        auto expected_genotype_16 = expected.GetAncestorGenotypes();
        DiploidProbsIndex10 expected_genotype;
        for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
            int index16 = LookupTable::index_converter_10_to_16[index10];
            expected_genotype[index10] = expected_genotype_16[index16] * expected_ancestor_prior[index10];
        }

        auto foctory_genotype = convertIndexToDiploidIndex10[factory.GetAncestorIndex()];
        ASSERT_GENOTYPES(expected_genotype, foctory_genotype);

        std::vector<uint32_t> const &descendantIndex = factory.GetDescendantIndex();
        for (int j = 0; j < expected.GetDescendantCount(); ++j) {
            auto foctory_genotype = convertIndexToHaploid[descendantIndex[j]];
            auto expected_genotype = expected.GetDescendantGenotypes(j);
            ASSERT_GENOTYPES(expected_genotype, foctory_genotype);
        }
    }


}
