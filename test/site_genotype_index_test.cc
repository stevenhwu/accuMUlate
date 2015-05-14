//
// Created by steven on 4/8/15.
//

//
// Created by steven on 4/2/15.
//

#include <mutations/sequencing_factory.h>
#include "gtest/gtest.h"
#include "unittest_utils.h"
class SiteGenotypeIndexTest: public ::testing::Test {

public:

protected:
    const std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};
    ModelParams params_not_equal;
    double mu_4 = 0.0001;

    virtual void SetUp() {
        params_not_equal = {0.01, freq_not_equal, mu_4, 0.01, 0.01, 0.01};
    }
};





TEST_F(SiteGenotypeIndexTest, TestInit){

    GenomeData genome_data;
    int descendant_count = 5;
    int site_count = 10;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0.5);
    GenomeData genome_data2 = genome_data;

    std::vector<SiteGenotypesIndex> sp;
    SequencingFactory sequencing_factory (params_not_equal);
    sequencing_factory.CreateSequenceProbsVector(sp, genome_data);

    SequencingFactory sequencing_factory2 (params_not_equal);
    sequencing_factory2.CreateSequenceProbsVector(genome_data2);
    std::vector<SiteGenotypesIndex> sp2 = sequencing_factory2.RemoveSiteGenotypeIndexVector();

    for (int i = 0; i < sp.size(); ++i) {
        ASSERT_EQ(sp[i].GetAncestorIndex(), sp2[i].GetAncestorIndex());
        ASSERT_VECTOR_LIKE(sp[i].GetDescendantIndex(), sp2[i].GetDescendantIndex());
    }

}

TEST_F(SiteGenotypeIndexTest, TestMove){

    GenomeData genome_data;
    int descendant_count = 5;
    int site_count = 10;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0.5);

    std::vector<SiteGenotypesIndex> sp;
    SequencingFactory sequencing_factory (params_not_equal);
    sequencing_factory.CreateSequenceProbsVector(sp, genome_data);

    std::vector<SiteGenotypesIndex> expected_spi;

    SiteGenotypesIndex test_assignment(0);
    std::vector<uint32_t> expected_empty_vector;
    for (int i = 0; i < sp.size(); ++i) {
        expected_spi.push_back(sp[i]);
        ASSERT_EQ(sp[i].GetAncestorIndex(), expected_spi[i].GetAncestorIndex());
        ASSERT_VECTOR_LIKE(sp[i].GetDescendantIndex(), expected_spi[i].GetDescendantIndex());

        SiteGenotypesIndex temp_index = std::move(sp[i]);
        ASSERT_EQ(0, sp[i].GetAncestorIndex());
        ASSERT_VECTOR_LIKE(expected_empty_vector, sp[i].GetDescendantIndex());

        ASSERT_EQ(expected_spi[i].GetAncestorIndex(), temp_index.GetAncestorIndex());
        ASSERT_VECTOR_LIKE(expected_spi[i].GetDescendantIndex(), temp_index.GetDescendantIndex());

        test_assignment = std::move(temp_index);
        ASSERT_EQ(0, temp_index.GetAncestorIndex());
        ASSERT_VECTOR_LIKE(expected_empty_vector, temp_index.GetDescendantIndex());

        ASSERT_EQ(expected_spi[i].GetAncestorIndex(), test_assignment.GetAncestorIndex());
        ASSERT_VECTOR_LIKE(expected_spi[i].GetDescendantIndex(), test_assignment.GetDescendantIndex());


    }

}
