//
// Created by steven on 4/8/15.
//
#include <mutations/model.h>
#include "gtest/gtest.h"
#include "unittest_utils.h"



class ModelTest : public ::testing::Test {

public:

protected:
    ModelTest() : base_custom1(0, 0), base_custom10(0, 0) {}

    const std::vector<double> freq_equal{0.25, 0.25, 0.25, 0.25};     //beta0= 1.333333
    const std::vector<double> freq_not_equal{0.1, 0.2, 0.3, 0.4};

    ModelParams params_equal;
    ModelParams params_not_equal;


    double mu = 0.1;
    double mu_4 = 0.0001;
    ModelInput base_custom3;
    //1 anc 3 des
    ModelInput base_custom1;
    //1 anc 1 des
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
        for (int i = 0; i < (3 + 1); ++i) {
            r.key = 0;
            for (int j = 0; j < 4; ++j) {
                r.reads[j] = (uint16_t) (j + i + 1);
            }
            base_custom3.all_reads.push_back(r);
        }


        base_custom10.reference = 0;
        for (int i = 0; i < (1000 + 1); ++i) {
            r.key = 0;
            for (int j = 0; j < 4; ++j) {
                r.reads[j] = (uint16_t) (j + i + 1);
            }
            base_custom10.all_reads.push_back(r);
        }
    }
};

TEST_F(ModelTest, TestModelInput) {

    GenomeData genome_data;
    int descendant_count = 5;
    int site_count = 10;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0.5);

    for (int i = 0; i < genome_data.size(); ++i) {

        ModelInput backup = genome_data[i];
        ModelInput m = genome_data[i];
        ASSERT_EQ(m.reference, genome_data[i].reference);
        for (int j = 0; j < m.all_reads.size(); ++j) {
            ASSERT_EQ(m.all_reads[j].key, genome_data[i].all_reads[j].key);
        }

        ModelInput moved_model = std::move(genome_data[i]);
        ASSERT_EQ(0, genome_data[i].reference);
        for (int j = 0; j < genome_data[i].all_reads.size(); ++j) {
            ASSERT_EQ(0, genome_data[i].all_reads[j].key);
        }
        ASSERT_EQ(moved_model.reference, m.reference);
        for (int j = 0; j < moved_model.all_reads.size(); ++j) {
            ASSERT_EQ(moved_model.all_reads[j].key, m.all_reads[j].key);
        }

        backup = std::move(moved_model);
        ASSERT_EQ(0, moved_model.reference);
        for (int j = 0; j < genome_data[i].all_reads.size(); ++j) {
            ASSERT_EQ(0, moved_model.all_reads[j].key);
        }
        ASSERT_EQ(backup.reference, m.reference);
        for (int j = 0; j < backup.all_reads.size(); ++j) {
            ASSERT_EQ(backup.all_reads[j].key, m.all_reads[j].key);
        }


    }


}

TEST_F(ModelTest, TestReadData) {

    ReadData rdata = {0};
    ASSERT_EQ(0, rdata.key);
    for (int i = 0; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }

    rdata = {10};
    ASSERT_EQ(10, rdata.key);
    for (int i = 1; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(10, rdata.reads[0]);

    rdata = {65535};
    ASSERT_EQ(65535, rdata.key);
    for (int i = 1; i < BASE_COUNT; ++i) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(65535, rdata.reads[0]);

    rdata = {65536};
    ASSERT_EQ(65536, rdata.key);
    for (int i : {0, 2, 3}) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(1, rdata.reads[1]);

    rdata = 0;
    rdata.reads[2] = 1;
    ASSERT_EQ(4294967296, rdata.key);

    rdata = 281474976710656;
    for (int i : {0, 1, 2}) {
        ASSERT_EQ(0, rdata.reads[i]);
    }
    ASSERT_EQ(1, rdata.reads[3]);

//    281474976710656 + 4294967296 + 65536 + 1
//    [1] 281479271743489
    rdata = 281479271743489;
    for (int i : {0, 1, 2, 3}) {
        ASSERT_EQ(1, rdata.reads[i]);
    }

}

