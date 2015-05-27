
#include <evolution_models/F81.h>

#include <mutations/sequence_prob_v1.h>
#include <mutations/site_prob.h>

#include "gtest/gtest.h"


class SiteProbTest : public ::testing::Test {
public:


protected:
    SiteProbTest() : base_custom1(0, 0) {}

    const std::vector<double> freq_equal {0.25,0.25,0.25,0.25};     //beta0= 1.333333
    const std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    ModelParams params_equal;
    ModelParams params_not_equal;

    double mu = 0.1;
    double mu_4 = 0.0001;

    ModelInput base_custom3;//1 anc 3 des
    ModelInput base_custom1;//1 anc 1 des

    virtual void SetUp() {
        params_equal = {0.01, freq_equal, mu, 0.01, 0.01, 0.01};
        params_not_equal = {0.01, freq_not_equal, mu_4, 0.01, 0.01, 0.01};

//        base_custom1(2);
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
    }
};


TEST_F(SiteProbTest, TestCalculateOneDescendantGivenAncestor){

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;

    double new_mu =  (-log(0.5)) / (1.0 + (1.0 / 3.0));//p = 1-exp(beta*mu) = 0.5, => (1-p) = 0.5
    params_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom1, params_equal);
    MutationProb mutation_prob = MutationProb(params_equal);
    F81 evo_model(mutation_prob);

//    SiteProb site(sp, mutation_prob, evo_model);
    SiteProb site(sp, evo_model);
    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    //Anc = AA
    int anc_index = 0; //AA
    site.CalculateOneDescendantGivenAncestor(anc_index, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    double expected_prob = 0;
    double expected_diff = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*des_probs[b]; // P_anc_to_D(AA) * P(R|D)
        expected_diff += freq_equal[b] * 0.5 * des_probs[b];
    }

    expected_diff/= expected_prob;
    ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);

    //Anc = AC
    anc_index = 1; //AC
    site.CalculateOneDescendantGivenAncestor(anc_index, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    expected_diff = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*des_probs[b]; // P_anc_to_D(AA) * P(R|D)
        expected_diff += freq_equal[b] * 0.5 * des_probs[b];
    }

    expected_diff/= expected_prob;
    ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);

    //Anc = AT
    anc_index = 3; //AT
    site.CalculateOneDescendantGivenAncestor(anc_index, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    expected_diff = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*des_probs[b]; // P_anc_to_D(AA) * P(R|D)
        expected_diff += freq_equal[b] * 0.5 * des_probs[b];
    }

    expected_diff/= expected_prob;
    ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);
    /*R Code:
    freq = rep(0.25*4)
    beta0 = 1.3333, mu = 1/1.3333, p = 1-exp(-beta0*mu) = 0 => 1-p = 1
    same = P(R|D) * (1-p) * lookup
    diff = P(R|D) * p * freq[i]
    */
};



TEST_F(SiteProbTest, TestCalculateOneDescendantGivenAncestor2)
{

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;


    double new_mu = (-log(0.8)) / (1.0 + (1.0 / 3.0));//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_equal);
    MutationProb mutation_prob = MutationProb(params_equal);
    F81 evo_model(mutation_prob);

    SiteProb site(sp, evo_model);
    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();


    //Anc = AA
    int anc_index = 0; //AA
    site.CalculateOneDescendantGivenAncestor(anc_index, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    double expected_prob = 0;
    double expected_diff = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*des_probs[b]; // P_anc_to_D(AA) * P(R|D)
        expected_diff += freq_equal[b] * 0.2 * des_probs[b];
    }

    expected_diff /= expected_prob;
    ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);

    //Anc = AG
    anc_index = 2; //AG
    site.CalculateOneDescendantGivenAncestor(anc_index, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    expected_diff = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*des_probs[b]; // P_anc_to_D(AA) * P(R|D)
        expected_diff += freq_equal[b] * 0.2 * des_probs[b];
    }

    expected_diff /= expected_prob;
    ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);

}


TEST_F(SiteProbTest, TestCalculateOneDescendantGivenAncestor3) {

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;

    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    SiteProb site(sp, evo_model);

    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    for (int anc_index = 0; anc_index < ANCESTOR_COUNT; ++anc_index) {
        int index16 = LookupTable::index_converter_10_to_16[anc_index];
        site.CalculateOneDescendantGivenAncestor(index16, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

        double expected_prob = 0;
        double expected_diff = 0;
        for (int b = 0; b < 4; ++b) {
            double transition_a_to_d = matrix(index16, b);
            expected_prob += transition_a_to_d * des_probs[b]; // P_anc_to_D(AA) * P(R|D)
            expected_diff += freq_not_equal[b] * 0.2 * des_probs[b];
        }

        expected_diff /= expected_prob;
        ASSERT_DOUBLE_EQ(expected_prob, prob_reads_d_given_a);
        ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff);

    }



}


TEST_F(SiteProbTest, TestCalculateAllDescendantGivenAncestor) {

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;



    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, evo_model);
//    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);
//    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    double summary_stat_diff_ancestor = 0;
    double summary_stat_same_ancestor = 0;
    double sum_prob_ancestor = 0;
    for (int anc_index = 0; anc_index < ANCESTOR_COUNT; ++anc_index) {
        int index16 = LookupTable::index_converter_10_to_16[anc_index];
        site.CalculateAllDescendantGivenAncestor(index16, sum_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);

        double expected_prob = 1;
        double expected_same = 0;
        double expected_diff = 0;
        for (int d = 0; d < sp.GetDescendantCount(); ++d) {
            site.CalculateOneDescendantGivenAncestor(index16, d, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);
            expected_prob *= prob_reads_d_given_a;
            expected_diff += summary_stat_diff;

        }

        ASSERT_DOUBLE_EQ(expected_prob, sum_prob_ancestor);
        ASSERT_DOUBLE_EQ(expected_diff, summary_stat_diff_ancestor);
    }
}



TEST_F(SiteProbTest, TestCalculateAncestorToDescendant) {

    
    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom1, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, evo_model);

    double summary_stat_diff_ancestor = 0;
    double summary_stat_same_ancestor = 0;
    double sum_prob_ancestor = 0;

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;

    double expected_prob = 0;
    double expected_same = 0;
    double expected_diff = 0;

    HaploidProbs prob_reads_given_descent = sp.GetDescendantGenotypes(0);
    Array10D ancestor_prior = mutation_prob.GetAncestorPrior();
    DiploidProbs ancestor_genotypes = sp.GetAncestorGenotypes();
    for (int anc_index = 0; anc_index < 10; ++anc_index) {
        int index16 = LookupTable::index_converter_10_to_16[anc_index];
        site.CalculateAllDescendantGivenAncestor(index16, sum_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);
        site.CalculateOneDescendantGivenAncestor(index16, 0, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);
        ASSERT_DOUBLE_EQ(prob_reads_d_given_a, sum_prob_ancestor);
        ASSERT_DOUBLE_EQ(summary_stat_diff, summary_stat_diff_ancestor);


        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] *  sum_prob_ancestor;
        expected_prob += prob_reads_given_a;
        expected_diff += summary_stat_diff*prob_reads_given_a;

    }
    expected_diff /= expected_prob;


    double overall_stat_same;
    double overall_stat_diff;
    double prob;
    site.CalculateAncestorToDescendant(prob, overall_stat_same, overall_stat_diff);
    ASSERT_DOUBLE_EQ(expected_prob, prob);
    ASSERT_DOUBLE_EQ(expected_diff, overall_stat_diff);

}



TEST_F(SiteProbTest, TestCalculateAncestorToDescendant3) {


    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, evo_model);

    double summary_stat_diff_ancestor = 0;
    double summary_stat_same_ancestor = 0;
    double prod_prob_ancestor = 0;

    double summary_stat_diff_descendant = 0;
    double summary_stat_same_descendant = 0;
    double sum_prob_descendant = 0;


    double expected_prob = 0;
    double expected_same = 0;
    double expected_diff = 0;


    Array10D ancestor_prior = mutation_prob.GetAncestorPrior();
    DiploidProbs ancestor_genotypes = sp.GetAncestorGenotypes();
    for (int anc_index = 0; anc_index < 10; ++anc_index) {
        int index16 = LookupTable::index_converter_10_to_16[anc_index];

        double prob_reads_d_given_a = 1;
        double summary_stat_diff = 0;
        double summary_stat_same = 0;
        for (int d = 0; d < sp.GetDescendantCount(); ++d) {
            site.CalculateOneDescendantGivenAncestor(index16, d, sum_prob_descendant, summary_stat_same_descendant, summary_stat_diff_descendant);
            prob_reads_d_given_a *= sum_prob_descendant;
            summary_stat_diff += summary_stat_diff_descendant;
            summary_stat_same += summary_stat_same_descendant;
        }

        site.CalculateAllDescendantGivenAncestor(index16, prod_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);
        ASSERT_DOUBLE_EQ(prob_reads_d_given_a, prod_prob_ancestor);
        ASSERT_DOUBLE_EQ(summary_stat_diff, summary_stat_diff_ancestor);
        ASSERT_DOUBLE_EQ(summary_stat_same, summary_stat_same_ancestor);


        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] * prod_prob_ancestor;
        expected_prob += prob_reads_given_a;
        expected_diff += summary_stat_diff*prob_reads_given_a;
        expected_same += summary_stat_same*prob_reads_given_a;
    }

    expected_diff /= expected_prob * sp.GetDescendantCount();
    expected_same = 1 - expected_diff;

    double prob;
    double overall_stat_diff;
    double overall_stat_same;

    site.CalculateAncestorToDescendant(prob, overall_stat_same, overall_stat_diff);
    ASSERT_DOUBLE_EQ(expected_prob, prob);
    ASSERT_DOUBLE_EQ(expected_diff, overall_stat_diff);
    ASSERT_DOUBLE_EQ(expected_same, overall_stat_same);


}


TEST_F(SiteProbTest, TestFullCalculation){

    double mutation_rate = 0.2;
    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log( 1-mutation_rate )) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, evo_model);

    double site_prob = 0;
    double site_same = 0;
    double site_diff = 0;

    Array10D ancestor_prior = mutation_prob.GetAncestorPrior();
    DiploidProbs ancestor_genotypes = sp.GetAncestorGenotypes();

    double expected_prob = 0;
    double expected_diff = 0;
    for (int anc_index = 0; anc_index < 10; ++anc_index) {
        int index16 = LookupTable::index_converter_10_to_16[anc_index];
        double prod_prob_ancestor = 1;
        double summary_stat_diff_ancestor = 0;

        for (int des_index = 0; des_index < sp.GetDescendantCount(); ++des_index) {
            HaploidProbs des_probs = sp.GetDescendantGenotypes(des_index);
            double sum_prob_descendant = 0;
            double summary_stat_diff_descendant = 0;

            for (int b = 0; b < 4; ++b) {
                double transition_a_to_d = matrix(index16, b);
                sum_prob_descendant += transition_a_to_d * des_probs[b]; // P_anc_to_D(AA) * P(R|D)
                summary_stat_diff_descendant += freq_not_equal[b] * mutation_rate  * des_probs[b];
            }
            summary_stat_diff_descendant /= sum_prob_descendant;

            site.CalculateOneDescendantGivenAncestor(index16, des_index, site_prob, site_same, site_diff);
            ASSERT_DOUBLE_EQ(sum_prob_descendant, site_prob);
            ASSERT_DOUBLE_EQ(summary_stat_diff_descendant, site_diff);

            prod_prob_ancestor *= sum_prob_descendant;
            summary_stat_diff_ancestor += summary_stat_diff_descendant;
        }

        site.CalculateAllDescendantGivenAncestor(index16, site_prob, site_same, site_diff);
        ASSERT_DOUBLE_EQ(prod_prob_ancestor, site_prob);
        ASSERT_DOUBLE_EQ(summary_stat_diff_ancestor, site_diff);


        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] * prod_prob_ancestor;
        expected_prob += prob_reads_given_a;
        expected_diff += summary_stat_diff_ancestor * prob_reads_given_a;
    }
    expected_diff /= expected_prob * sp.GetDescendantCount();

    site.CalculateAncestorToDescendant(site_prob, site_same, site_diff);
    ASSERT_DOUBLE_EQ(expected_prob, site_prob);
    ASSERT_DOUBLE_EQ(expected_diff, site_diff);


}


TEST_F(SiteProbTest, TestFullCalculationSimulatedGenomeDate){

    double mutation_rate = 0.2;
    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log( 1-mutation_rate )) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(mutation_prob);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    GenomeData genome_data;
    int descendant_count = 10;
    int site_count = 100;
    SimulateGenomeData(genome_data, descendant_count, site_count, 0);

    std::vector<SequenceProb> sp_vector;
    std::vector<SiteProb> site_prob_vector;
    sp_vector.reserve(genome_data.size());
    for (int i = 0; i < genome_data.size(); ++i) {
        sp_vector.emplace_back(genome_data[i], params_not_equal);
        site_prob_vector.emplace_back(sp_vector[i], evo_model);
    }


    double site_prob = 0;
    double site_same = 0;
    double site_diff = 0;
    Array10D ancestor_prior = mutation_prob.GetAncestorPrior();
    for (int i = 0; i < site_prob_vector.size(); ++i) {

        SequenceProb sp = sp_vector[i];
        SiteProb site = site_prob_vector[i];
        DiploidProbs ancestor_genotypes = sp.GetAncestorGenotypes();


        double expected_prob = 0;
        double expected_diff = 0;
        for (int anc_index = 0; anc_index < 10; ++anc_index) {
            int index16 = LookupTable::index_converter_10_to_16[anc_index];
            double prod_prob_ancestor = 1;
            double summary_stat_diff_ancestor = 0;

            for (int des_index = 0; des_index < sp.GetDescendantCount(); ++des_index) {
                HaploidProbs des_probs = sp.GetDescendantGenotypes(des_index);
                double sum_prob_descendant = 0;
                double summary_stat_diff_descendant = 0;

                for (int b = 0; b < 4; ++b) {
                    double transition_a_to_d = matrix(index16, b);
                    sum_prob_descendant += transition_a_to_d * des_probs[b]; // P_anc_to_D(AA) * P(R|D)
                    summary_stat_diff_descendant += freq_not_equal[b] * mutation_rate * des_probs[b];
                }
                summary_stat_diff_descendant /= sum_prob_descendant;

                site.CalculateOneDescendantGivenAncestor(index16, des_index, site_prob, site_same, site_diff);
                ASSERT_NEAR(sum_prob_descendant, site_prob, ERROR_THRESHOLD);
                ASSERT_NEAR(summary_stat_diff_descendant, site_diff, ERROR_THRESHOLD);

                prod_prob_ancestor *= sum_prob_descendant;
                summary_stat_diff_ancestor += summary_stat_diff_descendant;
            }

            site.CalculateAllDescendantGivenAncestor(index16, site_prob, site_same, site_diff);
            ASSERT_NEAR(prod_prob_ancestor, site_prob, ERROR_THRESHOLD);
//            ASSERT_DOUBLE_EQ(summary_stat_diff_ancestor, site_diff);
            ASSERT_NEAR(summary_stat_diff_ancestor, site_diff, ERROR_THRESHOLD);

            double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] * prod_prob_ancestor;
            expected_prob += prob_reads_given_a;
            expected_diff += summary_stat_diff_ancestor * prob_reads_given_a;
        }
        expected_diff /= expected_prob * sp.GetDescendantCount();

        site.CalculateAncestorToDescendant(site_prob, site_same, site_diff);
        ASSERT_NEAR(expected_prob, site_prob, ERROR_THRESHOLD);
        ASSERT_NEAR(expected_diff, site_diff, ERROR_THRESHOLD);

    }

}