#include <SiteProb.h>
#include <evolutionModels/F81.h>

#include "gtest/gtest.h"
#include "SequenceProb.h"

class SiteProbTest : public ::testing::Test {
public:
    const double ERROR_THRESHOLD = 1e-10;
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
        cout << "setup\n";
    }
};


TEST_F(SiteProbTest, TestCalculateOneDescendantGivenAncestor){

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;

    HaploidProbs prob_reads_given_descent = {0.25,0.25,0.25,0.25}; //Fixed value for now

    double new_mu =  (-log(0.5)) / (1.0 + (1.0 / 3.0));//p = 1-exp(beta*mu) = 0.5, => (1-p) = 0.5
    params_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_equal);
    MutationProb mutation_prob = MutationProb(params_equal);
    F81 evo_model(new_mu, params_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

    //Anc = AA
    int anc_index = 0; //AA
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    double expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*0.25; // P_anc_to_D(AA) * P(R|D)

    }

    double expected_same = prob_reads_given_descent[0] * 0.5 * 1; //anc=AA, lookup = 1
    double expected_diff = 4 * 0.25 * 0.5 * 0.25;
    expected_same/= expected_prob;
    expected_diff/= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_EQ(expected_diff, summary_stat_diff);

    //Anc = AC
    anc_index = 1; //AC
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*0.25; // P_anc_to_D(AA) * P(R|D)
    }

    expected_same = prob_reads_given_descent[0] * 0.5 * 0.5 + prob_reads_given_descent[1] * 0.5 * 0.5; //anc=AC, lookup = 0.5, 0.5, 0, 0
    expected_diff = 4 * 0.25 * 0.5 * 0.25;
    expected_same/= expected_prob;
    expected_diff/= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_EQ(expected_diff, summary_stat_diff);

    //Anc = AT
    anc_index = 3; //AT
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index,b);
        expected_prob += transition_a_to_d*0.25; // P_anc_to_D(AA) * P(R|D)
    }

    expected_same = prob_reads_given_descent[0] * 0.5 * 0.5 + prob_reads_given_descent[3] * 0.5 * 0.5; //anc=AT, lookup = 0.5, 0, 0, 0.5
    expected_diff = 4 * 0.25 * 0.5 * 0.25;
    expected_same/= expected_prob;
    expected_diff/= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_EQ(expected_diff, summary_stat_diff);
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

    HaploidProbs prob_reads_given_descent = {0.4, 0.3, 0.2, 0.1}; //Fixed value for now

    double new_mu = (-log(0.8)) / (1.0 + (1.0 / 3.0));//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_equal);
    MutationProb mutation_prob = MutationProb(params_equal);
    F81 evo_model(new_mu, params_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

    //Anc = AA
    int anc_index = 0; //AA
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    double expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index, b);
        expected_prob += transition_a_to_d * prob_reads_given_descent[b]; // P_anc_to_D(AA) * P(R|D)
    }

    double expected_same = 0.4 * 0.8 * 1; //anc=AA, lookup = 1
    double expected_diff = (0.4+0.3+0.2+0.1) * 0.2 * 0.25;
    expected_same /= expected_prob;
    expected_diff /= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_NEAR(expected_diff, summary_stat_diff, ERROR_THRESHOLD);

    //Anc = AG
    anc_index = 2; //AG
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index, b);
        expected_prob += transition_a_to_d * prob_reads_given_descent[b]; // P_anc_to_D(AA) * P(R|D)
    }

    expected_same = 0.4 * 0.8 * 0.5 + 0.2 * 0.8 * 0.5; //anc=AG, lookup = 0.5, 0, 0.5, 0
    expected_diff = (0.4+0.3+0.2+0.1) * 0.2 * 0.25;
    expected_same /= expected_prob;
    expected_diff /= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_NEAR(expected_diff, summary_stat_diff, ERROR_THRESHOLD);

//    expected_diff = 0.2*0.25;
////    expected_diff /= expected_prob;
//    ASSERT_NEAR(expected_diff, summary_stat_diff*prob_reads_d_given_a, 1e-20);
//    expected_diff = 0.2*0.25;
//    expected_diff /= expected_prob;
//    ASSERT_NEAR(expected_diff, summary_stat_diff, 1e-20);


}


TEST_F(SiteProbTest, TestCalculateOneDescendantGivenAncestor3) {

    double summary_stat_same = 0;
    double summary_stat_diff = 0;
    double prob_reads_d_given_a = 0;

    HaploidProbs prob_reads_given_descent = {0.4, 0.3, 0.2, 0.1}; //Fixed value for now

    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(new_mu, params_not_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

    //Anc = AA
    int anc_index = 0; //AA
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    double expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index, b);
        expected_prob += transition_a_to_d * prob_reads_given_descent[b]; // P_anc_to_D(AA) * P(R|D)
    }

    double expected_same = 0.4 * 0.8 * 1; //anc=AA, lookup = 1
    double expected_diff = 0.4 * 0.2 * 0.1 + 0.3 * 0.2 * 0.2 + 0.2 * 0.2 * 0.3 + 0.1 * 0.2 * 0.4;
    expected_same /= expected_prob;
    expected_diff /= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_NEAR(expected_diff, summary_stat_diff, ERROR_THRESHOLD);


    //Anc = AT
    anc_index = 3; //AT
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(anc_index, b);
        expected_prob += transition_a_to_d * prob_reads_given_descent[b]; // P_anc_to_D(AT) * P(R|D)
    }

    expected_same = 0.4 * 0.8 * 0.5 + 0.1 * 0.8 * 0.5; //anc=AT, lookup = 0.5, 0, 0, 0.5
    expected_diff = 0.4 * 0.2 * 0.1 + 0.3 * 0.2 * 0.2 + 0.2 * 0.2 * 0.3 + 0.1 * 0.2 * 0.4;
    expected_same /= expected_prob;
    expected_diff /= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_NEAR(expected_diff, summary_stat_diff, ERROR_THRESHOLD);

    //Anc = CG
    anc_index = 5; //CG index10 = 5, index16=6
    site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);

    expected_prob = 0;
    for (int b = 0; b < 4; ++b) {
        double transition_a_to_d = matrix(6, b);//index_16
        expected_prob += transition_a_to_d * prob_reads_given_descent[b]; // P_anc_to_D(AT) * P(R|D)
    }

    expected_same = 0.3 * 0.8 * 0.5 + 0.2 * 0.8 * 0.5; //anc=AT, lookup = 0, 0.5, 0.5, 0
    expected_diff = 0.4 * 0.2 * 0.1 + 0.3 * 0.2 * 0.2 + 0.2 * 0.2 * 0.3 + 0.1 * 0.2 * 0.4;
    expected_same /= expected_prob;
    expected_diff /= expected_prob;
    ASSERT_EQ(expected_prob, prob_reads_d_given_a);
    ASSERT_EQ(expected_same, summary_stat_same);
    ASSERT_NEAR(expected_diff, summary_stat_diff, ERROR_THRESHOLD);



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
    F81 evo_model(new_mu, params_not_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

    //Anc = AA
    int anc_index = 0; //AA

    double summary_stat_diff_ancestor = 0;
    double summary_stat_same_ancestor = 0;
    double sum_prob_ancestor = 0;
    site.CalculateAllDescendantGivenAncestor(anc_index, sum_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);


    double expected_prob = 1;
    double expected_same = 0;
    double expected_diff = 0;
    for (int i = 0; i < sp.GetDescendantCount(); ++i) {
        HaploidProbs prob_reads_given_descent = sp.GetDescendantGenotypes(i);
        site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);
        expected_prob *= prob_reads_d_given_a;
        expected_same += summary_stat_same;
        expected_diff += summary_stat_diff;
    }

    ASSERT_EQ(expected_prob, sum_prob_ancestor);
    ASSERT_EQ(expected_same, summary_stat_same_ancestor);
    ASSERT_EQ(expected_diff, summary_stat_diff_ancestor);

}



TEST_F(SiteProbTest, TestCalculateAncestorToDescendant) {

    
    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom1, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(new_mu, params_not_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

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
        site.CalculateAllDescendantGivenAncestor(anc_index, sum_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);
        site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, prob_reads_d_given_a, summary_stat_same, summary_stat_diff);
        ASSERT_EQ(prob_reads_d_given_a, sum_prob_ancestor);
        ASSERT_EQ(summary_stat_same, summary_stat_same_ancestor);
        ASSERT_EQ(summary_stat_diff, summary_stat_diff_ancestor);


        int index16 = LookupTable::index_converter_10_to_16[anc_index];
        
        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] *  sum_prob_ancestor;
        expected_prob += prob_reads_given_a;
        expected_same += summary_stat_same*prob_reads_given_a;
        expected_diff += summary_stat_diff*prob_reads_given_a;
    }
    expected_same /= expected_prob;
    expected_diff /= expected_prob;


    double overall_stat_same;
    double overall_stat_diff;
    site.CalculateAncestorToDescendant(<#initializer#>, overall_stat_same, overall_stat_diff);
//    ASSERT_EQ(expected_prob, sum_prob_ancestor);
    ASSERT_EQ(expected_same, overall_stat_same);
    ASSERT_EQ(expected_diff, overall_stat_diff);
    ASSERT_EQ(1, overall_stat_diff+overall_stat_same);

}



TEST_F(SiteProbTest, TestCalculateAncestorToDescendant3) {


    double beta0 = 1.0 / ( 1 - 0.1*0.1 - 0.2*0.2 - 0.3*0.3 - 0.4*0.4  );
    double new_mu = (-log(0.8)) / beta0;//p = 1-exp(beta*mu) = 0.2, => (1-p) = 0.8
    params_not_equal.mutation_rate = new_mu;
    SequenceProb sp(base_custom3, params_not_equal);
    MutationProb mutation_prob = MutationProb(params_not_equal);
    F81 evo_model(new_mu, params_not_equal.nuc_freq);
    MutationMatrix matrix = evo_model.GetTranstionMatirxAToD();

    SiteProb site(sp, mutation_prob, evo_model);

    double summary_stat_diff_ancestor = 0;
    double summary_stat_same_ancestor = 0;
    double sum_prob_ancestor = 0;

    double summary_stat_diff_descendant = 0;
    double summary_stat_same_descendant = 0;
    double sum_prob_descendant = 0;


    double expected_prob = 0;
    double expected_same = 0;
    double expected_diff = 0;


    Array10D ancestor_prior = mutation_prob.GetAncestorPrior();
    DiploidProbs ancestor_genotypes = sp.GetAncestorGenotypes();
    for (int anc_index = 0; anc_index < 10; ++anc_index) {

        double summary_stat_same = 0;
        double summary_stat_diff = 0;
        double prob_reads_d_given_a = 1;
        for (int d = 0; d < sp.GetDescendantCount(); ++d) {
            HaploidProbs prob_reads_given_descent = sp.GetDescendantGenotypes(d);
            site.CalculateOneDescendantGivenAncestor(anc_index, prob_reads_given_descent, sum_prob_descendant, summary_stat_same_descendant, summary_stat_diff_descendant);
            prob_reads_d_given_a *= sum_prob_descendant;
            summary_stat_same += summary_stat_same_descendant;
            summary_stat_diff += summary_stat_diff_descendant;
        }

        site.CalculateAllDescendantGivenAncestor(anc_index, sum_prob_ancestor, summary_stat_same_ancestor, summary_stat_diff_ancestor);

        ASSERT_EQ(prob_reads_d_given_a, sum_prob_ancestor);
        ASSERT_EQ(summary_stat_same, summary_stat_same_ancestor);
        ASSERT_EQ(summary_stat_diff, summary_stat_diff_ancestor);


        int index16 = LookupTable::index_converter_10_to_16[anc_index];

        double prob_reads_given_a = ancestor_genotypes[index16] * ancestor_prior[anc_index] *  sum_prob_ancestor;
        expected_prob += prob_reads_given_a;
        expected_same += summary_stat_same*prob_reads_given_a;
        expected_diff += summary_stat_diff*prob_reads_given_a;
    }
    expected_same /= expected_prob;
    expected_diff /= expected_prob;


    double overall_stat_same;
    double overall_stat_diff;
    site.CalculateAncestorToDescendant(<#initializer#>, overall_stat_same, overall_stat_diff);
//    ASSERT_EQ(expected_prob, sum_prob_ancestor);
    ASSERT_EQ(expected_same, overall_stat_same);
    ASSERT_EQ(expected_diff, overall_stat_diff);

    ASSERT_EQ(1*sp.GetDescendantCount(), overall_stat_diff+overall_stat_same);
//
}