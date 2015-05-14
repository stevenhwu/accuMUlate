#include <evolution_models/F81.h>
#include <mutations/sequence_prob_v1.h>
#include "gtest/gtest.h"
#include "unittest_utils.h"
//#include "site_genotypes.h"
//#include "sequence_prob_v1.h"

class SequenceProbTest : public ::testing::Test {
public:

protected:
    int foo;
    ModelParams params_equal;
    ModelParams params_not_equal;

    std::vector<double> freq_equal {0.25,0.25,0.25,0.25};     //beta0= 1.333333
    std::vector<double> freq_not_equal {0.1, 0.2, 0.3, 0.4};

    double mu = 0.1;
    double mu_4 = 0.0001;

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
    }
};


TEST_F(SequenceProbTest, EqualFreqsGenotypes) {


    int read_count = 2;
    ModelInput base_custom;
    base_custom.reference = 0;
    for (int i = 0; i < read_count; ++i) {
        ReadData r;

        for (int j = 0; j < 4; ++j) {
            r.reads[j] = (uint16_t) (j + 1);
        }
        base_custom.all_reads.push_back(r);
    }

    SequenceProb sp(base_custom, params_equal);


    DiploidProbs anc_probs = sp.GetAncestorGenotypes();
    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);

    double expected[] = {1.03172733389759e-05, 0.000768055062625644, 0.0329670329670333, 1};
    for (int i = 0; i < 4; ++i) { //ddm( c(1,2,3,4), )
        ASSERT_NEAR(expected[i], des_probs[i], ERROR_THRESHOLD);
    }

    double expected16[] = {1.24423781767733e-09, 3.50842296285468e-06, 7.70497866992586e-05, 0.00120711332495508, 3.50842296285468e-06, 9.26255536302597e-08, 0.00290645699150511, 0.0455344928669145, 7.70497866992586e-05, 0.00290645699150511, 3.97574318393215e-06, 1, 0.00120711332495508, 0.0455344928669145, 1, 0.000120597543245941};
    double expectedRef[] = {0.24721765553421, 0.000308636274075169, 0.000308636274075169, 0.000308636274075169, 0.000308636274075169, 0.000308636274075169, 7.69666518890692e-07, 7.69666518890692e-07, 0.000308636274075169, 7.69666518890692e-07, 0.000308636274075169, 7.69666518890692e-07, 0.000308636274075169, 7.69666518890692e-07, 7.69666518890692e-07, 0.000308636274075169};
    double expectedBoth[] = {3.07597556213192e-10, 1.08282659113523e-09, 2.37803590851457e-08, 3.72558959000623e-07, 1.08282659113523e-09, 2.85876057565931e-11, 2.23700263495725e-09, 3.50463746143312e-08, 2.37803590851457e-08, 2.23700263495725e-09, 1.22705856296857e-09, 7.69666518890692e-07, 3.72558959000623e-07, 3.50463746143312e-08, 7.69666518890692e-07, 3.72207764100462e-08};
    for (int i = 0; i < 16; ++i) { //ddm( c(1,2,3,4), )
        ASSERT_NEAR(expectedBoth[i], anc_probs[i], ERROR_THRESHOLD);
    }

/*
phi = 0.01
err = 0.01
alpha = (1-phi)/phi
AT = (1-phi)/phi
a_same = (1-err)*AT
a_diff = err/3*AT
data = c(1,2,3,4)
result = vector(length=4)
result2 = vector(length=4)

for(i in 1:4){
    aa = rep(a_diff, 4)
    aa[i] = a_same
    c1 = gamma(sum(aa)) / gamma(sum(data)+sum(aa))
    c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))

    x1 = prod( gamma(data+aa) / gamma(aa) )
    x2 = sum( lgamma(data+aa) - lgamma(aa) )
    result[i] = c1*x1
    result2[i] = c2+x2

}
r1 = result / max(result)
paste(r1,sep=",", collapse=", ")
#[1] "1.03172733389759e-05, 0.000768055062625644, 0.0329670329670333, 1

r2 = result2 - max(result2)
r2 = exp(r2)
paste(r2,sep=",", collapse=", ")
#[1] "1.03172733389759e-05, 0.000768055062625643, 0.0329670329670333, 1"


////////////////////////////////////////////////

phi = 0.01
err = 0.01
theta = 0.01
alpha = (1-phi)/phi
AT = (1-phi)/phi
a_same = (1-err)*AT
a_diff = err/3*AT
a_ij_same = (0.5-err/3)*AT

result16 = vector(length=16)
result_ref = vector(length=16)

## anc
    data<- c(1,2,3,4)
	for(i in 1:4){for(j in 1:4){
	    index16 = (i-1)*4+j
	    if(i==j){
            aa = rep(a_diff, 4)
            aa[i] = a_same
            c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
            x2 = sum( lgamma(data+aa) - lgamma(aa) )
            result16[index16] = c2+x2

        }
        else{
            aa = rep(a_diff, 4)
            aa[i] = a_ij_same
            aa[j] = a_ij_same

            c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
            x2 = sum( lgamma(data+aa) - lgamma(aa) )
            result16[index16] = c2+x2

        }

    }}
    r16 = result16 - max(result16)
    r16 = exp(r16)
    paste(r16,sep=",", collapse=", ")
    r16m<- matrix(r16,4,4); r16m

## ref
    aa = rep(theta * 0.25, 4)
    for(i in 1:4){for(j in 1:4){
        data <- rep(0,4)
	    index16 = (i-1)*4+j
	    data[i] = data[i]+1
	    data[j] = data[j]+1
	    data[1] = data[1]+1

        c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
        x2 = sum( lgamma(data+aa) - lgamma(aa) )
        result_ref[index16] = c2+x2

    }}
    r_ref = exp(result_ref)
    paste(r_ref,sep=",", collapse=", ")
    paste(r16*r_ref,sep=",", collapse=", ")

*/

}


TEST_F(SequenceProbTest, NotEqualFreqsGenotypes) {


    int read_count = 2;
    ModelInput base_custom;
    base_custom.reference = 0;
    for (int i = 0; i < read_count; ++i) {
        ReadData r;

        for (int j = 0; j < 4; ++j) {
            r.reads[j] = (uint16_t) (j * 2)+1 ;
        }
        base_custom.all_reads.push_back(r);
    }

    SequenceProb sp(base_custom, params_not_equal);


    DiploidProbs anc_probs = sp.GetAncestorGenotypes();
    HaploidProbs des_probs = sp.GetDescendantGenotypes(0);

    double expected[] = {1.37907066641515e-09, 4.40657784570141e-06, 0.00314902737102689, 1};
    for (int i = 0; i < 4; ++i) { //ddm( c(1,2,3,4), )
        ASSERT_NEAR(expected[i], des_probs[i], ERROR_THRESHOLD);
    }

    double expected16[] = {5.94618429635949e-15, 3.68219182217315e-10, 7.08372636004037e-08, 6.27469581357529e-06, 3.68219182217315e-10, 1.89999864582048e-11, 5.86831924857096e-05, 0.00519810003241312, 7.08372636004037e-08, 5.86831924857096e-05, 1.35777647646444e-08, 1, 6.27469581357529e-06, 0.00519810003241312, 1, 4.31173285109194e-06};
    double expectedRef[] = {0.0986651396482933, 9.86158317324269e-05, 0.000147923747598641, 0.000197231663464854, 9.86158317324269e-05, 9.8714349046845e-05, 2.95551943254027e-07, 3.94069257672036e-07, 0.000147923747598641, 2.95551943254027e-07, 0.000148219299541895, 5.91103886508055e-07, 0.000197231663464854, 3.94069257672036e-07, 5.91103886508055e-07, 0.000197822767351362};
    double expectedBoth[] = {5.86681103974797e-16, 3.63122409141946e-14, 1.04785135014045e-11, 1.23756869304741e-09, 3.63122409141946e-14, 1.87557129512056e-15, 1.73439315755016e-11, 2.04841142107802e-09, 1.04785135014045e-11, 1.73439315755016e-11, 2.01248678276021e-12, 5.91103886508055e-07, 1.23756869304741e-09, 2.04841142107802e-09, 5.91103886508055e-07, 8.52958924682786e-10};
    for (int i = 0; i < 16; ++i) { //ddm( c(1,2,3,4), )
        ASSERT_NEAR(expectedBoth[i], anc_probs[i], ERROR_THRESHOLD);
    }


/*
phi = 0.01
err = 0.01
alpha = (1-phi)/phi
AT = (1-phi)/phi
a_same = (1-err)*AT
a_diff = err/3*AT

data = c(1,3,5,7)
result = vector(length=4)
result2 = vector(length=4)

for(i in 1:4){
    aa = rep(a_diff, 4)
    aa[i] = a_same
    c1 = gamma(sum(aa)) / gamma(sum(data)+sum(aa))
    c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))

    x1 = prod( gamma(data+aa) / gamma(aa) )
    x2 = sum( lgamma(data+aa) - lgamma(aa) )
    result[i] = c1*x1
    result2[i] = c2+x2

}
r1 = result / max(result)
paste(r1,sep=",", collapse=", ")
#[1] "1.03172733389759e-05, 0.000768055062625644, 0.0329670329670333, 1

r2 = result2 - max(result2)
r2 = exp(r2)
paste(r2,sep=",", collapse=", ")
#[1] "1.03172733389759e-05, 0.000768055062625643, 0.0329670329670333, 1"


////////////////////////////////////////////////

phi = 0.01
err = 0.01
theta = 0.01
alpha = (1-phi)/phi
AT = (1-phi)/phi
a_same = (1-err)*AT
a_diff = err/3*AT
a_ij_same = (0.5-err/3)*AT

result16 = vector(length=16)
result_ref = vector(length=16)
data<- c(1,3,5,7)
## anc

	for(i in 1:4){for(j in 1:4){
	    index16 = (i-1)*4+j
	    if(i==j){
            aa = rep(a_diff, 4)
            aa[i] = a_same
            c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
            x2 = sum( lgamma(data+aa) - lgamma(aa) )
            result16[index16] = c2+x2

        }
        else{
            aa = rep(a_diff, 4)
            aa[i] = a_ij_same
            aa[j] = a_ij_same

            c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
            x2 = sum( lgamma(data+aa) - lgamma(aa) )
            result16[index16] = c2+x2

        }

    }}
    r16 = result16 - max(result16)
    r16 = exp(r16)
    paste(r16,sep=",", collapse=", ")
    r16m<- matrix(r16,4,4); r16m

## ref
    aa = theta * c(0.1, 0.2, 0.3, 0.4)
    for(i in 1:4){for(j in 1:4){
        data <- rep(0,4)
	    index16 = (i-1)*4+j
	    data[i] = data[i]+1
	    data[j] = data[j]+1
	    data[1] = data[1]+1

        c2 = lgamma(sum(aa)) - lgamma(sum(data)+sum(aa))
        x2 = sum( lgamma(data+aa) - lgamma(aa) )
        result_ref[index16] = c2+x2

    }}
    r_ref = exp(result_ref)
    paste(r_ref,sep=",", collapse=", ")
    paste(r16*r_ref,sep=",", collapse=", ")

*/

}

