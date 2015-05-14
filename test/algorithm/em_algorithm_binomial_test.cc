/*
 * em_algorithm_binomial_test.cc
 *
 *  Created on: 12/7/14
 *      Author: Steven Wu
 */



class EmAlgorithmBinomialTest : public ::testing::Test {
protected:
//
//    EmModelBinomial em_model_binomial;


};


TEST_F(EmAlgorithmBinomialTest, EmAlgorithmBinomialTest1) {


    std::vector<std::unique_ptr<EmData>> em_data_binomial;

    em_data_binomial.emplace_back(new EmDataBinomial(10, 5));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 4));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 7));


    EmModelBinomial em_model_binomial(10, 0.5);

    EmAlgorithmBinomial em_bin(2, em_data_binomial, em_model_binomial);

    em_bin.Run();

    std::vector<double> p1 = em_bin.GetParameters();
    std::vector<double> proportion = em_bin.GetProportion();
    for (auto item :p1) {
        printf("%.20f\t", item);
    }

    printf("\n");//0.80, 0.52
    for (auto item :proportion) {
        printf("%.20f\t", item);
    }
    ASSERT_NEAR(0.5227513173839890559, proportion[0], 1e-9);
    ASSERT_NEAR(0.79336764950715998879, p1[0], 1e-9);
    ASSERT_NEAR(0.51391659104406917091, p1[1], 1e-9);
/*
    [[1]]
    [1] 0.5227513173839890559

    [[2]]
    [1] 0.79336764950715998879

    [[3]]
    [1] 0.51391659104406917091
inits, 0.5, 0.6, 0.5, tol=1e-10
*/

    exit(35);

}



TEST_F(EmAlgorithmBinomialTest, EmAlgorithmBinomialTest2) {

    exit(35);
    std::vector<std::unique_ptr<EmData>> em_data_binomial;

    em_data_binomial.emplace_back(new EmDataBinomial(10, 1));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 2));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 1));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 2));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 1));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 2));

    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 9));
    em_data_binomial.emplace_back(new EmDataBinomial(10, 8));


    EmModelBinomial em_model_binomial(10, 0.5);

    EmAlgorithmBinomial em_bin(2, em_data_binomial, em_model_binomial);

//    em.Run();
    em_bin.Run();

    std::vector<double> p1 = em_bin.GetParameters();
    for (
        auto item :
            p1) {
        printf("%f\t", item);
    }
    printf("\n");//0.80, 0.52

    exit(35);

}


/*
N=10
## Expectation Step
estep <- function(obs,pi,p,q){
  pi_estep <- (pi*dbinom(obs,N,p)) / ( pi*dbinom(obs,N,p) + (1-pi)*dbinom(obs,N,q) )
	#print(pi_estep)

  return(pi_estep)
}


## Maximization Step
mstep <- function(obs,e.step){

  # estimate pi
  pi_temp <- mean(e.step)

  # estimate p,q
  p_temp <- sum(obs*e.step) / sum(e.step)/N
  q_temp <- sum(obs*(1-e.step)) / sum(1-e.step)/N

  list(pi_temp,p_temp,q_temp)
}


## EM Algorithm
em.algo <- function(obs,pi_inits,p_inits,q_inits,maxit=1000,tol=1e-10){
  # Initial parameter estimates
  flag <- 0
  pi_cur <- pi_inits; p_cur <- p_inits; q_cur <- q_inits

  # Iterate between expectation and maximization steps
  for(i in 1:maxit){

    cur <- c(pi_cur,p_cur,q_cur)
    new <- mstep(obs,estep(obs, pi_cur, p_cur, q_cur))
    pi_new <- new[[1]]; p_new <- new[[2]]; q_new <- new[[3]]
    new_step <- c(pi_new,p_new,q_new)
#print(new_step)
#print("==============")
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}


    # Otherwise continue iteration
    pi_cur <- pi_new; p_cur <- p_new; q_cur <- q_new
  }
  if(!flag) warning("Didn’t converge\n")

  list(pi_cur, p_cur, q_cur)
}



## Calculate Information matrix
Info.Mat.function <- function(obs, pi.est, p.est, q.est){
  estep.est <- estep(obs,pi.est, p.est, q.est)
  info.mat <- matrix(rep(0,9),3,3)
  info.mat[1,1] <- - sum(estep.est)/(pi.est^2)  - sum((1-estep.est))/((1-pi.est)^2)
  info.mat[2,2] <- - sum(estep.est*obs)/(p.est^2) - sum(estep.est*(1-obs))/((1-p.est)^2)
  info.mat[3,3] <- - sum((1-estep.est)*obs)/(q.est^2) - sum((1-estep.est)*(1-obs))/((1-q.est)^2)
  return(-info.mat)
}


## Generate sample data
n <- 5000
pi_true <- 0.50 # prob of using first coin
p_true <-  0.90 # the first coin has P(heads) = 0.70
q_true <-  0.10 # the second coin has P(heads) = 0.30
true <- c(pi_true,p_true,q_true)
u <- ifelse(runif(n)<pi_true, rbinom(N,10,p_true),rbinom(N,10,q_true))



## Set parameter estimates
pi_init = 0.50; p_init = 0.60; q_init = 0.50
pi_inits = 0.50; p_inits = 0.60; q_inits = 0.50


u <- c(5,9,8,4,7)
N <- 10
## Run EM Algorithm
output <- em.algo(u, pi_init, p_init, q_init)


## Calculate Confidence Intervals
sd.out <- sqrt(diag(solve(Info.Mat.function(u,output[[1]],output[[2]],output[[3]]))))
data.frame("Truth" = true, "EM Estimate" = unlist(output), "Lower CI" = unlist(output) - qnorm(.975)*sd.out, "Upper CI" = unlist(output) + qnorm(.975)*sd.out)

Example output:
  Truth EM.Estimate  Lower.CI  Upper.CI
1   0.9   0.8941314 0.8856034 0.9026594
2   0.6   0.6081113 0.5938014 0.6224211
3   0.5   0.4483729 0.4060064 0.4907393

 */