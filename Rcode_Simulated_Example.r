
## Rcode_Simulated_Example_Dec_29.r

## This file provides code and comments for a step by step illustrative example to accompany
## Online Appendix D of "The First-Order Markov Conditional Linear Expectation Approach for Analysis of Longitudinal Data"
## by Bender, Gamerman, Reese, Gray, Li** and Shults ** (** indicates co-senior authors, equal contribution).
	
## This example recreates Online Appendix D and expands on the analysis that is provided in Online Appendix D.	

## First, load the libraries that will be used for this illustrative example.
  library(xtable)
  library(MASS)
  library(geeM)


## Next, load the code for the analysis:
## (You will need to modify the path to indicate the location of 20160610_Project2_Code.R on your computer.)
  source("h:\\bender\\20160610_Project2_Code.R")

## Read the simulated data:
## (You will need to modify the path to indicate the location of OneSimulatedData.csv on your computer.)
data2 = read.csv("h:\\bender\\OneSimulatedData.csv")

## Next, summarize the variables in the simulated data set:
summary(data2)
## subject	= numeric subject identification variable
## y = outcome variable = number of seizures
## trt = variable that indicates treatment, takes value placebo or progabide (treatment)	
## base = number of seizures during 8 week period prior to baseline
## age = participant age at baseline	
## period = variable that indicates timing of measurement

## period = 1, 2, 4, 8 for each participant, so that measurements are unequally spaced in time.

## The data set OneSimulatedData.csv was simulated under the assumption of a Poisson distribution for
## the following true model for the expected seizure count on subject i at measurement occasion j

## E(count) = exp(beta_0 + beta_1 trt + beta_2 base + beta_3 age )
## where trt  = placebo or progabide (treatment)
##       base = number of seizures in the period just prior to baseline
##       age  = age at baseline.

## The assumed values for (beta_0, beta_1, beta_2, beta_3) = (0.66, -0.25, 0.02, 0.02, -0.03).

## In addition, each participant had measurements at period = (1,2,4,8).
   
## The assumed correlation structure for the simulated data set was Markov, with correlation
## parameter alpha = 0.50.

## The number of subjects in the simulated data set = 177.

## Next,let's run the analysis:

## We will first fit a model that correctly specifies the Poisson model for the distribution of the
## outcome variable and that correctly specifies the model for the expected seizure counts on each
## subject. However, we will specify an AR(1) correlation strucure, which is not equal to the
## Markov correlation structure for this example, due to the unequal spacing of measurements on
## each subject.

## Analysis of simulated epilepsy data set with Poisson and AR(1) assumptions
   PoisAR = EndResults(y ~ trt + base + age + period, "AR(1)", "Poisson", data2, data2$subject, data2$period, rep(0,6))
   PoisAR
## The AIC and BIC for this model are 3304.503569 and 3323.560467, respectively.

## Next, let's fit the same model for the expected seizure counts, but change the assumed
## working correlation structure to the structure that was assumed for the simulated data set.

## Analysis of epilepsy data with Poisson and Markov assumptions
  PoisMark = EndResults(y ~ trt + base + age + period, "Markov", "Poisson", data2, data2$subject, data2$period, rep(0,6))
  PoisMark
## The AIC and BIC for this model are 3290.940209 and 3309.997107, respectively.
## Both the AIC and BIC are lower for the model with the Markov structure, which suggests that the fit is
## better for the Markov structure. (This is as anticipated because the data was simulated for a model with true Markov structure.)

## Note also that we did not use the likelihood ratio test to compare the AR(1) and Markov structures for this example,
## since the AR(1) structure is not a special case of the Markov structure for this example.

## The analyses thus far were also provided in Online Appendix D. We now provide some additional analysis,
## to further demonstrate the application of the proposed method.

## Note that for the Markov structure, the estimated value of the correlation parameter was alpha = 0.5414033986.
## The timings of measurements on each subject are (1,2,4,8).
## The measurements on each subject are therefore unequally spaced in time, but the timings are
## are the same for each subject.

## The estimated correlation between consecutive measurements on each subject for the Markov structure is

## (0.5414033986)^(2-1) = 0.5414033986 between timing 1 and 2

(0.5414033986)^(4-2)
## (0.5414033986)^(4-2) = 0.29311764 between timing 2 and 4

(0.5414033986)^(8-4)
## (0.5414033986)^(8-4) = 0.08591795089 between timing 4 and 8

## Let's next fit the correct model for the expected seizure counts, but fit an AD(1) correlation structure.

## For the AD(1) structure, the assumed correlation between consecutive measurements on each subject is

## alpha1 between timing 1 and 2
## alpha2 between timing 2 and 4
## alpha3 between timing 4 and 8

## Because the timings are the same for each subject, the Markov structure is a special case of the AD(1)
## structure for this example.

  PoisAD1 = EndResults(y ~ trt + base + age + period, "AD(1)", "Poisson", data2, data2$subject, data2$period, rep(0,8))
  PoisAD1
## The AIC and BIC for this model are 3288.709272 and 3307.766170, which are both slightly lower than for the
## model with the Markov structure.

## The estimated correlation between consecutive measurements on each subject for the AD(1) structure is:

##  0.5037994848 between timing 1 and 2
##  0.2963333571 between timing 2 and 4 
##  0.1751822128 between timing 4 and 8

## Note that the estimated correlations are similar to the estimated correlations for the Markov structure
## that were provided above.

## Since the Markov structure is a special case of the AD(1) structure for this example, we can
## also apply the likelihood ratio test to compare the models for the Markov and AD(1) structures.

## Likelihood ratio test comparing Markov to AD(1) under Poisson assumption

## The degrees of freedom for the test statistic are 2 because the model with the
## AD(1) structure has 2 additional parameters than the model with the Markov structure
dchisq(2 * (PoisAD1[[1]][1] - PoisMark[[1]][1]), 2)
## Note that the p-value for this test is 0.1638808106.
## Therefore, although the AIC and BIC are slightly smaller for the model with the AD(1) structure,
## the likelihood ratio test (for a 0.05 level of significance) does not cause us to reject the
## Markov structure in favor of the AD(1) structure. 
## Since the AIC and BIC are very close and the p-value for the likelihood ratio test exceeds 0.05,
## we will select the model with the Markov structure as our final model.

## Let's next use the likelihood ratio test to compare the fit of the final model with the
## Markov structure and an assumed Poisson distribution, with the same model for an
## assumed negative Binomial distribution.

  NegBMark = EndResults(y ~ trt + base + age + period, "Markov", "Negative-Binomial", data2, data2$subject, data2$period, c(rep(0,6),1))
  NegBMark
## The AIC and BIC are 3290.536266 and 3309.593164 for this model.
## The values were 3290.940209 and 3309.997107 for the Markov structure and Poisson distribution.
## Therefore, the AIC and BIC are only very slightly smaller for the model with the Negative-Binomial distribution.

## Note also that the estimated value of the r parameter for the negative-Binomial distribution is very close to zero:
## The estimated value of r = 0.002191353167.

## Let's next apply the likelihood ratio test to compare the models for the Poisson and negative-Binomial distributions:
## (The degree of freedom for this test is 1 because the model for the negative-Binomial distribution only involves one 
## additional parameter (r).)
dchisq(2 * (NegBMark[[1]][1] - PoisMark[[1]][1]), 1)
## The p-value for this test is 0.5129022834.
## Note that since we are testing if r is equal to zero, which is on the boundary of its parameter space (since r must be positive)
## the p-value for this test should be divided by 2.
0.5129022834/2
## The p-value for the likelihood ratio test is 0.2564511417.
## Since the likelihood ratio test fails to reject the model with the Poisson distribution in favor of the model
## with the negative Binomial dstribution, and since the AIC and BIC values are so close in value, we base our
## final model on the Poisson distribution, with a Markov correlation structure.

## Suppose that we wanted to fit a different model to a different data set.

## The syntax for a model with 3 covariates, cov1, cov2, and cov3:

## model =  EndResults(y ~ cov1 + cov2 + cov3, "correlation", "distribution", data, data$subject, data$time, rep(0,Number))

## where
## correlation is AR(1) or Markov or AD(1)
## distribution is Binomial or Poisson or Negative-Binomial
## data is the name of the data set
## data$subject indicates that subject is the name of the subject identification variable in the data set
## data$time indicates that time is the timing variable in the data set
## Number is a numeric value that equals the number of regression parameters plus the number of correlation parameters in the model.
## (For example, Number = 5 for a model with three covariates and an AR(1) structure, 
##  because the regression parameter = (beta0,beta1,beta2,beta3) involves a regression parameter for the constant
##  and each covariate. In addition, the AR(1) structure involves one parameter alpha. In general,
##  Number = the number of covariates in the model + 1 + the number of correlation parameters in the 
##  assumed correlation structure.)
