# README 

GLMM since we use count data (discrete values + 0-infinite) and therefore not continuous values + not in the range ]-Inf;+Inf[

g(number of larvae) = genotype + 1|genotype:batch

# Statistical info 

Model to use for count data = Generalized Linear Model with Poisson as the link function.

## Dealing with batch effects (sessions: dates)

Since we have batch effects due to series of data acquired during different months, we code the month as a series of batch = batch_01, batch_02, etc.

```r
library(lme4)

my_model <- lmer(number_larvae ~ genotype + (1|genotype:batch))
```

## Equivalent of post-hoc tests with a GLMM

[Stat Exchange link](https://stats.stackexchange.com/questions/499089/emmeans-output-interpretation-of-a-glmer-fit-with-nesting)  

emmeans to estimate the 
