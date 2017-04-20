



T test


```r
> # set working directory
> setwd("~/Desktop/evol_sim_vs_rosetta//t_test/")
> 
> #install.packages("tidyr")
> #install.packages("dplyr")
> #install.packages("ggplot2")
> 
> # Load packages
> library(tidyr)
> library(ggplot2)
> library(dplyr)
```

DATA


```r
> # read data
> natural <- read.csv("graph_mean_data_natural.csv", header = TRUE, sep = "")
> evolved <- read.csv("graph_mean_data_evolved.csv", header = TRUE, sep = "")
> rosetta <- read.csv("graph_mean_data_rosetta.csv", header = TRUE, sep = "")
```

t-test of mean entropy


```r
> # t-test of mean entropy
> X <- c(natural$mean_entropy, evolved$mean_entropy, rosetta$mean_entropy)
> A <- factor(rep(1:3, each=38), labels = c("natural","evolved","rosetta"))
> entropy_data <- data.frame(X, A)
> 
> aov.entropy <- aov(X~A, data = entropy_data)
> summary(aov.entropy)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2  58.44  29.221   98.63 <2e-16 ***
Residuals   111  32.89   0.296                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> plot(entropy_data$X~entropy_data$A)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
> # compare
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved 1       -      
rosetta <2e-16  <2e-16 

P value adjustment method: bonferroni 
```

```r
> pairwise.t.test(X, A, p.adjust.method="none")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved 0.72    -      
rosetta <2e-16  <2e-16 

P value adjustment method: none 
```

```r
> TukeyHSD(aov(X~A, entropy_data))
```

```
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = X ~ A, data = entropy_data)

$A
                       diff        lwr        upr     p adj
evolved-natural -0.04546327 -0.3421059  0.2511793 0.9295905
rosetta-natural -1.54107619 -1.8377188 -1.2444336 0.0000000
rosetta-evolved -1.49561292 -1.7922555 -1.1989703 0.0000000
```

t-test of correlation of entropy


```r
> # t-test of correlation of entropy
> 
> cor_entropy <- read.csv("graph_entropy_corr.csv", header = TRUE, sep = "")
> 
> # paired t-test
> t.test(cor_entropy$natural_evolved_corr, cor_entropy$natural_rosetta_corr, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  cor_entropy$natural_evolved_corr and cor_entropy$natural_rosetta_corr
t = 14.906, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.2576425       Inf
sample estimates:
mean of the differences 
              0.2905249 
```

t-test of correlation of entropy and RSA


```r
> # paired t-test
> t.test(evolved$cor_entropy_RSA, natural$cor_entropy_RSA, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  evolved$cor_entropy_RSA and natural$cor_entropy_RSA
t = 9.6539, df = 37, p-value = 5.939e-12
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.1790246       Inf
sample estimates:
mean of the differences 
              0.2169357 
```

```r
> t.test(evolved$cor_entropy_RSA, rosetta$cor_entropy_RSA, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  evolved$cor_entropy_RSA and rosetta$cor_entropy_RSA
t = 26.556, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.5036904       Inf
sample estimates:
mean of the differences 
              0.5378606 
```

```r
> #paired t-test with bonferroni correction
> X <- c(natural$cor_entropy_RSA, evolved$cor_entropy_RSA, rosetta$cor_entropy_RSA)
> cor_entropy_RSA_data <- data.frame(X, A)
> 
> aov.rsa <- aov(X~A, data = cor_entropy_RSA_data)
> summary(aov.rsa)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2  5.565  2.7825   245.1 <2e-16 ***
Residuals   111  1.260  0.0114                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved 4.1e-14 -      
rosetta < 2e-16 < 2e-16

P value adjustment method: bonferroni 
```

t-test of correlation of entropy and icn


```r
> # t-test of correlation of entropy and icn
> 
> X <- c(natural$cor_entropy_icn, evolved$cor_entropy_icn, rosetta$cor_entropy_icn)
> cor_entropy_icn_data <- data.frame(X, A)
> 
> aov.icn <- aov(X~A, data = cor_entropy_icn_data)
> summary(aov.icn)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2  5.070  2.5350   237.8 <2e-16 ***
Residuals   111  1.183  0.0107                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved 0.00046 -      
rosetta < 2e-16 < 2e-16

P value adjustment method: bonferroni 
```

t-test of correlation of entropy and iwcn


```r
> # t-test of correlation of entropy and iwcn
> t.test(evolved$cor_entropy_iwcn, natural$cor_entropy_iwcn, paired = TRUE)
```

```

	Paired t-test

data:  evolved$cor_entropy_iwcn and natural$cor_entropy_iwcn
t = 2.8162, df = 37, p-value = 0.007748
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.01864326 0.11427761
sample estimates:
mean of the differences 
             0.06646044 
```

```r
> X <- c(natural$cor_entropy_iwcn, evolved$cor_entropy_iwcn, rosetta$cor_entropy_iwcn)
> cor_entropy_iwcn_data <- data.frame(X, A)
> 
> aov.iwcn <- aov(X~A, data = cor_entropy_iwcn_data)
> summary(aov.iwcn)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2  5.140  2.5702   225.4 <2e-16 ***
Residuals   111  1.266  0.0114                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved 0.023   -      
rosetta <2e-16  <2e-16 

P value adjustment method: bonferroni 
```

t-test of KL unorder


```r
> # t-test of KL unorder
> KL_unorder <- read.csv("graph_mean_KL_all_method_data.csv", header = TRUE, sep = "")
> 
> X <- c(KL_unorder$mean_KL_split_natural, KL_unorder$mean_KL_method_evolved, KL_unorder$mean_KL_method_rosetta)
> KL_unorder_data <- data.frame(X, A)
> 
> aov.kl <- aov(X~A, data = KL_unorder_data)
> summary(aov.kl)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2  224.3  112.14   864.4 <2e-16 ***
Residuals   111   14.4    0.13                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved <2e-16  -      
rosetta <2e-16  <2e-16 

P value adjustment method: bonferroni 
```

t-test of KL ordered


```r
> # t-test of KL ordered
> 
> KL_ordered <- read.csv("graph_mean_KL_all_method_data_ordered.csv", header = TRUE, sep = "")
> 
> X <- c(KL_ordered$mean_KL_split_natural, KL_ordered$mean_KL_method_evolved, KL_ordered$mean_KL_method_rosetta)
> KL_ordered_data <- data.frame(X, A)
> 
> aov.kl_ordered <- aov(X~A, data = KL_ordered_data)
> summary(aov.kl_ordered)
```

```
             Df Sum Sq Mean Sq F value Pr(>F)    
A             2 17.960   8.980   373.4 <2e-16 ***
Residuals   111  2.669   0.024                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
> pairwise.t.test(X, A, p.adjust.method="bonferroni")
```

```

	Pairwise comparisons using t tests with pooled SD 

data:  X and A 

        natural evolved
evolved < 2e-16 -      
rosetta < 2e-16 2.4e-16

P value adjustment method: bonferroni 
```

