



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
> t.test(natural$mean_entropy, evolved$mean_entropy, paired = T)
```

```

	Paired t-test

data:  natural$mean_entropy and evolved$mean_entropy
t = 0.4058, df = 37, p-value = 0.6872
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1815397  0.2724662
sample estimates:
mean of the differences 
             0.04546327 
```

```r
> t.test(rosetta$mean_entropy, natural$mean_entropy, alternative = "less", paired = T)
```

```

	Paired t-test

data:  rosetta$mean_entropy and natural$mean_entropy
t = -19.167, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
      -Inf -1.405433
sample estimates:
mean of the differences 
              -1.541076 
```

t-test of correlation of entropy


```r
> cor_entropy <- read.csv("graph_entropy_corr.csv", header = TRUE, sep = "")
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

t-test of correlation of entropy and icn


```r
> # paired t-test
> t.test(natural$cor_entropy_icn, evolved$cor_entropy_icn, paired = TRUE)
```

```

	Paired t-test

data:  natural$cor_entropy_icn and evolved$cor_entropy_icn
t = -4.1254, df = 37, p-value = 0.0002013
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.13849664 -0.04726081
sample estimates:
mean of the differences 
            -0.09287872 
```

```r
> t.test(natural$cor_entropy_icn, rosetta$cor_entropy_icn, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  natural$cor_entropy_icn and rosetta$cor_entropy_icn
t = 14.798, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.3487537       Inf
sample estimates:
mean of the differences 
              0.3936325 
```

```r
> t.test(evolved$cor_entropy_icn, rosetta$cor_entropy_icn, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  evolved$cor_entropy_icn and rosetta$cor_entropy_icn
t = 29.447, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.4586377       Inf
sample estimates:
mean of the differences 
              0.4865112 
```

t-test of correlation of entropy and iwcn


```r
> # paired t-test
> t.test(natural$cor_entropy_iwcn, evolved$cor_entropy_iwcn, paired = TRUE)
```

```

	Paired t-test

data:  natural$cor_entropy_iwcn and evolved$cor_entropy_iwcn
t = -2.8162, df = 37, p-value = 0.007748
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.11427761 -0.01864326
sample estimates:
mean of the differences 
            -0.06646044 
```

```r
> t.test(natural$cor_entropy_iwcn, rosetta$cor_entropy_iwcn, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  natural$cor_entropy_iwcn and rosetta$cor_entropy_iwcn
t = 15.68, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.369037      Inf
sample estimates:
mean of the differences 
              0.4135302 
```

```r
> t.test(evolved$cor_entropy_iwcn, rosetta$cor_entropy_iwcn, alternative = c("greater"),paired = TRUE)
```

```

	Paired t-test

data:  evolved$cor_entropy_iwcn and rosetta$cor_entropy_iwcn
t = 29.14, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is greater than 0
95 percent confidence interval:
 0.452201      Inf
sample estimates:
mean of the differences 
              0.4799906 
```

t-test of KL unorder


```r
> # t-test of KL unorder
> KL_unorder <- read.csv("graph_mean_KL_all_method_data.csv", header = TRUE, sep = "")
> 
> A <- factor(rep(1:3, each=38), labels = c("natural","evolved","rosetta"))
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

```r
> t.test(KL_unorder$mean_KL_method_evolved, KL_unorder$mean_KL_method_rosetta, alternative = c("less"),paired = TRUE)
```

```

	Paired t-test

data:  KL_unorder$mean_KL_method_evolved and KL_unorder$mean_KL_method_rosetta
t = -23.051, df = 37, p-value < 2.2e-16
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
      -Inf -1.483286
sample estimates:
mean of the differences 
              -1.600423 
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

```r
> t.test(KL_ordered$mean_KL_method_evolved, KL_ordered$mean_KL_method_rosetta, alternative = c("less"),paired = TRUE)
```

```

	Paired t-test

data:  KL_ordered$mean_KL_method_evolved and KL_ordered$mean_KL_method_rosetta
t = -8.4473, df = 37, p-value = 1.845e-10
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
       -Inf -0.2803324
sample estimates:
mean of the differences 
             -0.3502929 
```

