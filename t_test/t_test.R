#+ setup,include=FALSE
options (max.print=10000)

knitr::opts_chunk$set (prompt=TRUE, tidy=FALSE, comment="")

options (markdown.HTML.options = 
           setdiff(markdown::markdownHTMLOptions(default=TRUE),"highlight_code"))

options (warn=-1)
#+

#' T test
# set working directory
setwd("~/Desktop/evol_sim_vs_rosetta//t_test/")

#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("ggplot2")

# Load packages
library(tidyr)
library(ggplot2)
library(dplyr)

#' DATA
# read data
natural <- read.csv("graph_mean_data_natural.csv", header = TRUE, sep = "")
evolved <- read.csv("graph_mean_data_evolved.csv", header = TRUE, sep = "")
rosetta <- read.csv("graph_mean_data_rosetta.csv", header = TRUE, sep = "")


#' t-test of mean entropy
# t-test of mean entropy
t.test(natural$mean_entropy, evolved$mean_entropy, paired = T)
t.test(rosetta$mean_entropy, natural$mean_entropy, alternative = "less", paired = T)

#' t-test of correlation of entropy
cor_entropy <- read.csv("graph_entropy_corr.csv", header = TRUE, sep = "")
# paired t-test
t.test(cor_entropy$natural_evolved_corr, cor_entropy$natural_rosetta_corr, alternative = c("greater"),paired = TRUE)

#' t-test of correlation of entropy and RSA
# paired t-test
t.test(evolved$cor_entropy_RSA, natural$cor_entropy_RSA, alternative = c("greater"),paired = TRUE)
t.test(evolved$cor_entropy_RSA, rosetta$cor_entropy_RSA, alternative = c("greater"),paired = TRUE)

#' t-test of correlation of entropy and icn
# paired t-test
t.test(natural$cor_entropy_icn, evolved$cor_entropy_icn, paired = TRUE)
t.test(natural$cor_entropy_icn, rosetta$cor_entropy_icn, alternative = c("greater"),paired = TRUE)
t.test(evolved$cor_entropy_icn, rosetta$cor_entropy_icn, alternative = c("greater"),paired = TRUE)

#' t-test of correlation of entropy and iwcn
# paired t-test
t.test(natural$cor_entropy_iwcn, evolved$cor_entropy_iwcn, paired = TRUE)
t.test(natural$cor_entropy_iwcn, rosetta$cor_entropy_iwcn, alternative = c("greater"),paired = TRUE)
t.test(evolved$cor_entropy_iwcn, rosetta$cor_entropy_iwcn, alternative = c("greater"),paired = TRUE)

#' t-test of KL unorder
# t-test of KL unorder
KL_unorder <- read.csv("graph_mean_KL_all_method_data.csv", header = TRUE, sep = "")

A <- factor(rep(1:3, each=38), labels = c("natural","evolved","rosetta"))
X <- c(KL_unorder$mean_KL_split_natural, KL_unorder$mean_KL_method_evolved, KL_unorder$mean_KL_method_rosetta)
KL_unorder_data <- data.frame(X, A)

aov.kl <- aov(X~A, data = KL_unorder_data)
summary(aov.kl)

pairwise.t.test(X, A, p.adjust.method="bonferroni")

t.test(KL_unorder$mean_KL_method_evolved, KL_unorder$mean_KL_method_rosetta, alternative = c("less"),paired = TRUE)

#' t-test of KL ordered
# t-test of KL ordered

KL_ordered <- read.csv("graph_mean_KL_all_method_data_ordered.csv", header = TRUE, sep = "")

X <- c(KL_ordered$mean_KL_split_natural, KL_ordered$mean_KL_method_evolved, KL_ordered$mean_KL_method_rosetta)
KL_ordered_data <- data.frame(X, A)

aov.kl_ordered <- aov(X~A, data = KL_ordered_data)
summary(aov.kl_ordered)

pairwise.t.test(X, A, p.adjust.method="bonferroni")

t.test(KL_ordered$mean_KL_method_evolved, KL_ordered$mean_KL_method_rosetta, alternative = c("less"),paired = TRUE)


