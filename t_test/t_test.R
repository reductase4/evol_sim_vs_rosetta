#+ setup,include=FALSE
options (max.print=10000)

knitr::opts_chunk$set (prompt=TRUE, tidy=FALSE, comment="")

options (markdown.HTML.options = 
           setdiff(markdown::markdownHTMLOptions(default=TRUE),"highlight_code"))

options (warn=-1)
#+

#' T test
# set working directory
setwd("~/Desktop/evol_sim_vs_rosetta/t_test/")

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
t.test(natural$cor_entropy_RSA, rosetta$cor_entropy_RSA, alternative = c("greater"),paired = TRUE)


#' t-test of KL unorder
# t-test of KL unorder
KL_unorder <- read.csv("graph_mean_KL_all_method_data.csv", header = TRUE, sep = "")

t.test(KL_unorder$mean_KL_method_evolved, KL_unorder$mean_KL_split_natural, paired = TRUE)
t.test(KL_unorder$mean_KL_method_rosetta, KL_unorder$mean_KL_split_natural, paired = TRUE)
t.test(KL_unorder$mean_KL_method_rosetta, KL_unorder$mean_KL_method_evolved, paired = TRUE)

#' t-test of mean seq divergence
seq_divergence <- read.csv("mean_seq_divergence.csv", header = TRUE, sep = "")
# paired t-test
t.test(seq_divergence$rosetta, seq_divergence$natural, alternative = c("greater"), paired = TRUE)
t.test(seq_divergence$evolved, seq_divergence$natural, paired = TRUE)
t.test(seq_divergence$rosetta, seq_divergence$evolved, alternative = c("greater"), paired = TRUE)

#' t-test of mean seq divergence for ten proteins
seq_divergence_ten <- read.csv("mean_seq_divergence_ten.csv", header = TRUE, sep = "")
# paired t-test
t.test(seq_divergence_ten$rosetta, seq_divergence_ten$natural, alternative = c("greater"), paired = TRUE)
t.test(seq_divergence_ten$evolved_from_designed, seq_divergence_ten$natural, alternative = c("greater"), paired = TRUE)
#t.test(seq_divergence_ten$evolved, seq_divergence_ten$natural, paired = TRUE)
#t.test(seq_divergence_ten$rosetta, seq_divergence_ten$evolved, paired = TRUE)
t.test(seq_divergence_ten$rosetta, seq_divergence_ten$evolved_from_designed, alternative = c("greater"), paired = TRUE)
#t.test(seq_divergence_ten$evolved_from_designed, seq_divergence_ten$evolved, paired = TRUE)

#' t-test for sequences evolved from design
KL <- read.csv("evolved_from_designed/graph_mean_KL_all_method_data.csv", header = TRUE, sep = "")
rosetta_2 <- read.csv("evolved_from_designed/graph_mean_data_rosetta.csv", header = TRUE, sep = "")
evolved_2 <- read.csv("evolved_from_designed/graph_mean_data_evolved.csv", header = TRUE, sep = "")
natural_2 <- read.csv("evolved_from_designed/graph_mean_data_natural.csv", header = TRUE, sep = "")
cor_entropy_2 <- read.csv("evolved_from_designed/graph_entropy_corr.csv", header = TRUE, sep = "")

# KL
t.test(KL$mean_KL_method_rosetta, KL$mean_KL_method_evolved, alternative = c("greater"),paired = TRUE)

# effective number
t.test(evolved_2$mean_entropy, rosetta_2$mean_entropy, alternative = c("greater"), paired = TRUE)

# cor_effective_number
t.test(cor_entropy_2$natural_evolved_corr, cor_entropy_2$natural_rosetta_corr, alternative = c("greater"), paired = TRUE)

# cor_RSA_effective_number
t.test(evolved_2$cor_entropy_RSA, rosetta_2$cor_entropy_RSA, alternative = c("greater"), paired = TRUE)
t.test(evolved_2$cor_entropy_RSA, natural_2$cor_entropy_RSA, paired = TRUE)




