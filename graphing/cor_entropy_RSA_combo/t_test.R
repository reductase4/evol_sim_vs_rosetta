

# T test
# set working directory
setwd("~/Desktop/calc_score/matched_stability/graphing/cor_entropy_RSA/")

# Load packages
library(tidyr)
library(ggplot2)
library(dplyr)

# read data
evolved_matched <- read.csv("graph_mean_data_evolved_matched.csv", header = TRUE, sep = "")
rosetta_matched <- read.csv("graph_mean_data_rosetta_matched.csv", header = TRUE, sep = "")
evolved_max <- read.csv("graph_mean_data_evolved_max.csv", header = TRUE, sep = "")
evolved_mean <- read.csv("graph_mean_data_evolved_mean.csv", header = TRUE, sep = "")
rosetta <- read.csv("graph_mean_data_rosetta.csv", header = TRUE, sep = "")


# t-test 
t.test(evolved_matched$cor_entropy_RSA, rosetta_matched$cor_entropy_RSA, 
       alternative = "greater", paired = T)
t.test(evolved_mean$cor_entropy_RSA, rosetta$cor_entropy_RSA, 
       alternative = "greater", paired = T)
t.test(evolved_max$cor_entropy_RSA, rosetta$cor_entropy_RSA, 
       alternative = "greater", paired = T)
t.test(evolved_mean$cor_entropy_RSA, evolved_max$cor_entropy_RSA, 
       paired = T)


