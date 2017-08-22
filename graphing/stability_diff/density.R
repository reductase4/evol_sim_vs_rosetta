setwd("~/Desktop/calc_score/stability_diff/")

require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
require(ggplot2)

# get pdb names
PDBS <- c("1ci0A", "1g58B", "1hujA", "1ibsA", "1jlwA", "1kzlA", "1m3uA", "1mozA", "1pv1A", "1qmvA", "1riiA", "1v9sB", "1w7wB", "1x1oB", "1ypiA", "1znnA", "2a84A", "2bcgY", "2br9A", "2cjmC", "2esfA", "2fliA", "2gv5D", "1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A", "1okcA", "1r6mA", "1xtdA", "1ysbA", "1zwkA", "2aiuA", "2cfeA", "2cnvA", "2eu8A", "2g0nB")
pdb <- NULL
for (idx in 1:length(PDBS)) {
  protein <- PDBS[idx]
  pdb_id <- toupper(substr(protein, 1, 4))
  chain_id <- toupper(substr(protein, 5, 5))
  pdb[idx] <- paste(pdb_id, chain_id, sep = "_")
}

pdb_names <- data_frame(pdb = pdb)

data_path <- "~/Desktop/calc_score/scores"   # path to the data

# get evolved data
data_evolved <- pdb_names %>%
  # create filenames
  mutate(filename = paste(pdb, "_evolved_scores.csv", sep="")) %>%
  # read in data
  mutate(file_contents = map(filename,
                             ~ read_csv(file.path(data_path, .)))
  ) %>% 
  select(-filename) %>% # remove filenames, not needed anynmore
  unnest() # unnest
  
data_evolved <- data_evolved[-2]
data_evolved$method <- "evolved"
colnames(data_evolved) <- c("pdb", "score", "method")

# get designed data
data_designed <- pdb_names %>%
  # create filenames
  mutate(filename = paste(pdb, "_rosetta_scores", sep="")) %>%
  # read in data
  mutate(file_contents = map(filename,
                             ~ read_csv(file.path(data_path, .), col_names=F, n_max=500))
  ) %>% 
  select(-filename) %>% # remove filenames, not needed anynmore
  unnest() # unnest

data_designed$method <- "designed"
colnames(data_designed) <- c("pdb", "score", "method")

all_data <- rbind(data_evolved, data_designed)

p1 <- all_data %>%
  filter(pdb=="1B4T_A"|pdb=="1CI0_A") %>%
  ggplot(aes(x=score, fill=method)) + geom_density() + facet_wrap(~pdb, scales = "free") 
p1

p2 <- all_data %>%
  filter(pdb=="1EFV_B") %>%
  ggplot(aes(x=score, fill=method)) + geom_density() + facet_wrap(~pdb) 
p2

p <- ggplot(all_data, aes(x=score, fill=method)) + geom_density() + facet_wrap(~pdb, scales = "free")
p

equal_breaks <- function(n, s, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(round(min(x)+d, digits=1), round(max(x)-d, digits=1), length=n)
  }
}


p <- ggplot(all_data, aes(x=score, fill=method)) + geom_density(alpha=0.5) + 
  scale_fill_manual(values = c(designed="#E69F00", evolved="#56B4E9")) +
  facet_wrap(~pdb, scales = "free_x", ncol = 6) +
  xlab("Scores") + ylab("Density") +
  scale_x_continuous(breaks=equal_breaks(n=3, s=0.15), 
                     expand = c(0.15, 0)) +
  theme(text = element_text(size=10), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(1, "lines")) 
p

ggsave("scores_original.eps", plot = p, width = 15, height = 10, units = "cm", dpi=500)
