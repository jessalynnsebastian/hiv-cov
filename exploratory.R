source("data_setup_BEAST.R")
library(table1)
library(TransPhylo)
# In this file: table1 for metadata, some exploration of transphylo outputs

#### formatting metadata for table ####
tbl_dat <- meta[,c("new_hiv_stat", "cal_age", "gender", "cvid_syp", "vacc_sts", "travel", "risk_lvl")]
tbl_dat$new_hiv_stat <- as.factor(tbl_dat$new_hiv_stat)
tbl_dat$gender <- as.factor(tbl_dat$gender)
tbl_dat$cvid_syp <- as.factor(tbl_dat$cvid_syp)
tbl_dat$vacc_sts <- as.factor(tbl_dat$vacc_sts)
tbl_dat$travel <- as.factor(tbl_dat$travel)
tbl_dat$risk_lvl <- factor(tbl_dat$risk_lvl, labels = c("High", "Low", "Medium", "None"))
###############

# table1 comparing hiv group to no hiv group
table1( ~ cal_age + gender + cvid_syp + vacc_sts + travel + risk_lvl | new_hiv_stat, data = tbl_dat)

# distribution of snp distances
snp_dists <- dist.dna(seq_2023.08, model = "N")
hist(snp_dists)


## transphylo results
results_BA <- readRDS("./results/tp_BA.rds")
results_BE <- readRDS("./results/tp_BE.rds")
results_BQ <- readRDS("./results/tp_BQ.rds")
results_FN <- readRDS("./results/tp_FN.rds")
results_XBB <- readRDS("./results/tp_XBB.rds")

# traceplots (hard to look at all of these since there're a lot)
# convert each TP to coda mcmc object
# then use coda::traceplot


# analyses by phylogenetic tree group: how much uncertainty comes from the phylogeny vs. from the transmission tree?

### BA trees
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_BA <- lapply(results_BA, function(x) getOffspringDist(x, k = results_BA[[1]][[1]]$ctree$nam ) ) # using default burnin of 0.5, get offspring dist of all BA seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_BA <- lapply(offspring_BA, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# for each of the 100 runs, find how many of the 17 BA individuals were labeled infectors (say prob >= 0.6)
num_infector <- sapply(prob_source_BA, function(x) sum( x >= 0.6 ) )
hist(num_infector, breaks = c(-.5,0.5,1.5,2.5,3.5))

### BE trees
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_BE <- lapply(results_BE, function(x) getOffspringDist(x, k = results_BE[[1]][[1]]$ctree$nam ) ) # using default burnin of 0.5, get offspring dist of all BE seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_BE <- lapply(offspring_BE, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# for each of the 100 runs, find how many of the 10 BE individuals were labeled infectors (say prob >= 0.6)
num_infector <- sapply(prob_source_BE, function(x) sum( x >= 0.6 ) )
hist(num_infector, breaks = c(-.5,0.5,1.5,2.5,3.5))

### BQ - do these later on cluster because small RAM :(

### FN trees
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_FN <- lapply(results_FN, function(x) getOffspringDist(x, k = results_FN[[1]][[1]]$ctree$nam ) ) # using default burnin of 0.5, get offspring dist of all FN seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_FN <- lapply(offspring_FN, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# for each of the 100 runs, find how many of the 10 FN individuals were labeled infectors (say prob >= 0.6)
num_infector <- sapply(prob_source_FN, function(x) sum( x >= 0.6 ) )
hist(num_infector, breaks = c(-.5,0.5,1.5,2.5,3.5))

### XBB trees
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_XBB <- lapply(results_XBB, function(x) getOffspringDist(x, k = results_XBB[[1]][[1]]$ctree$nam ) ) # using default burnin of 0.5, get offspring dist of all XBB seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_XBB <- lapply(offspring_XBB, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# for each of the 100 runs, find how many of the 7 XBB individuals were labeled infectors (say prob >= 0.6)
num_infector <- sapply(prob_source_XBB, function(x) sum( x >= 0.6 ) )
hist(num_infector, breaks = c(-.5,0.5,1.5,2.5,3.5))
