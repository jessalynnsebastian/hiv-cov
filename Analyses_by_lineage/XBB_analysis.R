source("./data_setup_BEAST.R")
# set seed
set.seed(4935)

# get BEAST results, sample 100, convert to TP ptrees
tree_XBB <- ape::read.nexus(file = "../XBB_beast2/XBB-XBB_sequences.trees")
# discard 10% burnin
tree_XBB <- tree_XBB[2001:20001]
# take sample of 100
tree_XBB <- sample(tree_XBB, size = 100)
# convert to TP format
tree_XBB <- lapply(tree_XBB, function(t) toTransPhylo(t) )
# # and run TP
if (!file.exists("../results/tp_BA.rds")){
  tp_res_XBB <- infer_multittree_share_param( ptree_lst = tree_XBB,
                                             w.shape = (9/2), # shape param of gamma for generation time https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/#:~:text=When%20allowing%20only%20positive%20serial,estimated%20to%20be%202.57%20days.
                                             w.scale = (2/3),
                                             ws.shape = (2), # shape param of gamma for sampling time
                                             ws.scale = (3/2),
                                             mcmcIterations = 10000,
                                             thinning = 10,
                                             startNeg = (3), # starting values for: within host coalescent param/effective population size
                                             startOff.r = 1, # r in negative binomial offspring dist
                                             startOff.p = 0.5, # p in above
                                             startPi = 0.1, # sampling proportion pi
                                             optiStart = 2, # type of optimisation to apply to MCMC (0 = none, 1 = slow, 2 = fast)
                                             dateT = dateLastSample, # date process stops - should only be inf for outbreak that's over
                                             delta_t = 0.5,
                                             verbose = F
  )
  # save TP results as R object
  saveRDS(tp_res_XBB, "../results/tp_XBB.rds")
}

## READ IN RESULTS FROM HPC3 & LABEL INFECTOR VS. NOT
tp_res_XBB <- readRDS("../results/tp_XBB.rds")
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_XBB <- lapply(tp_res_XBB, function(x) getOffspringDist(x, k = tree_XBB[[1]]$nam ) ) # using default burnin of 0.5, get offspring dist of all XBB seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_XBB <- lapply(offspring_XBB, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# each run had the same number of posterior samples so the proXBBbility for each individual is just the arithmetic mean of the proXBBbilities in each run
prob_source_XBB <- rowMeans(simplify2array(prob_source_XBB))
names(prob_source_XBB) <- tree_XBB[[1]]$nam
saveRDS(prob_source_XBB, "../prob_source/prob_source_XBB.rds")
