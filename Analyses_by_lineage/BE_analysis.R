source("./data_setup_BEAST.R")
# set seed
set.seed(4935)

# get BEAST results, sample 100, convert to TP ptrees
tree_BE <- ape::read.nexus(file = "../BE_beast2/BE-BE_sequences.trees")
# discard 10% burnin
tree_BE <- tree_BE[2001:20001]
# take sample of 100
tree_BE <- sample(tree_BE, size = 100)
# convert to TP format
tree_BE <- lapply(tree_BE, function(t) toTransPhylo(t) )
# and run TP
if (!file.exists("../results/tp_BE.rds")){
  tp_res_BE <- infer_multittree_share_param( ptree_lst = tree_BE,
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
  saveRDS(tp_res_BE, "../results/tp_BE.rds")
}

## Read in results and label infectors vs. non-infectors
tp_res_BE <- readRDS("../results/tp_BE.rds")
# in each of the 100 runs, get offspring distribution for all sampled individuals
offspring_BE <- lapply(tp_res_BE, function(x) getOffspringDist(x, k = tree_BE[[1]]$nam ) ) # using default burnin of 0.5, get offspring dist of all BE seqs
# then get the proportion of samps for each host where they infected at least one other individual
prob_source_BE <- lapply(offspring_BE, function(x) {1 - rowSums(x == 0)/ncol(x)} )
# each run had the same number of posterior samples so the proBEbility for each individual is just the arithmetic mean of the proBEbilities in each run
prob_source_BE <- rowMeans(simplify2array(prob_source_BE))
saveRDS(prob_source_BE, "../prob_source/prob_source_BE.rds")
