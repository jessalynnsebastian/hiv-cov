source("data_setup_BEAST.R")

prob_source_BA <- readRDS("./prob_source/prob_source_BA.rds")
prob_source_BE <- readRDS("./prob_source/prob_source_BE.rds")
prob_source_BQ <- readRDS("./prob_source/prob_source_BQ.rds")
prob_source_FN <- readRDS("./prob_source/prob_source_FN.rds")
prob_source_XBB <- readRDS("./prob_source/prob_source_XBB.rds")

### Mark individuals as infectors/non-infectors
# concatenate groups and (for now) make the cutoff a probability of 0.1 just so i can run the model
prob_infector <- c(prob_source_BA, prob_source_BE, prob_source_BQ, prob_source_FN, prob_source_XBB)
ind_infector <- as.numeric(prob_infector >= 0.1)
sum(ind_infector)
### Get columns of metadata to include as adjustment covariates
# below are the ones without too much missingness
covariates <- meta[,c("new_hiv_stat", "cal_age", "gender", "cvid_syp", "vacc_sts", "travel", "risk_lvl")]
dat <- cbind(ind_infector, covariates)
dat$ind_infector <- as.factor(ind_infector)
dat$new_hiv_stat <- as.factor(dat$new_hiv_stat)
dat <- dat[dat$new_hiv_stat != "Unknown",] # remove unknown hiv stat

# for now not including covariates because there aren't very many cases
test_model <- glm(ind_infector ~ new_hiv_stat, data = dat, family = "binomial") # can't use cvd syp
summary(test_model)

