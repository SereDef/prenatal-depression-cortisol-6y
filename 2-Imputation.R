
# Hi again, 
# the following code runs the imputation of missing data, necessary for 
# the analysis of the association between prenatal depression and hair cortisol
# in children. 

# All you need it the file we created in the 'Sample_selection.R' script: the 
# "dataset_raw.rds" 
# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

# Point to the necessary libraries
library(mice);
library(miceadds);

# check if the path to the datasets is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

################################################################################
# Load datasets
depcor_raw <- readRDS(paste(pathtodata, 'dataset_raw.rds', sep = ""))

#-------------------------------------------------------------------------------
# Calculate the percentage missing data for every column
percent_missing <- function(var) { sum(is.na(var)) / length(var) * 100 }
lapply(depcor_raw, percent_missing)

# Let's select only the variables we need and also determine their order that is important
# for mice to work 
depcor <- depcor_raw[, c('idc', 'bsiScore', 'bsiScoreFlg', 'cortisol', 'lncortisol', # main
                         'sex', 'childAge', # minimally adjusted
                         'prePregBMI', 'maternalEdu', 'maternalIncome', # fully adjusted 
                         'familySS', 'maternalMaritalStatus', 'maternalSmoking', 'parity', # confounders ?
                         'maternalAge', 'gestAge', 'childRaceEth')] 

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We started with a dry run to specify the default arguments.
imp0 <- mice(depcor, maxit = 0, 
             defaultMethod = rep('pmm',4)) # set the imputation method to predictive mean matching (PMM)* 

# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- make.method(depcor)

# We use passive imputation for dichotomous depression. This means that the continuous 
# score is imputed first, and then, using complete version, dichotomized values are 
# derived by the formula specified below.
meth['bsiScoreFlg'] <- "~I( ifelse( bsiScore > 0.75, yes = 1, no = 0) )" 

# We are going to need a different set of predictors for the different variables 
# so let's define them using the predictormatrix, that gives instructions to mice
predictormatrix <- imp0$predictorMatrix

# Do not use IDC as predictor:
predictormatrix[, "idc"] <- predictormatrix["idc",] <- 0
predictormatrix[, "bsiScoreFlg"] <- predictormatrix["bsiScoreFlg",] <- 0
# maternalAge, childAge and sex are complete

### Impute auxiliary variables and covariates ###

# Auxiliary variables were selected because they are either related to missingness 
# or to the variables themselves. 
# To prevent multicollinearity, we adjust the predictor matrix such that the 
# not all available items are used, thereby also reducing computational load.
# The technique has been found to reduce standard error substantially compared to 
# complete-case analysis (Plumpton et al., 2010), and it outperforms other existing 
# techniques (Eekhout et al., 2018). 

# So, let's adjust the predictormatrices such that ‘redundant’ items were not used as a predictor.

predictormatrix[c('prePregBMI', 'maternalEdu', 'maternalIncome', 'familySS', 'maternalMaritalStatus', 'maternalSmoking'), 
                     # USE: each other (except for familySS) & maternalAge
                c('bsiScore', 'bsiScoreFlg', 'cortisol', 'lncortisol', 'sex', 'childAge', 
                  'familySS', 'parity', 'gestAge', 'childRaceEth')] <- 0

predictormatrix[c('bsiScore'), # USE: maternalAge, prePregBMI, maternalEdu, maternalSmoking, 
                               #      maternalMaritalStatus, maternalIncome, familySS
                c('bsiScoreFlg', 'cortisol', 'lncortisol', 'sex', 'childAge',
                  'parity', 'gestAge', 'childRaceEth')] <- 0

predictormatrix[c('gestAge'), # USE: childRaceEth, sex
                c('bsiScore', 'bsiScoreFlg', 'cortisol', 'lncortisol', 'childAge', 
                  'prePregBMI', 'maternalEdu', 'maternalIncome', 
                  'familySS', 'maternalMaritalStatus', 'maternalSmoking', 'parity', 
                  'maternalAge')] <- 0

predictormatrix[c('parity'), # USE: maternalAge
                c('bsiScore', 'bsiScoreFlg', 'cortisol', 'lncortisol', 'sex', 'childAge', 
                  'prePregBMI', 'maternalEdu', 'maternalIncome', 
                  'familySS', 'maternalMaritalStatus', 'maternalSmoking', 'gestAge', 'childRaceEth')] <- 0

# OPTIONAL :Quickly check the matrix to make sure it looks legit
# library(pheatmap)
# pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain
# low computational load, we do 20 iterations using 5 imputed datasets

imputation <- mice(depcor, m = 5, # nr of imputed datasets
                   maxit = 20, # nr of iteration taken to impute missing values
                   seed = 373844, # set a seed for the random number generation in case i need to generate the same dataset again
                   method = meth,
                   visitSequence = VisSeq, 
                   predictorMatrix = predictormatrix)

### There were no logged events in the imputation. 
### Visual inspection of the convergence graphs showed convergence after 20 iterations.
# plot(imputation)

################### OPTIONAL CHECKS (beware: it takes time) ####################
# # Inspecting the distribution of observed and imputed values
# stripplot(imputation, pch = 20, cex = 1.2) # red dots are imputed values
# # A scatter plot is also useful to spot unusual patterns in two vars
# xyplot(imputation, bsiScore ~ familySS | .imp, pch = 20, cex = 1.4)

#------------------------------------------------------------------------------#
##------------------------ Two samples (ethnicity) -------------------------- ##
#------------------------------------------------------------------------------#

# Select only white participants. Note because ethnicity is imputed, each dataset 
# will have different number of obs. Also note this is not a mids but a datlist object. 
after_racism <- miceadds::subset_datlist(imputation, expr_subset = expression(childRaceEth == 0))

################################################################################
#### ------------------------- complete and save -------------------------- ####
################################################################################

# Let's save the mids object (i.e. list of imputed datasets)
saveRDS(after_racism, paste(pathtodata,'sampleA.rds', sep = "")) # white only 
saveRDS(imputation, paste(pathtodata,'sampleB.rds', sep = "")) # complete set

# I also save the last imputed dataset for sanity checks
depcor_imputed <- complete(imputation, 5) 
saveRDS(depcor_imputed, paste(pathtodata,'B_imputed5.rds', sep = ""))


