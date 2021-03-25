
# Final step in our journey,
# The following code runs a set of analyses for the project "Association of Prenatal 
# Depression with Cord Blood Glucocorticoids and Child Cortisol Levels" in generation R 

# All it requires is just the two imputed dataset we built using the "Imputation.R"
# script, that you can find in this repository. 
# In the previous step, we created several complete versions of the data by replacing
# the missing values with plausible data values. Now we need to estimate the 
# parameters of interest from each imputed dataset and pool the estimates into one.

# Ok, let's get started!

# As usual, here are the packages we need 
library(mice);
library(openxlsx)

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#------------------------------------------------------------------------------#

# Load the list of imputed dataset
impA <- readRDS(paste(pathtodata,'sampleA.rds', sep = ""))
impB <- readRDS(paste(pathtodata,'sampleB.rds', sep = ""))

################################################################################
############################## DICH DEPRESSION #################################
################################################################################

#                           SAMPLE A - WHITE ONLY                               
#------------------------------------------------------------------------------#
# First, we assess the effect of prenatal stress (i.e. maternal depression) on 
# child cortisol levels at age 6. Then we perform two levels of adjustment: a 
# minimally adjusted (sex and age) and fully adjusted version (sex, age, maternal 
# pre pregnancy BMI, education, and household income). 
# In a last version we also adjust for other potential confounders: mid-pregnancy 
# social support, marital status, smoking status, and parity.

# There we go:

# base model
fit1 <- with(impA, lm(lncortisol ~  bsiScoreFlg))
mod1 <- summary(pool(fit1))

# minimally adjusted model
fit2 <- with(impA, lm(lncortisol ~  bsiScoreFlg + sex + childAge))
mod2 <- summary(pool(fit2))

# fully adjusted model
fit3 <- with(impA, lm(lncortisol ~  bsiScoreFlg + sex + childAge + 
                        prePregBMI + maternalEdu + maternalIncome))
mod3 <- summary(pool(fit3))

# additionally adjusted model
fit4 <- with(impA, lm(lncortisol ~  bsiScoreFlg + sex + childAge + 
                        prePregBMI + maternalEdu + maternalIncome + 
                        familySS + maternalMaritalStatus + maternalSmoking))
mod4 <- summary(pool(fit4))

################################################################################

#                           SAMPLE B - ALL ETHNICITIES                               
#------------------------------------------------------------------------------#
# Now let's examine the full dataset, including all ethnicities.

# base model
fit5 <- with(impB, lm(lncortisol ~  bsiScoreFlg))
mod5 <- summary(pool(fit5))

# minimally adjusted model
fit6 <- with(impB, lm(lncortisol ~  bsiScoreFlg + sex + childAge + childRaceEth))
mod6 <- summary(pool(fit6))

# fully adjusted model
fit7 <- with(impB, lm(lncortisol ~  bsiScoreFlg + sex + childAge + childRaceEth +
                        prePregBMI + maternalEdu + maternalIncome))
mod7 <- summary(pool(fit7))

# additionally adjusted model
fit8 <- with(impB, lm(lncortisol ~  bsiScoreFlg + sex + childAge + childRaceEth +
                        prePregBMI + maternalEdu + maternalIncome + 
                        familySS + maternalMaritalStatus + maternalSmoking))
mod8 <- summary(pool(fit8))


################################################################################
############################## CONT DEPRESSION #################################
################################################################################

#                           SAMPLE A - WHITE ONLY                               
#------------------------------------------------------------------------------#
# Second, we assess the effect of maternal depressive symptomatology  on 
# child cortisol levels at age 6. Then we perform two levels of adjustment: a 
# minimally adjusted (sex and age) and fully adjusted version (sex, age, maternal 
# pre pregnancy BMI, education, and household income). 
# In a last version we also adjust for other potential confounders: mid-pregnancy 
# social support, marital status, smoking status, and parity.

# There we go:

# base model
fit9 <- with(impA, lm(lncortisol ~  bsiScore))
mod9 <- summary(pool(fit9))

# minimally adjusted model
fit10 <- with(impA, lm(lncortisol ~  bsiScore + sex + childAge))
mod10 <- summary(pool(fit10))

# fully adjusted model
fit11 <- with(impA, lm(lncortisol ~  bsiScore + sex + childAge + 
                        prePregBMI + maternalEdu + maternalIncome))
mod11 <- summary(pool(fit11))

# additionally adjusted model
fit12 <- with(impA, lm(lncortisol ~  bsiScore + sex + childAge + 
                        prePregBMI + maternalEdu + maternalIncome + 
                        familySS + maternalMaritalStatus + maternalSmoking))
mod12 <- summary(pool(fit12))

################################################################################

#                           SAMPLE B - ALL ETHNICITIES                               
#------------------------------------------------------------------------------#

# Now let's examine the full dataset, including all ethnicities. 
# base model
fit13 <- with(impB, lm(lncortisol ~  bsiScore))
mod13 <- summary(pool(fit13))

# minimally adjusted model
fit14 <- with(impB, lm(lncortisol ~  bsiScore + sex + childAge + childRaceEth))
mod14 <- summary(pool(fit14))

# fully adjusted model
fit15 <- with(impB, lm(lncortisol ~  bsiScore + sex + childAge + childRaceEth +
                        prePregBMI + maternalEdu + maternalIncome))
mod15 <- summary(pool(fit15))

# additionally adjusted model
fit16 <- with(impB, lm(lncortisol ~  bsiScore + sex + childAge + childRaceEth +
                        prePregBMI + maternalEdu + maternalIncome + 
                        familySS + maternalMaritalStatus + maternalSmoking))
mod16 <- summary(pool(fit16))

################################################################################
################################################################################

# Export the outputs of summary statistics into an xlsx file 

modls <- list("1.base_white" = mod1, "2.min_white" = mod2, 
              "3.full_white" = mod3, "4.addit_white" = mod4, 
              "5.base_all" = mod5, "6.min_all" = mod6, 
              "7.full_all" = mod7, "8.addit_all" = mod8, 
              "9.cont_base_white" = mod9, "10.cont_min_white" = mod10, 
              "11.cont_full_white" = mod11, "12.cont_addit_white" = mod12, 
              "13.cont_base_all" = mod13, "14.cont_min_all" = mod14, 
              "15.cont_full_all" = mod15, "16.cont_addit_all" = mod16)

write.xlsx(modls, file = paste0(pathtodata, "Results.xlsx"))


################################################################################
############################### THE END ########################################
################################################################################