
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
library(mice)
library(miceadds)
library(openxlsx)

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#------------------------------------------------------------------------------#

# Load the list of imputed dataset
impA_orig <- readRDS(paste(pathtodata,'sampleA.rds', sep = ""))
impB_orig <- readRDS(paste(pathtodata,'sampleB.rds', sep = ""))
#------------------------------------------------------------------------------#
# Before running the models, we need to standardize all numeric variables so 
# that the betas obtained are scaled. 

# Transform mids into datalist
datlistA <- miceadds::mids2datlist(impA_orig)
datlistB <- miceadds::mids2datlist(impB_orig)

varnames <- c('bsiScore', 'lncortisol', 'childAge', 'prePregBMI', 
              'familySS', 'maternalSmoking', 'maternalAge', 'gestAge')

# Scale
newA = miceadds::scale_datlist(datlistA, varnames, 
                               paste0("z", varnames), weights=NULL, M=0, SD=1,
                               digits=NULL)
newB = miceadds::scale_datlist(datlistB, varnames, 
                               paste0("z", varnames), weights=NULL, M=0, SD=1,
                               digits=NULL)

# Transform back into mids objects
impA <- miceadds::datlist2mids(newA)
impB <- miceadds::datlist2mids(newB)

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
# social support, marital status and smoking status.
# For all adjusted models we also run a version including sex as an effect modifier.

# There we go:

# base model
fit1 <- with(impA, lm(zlncortisol ~ bsiScoreFlg))
mod1 <- summary(pool(fit1)); mod1$sign <- ifelse(mod1$p.value < 0.05, "*", " ")

# minimally adjusted model
fit2 <- with(impA, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge))
mod2 <- summary(pool(fit2)); mod2$sign <- ifelse(mod2$p.value < 0.05, "*", " ")

# fully adjusted model
fit3 <- with(impA, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome))
mod3 <- summary(pool(fit3)); mod3$sign <- ifelse(mod3$p.value < 0.05, "*", " ")

# additionally adjusted model
fit4 <- with(impA, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod4 <- summary(pool(fit4)); mod4$sign <- ifelse(mod4$p.value < 0.05, "*", " ")

# minimally adjusted model with interaction
fit2int <- with(impA, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge))
mod2int <- summary(pool(fit2int)); mod2int$sign <- ifelse(mod2int$p.value < 0.05, "*", " ")

# fully adjusted model with interaction
fit3int <- with(impA, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome))
mod3int <- summary(pool(fit3int)); mod3int$sign <- ifelse(mod3int$p.value < 0.05, "*", " ")

# additionally adjusted model with interaction
fit4int <- with(impA, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod4int <- summary(pool(fit4int)); mod4int$sign <- ifelse(mod4int$p.value < 0.05, "*", " ")

################################################################################

#                           SAMPLE B - ALL ETHNICITIES                               
#------------------------------------------------------------------------------#
# Now let's examine the full dataset, including all ethnicities.

# base model
fit5 <- with(impB, lm(zlncortisol ~  bsiScoreFlg))
mod5 <- summary(pool(fit5)); mod5$sign <- ifelse(mod5$p.value < 0.05, "*", " ")

# minimally adjusted model
fit6 <- with(impB, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge + childRaceEth))
mod6 <- summary(pool(fit6)); mod6$sign <- ifelse(mod6$p.value < 0.05, "*", " ")

# fully adjusted model
fit7 <- with(impB, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome))
mod7 <- summary(pool(fit7)); mod7$sign <- ifelse(mod7$p.value < 0.05, "*", " ")

# additionally adjusted model
fit8 <- with(impB, lm(zlncortisol ~  bsiScoreFlg + sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod8 <- summary(pool(fit8)); mod8$sign <- ifelse(mod8$p.value < 0.05, "*", " ")

# minimally adjusted model with interaction 
fit6int <- with(impB, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge + childRaceEth))
mod6int <- summary(pool(fit6int)); mod6int$sign <- ifelse(mod6int$p.value < 0.05, "*", " ")

# fully adjusted model with interaction 
fit7int <- with(impB, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome))
mod7int <- summary(pool(fit7int)); mod7int$sign <- ifelse(mod7int$p.value < 0.05, "*", " ")

# additionally adjusted model with interaction 
fit8int <- with(impB, lm(zlncortisol ~  bsiScoreFlg * sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod8int <- summary(pool(fit8int)); mod8int$sign <- ifelse(mod8int$p.value < 0.05, "*", " ")

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
fit9 <- with(impA, lm(zlncortisol ~  zbsiScore))
mod9 <- summary(pool(fit9)); mod9$sign <- ifelse(mod9$p.value < 0.05, "*", " ")

# minimally adjusted model
fit10 <- with(impA, lm(zlncortisol ~  zbsiScore + sex + zchildAge))
mod10 <- summary(pool(fit10)); mod10$sign <- ifelse(mod10$p.value < 0.05, "*", " ")

# fully adjusted model
fit11 <- with(impA, lm(zlncortisol ~  zbsiScore + sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome))
mod11 <- summary(pool(fit11)); mod11$sign <- ifelse(mod11$p.value < 0.05, "*", " ")

# additionally adjusted model
fit12 <- with(impA, lm(zlncortisol ~  zbsiScore + sex + zchildAge + 
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod12 <- summary(pool(fit12)); mod12$sign <- ifelse(mod12$p.value < 0.05, "*", " ")

# minimally adjusted model with interaction
fit10int <- with(impA, lm(zlncortisol ~  zbsiScore * sex + zchildAge))
mod10int <- summary(pool(fit10int)); mod10int$sign <- ifelse(mod10int$p.value < 0.05, "*", " ")

# fully adjusted model with interaction
fit11int <- with(impA, lm(zlncortisol ~  zbsiScore * sex + zchildAge + 
                         zprePregBMI + maternalEdu + maternalIncome))
mod11int <- summary(pool(fit11int)); mod11int$sign <- ifelse(mod11int$p.value < 0.05, "*", " ")

# additionally adjusted model with interaction
fit12int <- with(impA, lm(zlncortisol ~  zbsiScore * sex + zchildAge + 
                         zprePregBMI + maternalEdu + maternalIncome + 
                         zfamilySS + maternalMaritalStatus + maternalSmoking))
mod12int <- summary(pool(fit12int)); mod12int$sign <- ifelse(mod12int$p.value < 0.05, "*", " ")


################################################################################

#                           SAMPLE B - ALL ETHNICITIES                               
#------------------------------------------------------------------------------#

# Now let's examine the full dataset, including all ethnicities. 
# base model
fit13 <- with(impB, lm(zlncortisol ~  zbsiScore))
mod13 <- summary(pool(fit13)); mod13$sign <- ifelse(mod13$p.value < 0.05, "*", " ")

# minimally adjusted model
fit14 <- with(impB, lm(zlncortisol ~  zbsiScore + sex + zchildAge + childRaceEth))
mod14 <- summary(pool(fit14)); mod14$sign <- ifelse(mod14$p.value < 0.05, "*", " ")

# fully adjusted model
fit15 <- with(impB, lm(zlncortisol ~  zbsiScore + sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome))
mod15 <- summary(pool(fit15)); mod15$sign <- ifelse(mod15$p.value < 0.05, "*", " ")

# additionally adjusted model
fit16 <- with(impB, lm(zlncortisol ~  zbsiScore + sex + zchildAge + childRaceEth +
                        zprePregBMI + maternalEdu + maternalIncome + 
                        zfamilySS + maternalMaritalStatus + maternalSmoking))
mod16 <- summary(pool(fit16)); mod16$sign <- ifelse(mod16$p.value < 0.05, "*", " ")

# minimally adjusted model with interaction
fit14int <- with(impB, lm(zlncortisol ~  zbsiScore * sex + zchildAge + childRaceEth))
mod14int <- summary(pool(fit14int)); mod14int$sign <- ifelse(mod14int$p.value < 0.05, "*", " ")

# fully adjusted model with interaction
fit15int <- with(impB, lm(zlncortisol ~  zbsiScore * sex + zchildAge + childRaceEth +
                         zprePregBMI + maternalEdu + maternalIncome))
mod15int <- summary(pool(fit15int)); mod15int$sign <- ifelse(mod15int$p.value < 0.05, "*", " ")

# additionally adjusted model with interaction
fit16int <- with(impB, lm(zlncortisol ~  zbsiScore * sex + zchildAge + childRaceEth +
                         zprePregBMI + maternalEdu + maternalIncome + 
                         zfamilySS + maternalMaritalStatus + maternalSmoking))
mod16int <- summary(pool(fit16int)); mod16int$sign <- ifelse(mod16int$p.value < 0.05, "*", " ")

################################################################################
################################################################################

# Export the outputs of summary statistics into an xlsx file 

modls <- list("1.base_white" = mod1, "2.min_white" = mod2, 
              "3.full_white" = mod3, "4.addit_white" = mod4,
              "2int.min_white" = mod2int, "3int.full_white" = mod3int, "4int.addit_white" = mod4int,
              "5.base_all" = mod5, "6.min_all" = mod6, 
              "7.full_all" = mod7, "8.addit_all" = mod8, 
              "6int.min_all" = mod6int, "7int.full_all" = mod7int, "8int.addit_all" = mod8int, 
              "9.cont_base_white" = mod9, "10.cont_min_white" = mod10, 
              "11.cont_full_white" = mod11, "12.cont_addit_white" = mod12,
              "10int.cont_min_white" = mod10int, "11int.cont_full_white" = mod11int, "12int.cont_addit_white" = mod12int,
              "13.cont_base_all" = mod13, "14.cont_min_all" = mod14, 
              "15.cont_full_all" = mod15, "16.cont_addit_all" = mod16, 
              "14int.cont_min_all" = mod14int, "15int.cont_full_all" = mod15int, "16int.cont_addit_all" = mod16int)

write.xlsx(modls, file = paste0(pathtodata, "Results.xlsx"))


################################################################################
############################### THE END ########################################
################################################################################