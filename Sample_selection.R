# Hi, 
# the following code is building a dataset with all necessary variables for the analysis
# of the association between prenatal depression and hair cortisol at age 6. 
# This includes exposures and outcome of interest, covariates, as well as the auxiliary  
# variables used in the imputation of the final dataset. 

## For this script you are going to need the following datasets from datamanagement
# GR1003-BSI D1_22112016.sav
# CORTISOLHAIRF5_08062015.sav
# CHILD-ALLGENERALDATA_07072020.sav
# GR1003-C_08042016.sav, (GR1003-Family Assessment Device J1-J12_22112016.sav,)
# MATERNALSMOKING_22112016.sav, GEDRAGSGROEP_MaternalDrinking_22112016.sav, 
# (GR1019-E_01072012.sav and BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav)

# Ok let's get started

# Point to the necessary libraries
library(foreign)
library(car)

# check if the path to the data is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") } # -> CUSTOMIZE TO APPROPRIATE PATH /Users/Serena/Desktop/CortisolData/
# The code assumes that all raw data is stored in ONE folder.

#### ----------------------------- FUNCTIONS ----------------------------- ####

# read in data quickly
readquick <- function(filename, rootdir = pathtodata, exclude_col = "") { # only works for SPSS files
  dat <- read.spss(paste(rootdir, filename, sep=""), 
                   use.value.labels = F, to.data.frame = T)
  # Get rid of all capital letters in column names (so you don't have to worry)
  names(dat) <- tolower(names(dat))
  # Replace values of 777, 888 or 999 with NAs unless they are IDCs or IDMs 
  # If you do not want this to happen for any other column use the exclude_col argument. 
  for (i in 1:length(dat)) {
    if (colnames(dat)[i] == "idm" | colnames(dat)[i] == "idc" | colnames(dat)[i] == exclude_col) {
      dat[,i] <- dat[,i]
    } else {
      dat[,i] <- ifelse(dat[,i] == 777 | dat[,i] == 888 | dat[,i] == 999, NA, dat[,i]) }
  } 
  return(dat)
}
#-------------------------------------------------------------------------------
# Calculate the percentage missing data
percent_missing <- function(var) { sum(is.na(var)) / length(var) * 100 }

#### ---------------------------------------------------------------------- ####

# Read in the datasets
pre_bsi <- readquick("GR1003-BSI D1_22112016.sav") # 9778 obs. of 261 vars
hair <- readquick("CORTISOLHAIRF5_08062015.sav") # 6690 obs. of 35 vars

demogr <- readquick("CHILD-ALLGENERALDATA_07072020.sav") # 9901 obs. of 112 vars
social <- readquick("GR1003-C_08042016.sav")  #  9778 obs of 90 vars
smokingv1 <- readquick("MATERNALSMOKING_22112016.sav") #  9778 obs of 11 vars
drinkingv1 <- readquick("GEDRAGSGROEP_MaternalDrinking_22112016.sav") #  9778 obs of 18 vars

################################################################################
############################# PRENATAL DEPRESSION ##############################
################################################################################

# Depressive symptoms were assessed at about 20.6 weeks of gestation with the 
# Brief Symptom Inventory (BSI); a validated self‐report questionnaire defining a 
# spectrum of psychiatric symptoms. We use the 6‐item depressive symptoms scale (total scores). 
# No more than one missing item was allowed to minimize selective non-response. 
# Internal consistency for the depression scale: α = .80. 
# Mothers with a score higher than 0.75 have clinically relevant depressive symptoms 
# according to the Dutch manual. 
predep <- data.frame(pre_bsi$idm, pre_bsi$dep, # min: 0 ; max: 4
                     ifelse(pre_bsi$dep > 0.75, yes = 1, no = 0)) # 1 = depression; 0 = no depression
colnames(predep) <- c("IDM", "cont_prenatal_depression", "prenatal_depression")


################################################################################
########################## CHILD HAIR CORTISOL @6y #############################
################################################################################

haircort <- data.frame(hair$idc, hair$agechildyf5, 
                       hair$hair_cortisol_pgmg) # min: 0.076 ; max: 292.736
colnames(haircort) <- c("IDC", "child_age", "cortisol")

################################################################################
################################ COVARIATES ####################################
################################################################################

child_general = demogr[,c( 'idc', 'idm', 
                     'gender',    # child sex
                     'ethnfv2', # ethnicity original variable
                     'educm',     # highest education finished (mother, during pregnancy)
                     'income5',   # net income of household when the child is 5. 
                     'twin',      # used for exclusion criteria 
                     'mother',    # mother id used to identify siblings (for exclusion)
                     'parity',    # parity (used for imputation)
                     'gestbir',   # gestational age at birth (used for imputation)
                     'weight',    # gestational weight (used for imputation)
                     'bmi_0',     # maternal BMI (self-reported), before pregnancy
                     'bmi_1',     # maternal BMI during pregnancy (used for imputation)
                     'age_m_v2')] # maternal age at intake (used for imputation) 
# Again, let's try to keep it user friendly 
colnames(child_general) = c("IDC", "IDM", "sex", "ethnicity_orig", "m_education_pregnancy", 
                            "income_5y", "twin", "mother_id", "parity", "gest_age_birth", 
                           "gest_weight", "m_bmi_berore_pregnancy", "m_bmi_pregnancy", "m_age_cont")

# Ethnicity recode – dichotomized into: dutch, western and non-western;
for (i in 1:9901) {
  if (is.na(child_general$ethnicity_orig[i])) { child_general$ethnicity[i] <- NA
  } else if (child_general$ethnicity_orig[i] == 1 | child_general$ethnicity_orig[i] == 300 | 
             child_general$ethnicity_orig[i] == 500 | child_general$ethnicity_orig[i] >= 700) { 
    # Dutch (1), American, western (300) Asian, western (500) European (700), Oceanie (800)
    child_general$ethnicity[i] <- 1 
  } else { 
    child_general$ethnicity[i] <- 2 } 
    # Indonesian (2), Cape Verdian (3), Maroccan (4) Dutch Antilles (5) Surinamese 
    # (6) Turkish (7) African (200), American, non western (400), Asian, non western (600)
} 
#-------------------------------------------------------------------------------
### MATERNAL SMOKING during pregnancy 

smoking = smokingv1[,c('idm', 'smoke_all')] # (1) never a smoker; 
# (2) smoked until pregnancy was known (i.e., first trimester only); 
# (3) continued smoking during pregnancy.
colnames(smoking) = c("IDM", "m_smoking")

#-------------------------------------------------------------------------------
### MATERNAL ALCOHOL CONSUMPTION during pregnancy

drinking = drinkingv1[,c('idm', 'mdrink_updated')] # (0) never; 
# (1) until pregnancy was known (i.e., first trimester only); 
# (2) continued during pregnancy occasionally;
# (3) continued during pregnancy frequently.
colnames(drinking) = c("IDM", "m_drinking")

#-------------------------------------------------------------------------------
#  social support disuring pregnancy 

# First I need to recode some "Not applicable" answers into NAs.
social$c0700103[social$c0700103 == 5] <- NA

support = social[, c('idm', 'c0700103')] # Difficulties between you and your partner? "No"/ "Slight"/"Moderate"/"Serious"
colnames(support) = c("IDM", "support_partner")

################################################################################
# Merge all variables together
depre_covs <- Reduce(function(x,y) merge(x = x, y = y, by = 'IDM',  all.x = TRUE),
                                   list(predep, child_general, smoking, drinking, support))

full_dataset = merge(depre_covs, haircort, by = "IDC" , all.x = TRUE)
################################################################################
##----------------------------------------------------------------------------##
## -------------------- Exclude participants (flowchart) -------------------- ##
##----------------------------------------------------------------------------##

initial_sample <- 9901

## First exclusion step:

# Exclude children with missing internalizing score
outcome_meas <- full_dataset[!is.na(full_dataset$cortisol),] 
after_cort <- nrow(outcome_meas)
suba <- after_cort - initial_sample

# Exclude all non-white participants 
white <- outcome_meas[outcome_meas$ethnicity == 1, ]
after_racism <- nrow(white)
subb <- after_racism - after_cort

# Exclude twins
no_twins <- white[white$twin == 0, ]
after_twins <- nrow(no_twins)
subc <- after_twins - after_racism

# Flowchart
flowchart <- list(initial_sample, suba, after_cort, subb, after_racism, subc, 
                  after_twins)

# Rename final dataset:
dataset <- no_twins
cat(paste("Well, congrats! Your final dataset includes", after_twins ,"participants.")) # 1886
