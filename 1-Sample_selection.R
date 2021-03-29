# Hi, 
# the following code is building a dataset with all necessary variables for the analysis
# of the association between prenatal depression and hair cortisol at age 6. 
# This includes exposures and outcome of interest, covariates, as well as the auxiliary  
# variables used in the imputation of the final dataset. 

## For this script you are going to need the following datasets from datamanagement
# GR1003-BSI D1_22112016.sav
# CORTISOLHAIRF5_08062015.sav
# CHILD-ALLGENERALDATA_07072020.sav
# GR1003-Family Assessment Device J1-J12_22112016.sav,
# MATERNALSMOKING_22112016.sav,
# (GR1019-E_01072012.sav and BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav)

# Ok let's get started

# Point to the necessary libraries
library(foreign)
library(car)

# check if the path to the data is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") } # -> CUSTOMIZE TO APPROPRIATE PATH 
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
# General functioning of the family was measured with the family assessment device (FAD).
# This function calculates FAD scores accounting for (max 25%) missing values and 
# dichotomizes the scores according to the cutoff for unhealthy family functioning 
# provided by the manual (Byles et al., 1988), = 2.17. 
fad_scores <- function(set, dich = FALSE){
  fadtotal <- rowSums(abs(set)) # Sum items
  no_na_fad <- rowSums(!is.na(set)) # Number of endorsed items
  # Compute mean score + do not calculate when more than 25% of items are missing
  fad_score <- ifelse(no_na_fad >= 9, yes = fadtotal/no_na_fad, no = NA) 
  if (dich == T) { fad_score <- ifelse(fad_score > 2.17, yes = 1, no = 0) }
  return(fad_score) 
}

#-------------------------------------------------------------------------------
# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

#### ---------------------------------------------------------------------- ####

# Read in the datasets
pre_bsi <- readquick("GR1003-BSI D1_22112016.sav") # 9778 obs. of 261 vars
hair <- readquick("CORTISOLHAIRF5_08062015.sav") # 6690 obs. of 35 vars

demogr <- readquick("CHILD-ALLGENERALDATA_07072020.sav") # 9901 obs. of 112 vars
smokingv1 <- readquick("MATERNALSMOKING_22112016.sav") #  9778 obs of 11 vars
# drinkingv1 <- readquick("GEDRAGSGROEP_MaternalDrinking_22112016.sav") # excluded from analysis because not available in PV
social <- readquick("GR1003-Family Assessment Device J1-J12_22112016.sav")

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
predep <- data.frame(pre_bsi$idm, 
                     pre_bsi$dep, # min: 0 ; max: 4 
                     ifelse(pre_bsi$dep > 0.75, yes = 1, no = 0)) # 1 = depression; 0 = no depression
colnames(predep) <- c("idm", "bsiScore", "bsiScoreFlg")

################################################################################
########################## CHILD HAIR CORTISOL @6y #############################
################################################################################

haircort <- data.frame(hair$idc, hair$agechildyf5, 
                       hair$hair_cortisol_pgmg, 
                       log(hair$hair_cortisol_pgmg))
colnames(haircort) <- c("idc", "childAge", "cortisol", "lncortisol")

################################################################################
################################ COVARIATES ####################################
################################################################################

child_general = demogr[,c( 'idc', 'idm', 
                     'gender',    # child sex
                     'ethnfv2',   # ethnicity original variable
                     'educm',     # highest education finished (mother, during pregnancy)
                     'income5',   # net income of household when the child is 5. 
                     'twin',      # used for exclusion criteria 
                     'mother',    # mother id used to identify siblings (for exclusion)
                     'parity',    # parity (used for imputation)
                     'gestbir',   # gestational age at birth (used for imputation)
                     'bmi_0',     # maternal BMI (self-reported), before pregnancy
                     'mardich',   # marital status during pregnancy
                     'age_m_v2')] # maternal age at intake (used for imputation) 
# Again, let's try to keep it user friendly 
colnames(child_general)[c(10:11, 13)] = c("gestAge", "prePregBMI", "maternalAge")

# Ethnicity recode – dichotomized into: western and non-western;
child_general$childRaceEth <- ifelse(child_general$ethnfv2 == 1 | child_general$ethnfv2 == 300 | 
                                       child_general$ethnfv2 == 500 | child_general$ethnfv2 >= 700, 
                                    # Dutch (1), American, western (300) Asian, western (500) European (700), Oceanie (800)
                                    yes = 0, # White 
                                    no = 1) # Non-white
                                    # Indonesian (2), Cape Verdian (3), Maroccan (4) Dutch Antilles (5) Surinamese 
                                    # (6) Turkish (7) African (200), American, non western (400), Asian, non western (600)
# Sex recode 0 = boy, 1 = girl
child_general$sex <- child_general$gender - 1 

# Marital status recode 0 = no partner 1 = married or cohabiting
child_general$maternalMaritalStatus <- recode(child_general$mardich, '1=1; 2=0')

# Maternal education recode 
# DICH: "No education"/"Primary"/"Secondary-phase 1"/"Secondary-phase 2" = 0 "Higher-phase 1"/"Higher-phase 2" = 1.
# based on Centraal Bureau voor de Statistiek (2016).
child_general$maternalEdu <- ifelse(child_general$educm <= 3, yes = 0, no = 1) 

# Income recode 
# DICH: according to the Central Statistic Netherlands (2013). 
# Net household income below 1600 €/month (basic needs level) = risk.
child_general$maternalIncome <- ifelse(child_general$income5 < 4, yes = 0, no = 1)

#-------------------------------------------------------------------------------
# MATERNAL SMOKING during pregnancy 

smoking <- data.frame(smokingv1$idm, 
                      smokingv1$smoke_all - 1) # recode: 0 = never a smoker;
            # (1) smoked until pregnancy was known (i.e., first trimester only); 
            # (2) continued smoking during pregnancy.
colnames(smoking) <- c('idm', "maternalSmoking")

#-------------------------------------------------------------------------------
# MATERNAL ALCOHOL CONSUMPTION during pregnancy # not used
#-------------------------------------------------------------------------------
# SOCIAL SUPPORT during pregnancy 

# Recode items so higher scores reflect more social support 
fad <- data.frame(social$j0100103, 
                  social$j0300103,
                  social$j0500103,
                  social$j0700103,
                  social$j0900103,
                  social$j1100103,
                  5 - social$j0200103, # Recode inverse item
                  5 - social$j0400103, # Recode inverse item
                  5 - social$j0600103, # Recode inverse item
                  5 - social$j0800103, # Recode inverse item
                  5 - social$j1000103, # Recode inverse item
                  5 - social$j1200103) # Recode inverse item

social$familySS <- fad_scores(fad)

################################################################################
# Merge all variables together
depre_covs <- Reduce(function(x,y) merge(x = x, y = y, by = "idm",  all.x = TRUE),
                                   list(predep, child_general, smoking, social))

full_dataset = merge(depre_covs, haircort, by = "idc" , all.x = TRUE)

################################################################################
##----------------------------------------------------------------------------##
## -------------------- Exclude participants (flowchart) -------------------- ##
##----------------------------------------------------------------------------##

initial_sample <- 9901

## First exclusion step: Exclude children with missing cortisol
outcome_meas <- full_dataset[!is.na(full_dataset$cortisol),] 
after_cort <- nrow(outcome_meas)
sub1 <- after_cort - initial_sample

# Exclude participants with missing ethnicity values
ethn_present <- outcome_meas[!is.na(outcome_meas$childRaceEth),] 
after_ethn <- nrow(ethn_present)
sub2 <- after_ethn - after_cort

# Exclude twins 
no_twins <- ethn_present[ethn_present$twin == 0, ]
after_twins <- nrow(no_twins)
sub3 <- after_twins - after_ethn

# Select only one sibling (based on data availability or randomly).
# First, I determine a list of mothers that have more than one child in the set.
# NOTE: duplicated() is the best option I could find to determine which mother IDs
# recur more than once (very non-elegant, tbh, but using table() gets even uglier)
siblings_id = data.frame(no_twins$mother[duplicated(no_twins$mother)])
# duplicated() funtion does not allow to ignore NAs so I remove them manually.
# I also transform the numeric vector into a dataframe because of indexing problems.
siblings_id = siblings_id[!is.na(siblings_id),]; siblings_id = data.frame(siblings_id);
# Second, I create an empty vector to fill with the IDC of the sibling(s) with more 
# missing items or with a randomly picked sibling in case they have the same nr of missing. 
worse_sibling = rep(NA, dim(siblings_id)[1] + 1) # I will need the "+1" for triplets! 
# Loop through the mother IDs I previously identified and link them to IDCs
for (i in 1:dim(siblings_id)[1]) {
  siblings = no_twins[no_twins$mother == siblings_id[i,1], ] # identify the couples of siblings
  # For some reason when I run the line above 2 rows of NAs are created too, go figure. 
  # Let's get rid of them:
  siblings = siblings[rowSums(is.na(siblings)) != ncol(siblings), ]
  # There is one mother with 3 siblings, let's select the "best" one and get rid of the other two
  if (dim(siblings)[1] > 2) {
    nmiss = c( sum(is.na(siblings[1,])), sum(is.na(siblings[2,])), sum(is.na(siblings[3,])) )
    if (which.min(nmiss) == 1) { worse_sibling[i] = siblings[2,'idc']
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[3,'idc']
    } else if (which.min(nmiss) == 2) { worse_sibling[i] = siblings[1,'idc']
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[3,'idc']
    } else { worse_sibling[i] = siblings[1,'idc'] 
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[2,'idc'] }
  }
  # otherwise, select the "worse" sibling (with more missing) and add to it the the black list
  if ( sum(is.na(siblings[1,])) > sum(is.na(siblings[2,])) ) {
    worse_sibling[i] = siblings[1,'idc']
  } else if ( sum(is.na(siblings[1,])) == sum(is.na(siblings[2,])) ) {
    worse_sibling[i] = siblings[sample(1:2, 1),'idc']
  } else { worse_sibling[i] = siblings[2,'idc'] }
}

# Now we are finally ready to exclude siblings
final <- no_twins[no_twins$idc %notin% worse_sibling, ]
after_siblings <- nrow(final)
sub4 <- after_siblings - after_twins

# Flowchart
flowchart <- list(initial_sample, sub1, after_cort, sub2, after_ethn, sub3, after_twins, sub4, after_siblings)

cat(paste("Well, congrats! Your final dataset includes", after_siblings ,"participants.")) # 2877

# ------------------------------------------------------------------------------
# Save the raw dataset
saveRDS(final, paste(pathtodata,'dataset_raw.rds', sep = ""))
