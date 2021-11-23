# plotting and descriptives 

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }
################################################################################
# Load datasets
depcor <- readRDS(paste(pathtodata, 'dataset_raw.rds', sep = ""))

# white participants only 
white <- depcor[depcor$childRaceEth == 0, ]
################################################################################
#                             BASIC DESCRIPTIVES

N <- nrow(depcor) # total sample N
NODEP <- summary(as.factor(depcor$bsiScoreFlg))[1] # N controls
DEP <- summary(as.factor(depcor$bsiScoreFlg))[2] # N cases

wN <- nrow(white) # total sample N
wNODEP <- summary(as.factor(white$bsiScoreFlg))[1] # N controls
wDEP <- summary(as.factor(white$bsiScoreFlg))[2] # N cases

cat(paste("FULL\nOverall:   ", N, 
          "\nBSI < 0.75:", NODEP, 
          "\nBSI > 0.75:", DEP ))

cat(paste("WHITE\nOverall:   ", wN, 
          "\nBSI < 0.75:", wNODEP, 
          "\nBSI > 0.75:", wDEP ))
# ------------------------------------------------------------------------------
# compute mean and standard deviation of entire sample and in cases vs controls
meansd <- function(var, splitvar = depcor$bsiScoreFlg) {
  totm <- format(round(mean(var, na.rm = T), 1), nsmall = 1)
  tots <- format(round(sd(var, na.rm = T), 2), nsmall = 2)
  splitm <- format(round(tapply(var, splitvar, mean, na.rm = T), 1), nsmall = 1)
  splits <- format(round(tapply(var, splitvar, sd, na.rm = T), 2), nsmall = 2)
  
  cat(paste("Overall:    ", totm, " (± ", tots, ")",
            "\nBSI < 0.75: ", splitm[1], " (± ", splits[1], ")",
            "\nBSI > 0.75: ", splitm[2], " (± ", splits[2], ")", sep = ""))
}
# ------------------------------------------------------------------------------
# compute n and percent of entire sample and in cases vs controls for binary variables
npercent <- function(var, splitvar = depcor$bsiScoreFlg) {
  ns <- summary(as.factor(var)); 
  prc0 = round((ns[1]/N)*100, 1); prc1 = round((ns[2]/N)*100, 1); 
  
  spt <- split(var, splitvar)
  depns <- summary(as.factor(spt$`1`))
  depprc0 = round((depns[1]/DEP)*100, 1); depprc1 = round((depns[2]/DEP)*100, 1); 
  nodns <- summary(as.factor(spt$`0`))
  nodepprc0 = round((nodns[1]/NODEP)*100, 1); nodepprc1 = round((nodns[2]/NODEP)*100, 1); 
  
  cat(paste("Overall:    ", ns[1], " (", prc0, "%)  ",  ns[2], " (", prc1, "%)", 
            "\nBSI < 0.75: ", depns[1], " (", depprc0, "%)  ",  depns[2], " (", depprc1, "%)",
            "\nBSI > 0.75: ", nodns[1], " (", nodepprc0, "%)  ",  nodns[2], " (", nodepprc1, "%)",
            sep = ""))
}
# ------------------------------------------------------------------------------
# compute n and percent of entire sample and in cases vs controls for variables with 3 lvls
npercent3 <- function(var, splitvar = depcor$bsiScoreFlg) {
  ns <- summary(as.factor(var)); 
  prc0 = round((ns[1]/N)*100, 1); prc1 = round((ns[2]/N)*100, 1); prc2 = round((ns[3]/N)*100, 1); 
  
  spt <- split(var, splitvar)
  depns <- summary(as.factor(spt$`1`))
  depprc0 = round((depns[1]/DEP)*100, 1); depprc1 = round((depns[2]/DEP)*100, 1); depprc2 = round((depns[3]/DEP)*100, 1); 
  nodns <- summary(as.factor(spt$`0`))
  nodepprc0 = round((nodns[1]/NODEP)*100, 1); nodepprc1 = round((nodns[2]/NODEP)*100, 1); nodepprc2 = round((nodns[3]/NODEP)*100, 1); 
  
  cat(paste("Overall:    ", ns[1], " (", prc0, "%)  ",  ns[2], " (", prc1, "%)  ", ns[3], " (", prc2, "%)",
            "\nBSI < 0.75: ", depns[1], " (", depprc0, "%)  ",  depns[2], " (", depprc1, "%)  ",  depns[3], " (", depprc2, "%)",
            "\nBSI > 0.75: ", nodns[1], " (", nodepprc0, "%)  ",  nodns[2], " (", nodepprc1, "%)  ", nodns[3], " (", nodepprc2, "%)",
            sep = ""))
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

meansd(depcor$maternalAge)
meansd(depcor$prePregBMI)
npercent(depcor$maternalEdu)
npercent3(depcor$maternalSmoking)
npercent(depcor$maternalMaritalStatus)
npercent(depcor$maternalIncome)
npercent(depcor$sex)
npercent(depcor$childRaceEth)
meansd(depcor$gestAge)
meansd(depcor$familySS)
meansd(depcor$childAge)
meansd(depcor$cortisol)

N <- wN
DEP <- wDEP
NODEP <- wNODEP
meansd(white$maternalAge, splitvar = white$bsiScoreFlg); # white only
meansd(white$prePregBMI, splitvar = white$bsiScoreFlg); # white only
npercent(white$maternalEdu, splitvar = white$bsiScoreFlg); # white only
npercent3(white$maternalSmoking, splitvar = white$bsiScoreFlg); # white only
npercent(white$maternalMaritalStatus, splitvar = white$bsiScoreFlg); # white only
npercent(white$maternalIncome, splitvar = white$bsiScoreFlg); # white only
npercent(white$sex, splitvar = white$bsiScoreFlg); # white only
meansd(white$gestAge, splitvar = white$bsiScoreFlg); # white only
meansd(white$familySS, splitvar = white$bsiScoreFlg); # white only
meansd(white$childAge, splitvar = white$bsiScoreFlg); # white only
meansd(white$lncortisol, splitvar = white$bsiScoreFlg); # white only
################################################################################
################################################################################
################################################################################

library(ggplot2)
library(cowplot)
library(raincloudplots)


# width and height variables for saved plots
w = 6
h = 3


#Rainclouds with boxplots
p6 <- ggplot(depcor, aes(x=bsiScoreFlg, y=lncortisol, fill = bsiScoreFlg, colour = bsiScoreFlg))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  geom_boxplot(aes(x = as.numeric(bsiScoreFlg)+0.25, y = lncortisol),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  
  ylab('ln(cortisol)')+xlab('Maternal Depression')+coord_flip()+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Figure")
ggsave('boxplots.png', width = w, height = h)

p6