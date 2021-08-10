# SKRIPT for perceptual/conceptual analysis with VRAT and CRAT data


rm(list= ls())

## load libraries ##############################################################
library(lme4)
library(lmerTest)
library(corrplot)
library(sjPlot)
library(sjmisc)
library(knitr)
library(magrittr)
library(sjlabelled)      
library(sjmisc)                                                                                    
library(sjstats) 
library(ggeffects)
library(performance)
library(parameters)
library(betareg)
library(tidyverse)
library(glmmTMB)
library(MKinfer)

#load datafile
data <- read.table("BeckerDavisCabeza_2021_data.csv", sep = ",", dec = ".", header=TRUE, na.strings=c("", " " , "NA", "NAN" )) 

data$AHA_fac = as.factor(data$AHA_fac)

########################################################################### ##

### Descriptives C-RAT & LI-RAT ##############################################

b <- data %>%
  count( ID,sex)
sum(b$sex)

mean(data[data$sample == "english" & data$sex == 0,]$age, na.rm =T)
mean(data[data$sample == "english"& data$sex == 1,]$age, na.rm =T)

min(data[data$sample == "english",]$age, na.rm =T)
max(data[data$sample == "english",]$age, na.rm =T)


### C-RAT ### ---------------------------------------------------------------------------- -- 
mean(data[data$task == "CRA" & data$sample == "english",]$cor2, na.rm = T)
sd(data[data$task == "CRA" & data$sample == "english",]$cor2, na.rm = T)

mean(data[data$task == "CRA" & data$sample == "english",]$RT, na.rm = T)
sd(data[data$task == "CRA" & data$sample == "english",]$RT, na.rm = T)

mean(data[data$task == "CRA" & data$sample == "english",]$AHA_num, na.rm = T)
sd(data[data$task == "CRA" & data$sample == "english",]$AHA_num, na.rm = T)

mean(data[data$task == "CRA" & data$sample == "english",]$surprise, na.rm = T) # surprise = suddenness
sd(data[data$task == "CRA" & data$sample == "english",]$surprise, na.rm = T)

### LI-RAT ----------------------------------------------------------------------- --
mean(data[data$task == "VRAT" & data$sample == "english",]$cor2, na.rm = T)
sd(data[data$task == "VRAT" & data$sample == "english",]$cor2, na.rm = T)

mean(data[data$task == "VRAT" & data$sample == "english",]$RT, na.rm = T)
sd(data[data$task == "VRAT" & data$sample == "english",]$RT, na.rm = T)

mean(data[data$task == "VRAT" & data$sample == "english",]$AHA_num, na.rm = T)
sd(data[data$task == "VRAT" & data$sample == "english",]$AHA_num, na.rm = T)

mean(data[data$task == "VRAT" & data$sample == "english",]$surprise, na.rm = T) # surprise = suddenness
sd(data[data$task == "VRAT" & data$sample == "english",]$surprise, na.rm = T)


#### General Linear Mixed Models #####################################
### C-RAT: cor2  #######################################################

### cor2 - C-RAT
enCRA_acc <- glmer(cor2~ scale(T123Solve_lch)+scale(T123_lch) + (1|ID),data= data,family = binomial(link ="logit"), na.action  = na.omit)
enCRA_acc_IA <-glmer(cor2~ scale(T123Solve_lch)*scale(T123_lch)+(1|ID),data= data,family = binomial(link ="logit"), na.action  = na.omit)
#plot(check_distribution(enCRA_acc))
anova(enCRA_acc, enCRA_acc_IA)
summary(enCRA_acc_IA)
CRAT_cor2 <- ggpredict(enCRA_acc_IA , c( 'T123Solve_lch', 'T123_lch')) %>% plot() + ggplot2::theme_classic()
#confint(enCRA_acc_IA)

### RT - C-RAT
enCRA_RT <- lmer(log(RT)~ scale(T123Solve_lch)+scale(T123_lch) + (1|ID),data= data, na.action  = na.omit)
enCRA_RT_IA <-lmer(log(RT)~ scale(T123Solve_lch)*scale(T123_lch)+(1|ID),data= data, na.action  = na.omit)
#plot(check_distribution(enCRA_RT))
anova(enCRA_RT, enCRA_RT_IA)
summary(enCRA_RT_IA)
CRAT_RT <-ggpredict(enCRA_RT_IA , c(  'T123Solve_lch','T123_lch')) %>% plot() + ggplot2::theme_classic()

### aha (pleasure) - C-RAT
enCRA_aha <- glmer(AHA_fac~ cor2 + scale(T123Solve_lch)+scale(T123_lch) + (1|ID),data= data,family = binomial(link ="logit"), na.action  = na.omit)
enCRA_aha_IA <-glmer(AHA_fac~ cor2 + scale(T123Solve_lch)*scale(T123_lch)+(1|ID),data= data,family = binomial(link ="logit"), na.action  = na.omit)
#plot(check_distribution(enCRA_aha_IA))
anova(enCRA_aha, enCRA_aha_IA)
summary(enCRA_aha)
CRAT_aha <-ggpredict(enCRA_aha_IA , c( 'T123Solve_lch', 'T123_lch')) %>% plot() + ggplot2::theme_classic()

### aha (suddenness) - C-RAT
enCRA_sud<- lmer(suddenness~ cor2 + scale(T123Solve_lch)+scale(T123_lch) + (1|ID),data= data, na.action  = na.omit)
enCRA_sud_IA <-lmer(suddenness~ cor2 + scale(T123Solve_lch)*scale(T123_lch)+(1|ID),data= data, na.action  = na.omit)
#plot(check_distribution(enCRA_sud))
anova(enCRA_sud, enCRA_sud_IA)
summary(enCRA_sud_IA)
CRAT_sud <- ggpredict(enCRA_sud_IA , c( 'T123Solve_lch','T123_lch')) %>% plot() + ggplot2::theme_classic()

############################################################################## ##
### LI-RAT  ################################################################

### accuracy
#data(cbpp, package="lme4")
LIRAT_cor2_WN1 <- glmmTMB(cor2 ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) +scale(PC_LCH) + scale(PS_LCH) +(1|ID),
                   family=binomial, data=data)

LIRAT_cor2_WN2 <- glmmTMB(cor2 ~ scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +(1|ID), # +scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR)  
                           family=binomial(link="logit"), data=data)

LIRAT_cor2_WN3 <- glmmTMB(cor2 ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) +scale(PS_LCH) +(1|ID),
                           family=binomial, data=data)

anova(LIRAT_cor2_WN1,LIRAT_cor2_WN2, LIRAT_cor2_WN3)
#plot(check_distribution(LIRAT_cor2_WN1))
LIRAT_WN_cor2 <- ggpredict(LIRAT_cor2_WN2, c("CS_LCH", "PC_LCH")) %>% plot + ggplot2::theme_classic()
summary(LIRAT_cor2_WN2)

### RT - LiRAT-LCH
LIRAT_RT_WN1 <- lmer(log(RT) ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) +scale(PC_LCH) + scale(PS_LCH) +(1|ID),data= data,na.action  = na.omit)
LIRAT_RT_WN2 <- lmer(log(RT) ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +(1|ID),data= data,na.action  = na.omit)
LIRAT_RT_WN3 <- lmer(log(RT) ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH)+(1|ID),data= data,na.action  = na.omit)
anova(LIRAT_RT_WN1,LIRAT_RT_WN2,LIRAT_RT_WN3 )
ggpredict(LIRAT_RT_WN3, c("PS_Vis_HR","PC_Vis_HR")) %>%plot + ggplot2::theme_classic()
summary(LIRAT_RT_WN1)

### AHA - pleasure - LiRAT-LCH
LIRAT_aha_WN1 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) +scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data, family = binomial(link = "logit" ),na.action  = na.omit)
LIRAT_aha_WN2 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data, family = binomial(link = "logit" ),na.action  = na.omit)
LIRAT_aha_WN3 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data, family = binomial(link = "logit" ),na.action  = na.omit)
anova(LIRAT_aha_WN1,LIRAT_aha_WN2, LIRAT_aha_WN3)
summary(LIRAT_aha_WN1)

### AHA - Suddenness - LiRAT-LCH
LIRAT_sud_WN1 <- lmer(suddenness ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) +scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data,na.action  = na.omit)
LIRAT_sud_WN2 <- lmer(suddenness ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data,na.action  = na.omit)
LIRAT_sud_WN3 <- lmer(suddenness ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_LCH) *scale(PC_LCH) + scale(PS_LCH) +cor2+(1|ID),data= data,na.action  = na.omit)
anova(LIRAT_sud_WN1,LIRAT_sud_WN2, LIRAT_sud_WN3)
summary(LIRAT_sud_WN1)

###########################################################################################################################################################
############ LIRAT with cosine similarity

### cor2 -  LIRAT - Cosine
LIRAT_cor2_CD1 <- glmmTMB(cor2 ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) +scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),
                           family=binomial, data=data)
LIRAT_cor2_CD2 <- glmmTMB(cor2 ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),
                           family=binomial, data=data)
LIRAT_cor2_CD3 <- glmmTMB(cor2 ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),
                           family=binomial, data=data)
#plot(check_distribution(LIRAT_cor2_CD2))
anova(LIRAT_cor2_CD1,LIRAT_cor2_CD2,LIRAT_cor2_CD3)
summary(LIRAT_cor2_CD2)
LIRAT_CD_cor2 <- ggpredict(LIRAT_cor2_CD2, c("CS_CD_conc","PC_CD_conc")) %>%plot + ggplot2::theme_classic()

### RT - LIRAT - Cosine
LIRAT_RT_CD1  <- lmer(log(RT) ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) +scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
LIRAT_RT_CD2  <- lmer(log(RT) ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
LIRAT_RT_CD3  <- lmer(log(RT) ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
#plot(check_distribution(LIRAT_RT_CD3))
anova(LIRAT_RT_CD1,LIRAT_RT_CD2,LIRAT_RT_CD3)
summary(LIRAT_RT_CD1)


### AHA (pleasure) - LIRAT - Cosine
LIRAT_aha_CD1 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) +  scale(CS_CD_conc) +scale(PC_CD_conc) + scale(PS_CD_conc) + cor2 + (1|ID),
                           family=binomial, data=data)
LIRAT_aha_CD2 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) +  scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +cor2 +(1|ID),
                           family=binomial, data=data)
LIRAT_aha_CD3 <- glmmTMB(AHA_fac ~  scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) +  scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +cor2 +(1|ID),
                           family=binomial, data=data)
#plot(check_distribution(LIRAT_aha_CD3))
anova(LIRAT_aha_CD1,LIRAT_aha_CD2,LIRAT_aha_CD3)
LIRAT_CD_aha <-ggpredict(LIRAT_aha_CD2, c("CS_CD_conc","PC_CD_conc")) %>%plot + ggplot2::theme_classic()
summary(LIRAT_aha_CD3)


### Suddenness - LIRAT Cosine
LIRAT_sud_CD1  <- lmer(suddenness ~  cor2 + scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) +scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
LIRAT_sud_CD2  <- lmer(suddenness ~  cor2 + scale(PS_Vis_HR)+scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
LIRAT_sud_CD3  <- lmer(suddenness ~  cor2 + scale(PS_Vis_HR)*scale(PC_Vis_HR)+scale(CS_Vis_HR) + scale(CS_CD_conc) *scale(PC_CD_conc) + scale(PS_CD_conc) +(1|ID),data= data,na.action  = na.omit)
#plot(check_distribution(LIRAT_sud_CD3))
anova( LIRAT_sud_CD1,LIRAT_sud_CD2,LIRAT_sud_CD3)
summary(LIRAT_sud_CD2)
LIRAT_CD_sud <-ggpredict(LIRAT_sud_CD2, c("CS_CD_conc","PC_CD_conc")) %>%plot + ggplot2::theme_classic()

############ make all Tables ###################################################

# C-RAT
tab_model(enCRA_acc_IA, enCRA_RT_IA )
tab_model(enCRA_aha, enCRA_sud_IA)

# LI-RAT (path length)
tab_model(LIRAT_cor2_WN2, LIRAT_RT_WN1) 
tab_model(LIRAT_aha_WN1, LIRAT_sud_WN1)

# LI-RAT (cosine similarity)
tab_model(LIRAT_cor2_CD2, LIRAT_RT_CD1)
tab_model(LIRAT_aha_CD3, LIRAT_sud_CD2)

######### Plot results ######################################################################
library(ggpubr)
fig3_acc <- ggarrange(CRAT_cor2, LIRAT_WN_cor2, LIRAT_CD_cor2,
                      common.legend = FALSE, legend = "bottom",
                      #labels = c("A", "B", "C", ),
                      ncol = 3, nrow = 1)

###### visual rating - sanity check #########
#load('VRAT_relational_properties.Rdata')
#perm.t.test(allSim$VisHu_PS, allSim$VisHu_CS, paired = T )

