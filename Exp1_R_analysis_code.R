rm(list=ls())
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyr)
library(ggplot2)
library(grid)
library(multcomp)


#rm(list=ls())
source("/Users/lewisforder/Dropbox/SharedFolders/LupyanExps/Psychopy\ Files/helpers.R")

#load data:
x1 <- read.table("~/Dropbox/SharedFolders/LupyanExps/LCD/analysis/Exp1_unblocked/Exp1_unblocked_Task3_allData.txt", sep="\t", header = TRUE)

#load euclidean distances (focal target to transitionary target and boundary target):
#(note that raw euclid dist values for prototype targets will of course be zero)
x2 <- read.table("~/Dropbox/SharedFolders/LupyanExps/LCD/analysis/Exp1_unblocked/LCD_Exp1_EucDists_Protos2Targets.csv", sep="\t", header = TRUE)

#concatenate
dataUnblocked <- cbind(x1,x2)

#tidy
rm(x1, x2)

#check for missing values
which(!complete.cases(dataUnblocked))

#toss low-performing subjects
dataUnblocked_agg = dataUnblocked %>% group_by(subjCode) %>% summarize(accuracy=mean(accuracy),RT=median(RT),numTrials=n()) %>% transform(normedAcc=scale(accuracy))

lowAccSubjs_unblocked = dataUnblocked_agg %>% filter(accuracy<.50)
# lowAccSubjs_unblocked = dataUnblocked_agg %>% filter(normedAcc< -2.5 | normedAcc>2.5)
dataUnblocked = dataUnblocked %>% filter(!subjCode %in% lowAccSubjs_unblocked$subjCode) #remove these subjects

#fix up the latency column
dataUnblocked = dataUnblocked %>% transform(latency = RT*accuracy*1000) %>% transform(latency=ifelse((accuracy=0 | latency<200 | latency>4000),NA,latency))

#split condition code into more informative columns. This is extra work due to inconvenient condition coding...
dataUnblocked = dataUnblocked %>% separate(condition_Code,c("condition_num2", "conditionSuffix"), sep='\\.') %>% separate(conditionSuffix,c("color", "boundary"),sep="2")
dataUnblocked$color = as.factor(dataUnblocked$color)
dataUnblocked$boundary = as.factor(dataUnblocked$boundary)

#create column to denote basic (Red-Orange / Yellow-Orange / Red-Blue / Blue-Purple)
# vs. non-basic conditions (Yellow-Lime / Green-Lime / Green-Teal / Blue-Teal)
dataUnblocked$isBasic <- ifelse(dataUnblocked$boundary %in% c("RYboundary", "BRboundary"), 1, 0)

#create column to denote if target was left/right
# vs. non-basic conditions (Yellow-Lime / Green-Lime / Green-Teal / Blue-Teal)
dataUnblocked$isLeft <- ifelse(dataUnblocked$position0TL3BR %in% c(0, 2), 1, 0)

#center variables
dataUnblocked$isLabel_c = myCenter(dataUnblocked$word1_vs_noise0)
dataUnblocked$comparisonType_c = myCenter(dataUnblocked$stim_pair_0focal1trans2bound)
dataUnblocked$color_c = myCenter(dataUnblocked$color)
dataUnblocked$eucDists_c = myCenter(dataUnblocked$EucDists)
dataUnblocked$isBasic_c = myCenter(dataUnblocked$isBasic)
dataUnblocked$isLeft_c = myCenter(dataUnblocked$isLeft)


######################
# Get descriptives
######################
dataUnblocked %>% group_by(word1_vs_noise0) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(stim_pair_0focal1trans2bound) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(stim_pair_0focal1trans2bound,word1_vs_noise0) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(color,stim_pair_0focal1trans2bound,word1_vs_noise0) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(color,stim_pair_0focal1trans2bound,word1_vs_noise0,boundary) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame

dataUnblocked %>% group_by(color) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(color) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame
dataUnblocked %>% group_by(isBasic,stim_pair_0focal1trans2bound,word1_vs_noise0) %>% summarize(isRight=mean(accuracy),latency=mean(latency,na.rm=TRUE)) %>% data.frame

######################
# Check and visualize data:
######################
#check to make sure things are ok
summary(dataUnblocked)
dataUnblocked_agg = dataUnblocked %>% group_by(subjCode) %>% summarize(accuracy=mean(accuracy),latency=mean(latency,na.rm=TRUE),numTrials=n()) %>% transform(normedLatency = scale(latency), normedAcc=scale(accuracy))
ggplot(dataUnblocked_agg, aes(x=latency)) + geom_dotplot() + ggtitle("1a Latency histogram")
ggplot(dataUnblocked_agg, aes(x=accuracy)) + geom_dotplot() + ggtitle("1b Accuracy histogram")


#let's visualize accuracy by condition
dataUnblocked_for_graphing = summarySEwithin(data=dataUnblocked,measurevar="accuracy",withinvars=c("color","stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=dataUnblocked_for_graphing, aes(y=accuracy,group=stim_pair_0focal1trans2bound,x=word1_vs_noise0,color=stim_pair_0focal1trans2bound))+
  facet_grid(.~color)+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.25)+
  myThemeBasic+
  ggtitle("1c Accuracy by condition")

#collapsing accuracy across colors
dataUnblocked_for_graphing_agg = summarySEwithin(data=dataUnblocked,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=dataUnblocked_for_graphing_agg, aes(y=accuracy,group=stim_pair_0focal1trans2bound,x=word1_vs_noise0,color=stim_pair_0focal1trans2bound))+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.25)+
  myThemeBasic+
  ggtitle("1d Accuracy collapsed across colors")

#now visualize latency
dataUnblocked_for_graphing_agg_latency = summarySEwithin(data=dataUnblocked,measurevar="latency",na.rm=TRUE,withinvars=c("color","stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=dataUnblocked_for_graphing_agg_latency, aes(y=latency,group=stim_pair_0focal1trans2bound,x=word1_vs_noise0,color=stim_pair_0focal1trans2bound))+
  facet_grid(.~color)+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymax=latency+se,ymin=latency-se),width=0.25)+
  myThemeBasic+
  ggtitle("1e Latency by condition")

#now visualize collapsed latency
dataUnblocked_for_graphing_agg_latency = summarySEwithin(data=dataUnblocked,measurevar="latency",na.rm=TRUE,withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=dataUnblocked_for_graphing_agg_latency, aes(y=latency,group=stim_pair_0focal1trans2bound,x=word1_vs_noise0,color=stim_pair_0focal1trans2bound))+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymax=latency+se,ymin=latency-se),width=0.25)+
  myThemeBasic+
  ggtitle("1f Latency collapsed across colors")


#now look at whether 'basic' colors produce different performance vs. non-basic:
dataUnblocked_for_graphing = summarySEwithin(data=dataUnblocked,measurevar="accuracy",withinvars=c("isBasic_c","stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=dataUnblocked_for_graphing, aes(y=accuracy,group=stim_pair_0focal1trans2bound,x=word1_vs_noise0,color=stim_pair_0focal1trans2bound))+
  facet_grid(.~isBasic_c)+
  geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.25)+
  myThemeBasic+
  ggtitle("Accuracy by condition for basic vs non-basic")


#########################
# GL model building: Accuracy
#########################
a=glmer(accuracy~1+(1|subjCode),data=dataUnblocked,family=binomial)
b=glmer(accuracy~isLabel_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
c=glmer(accuracy~isLabel_c+comparisonType_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
d=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
e=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,family=binomial)
# show the full model with coefficients / Z scores / p values/ CIs
anova(a,b,c,d,e)
summary(e)
e %>% confint(method="Wald")

# Next try adding color category as a fixed instead of random effect
f=glmer(accuracy~isLabel_c*comparisonType_c+color_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
g=glmer(accuracy~isLabel_c*comparisonType_c+comparisonType_c:color_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
h=glmer(accuracy~isLabel_c*comparisonType_c*comparisonType_c*color_c+(isLabel_c|subjCode),data=dataUnblocked,family=binomial)
summary(f)
summary(g)
summary(h)
h %>% confint(method="Wald")
anova(a,b,c,d,e,f)
anova(a,b,c,d,e,g)
anova(a,b,c,d,e,h)
# Moving color category from a random to a fixed effect reduces model accuracy


# Model the effect of cues separately for each of the 4 color categories
red_data    <- filter(dataUnblocked, condition_Num == 0 | condition_Num == 7)
yellow_data <- filter(dataUnblocked, condition_Num == 1 | condition_Num == 2)
green_data  <- filter(dataUnblocked, condition_Num == 3 | condition_Num == 4)
blue_data   <- filter(dataUnblocked, condition_Num == 5 | condition_Num == 6)

e_red=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=red_data,family=binomial)
summary(e_red)
e_red %>% confint(method="Wald")

e_yellow=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=yellow_data,family=binomial)
summary(e_yellow)
e_yellow %>% confint(method="Wald")

e_green=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=green_data,family=binomial)
summary(e_green)
e_green %>% confint(method="Wald")

e_blue=glmer(accuracy~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=blue_data,family=binomial)
summary(e_blue)
e_blue %>% confint(method="Wald")


# Next try adding color basicness to the model
h=glmer(accuracy~isLabel_c*comparisonType_c+isBasic_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,family=binomial)
i=glmer(accuracy~isLabel_c*comparisonType_c*isBasic_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,family=binomial)
anova(a,b,c,d,e,h)
anova(a,b,c,d,e,h,i)
summary(i)

j=glmer(accuracy~isLabel_c*comparisonType_c+isLabel_c:isBasic_c+comparisonType_c:isBasic_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,family=binomial)
summary(j)

k=glmer(accuracy~isLabel_c*comparisonType_c+comparisonType_c:isBasic_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,family=binomial)
# This (k) is the best-fitting model that includes basicness
summary(k)

anova(e,k) # including basicness is a significantly better fit (but the trend is ultimately the same for basic / non-basic colors after hearing a label vs. nosise).

#########################
#analyze the isLabel_c:comparisonType_c interaction:
# Label trials:
data_labels <- filter(dataUnblocked, word1_vs_noise0 ==1)
A = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels,family=binomial)
summary(A)
A %>% confint(method="Wald")

#labels proto vs trans:
data_labels_protoVsTrans <- subset(dataUnblocked, word1_vs_noise0 ==1 & stim_pair_0focal1trans2bound < 2)
summary(data_labels_protoVsTrans)
B = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_protoVsTrans,family=binomial)
summary(B)
B %>% confint(method="Wald")


#labels proto vs bound:
data_labels_protoVsBound <- filter(dataUnblocked, word1_vs_noise0 ==1)
data_labels_protoVsBound <- filter(data_labels_protoVsBound,stim_pair_0focal1trans2bound == 0 | stim_pair_0focal1trans2bound == 2)
summary(data_labels_protoVsBound)
C = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_protoVsBound,family=binomial)
summary(C)
C %>% confint(method="Wald")

#labels trans vs bound:
data_labels_TransVsBound <- filter(dataUnblocked, word1_vs_noise0 ==1)
data_labels_TransVsBound <- filter(data_labels_TransVsBound,stim_pair_0focal1trans2bound > 0)
summary(data_labels_TransVsBound)
D = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_TransVsBound,family=binomial)
summary(D)
D %>% confint(method="Wald")


#######################
# Noise trials:
data_labels <- filter(dataUnblocked, word1_vs_noise0 ==0)
E = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels,family=binomial)
summary(E)
E %>% confint(method="Wald")

#noise proto vs trans:
data_labels_protoVsTrans <- subset(dataUnblocked, word1_vs_noise0 ==0 & stim_pair_0focal1trans2bound < 2)
summary(data_labels_protoVsTrans)
G = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_protoVsTrans,family=binomial)
summary(G)
G %>% confint(method="Wald")


#noise proto vs bound:
data_labels_protoVsBound <- filter(dataUnblocked, word1_vs_noise0 ==0)
data_labels_protoVsBound <- filter(data_labels_protoVsBound,stim_pair_0focal1trans2bound == 0 | stim_pair_0focal1trans2bound == 2)
summary(data_labels_protoVsBound)
H = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_protoVsBound,family=binomial)
summary(H)
H %>% confint(method="Wald")

#noise trans vs bound:
data_labels_TransVsBound <- filter(dataUnblocked, word1_vs_noise0 ==0)
data_labels_TransVsBound <- filter(data_labels_TransVsBound,stim_pair_0focal1trans2bound > 0)
summary(data_labels_TransVsBound)
I = glmer(accuracy~comparisonType_c+(1|subjCode)+(1|color),data=data_labels_TransVsBound,family=binomial)
summary(I)
I %>% confint(method="Wald")


# Prototypes:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==0) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==0) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% confint(method="Wald")
# Transitionals:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==1) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==1) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% confint(method="Wald")
# Boundaries:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==2) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==2) %>% glmer(accuracy~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,family=binomial) %>% confint(method="Wald")


######################
# Model building: LATENCY:
######################

La=lmer(latency~1+(1|subjCode),data=dataUnblocked,REML=FALSE)
Lb=lmer(latency~isLabel_c+(isLabel_c|subjCode),data=dataUnblocked,REML=FALSE)
Lc=lmer(latency~isLabel_c+comparisonType_c+(isLabel_c|subjCode),data=dataUnblocked,REML=FALSE)
Ld=lmer(latency~isLabel_c*comparisonType_c+(isLabel_c|subjCode),data=dataUnblocked,REML=FALSE)
Le=lmer(latency~isLabel_c*comparisonType_c+(isLabel_c|subjCode)+(1|color),data=dataUnblocked,REML=FALSE)
# Le=lmer(latency~isLabel_c*comparisonType_c+color_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,REML=FALSE)
# Lf=lmer(latency~isLabel_c*comparisonType_c+isLabel_c:color_c+comparisonType_c:color_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,REML=FALSE)
# Lg=lmer(latency~isLabel_c*comparisonType_c*color_c+(isLabel_c|subjCode)+(1|color_c),data=dataUnblocked,REML=FALSE)
anova(La,Lb,Lc,Ld,Le)


# show the full model with coefficients / Z scores / p values/ CIs
summary(Le)
Le %>% confint(method="Wald")

#########################
#analyze the isLabel_c:comparison_c interaction:
# Prototypes:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==0) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==0) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% confint(method="Wald")
# Transitionals:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==1) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==1) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% confint(method="Wald")
# Boundaries:
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==2) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% summary
dataUnblocked %>% filter(stim_pair_0focal1trans2bound ==2) %>% lmer(latency~isLabel_c+(isLabel_c|subjCode)+(1|color),data=.,REML=FALSE) %>% confint(method="Wald")


#trying to figure post hoc tests (¿glht?):
dataUnblocked_red$condString <- as.factor(dataUnblocked_red$condString)
dataUnblocked_red$condString2 <- NA
dataUnblocked_red$condString2 <- dataUnblocked_red$stim_pair_0focal1trans2bound
dataUnblocked_red$condString2 <- as.factor(dataUnblocked_red$condString2)
Lf.z = lmer(latency~condString2+(isLabel_c|subjCode),data=dataUnblocked_red,REML=FALSE)
summary(Lf.z) # 
summary(glht(Lf.z, linfct = mcp(condString = "Tukey")), test = adjusted("holm"))


lsmeans(f, pairwise ~ isLabel_c | comparisonType_c)

#these should produce same output:
summary(glht(myModel, lsm(pairwise ~ factorC)))
summary(glht(myModel, mcp(factorC="Tukey")))

######################
# Make final plots for the MS:
######################
accuracy_for_graphing_agg = summarySEwithin(data=dataUnblocked,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=accuracy_for_graphing_agg, aes(y=accuracy,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(0.5, 1))+
  annotate("text", x = 2.5, y = 0.94, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 0.75, label = "no cue", color = "gray50")+
  ylab("Accuracy")+
  theme(axis.text = element_text(colour = "black"))
  ggsave("LCD_Exp1_accuracy.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)


latency_for_graphing_agg = summarySEwithin(data=dataUnblocked,measurevar="latency",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode", na.rm=TRUE)
ggplot(data=latency_for_graphing_agg, aes(y=latency,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=latency+se,ymin=latency-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(500, 1750))+
  ylab("Reaction Time (ms)")+
  theme(axis.text = element_text(colour = "black"))+
  annotate("text", x = 2.5, y = 1020, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 1450, label = "no cue", color = "gray50")
ggsave("LCD_Exp1_latency.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)



accuracy_for_graphing_agg = summarySEwithin(data=red_data,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=accuracy_for_graphing_agg, aes(y=accuracy,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(0.5, 1))+
  annotate("text", x = 2.5, y = 1.0, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 0.75, label = "no cue", color = "gray50")+
  ylab("Accuracy")+
  theme(axis.text = element_text(colour = "black"))
ggsave("LCD_Exp1_accuracy_red.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)



accuracy_for_graphing_agg = summarySEwithin(data=yellow_data,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=accuracy_for_graphing_agg, aes(y=accuracy,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(0.5, 1))+
  annotate("text", x = 2.5, y = 1.0, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 0.75, label = "no cue", color = "gray50")+
  ylab("Accuracy")+
  theme(axis.text = element_text(colour = "black"))
ggsave("LCD_Exp1_accuracy_yellow.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)

accuracy_for_graphing_agg = summarySEwithin(data=green_data,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=accuracy_for_graphing_agg, aes(y=accuracy,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(0.5, 1))+
  annotate("text", x = 2.5, y = 1.0, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 0.7, label = "no cue", color = "gray50")+
  ylab("Accuracy")+
  theme(axis.text = element_text(colour = "black"))
ggsave("LCD_Exp1_accuracy_green.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)


accuracy_for_graphing_agg = summarySEwithin(data=blue_data,measurevar="accuracy",withinvars=c("stim_pair_0focal1trans2bound","word1_vs_noise0"),idvar="subjCode")
ggplot(data=accuracy_for_graphing_agg, aes(y=accuracy,group=word1_vs_noise0,x=stim_pair_0focal1trans2bound,color=word1_vs_noise0))+  geom_line()+
  geom_point(size=1)+
  geom_errorbar(aes(ymax=accuracy+se,ymin=accuracy-se),width=0.10)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), panel.border=element_blank())+
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=10), 
        axis.title.y=element_text(size=10))+
  theme(axis.title.y=element_text(vjust=0.9,size=10))+
  theme(axis.title.x=element_text(vjust=0.9,size=10))+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(axis.title.x=element_text(vjust=0))+
  theme(axis.title.x=element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(legend.position="none")+
  scale_x_discrete(labels=c("0" = "Prototype", "1" = "Transitional",
                            "2" = "Boundary"))+
  scale_colour_manual(values = c("firebrick2", "dodgerblue4"))+
  coord_cartesian(ylim = c(0.5, 1))+
  annotate("text", x = 2.5, y = 1.0, label = "verbal color cue", color = "gray50")+
  annotate("text", x = 2.5, y = 0.75, label = "no cue", color = "gray50")+
  ylab("Accuracy")+
  theme(axis.text = element_text(colour = "black"))
ggsave("LCD_Exp1_accuracy_blue.tiff", plot = last_plot(), device = "tiff", path = '/Users/lewisforder/Desktop/', scale = 1, width = 3, height = 2.34, units = "in", dpi = 300, limitsize = TRUE)
