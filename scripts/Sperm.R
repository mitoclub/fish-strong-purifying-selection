rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(xlsx)
library(vegan)


path_to_sperm = '../data/SpermAnalysis/опыт 30.06.2023/'

sperm_cont = read.xlsx(paste(path_to_sperm, 'Контроль.xlsx', sep = ''), sheetIndex=1, startRow=6)
sperm_mut = read.xlsx(paste(path_to_sperm, 'Мутантный.xlsx', sep = ''), sheetIndex=1, startRow=6)


sperm_cont$Condition = 'Контроль'
sperm_mut$Condition = 'Мутант'

sperm_full = rbind(sperm_cont, sperm_mut)
sperm_full$Condition = as.factor(sperm_full$Condition)

group_count = as.data.frame(table(sperm_full$Condition))
names(group_count) = c('Condition', 'n')



### PLOT VCL, VSL, VAP
ggplot(data=sperm_full, aes(x=Condition, y=VSL, col=Condition))+
  geom_boxplot(show.legend = FALSE)+
  geom_text(data = group_count, aes(x = Condition, y = max(sperm_full$VSL)+1, label = paste("N =", n)), 
            position = position_dodge(width =0.25), vjust = -0.5, col='black')+
  theme_classic()

### Not significant
wilcox.test(sperm_full[sperm_full$Condition == 'Контроль',]$VSL, sperm_full[sperm_full$Condition == 'Мутант',]$VSL, paired=FALSE)  

ggplot(data=sperm_full, aes(x=Condition, y=VCL, col=Condition))+
  geom_boxplot(show.legend = FALSE)+
  geom_text(data = group_count, aes(x = Condition, y = max(sperm_full$VCL)+1, label = paste("N =", n)), 
            position = position_dodge(width =0.25), vjust = -0.5, col='black')+
  theme_classic()

###Significant
wilcox.test(sperm_full[sperm_full$Condition == 'Контроль',]$VCL, sperm_full[sperm_full$Condition == 'Мутант',]$VCL, paired=FALSE)  

ggplot(data=sperm_full, aes(x=Condition, y=VAP, col=Condition))+
  geom_boxplot(show.legend = FALSE)+
  geom_text(data = group_count, aes(x = Condition, y = max(sperm_full$VAP)+1, label = paste("N =", n)), 
            position = position_dodge(width =0.25), vjust = -0.5, col='black')+
  theme_classic()

###Significant
wilcox.test(sperm_full[sperm_full$Condition == 'Контроль',]$VAP, sperm_full[sperm_full$Condition == 'Мутант',]$VAP, paired=FALSE)  



### How many zeros in each feature

nrow(sperm_full[sperm_full$VAP==0,]) ### 464
nrow(sperm_full[sperm_full$VSL == 0,]) ### 87
nrow(sperm_full[sperm_full$VCL ==0,]) ### 42

### <----- Question: is it possible to have speed on one type of line but not in another
### If we have one zero on these features is it dead or what?

### Feature zero in each group

nrow(sperm_full[sperm_full$VAP==0 & sperm_full$Condition=='Контроль',] ) ### 72
nrow(sperm_full[sperm_full$VAP==0 & sperm_full$Condition=='Мутант',]) ### 392

nrow(sperm_full[sperm_full$VSL==0 & sperm_full$Condition=='Контроль',] ) ### 23
nrow(sperm_full[sperm_full$VSL==0 & sperm_full$Condition=='Мутант',]) ### 64

### So mutant sperm has more zeros - KEEP IN MIND mb check the concentration



### MANOVA Check

### All of them do not have normal distribution
shapiro.test(sperm_full$VCL)
#hist(sperm_full$VCL)
shapiro.test(sperm_full$VSL)
#hist(sperm_full$VSL)
shapiro.test(sperm_full$VAP)
#hist(sperm_full$VAP)


#If we still want it!
### Significant
res.man = manova(cbind(VCL, VSL, VAP) ~ Condition, data=sperm_full)
summary(res.man)

### Features that differs are the same compared to previous test (VSL is not significant)
summary.aov(res.man)

### Try another lib to perform perMANOVA (doesn't work because too much zeros)
# adonis2(sperm_full[sperm_full$VCL < 0,c('VCL', 'VSL', 'VAP')] ~ sperm_full$Condition, data=sperm_full, method = "bray")



