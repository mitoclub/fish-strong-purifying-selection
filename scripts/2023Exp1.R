rm(list=ls(all=TRUE))
library(ggplot2)

getwd()
setwd('C:/fish')
table = data.frame(read.csv('../data/2023_FISH_VNIIPRH - EXPERIMENT_1_DATE.csv', header=TRUE, sep = ","))
str(table)
View(table)
table$FERTperc = table$FERT/table$TOTAL
table$SHOCKperc = table$SHOCK/table$TOTAL
table$HATCHperc = table$HATCH/table$TOTAL
table$SWIMperc = table$SWIM/table$TOTAL
table$MODE = paste(as.character(table$TEMP), table$TIME, sep = "_")
table[table$MODE == "20_0",]$MODE = "C"
table$MutGroup = paste(table$MUT, table$KIND, sep="_")
#table[table$MutGroup == "C_S",]$MutGroup = "C"
#table[table$MutGroup == "C_R",]$MutGroup = "C"
View(table)

tableV = table[table$STAGE == "V",]
tableP = table[table$STAGE == "P",]

pdf('../figures/All_groups_fert_shock_hatch.pdf')

ggplot(table, aes(x=MUT, y=FERTperc, fill=KIND)) + 
  geom_boxplot()+
scale_fill_manual(values = c("#abb9c2", "#9c0633"))

ggplot(tableV, aes(x=MODE, y=SHOCKperc, fill=MutGroup)) + 
  geom_boxplot() +
scale_fill_manual(values = c("#70b9e6",
                              "#8fcf9b",
                              "#065f96",
                              "#139c2e"))+ggtitle("V stage")

ggplot(tableP, aes(x=MODE, y=SHOCKperc, fill=MutGroup)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#70b9e6",
                               "#8fcf9b",
                               "#065f96",
                               "#139c2e"))+ggtitle("P stage")

ggplot(tableV, aes(x=MODE, y=HATCHperc, fill=MutGroup)) + 
  geom_boxplot()+
  scale_fill_manual(values = c("#70b9e6",
                               "#8fcf9b",
                               "#065f96",
                               "#139c2e"))+ggtitle("V stage")

ggplot(tableP, aes(x=MODE, y=HATCHperc, fill=MutGroup)) + 
  geom_boxplot()+
  scale_fill_manual(values = c("#70b9e6",
                               "#8fcf9b",
                               "#065f96",
                               "#139c2e"))+ggtitle("P stage")
                               
ggplot(tableV, aes(x=MODE, y=SWIMperc, fill=MutGroup)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#70b9e6",
                                        "#8fcf9b",
                                        "#065f96",
                                        "#139c2e"))+ggtitle("V stage")
                                        
ggplot(tableP, aes(x=MODE, y=SWIMperc, fill=MutGroup)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#70b9e6",
                                        "#8fcf9b",
                                        "#065f96",
                                        "#139c2e"))+ggtitle("P stage")
                                        
                                        

dev.off()
