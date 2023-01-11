rm(list=ls(all=TRUE))
install.packages("ggplot2")
install.packages("gridExtra")
library(ggplot2)


#my currently directory (change directory)
getwd()
setwd('C:/Users/sozve/Desktop/fish-strong-purifying-selection/')

###ENU###
#Upload table
Fish = read.table('data/fishes_shok.txt', header=TRUE)

### Change numbers to factors
Fish$family = as.factor(Fish$family)
Fish$mode = as.factor(Fish$mode)
Fish$mode_scr = as.factor(Fish$mode_scr)
Fish$time = as.factor(Fish$time)
Fish$enu = as.factor(Fish$enu)


#create table with only ENU induced fishes/controls without temperature
enu_n_0 = paste('N=', nrow(Fish[Fish$enu == 0 & Fish$temperature == 20, ]), sep='')
enu_n_1 = paste('N=', nrow(Fish[Fish$enu == 1 & Fish$temperature == 20, ]), sep='')
enu_n_2 = paste('N=', nrow(Fish[Fish$enu == 2 & Fish$temperature == 20, ]), sep='')


### make labels with number for fertilization for different concentration ENU

enu_n_0_unclear = paste('N=', nrow(Fish[Fish$enu == 0,]), sep='')
enu_n_1_unclear = paste('N=', nrow(Fish[Fish$enu == 1,]), sep='')
enu_n_2_unclear = paste('N=', nrow(Fish[Fish$enu == 2,]), sep='')

#Starting PDF
pdf(file='figures/ENU_10.01.23.pdf', width=20, height=20)

ggtitle("Fish without temperature")

#GGplot fertilization
f0 <- ggplot(Fish[Fish$temperature==20,], (aes(enu, fertilization_per)))+
  geom_boxplot(colour = "black", notch = TRUE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", enu_n_0,sep='\n'), paste("1", enu_n_1,sep='\n'), paste("2", enu_n_2,sep='\n')))
f <- f0 + geom_point(aes(col=family))+ggtitle("Fertilization percent \n vs ENU")+
  xlab("Dose (ENU)") + ylab("Fertilization percent")

#Statistics
ANOVA_TEST_ENU_VS_fertilization = aov(Fish$fertilization_per ~ Fish$enu)
summary(ANOVA_TEST_ENU_VS_fertilization)
TukeyHSD(ANOVA_TEST_ENU_VS_fertilization, conf.level=.95)

#Hatched vs ENU
h0 <- ggplot(Fish[Fish$temperature==20,], (aes(enu, hatched_per)))+
  geom_boxplot(colour = "gray", notch = FALSE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", enu_n_0,sep='\n'), paste("1", enu_n_1,sep='\n'), paste("2", enu_n_2,sep='\n')))
h <- h0 + geom_point(aes(col=family))+ggtitle("Hatched percent \n vs ENU") +
  xlab("Dose (ENU)") + ylab("Hatched percent")

#Statistics
Fish_wth_tm = Fish[Fish$temperature ==20,]
ANOVA_TEST_ENU_VS_hatched = aov(Fish_wth_tm$hatched_per ~ Fish_wth_tm$enu)
summary(ANOVA_TEST_ENU_VS_hatched)
TukeyHSD(ANOVA_TEST_ENU_VS_hatched, conf.level=.95)

#Swimming vs ENU
s0 <- ggplot(Fish[Fish$temperature == 20, ], (aes(enu, swim_per)))+
  geom_boxplot(colour = "orange", notch = TRUE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", enu_n_0,sep='\n'), paste("1", enu_n_1,sep='\n'), paste("2", enu_n_2,sep='\n')))
s <- s0 + geom_point(aes(col=family))+ggtitle("Swimming percent \n vs ENU")+
  xlab("Dose (ENU)") + ylab("Swimming percent")

#ANOVA(является ли ену значимым в проценте вылупившимся?), Tukey(какие именно различаются?)
ANOVA_TEST_ENU_VS_swim = aov(Fish$swim_per ~ Fish$enu)
summary(ANOVA_TEST_ENU_VS_swim)
TukeyHSD(ANOVA_TEST_ENU_VS_swim, conf.level=.95)

#GGplot with freaks_per
fr0 <- ggplot(Fish[Fish$temperature == 20, ], (aes(enu, freaks_per)))+
  geom_boxplot(colour = "green", notch = TRUE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", enu_n_0,sep='\n'), paste("1", enu_n_1,sep='\n'), paste("2", enu_n_2,sep='\n')))
fr <- fr0 + geom_point(aes(col=family))+ggtitle("Freaks percent \n vs ENU") +
  xlab("Dose (ENU)") + ylab("Freaks percent")

#ANOVA(является ли ену значимым в проценте вылупившимся?), Tukey(какие именно различаются?)
ANOVA_TEST_ENU_VS_freaks = aov(Fish$freaks_per ~ Fish$enu)
summary(ANOVA_TEST_ENU_VS_freaks)
TukeyHSD(ANOVA_TEST_ENU_VS_freaks, conf.level=.95)

#sum graphics on 1 page
library(gridExtra)
grid.arrange(f, h, s, fr, ncol=2)

#end PDF
dev.off()

###TEMPERATURE###
#Upload table
Fish = read.table('data/fishes_shok.txt', header=TRUE)

### Change numbers to factors
Fish$family = as.factor(Fish$family)
Fish$mode = as.factor(Fish$mode)
Fish$mode_scr = as.factor(Fish$mode_scr)
Fish$time = as.factor(Fish$time)
Fish$enu = as.factor(Fish$enu)

#prepare table with ONLY Temperature effect
Fish_temp_intact <- Fish[Fish$enu == 0, ]

#Starting PDF
pdf('figures/Temperature__10.01.23.pdf', width=20, height=20)
ggtitle("Fish")

#Step sep our exposure time by 0 30 40 50 for the control controls + native
Fish_temp_intact$time = as.factor(Fish_temp_intact$time)
T_0 = paste('N=',nrow(Fish_temp_intact[Fish_temp_intact$time == 0,]),sep='')
T_30 = paste('N=',nrow(Fish_temp_intact[Fish_temp_intact$time == 30,]),sep='')
T_40 = paste('N=',nrow(Fish_temp_intact[Fish_temp_intact$time == 40,]),sep='')
T_50 = paste('N=',nrow(Fish_temp_intact[Fish_temp_intact$time == 50,]),sep='')

#Fertilization(intact fishes) doesn't depend on temperature

#Hatched(intact fishes)
th0 <- ggplot(Fish_temp_intact, (aes(time, hatched_per)))+
  geom_boxplot(notch = FALSE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", T_0,sep='\n'), paste("30", T_30,sep='\n'), paste("40", T_40,sep='\n'), paste("50", T_50,sep='\n')))
th <- th0 + geom_point(aes(col=family))+ggtitle("Hatched percent \n Without temperature")+
  xlab("Time NOT exposure") + ylab("Hatched percent")
th = th + theme_bw()

#Statistics
ANOVA_TEST_Temp_vs_hatched_per = aov(Fish_temp_intact$hatched_per ~ Fish_temp_intact$time)
summary(ANOVA_TEST_Temp_vs_hatched_per)

#Swimming(intact fishes)
ts0 <- ggplot(Fish_temp_intact, (aes(time, swim_per)))+
  geom_boxplot(notch = FALSE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", T_0,sep='\n'), paste("30", T_30,sep='\n'), paste("40", T_40,sep='\n'), paste("50", T_50,sep='\n')))
ts <- ts0 + geom_point(aes(col=family))+ggtitle("Swimming percent \n Without temperature")+
  xlab("Time NOT exposure") + ylab("Swimming percent")

#Statistics
ANOVA_TEST_Temp_VS_swim_per = aov(Fish_temp_intact$swim_per ~ Fish_temp_intact$time)
summary(ANOVA_TEST_Temp_VS_swim_per)


#Freaks(intact fishes)
tfr0 <- ggplot(Fish_temp_intact, (aes(time, freaks_per)))+
  geom_boxplot(notch = FALSE)+
  theme_classic()+
  scale_x_discrete(labels=c(paste("0", T_0,sep='\n'), paste("30", T_30,sep='\n'), paste("40", T_40,sep='\n'), paste("50", T_50,sep='\n')))
tfr <- tfr0 + geom_point(aes(col=family))+ggtitle("Freaks percent \n Without temperature")+
  xlab("Time NOT exposure") + ylab("Freaks percent")

#Statistics
ANOVA_TEST_Temp_VS_freaks_per = aov(Fish_temp_intact$freaks_per ~ Fish_temp_intact$time)
summary(ANOVA_TEST_Temp_VS_freaks_per)

#sum graphics on 1 page
library(gridExtra)
grid.arrange(th, ts, tfr, ncol=2)

#stop and save PDF file 
dev.off()

###ENU + TIME###
Fish = read.table('data/fishes_shok.txt', header=TRUE)

#Change numbers to factors
Fish$family = as.factor(Fish$family)
Fish$mode = as.factor(Fish$mode)
Fish$mode_scr = as.factor(Fish$mode_scr)
Fish$time = as.factor(Fish$time)
Fish$enu = as.factor(Fish$enu)

#Starting PDF #change directory
pdf(file='figures/ENU_Temperature_10.01.23.pdf', width=20, height=20)

#ENU VS Fertilization percent

ttf01 <- ggplot(Fish, (aes(mode_scr, fertilization_per, fill = enu)))+
  geom_boxplot(aes(col=enu))+
  theme_classic()
ttf011 <- ttf01+ geom_point(aes(pch=family))+ggtitle("Fertilization percent \n vs ENU")+
  xlab("ENU") + ylab("Fertilization percent")
ttf011 = ttf011 + theme_bw()

#Statistics
MANOVA_TEST_ENU_Temp_VS_Fertilization = manova(cbind(temperature, enu) ~ fertilization_per, data = Fish)
summary(MANOVA_TEST_ENU_Temp_VS_Fertilization)

#ENU VS Hatched percent + Temperature
tth01 <- ggplot(Fish, (aes(mode_scr, hatched_per, fill=enu)))+
  geom_boxplot()+
  theme_classic()
tth011 <- tth01+ geom_point(aes(pch=family))+ggtitle("Hatched percent \n vs ENU + Temperature")+
  xlab("ENU") + ylab("Hatched percent")
tth011 + theme_bw()

#Statistics
MANOVA_TEST_ENU_Temp_VS_Hatched = manova(cbind(temperature, enu) ~ hatched_per, data = Fish)
summary(MANOVA_TEST_ENU_Temp_VS_Hatched)

#TIME + ENU VS Swimming percent
tts01 <- ggplot(Fish, (aes(mode_scr, swim_per, fill=enu)))+
  geom_boxplot()+
  theme_classic()
tts011 <- tts01+ geom_point(aes(pch=family))+ggtitle("Swimming percent \n vs ENU + Temperature")+
  xlab("ENU + Temperature") + ylab("Swimming percent")
tts011 + theme_bw()

#Statistics
MANOVA_TEST_ENU_Temp_VS_SWimming = manova(cbind(temperature, enu) ~ swimming_per, data = Fish)
summary(MANOVA_TEST_ENU_Temp_VS_Swimming)

#TIME + ENU VS Freaks percent
ttfr01 <- ggplot(Fish, (aes(mode_scr, freaks_per, fill = enu)))+
  geom_boxplot()+
  theme_classic()
ttfr011 <- ttfr01+ geom_point(aes(shape=family))+ggtitle("Freaks percent \n vs ENU + Temperature")+
  xlab("ENU + Temperature") + ylab("Freaks percent")
ttfr011 + theme_bw()

#Statistics
MANOVA_TEST_ENU_Temp_VS_Freaks = manova(cbind(temperature, enu) ~ freaks_per, data = Fish)
summary(MANOVA_TEST_ENU_Temp_VS_Freaks)

#multiple graphics
library(gridExtra)
grid.arrange(ttf011, tth011, tts011, ttfr011, ncol=2)

#stop and save PDF file 
dev.off()

###ENU + TEMPERATURE IN FAMILIES###

#Upload table
Fish = read.table('data/fishes_shok.txt', header=TRUE)
Fish$family = as.factor(Fish$family)
Fish$mode = as.factor(Fish$mode)
Fish$mode_scr = as.factor(Fish$mode_scr)
Fish$time = as.factor(Fish$time)
Fish$enu = as.factor(Fish$enu)

View(Fish)
#Separate families
Fish_1x1 = Fish[Fish$family == '1x1',]
Fish_1x2 = Fish[Fish$family == '1x2',]
Fish_2x1 = Fish[Fish$family == '2x1',]
Fish_2x2 = Fish[Fish$family == '2x2',]

#Starting PDF
pdf(file='figures/ENU_Temperature_Fam_10.01.23.pdf', width=20, height=20)

#ENU VS Fertilization percent 1x1
ttf01_1x1 <- ggplot(Fish_1x1, (aes(mode_scr, fertilization_per)))+
  geom_boxplot(aes(col=enu))+
  theme_classic()

ttf011_1x1 <- ttf01_1x1+ geom_point(aes(col=mode_scr))+ggtitle("Fertilization percent 1x1\n vs ENU")+
  xlab("ENU") + ylab("Fertilization percent")

#Statistics
MANOVA_TEST_ENU_Temp_VS_Fertilization_1x1 = manova(cbind(temperature, enu) ~ fertilization_per, data = Fish_1x1)
summary(MANOVA_TEST_ENU_Temp_VS_Fertilization_1x1)

#ENU VS Fertilization percent family 1x2
ttf01_1x2 <- ggplot(Fish_1x2, (aes(mode_scr, fertilization_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttf011_1x2 <- ttf01_1x2+ geom_point(aes(col=mode_scr))+ggtitle("Fertilization 1x2 percent \n vs ENU")+
  xlab("ENU") + ylab("Fertilization percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Fertilization_1x2 = manova(cbind(temperature, enu) ~ fertilization_per, data = Fish_1x2)
summary(MANOVA_TEST_ENU_Temp_VS_Fertilization_1x2)

#ENU VS Fertilization percent family 2x1
Fish_2x1$enu=as.factor(Fish_2x1$enu)
ttf01_2x1 <- ggplot(Fish_2x1, (aes(mode_scr, fertilization_per, fill=enu)))+
  geom_boxplot(aes(col=enu))+
  theme_classic()
ttf011_2x1 <- ttf01_2x1+ geom_point(aes(col=enu))+ggtitle("Fertilization percent 2x1 \n vs ENU")+
  xlab("ENU") + ylab("Fertilization percent")
ttf011_2x1 + theme_bw()

#Statistics
MANOVA_TEST_ENU_Temp_VS_Fertilization_2x1 = manova(cbind(temperature, enu) ~ fertilization_per, data = Fish_2x1)
summary(MANOVA_TEST_ENU_Temp_VS_Fertilization_2x1)

#ENU VS Fertilization percent family 2x2
ttf01_2x2 <- ggplot(Fish_2x2, (aes(mode_scr, fertilization_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttf011_2x2 <- ttf01_2x2+ geom_point(aes(col=mode_scr))+ggtitle("Fertilization percent 2x2 \n vs ENU")+
  xlab("ENU") + ylab("Fertilization percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Fertilization_2x2 = manova(cbind(temperature, enu) ~ fertilization_per, data = Fish_2x2)
summary(MANOVA_TEST_ENU_Temp_VS_Fertilization_2x2)

#Time+ENU VS Hatched percent 1x1
tth01_1x1 <- ggplot(Fish_1x1, (aes(mode_scr, hatched_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tth011_1x1 <- ttf01_1x1+ geom_point(aes(col=mode_scr))+ggtitle("Hatched percent 1x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Hatched percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Hatched_1x1 = manova(cbind(temperature, enu) ~ hatched_per, data = Fish_1x1)
summary(MANOVA_TEST_ENU_Temp_VS_Hatched_1x1)

#Time+ENU VS Hatched percent 1x2
tth01_1x2 <- ggplot(Fish_1x2, (aes(mode_scr, hatched_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tth011_1x2 <- ttf01_1x2+ geom_point(aes(col=mode_scr))+ggtitle("Hatched percent 1x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Hatched percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Hatched_1x2 = manova(cbind(temperature, enu) ~ hatched_per, data = Fish_1x2)
summary(MANOVA_TEST_ENU_Temp_VS_Hatched_1x2)

#Time+ENU VS Hatched percent 2x1
tth01_2x1 <- ggplot(Fish_2x1, (aes(mode_scr, hatched_per, fill = enu)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tth011_2x1 <- ttf01_2x1+ geom_point(aes(col=enu))+ggtitle("Hatched percent 2x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Hatched percent")
tth011_2x1 + theme_bw()
#Statistics
MANOVA_TEST_ENU_Temp_VS_Hatched_2x1 = manova(cbind(temperature, enu) ~ hatched_per, data = Fish_2x1)
summary(MANOVA_TEST_ENU_Temp_VS_Hatched_2x1)

#Time+ENU VS Hatched percent 2x2
tth01_2x2 <- ggplot(Fish_2x2, (aes(mode_scr, hatched_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tth011_2x2 <- ttf01_2x2+ geom_point(aes(col=mode_scr))+ggtitle("Hatched percent 2x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Hatched percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Hatched_2x2 = manova(cbind(temperature, enu) ~ hatched_per, data = Fish_2x2)
summary(MANOVA_TEST_ENU_Temp_VS_Hatched_2x2)

#multiple 'fertilization' graphics
library(gridExtra)
grid.arrange(tth011_1x1, tth011_1x2, tth011_2x1, tth011_2x2, ncol=2)

#ENU VS Swimming percent 1x1
tts01_1x1 <- ggplot(Fish_1x1, (aes(mode_scr, swim_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tts011_1x1 <- tts01_1x1+ geom_point(aes(col=mode_scr))+ggtitle("Swimming percent 1x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Swimming percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_SWimming_1x1 = manova(cbind(temperature, enu) ~ swim_per, data = Fish_1x1)
summary(MANOVA_TEST_ENU_Temp_VS_Swimming_1x1)

#ENU VS Swimming percent 1x2
tts01_1x2 <- ggplot(Fish_1x2, (aes(mode_scr, swim_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tts011_1x2 <- tts01_1x2+ geom_point(aes(col=mode_scr))+ggtitle("Swimming percent 1x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Swimming percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_SWimming_1x2 = manova(cbind(temperature, enu) ~ swim_per, data = Fish_1x2)
summary(MANOVA_TEST_ENU_Temp_VS_SWimming_1x2)

#ENU VS Swimming percent 2x1
tts01_2x1 <- ggplot(Fish_2x1, (aes(mode_scr, swim_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tts011_2x1 <- tts01_2x1+ geom_point(aes(col=mode_scr))+ggtitle("Swimming percent 2x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Swimming percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_SWimming_2x1 = manova(cbind(temperature, enu) ~ swim_per, data = Fish_2x1)
summary(MANOVA_TEST_ENU_Temp_VS_Swimming_2x1)

#ENU VS Swimming percent 2x2
tts01_2x2 <- ggplot(Fish_2x2, (aes(mode_scr, swim_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
tts011_2x2 <- tts01_2x2+ geom_point(aes(col=mode_scr))+ggtitle("Swimming percent 2x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Swimming percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_SWimming_2x2 = manova(cbind(temperature, enu) ~ swim_per, data = Fish_2x2)
summary(MANOVA_TEST_ENU_Temp_VS_Swimming_2x2)

#multiple 'swim' graphics
library(gridExtra)
grid.arrange(tts011_1x1, tts011_1x2, tts011_2x1, tts011_2x2, ncol=2)

#ENU VS Freaks percent 1x1
ttfr01_1x1 <- ggplot(Fish_1x1, (aes(mode_scr, freaks_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttfr011_1x1 <- ttfr01_1x1+ geom_point(aes(col=mode_scr))+ggtitle("Freaks percent 1x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Freaks percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Freaks_1x1 = manova(cbind(temperature, enu) ~ freak_per, data = Fish_1x1)
summary(MANOVA_TEST_ENU_Temp_VS_Freaks_1x1)

#ENU VS Freaks percent 1x2
ttfr01_1x2 <- ggplot(Fish_1x2, (aes(mode_scr, freaks_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttfr011_1x2 <- ttfr01_1x2+ geom_point(aes(col=mode_scr))+ggtitle("Freaks percent 1x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Freaks percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Freaks_1x2 = manova(cbind(temperature, enu) ~ freak_per, data = Fish_1x2)
summary(MANOVA_TEST_ENU_Temp_VS_Freaks_1x2)

#ENU VS Freaks percent 2x1
ttfr01_2x1 <- ggplot(Fish_2x1, (aes(mode_scr, freaks_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttfr011_2x1 <- ttfr01_2x1+ geom_point(aes(col=mode_scr))+ggtitle("Freaks percent 2x1 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Freaks percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Freaks_2x1 = manova(cbind(temperature, enu) ~ freak_per, data = Fish_2x1)
summary(MANOVA_TEST_ENU_Temp_VS_Freaks_2x1)

#ENU VS Freaks percent 2x2
ttfr01_2x2 <- ggplot(Fish_2x2, (aes(mode_scr, freaks_per)))+
  geom_boxplot(aes(col=mode_scr))+
  theme_classic()
ttfr011_2x2 <- ttfr01_2x2+ geom_point(aes(col=mode_scr))+ggtitle("Freaks percent 2x2 \n vs Time+ENU")+
  xlab("Time+ENU") + ylab("Freaks percent")
#Statistics
MANOVA_TEST_ENU_Temp_VS_Freaks_2x2 = manova(cbind(temperature, enu) ~ freak_per, data = Fish_2x2)
summary(MANOVA_TEST_ENU_Temp_VS_Freaks_2x2)

#multiple 'fertilization' graphics
library(gridExtra)
grid.arrange(ttfr011_1x1, ttfr011_1x2, ttfr011_2x1, ttfr011_2x2, ncol=2)

#stop and save PDF file 
dev.off()

##old statistics example##
#ANOVA_TEST_ENU_Temp_VS_Hatched_1x1 = aov(Fish_1x1$hatched_per ~ Fish_1x1$mode_scr)
#summary(ANOVA_TEST_ENU_Temp_VS_Hatched_1x1)
#TukeyHSD(ANOVA_TEST_ENU_Temp_VS_Hatched_1x1, conf.level=.95)





