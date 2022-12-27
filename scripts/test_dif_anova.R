rm(list=ls(all=TRUE))

library(ggplot2)
library(dplyr)
library(car)
library(Hmisc)
library(multcomp)
fishes = read.table('../data/fishes_shok.txt', header = T)


### Change numbers to factors
fishes$family = as.factor(fishes$family)
fishes$mode = as.factor(fishes$mode)
fishes$mode_scr = as.factor(fishes$mode_scr)
fishes$time = as.factor(fishes$time)
fishes$enu = as.factor(fishes$enu)



### Check ANOVA for fertilization with family and ENU
ggplot(data = fishes, aes(x = enu , y = fertilization_per, color = family)) +
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, position = position_dodge(width = 0.2)) 


s_model <- lm(fertilization_per ~ enu * time, data=fishes)
simple<-Anova(mod = s_model)

simple_diag <- fortify(s_model) # функция из пакета ggplot2

ggplot(simple_diag, aes(x = 1:nrow(simple_diag), y = .cooksd)) +
  geom_bar(stat = 'identity')

ggplot(data = simple_diag, aes(x = family, y = .stdresid, colour = enu)) +
  geom_boxplot() + geom_hline(yintercept = 0)


qqPlot(s_model, id = FALSE) 

summary(simple)


### Test LM on hathced per with time of hit and enu
ggplot(data = fishes, aes(x = time , y = hatched_per, color = enu)) +
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, position = position_dodge(width = 0.2)) 


s_model <- lm(hatched_per ~ enu + time , data=fishes)
simple<-Anova(mod = s_model)

simple_diag <- fortify(s_model) # функция из пакета ggplot2

ggplot(simple_diag, aes(x = 1:nrow(simple_diag), y = .cooksd)) +
  geom_bar(stat = 'identity')

ggplot(data = simple_diag, aes(x = time, y = .stdresid, colour = enu)) +
  geom_boxplot() + geom_hline(yintercept = 0)


qqPlot(s_model, id = FALSE) 

summary(simple)



### Test LM on swimming per with time of hit and enu
ggplot(data = fishes, aes(x = enu , y = swim_per, color = family)) +
  stat_summary(geom = 'pointrange', fun.data = mean_cl_normal, position = position_dodge(width = 0.2))


s_model <- lm(swim_per ~ enu + time , data=fishes)
simple<-Anova(mod = s_model)

simple_diag <- fortify(s_model) # функция из пакета ggplot2

ggplot(simple_diag, aes(x = 1:nrow(simple_diag), y = .cooksd)) +
  geom_bar(stat = 'identity')

ggplot(data = simple_diag, aes(x = time, y = .stdresid, colour = enu)) +
  geom_boxplot() + geom_hline(yintercept = 0)


qqPlot(s_model, id = FALSE) 

summary(simple)
