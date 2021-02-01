#file originally written by Hannah Rubin June 2019
#edited for reproducibility and train/test update Nov/Dec 2020/Jan2021 by David Lutz and Bethel Steele

library(tidyverse)
# library(dplyr)
# library(tidyr)
# library(lubridate)
# library(ggplot2)
library(data.table)
library(Metrics)
library(gridExtra)
library(randomForest)
# library(writexl)
library(broom)

#add final theme
final_theme=theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=16, face='bold', hjust=0.5)) #save as a grom


# create directory paths, commented out for us on Dave's machine
datadir <- '5.a.90.10.datasets/'
dumpdir <- 'Rubin.MS/summary_data/'

# #local directory for Dave testing
# setwd("D:/NASA IDS/Manuscript 3 Rubin/DaveTestNov")

# Read in data
training = read_csv(paste0(datadir, "90percent_7.23.20.csv"),
                    col_types = cols(.default = col_character())) %>% 
  mutate(secchi_depth_m = as.numeric(secchi_depth_m)) %>% 
  mutate_at(vars(c(secchi_depth_m, B1:B4, B1B3:B2B1)),
            ~as.numeric(.))
testing = read_csv(paste0(datadir, "OOB_7.23.20.csv"),
                   col_types = cols(.default = col_character())) %>% 
  mutate(secchi_depth_m = as.numeric(secchi_depth_m)) %>% 
  mutate_at(vars(c(secchi_depth_m, B1:B4, B1B3:B2B1)),
            ~as.numeric(.))

# rsq <- function(x,y) summary(lm(y~x))$r.squared

######################################################################
#remove single sample from training set with very low B3 value (possibly erroneous but technically within confidence intervals)

training <- training %>% 
  filter(B3 >= 3)


##### look at histograms #####
hist(training$secchi_depth_m, xlim = c(0,19),
     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
     xlab = 'Secchi depth (m)', 
     main = NULL)
summary(training$secchi_depth_m)
hist(testing$secchi_depth_m,
     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
     xlab = 'Secchi depth (m)', 
     main = NULL)
summary(testing$secchi_depth_m)

train_test <- training %>% 
  mutate(dataset = 'training') %>% 
  full_join(., testing) %>% 
  mutate(dataset = case_when(is.na(dataset) ~ 'testing',
                             TRUE ~ dataset))

histogram <- ggplot(train_test, aes(x = secchi_depth_m)) +
  geom_histogram(bins = 30) +
  facet_grid(dataset ~ ., scales = 'free_y') +
  labs(x = 'Secchi depth (m)',
       y = 'frequency')  +
  theme_bw() +
  final_theme
histogram
# ggsave('Rubin.MS/manuscript/train_test_histogram.jpg', units = 'in', height = 5, width = 5, dpi = 250)


##########################Algorithms work flow, set up table##########
Table3 <- data.frame(matrix(nrow = 1, ncol = 10))
colnames(Table3)[1] <- "Model Name"
colnames(Table3)[2] <- "Training MAE (m)"
colnames(Table3)[3] <- "Training RMSE (m)"
colnames(Table3)[4] <- "Training pseudo-R2"
colnames(Table3)[5] <- "Testing MAE (m)"
colnames(Table3)[6] <- "Testing RMSE (m)"
colnames(Table3)[7] <- "Testing pseudo-R2"
colnames(Table3)[8] <- "Training p-value"
colnames(Table3)[9] <- "Training P-O slope"
colnames(Table3)[10] <- "Testing P-O slope"

###########################Allee AND JOHNSON ##############################
#Allee and Johnson, 1999

#create data frame and calculate needed values for algorithm
AlleeTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3)
AlleeTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3)
B3mean_Training <- mean(AlleeTrain[["B3"]])
AlleeTrain$B3dev <- (AlleeTrain$B3 - B3mean_Training)
AlleeTest$B3dev <- (AlleeTest$B3 - B3mean_Training)  
#check work
AlleeTrain
AlleeTest

#calculate constants for algorithm based on training data
AlleeModel <- lm(secchi_depth_m ~ poly(B3dev, 3), data = AlleeTrain)
AlleeTrainTable <- data.table(AlleeModel$model$secchi_depth_m, AlleeModel$fitted.values)
colnames(AlleeTrainTable) = c('secchi_depth_m', 'fitted.values')
AlleeTrainTable

#record training statistics
Table3[1,1] <- "Allee and Johnson"
Table3[1,2] <- mae(AlleeModel$model$secchi_depth_m, AlleeModel$fitted.values)
Table3[1,3] <- rmse(AlleeModel$model$secchi_depth_m, AlleeModel$fitted.values)
Table3[1,4] <-  1 - (sum((AlleeTrainTable$fitted.values - AlleeTrainTable$secchi_depth_m) ^ 2)/sum((AlleeTrainTable$secchi_depth_m - mean(AlleeTrainTable$secchi_depth_m)) ^ 2))
Table3[1,8] = glance(AlleeModel)$p.value

#apply model to test data
AlleeTestModel <- as.numeric(predict(AlleeModel, newdata = AlleeTest))
AlleeTestTable <- data.table(AlleeTestModel)
AlleeTestTable$Secchi <- as.numeric(AlleeTest$secchi_depth_m)
AlleeTestTable

#record test statitstics
Table3[1,5] = mae(AlleeTestTable$Secchi, AlleeTestTable$AlleeTestModel)
Table3[1,6] = rmse(AlleeTestTable$Secchi, AlleeTestTable$AlleeTestModel)
Table3[1,7] = 1 - (sum((AlleeTestTable$AlleeTestModel - AlleeTestTable$Secchi) ^ 2)/sum((AlleeTestTable$Secchi - mean(AlleeTestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
AlleeFigureTrain <- ggplot(AlleeTrainTable, aes(secchi_depth_m, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Allee/Johnson *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

AlleeFigureTest <- ggplot(AlleeTestTable, aes(Secchi, AlleeTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Allee/Johnson", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

AlleeFigureTrain
AlleeFigureTest

#Note, allee training creates negative secchi values

#Add columns for predicted/observed summary
Table3[1,9] = summary(lm(AlleeTrainTable$fitted.values~AlleeTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[1,10] = summary(lm(AlleeTestTable$AlleeTestModel~AlleeTestTable$Secchi))$coefficients[2,1]
Table3

##########################BABAN###########################################

#create data frame and calculate needed values for algorithm
BabanTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B1)
BabanTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B1)
BabanTrain
BabanTest

#calculate constants for algorithm based on training data
BabanModel <- lm(secchi_depth_m ~ B1, data = BabanTrain )
BabanTrainTable <- data.table(BabanModel$model$secchi_depth_m, BabanModel$fitted.values) 
colnames(BabanTrainTable) = c('secchi_depth_m', 'fitted.values')
BabanTrainTable

#record training statistics
Table3[2,1] <- "Baban"
Table3[2,2] <- mae(BabanModel$model$secchi_depth_m, BabanModel$fitted.values)
Table3[2,3] <- rmse(BabanModel$model$secchi_depth_m, BabanModel$fitted.values)
Table3[2,4] <-  1 - (sum((BabanTrainTable$fitted.values - BabanTrainTable$secchi_depth_m) ^ 2)/sum((BabanTrainTable$secchi_depth_m - mean(BabanTrainTable$secchi_depth_m)) ^ 2))

#apply model to test data
BabanTestModel <- predict(BabanModel, newdata = BabanTest)
BabanTestTable <- data.table(BabanTestModel)
BabanTestTable$Secchi <- BabanTest$secchi_depth_m
BabanTestTable

#record test statistics
Table3[2,5] = mae(BabanTestTable$Secchi, BabanTestTable$BabanTestModel)
Table3[2,6] = rmse(BabanTestTable$Secchi, BabanTestTable$BabanTestModel)
Table3[2,7] = 1 - (sum((BabanTestTable$BabanTestModel - BabanTestTable$Secchi) ^ 2)/sum((BabanTestTable$Secchi - mean(BabanTestTable$Secchi)) ^ 2))
Table3[2,8] = glance(BabanModel)$p.value
Table3

#create and save plots of observed-predicted
BabanFigureTrain <- ggplot(BabanTrainTable, aes(secchi_depth_m, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0, 1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Baban", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)
BabanFigureTest <- ggplot(BabanTestTable, aes(Secchi, BabanTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Baban", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

BabanFigureTrain
BabanFigureTest

#Add columns for predicted/observed summary
Table3[2,9] = summary(lm(BabanTrainTable$fitted.values~BabanTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[2,10] = summary(lm(BabanTestTable$BabanTestModel~BabanTestTable$Secchi))$coefficients[2,1]
Table3

#################Chipman et al.######################################

#create data frame and calculate needed values for algorithm
ChipmanTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B1B3)
ChipmanTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B1B3)
ChipmanTrain$lnSecchi <- log(ChipmanTrain$secchi_depth_m)
ChipmanTest$lnSecchi <- log(ChipmanTest$secchi_depth_m)
ChipmanTrain
ChipmanTest

#calculate constants for algorithm based on training data
ChipmanModel <- lm(lnSecchi ~ B1B3, data = ChipmanTrain)
ChipmanTrainTable <- data.table(ChipmanModel$model$lnSecchi, ChipmanModel$fitted.values)
colnames(ChipmanTrainTable) = c('lnSecchi', 'ln.fitted.values')
ChipmanTrainTable

#record training statistics
Table3[3,1] <- "Chipman et al."
Table3[3,2] <- mae(exp(ChipmanModel$model$lnSecchi), exp(ChipmanModel$fitted.values))
Table3[3,3] <- rmse(exp(ChipmanModel$model$lnSecchi), exp(ChipmanModel$fitted.values))
Table3[3,4] <-  1 - (sum((exp(ChipmanTrainTable$ln.fitted.values) - exp(ChipmanTrainTable$lnSecchi)) ^ 2)/sum((exp(ChipmanTrainTable$lnSecchi) - mean(exp(ChipmanTrainTable$lnSecchi))) ^ 2))
Table3[3,8] = glance(ChipmanModel)$p.value

#apply model to test data
ChipmanTestModel <- predict(ChipmanModel, newdata = ChipmanTest)
ChipmanTestTable <- data.table(ChipmanTestModel)
colnames(ChipmanTestTable) = 'lnChipmanTestModel'
ChipmanTestTable$lnSecchi <- log(ChipmanTest$secchi_depth_m)
ChipmanTestTable

#record test statistics
Table3[3,5] = mae((exp(ChipmanTestTable$lnSecchi)), (exp(ChipmanTestTable$lnChipmanTestModel)))
Table3[3,6] = rmse((exp(ChipmanTestTable$lnSecchi)), (exp(ChipmanTestTable$lnChipmanTestModel)))
Table3[3,7] <-  1 - (sum((exp(ChipmanTestTable$lnChipmanTestModel) - exp(ChipmanTestTable$lnSecchi)) ^ 2)/sum((exp(ChipmanTestTable$lnSecchi) - mean(exp(ChipmanTestTable$lnSecchi))) ^ 2))
Table3

#create and save plots of observed-predicted
ChipmanFigureTrain <- ggplot(ChipmanTrainTable, aes((exp(lnSecchi)), exp(ln.fitted.values))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Chipman et al. **", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

ChipmanFigureTest <- ggplot(ChipmanTestTable, aes(exp(lnSecchi), exp(lnChipmanTestModel))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Chipman et al. **", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

ChipmanFigureTrain
ChipmanFigureTest

#need to figure out what the error is about, because it's not negative values:

#reverse log transform to double check that's what's happening
ChipmanTrainTable <- ChipmanTrainTable %>% 
  mutate(Secchi = exp(lnSecchi),
         predicted = exp(ln.fitted.values))

#turns out it's astronomically high fitted values:
ChipmanTrainTable %>% 
  arrange(-predicted)

#Add columns for predicted/observed summary
Table3[3,9] = summary(lm(exp(ChipmanTrainTable$ln.fitted.values)~exp(ChipmanTrainTable$lnSecchi)))$coefficients[2,1]
Table3[3,10] = summary(lm(exp(ChipmanTestTable$lnChipmanTestModel)~exp(ChipmanTestTable$lnSecchi)))$coefficients[2,1]
Table3

#################DEKKER and PETERS, ln(B3)#################################################

#create data frame and calculate needed values for algorithm
DekkerTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3)
DekkerTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3)
DekkerTrain$lnSecchi <- log(DekkerTrain$secchi_depth_m)
DekkerTest$lnSecchi <- log(DekkerTest$secchi_depth_m)
DekkerTrain
DekkerTest

#calculate constants for algorithm based on training data
DekkerModel <- lm(lnSecchi ~ log(B3), data = DekkerTrain)
DekkerTrainTable <- data.table(DekkerModel$model$lnSecchi, DekkerModel$fitted.values)
colnames(DekkerTrainTable) = c('lnSecchi', 'ln.fitted.values')
DekkerTrainTable

#record training statistics
Table3[4,1] <- "Dekker and Peters"
Table3[4,2] <- mae(exp(DekkerModel$model$lnSecchi), exp(DekkerModel$fitted.values))
Table3[4,3] <- rmse(exp(DekkerModel$model$lnSecchi), exp(DekkerModel$fitted.values))
Table3[4,4] <-  1 - (sum((exp(DekkerTrainTable$ln.fitted.values) - exp(DekkerTrainTable$lnSecchi)) ^ 2)/sum((exp(DekkerTrainTable$lnSecchi) - mean(exp(DekkerTrainTable$lnSecchi))) ^ 2))
Table3[4,8] = glance(DekkerModel)$p.value

#apply model to test data
DekkerTestModel <- predict(DekkerModel, newdata = DekkerTest)
DekkerTestTable <- data.table(DekkerTestModel)
colnames(DekkerTestTable) = 'lnDekkerTestModel'
DekkerTestTable$lnSecchi <- log(DekkerTest$secchi_depth_m)
DekkerTestTable

#record test statistics
Table3[4,5] = mae(exp(DekkerTestTable$lnSecchi), exp(DekkerTestTable$lnDekkerTestModel))
Table3[4,6] = rmse(exp(DekkerTestTable$lnSecchi), exp(DekkerTestTable$lnDekkerTestModel))
Table3[4,7] <-  1 - (sum((exp(DekkerTestTable$lnDekkerTestModel) - exp(DekkerTestTable$lnSecchi)) ^ 2)/sum((exp(DekkerTestTable$lnSecchi) - mean(exp(DekkerTestTable$lnSecchi))) ^ 2))
Table3

#create and save plots of observed-predicted
DekkerFigureTrain <- ggplot(DekkerTrainTable, aes((exp(lnSecchi)), exp(ln.fitted.values))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dekker/Peters", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DekkerFigureTest <- ggplot(DekkerTestTable, aes((exp(lnSecchi)), exp(lnDekkerTestModel))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dekker/Peters", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DekkerFigureTrain
DekkerFigureTest

#Add columns for predicted/observed summary
Table3[4,9] = summary(lm(exp(DekkerTrainTable$ln.fitted.values)~exp(DekkerTrainTable$lnSecchi)))$coefficients[2,1]
Table3[4,10] = summary(lm(exp(DekkerTestTable$lnDekkerTestModel)~exp(DekkerTestTable$lnSecchi)))$coefficients[2,1]
Table3

#################DEKKER and PETERS, (B3)#################################################

#create data frame and calculate needed values for algorithm
DekkerB3Train <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3)
DekkerB3Test <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3)

#calculate constants for algorithm based on training data
DekkerB3Model <- lm(secchi_depth_m ~ B3, data = DekkerB3Train )
DekkerB3TrainTable <- data.table(DekkerB3Model$model$secchi_depth_m, DekkerB3Model$fitted.values)
colnames(DekkerB3TrainTable) = c('Secchi', 'fitted.values')
DekkerB3TrainTable

#record training statistics
Table3[5,1] <- "Dekker and Peters 2"
Table3[5,2] <- mae(DekkerB3Model$model$secchi_depth_m, DekkerB3Model$fitted.values)
Table3[5,3] <- rmse(DekkerB3Model$model$secchi_depth_m, DekkerB3Model$fitted.values)
Table3[5,4] <-  1 - (sum((DekkerB3TrainTable$fitted.values - DekkerB3TrainTable$Secchi) ^ 2)/sum((DekkerB3TrainTable$Secchi - mean(DekkerB3TrainTable$Secchi)) ^ 2))
Table3[5,8] = glance(DekkerB3Model)$p.value

#apply model to test data
DekkerB3TestModel <- predict(DekkerB3Model, newdata = DekkerB3Test)
DekkerB3TestTable <- data.table(DekkerB3TestModel)
DekkerB3TestTable$Secchi <- DekkerB3Test$secchi_depth_m
DekkerB3TestTable

#record test statistics
Table3[5,5] = mae(DekkerB3TestTable$Secchi, DekkerB3TestTable$DekkerB3TestModel)
Table3[5,6] = rmse(DekkerB3TestTable$Secchi, DekkerB3TestTable$DekkerB3TestModel)
Table3[5,7] = 1 - (sum((DekkerB3TestTable$DekkerB3TestModel - DekkerB3TestTable$Secchi) ^ 2)/sum((DekkerB3TestTable$Secchi - mean(DekkerB3TestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
DekkerB3FigureTrain <- ggplot(DekkerB3TrainTable, aes(Secchi, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dekker/Peters 2", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DekkerB3FigureTest <- ggplot(DekkerB3TestTable, aes(Secchi, DekkerB3TestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dekker/Peters 2", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DekkerB3FigureTrain
DekkerB3FigureTest

#Add columns for predicted/observed summary
Table3[5,9] = summary(lm(DekkerB3TrainTable$fitted.values~DekkerB3TrainTable$Secchi))$coefficients[2,1]
Table3[5,10] = summary(lm(DekkerB3TestTable$DekkerB3TestModel~DekkerB3TestTable$Secchi))$coefficients[2,1]
Table3

#################Dominguez Gomez et al.###########################################
#create data frame and calculate needed values for algorithm
DominguezGomezTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B2 )
DominguezGomezTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B2)

#calculate constants for algorithm based on training data
DominguezGomezModel <- nls(secchi_depth_m ~ b*(B2^z), data = DominguezGomezTrain, start = list(b = 1, z = 0))
DominguezGomezTrainTable <- data.table(DominguezGomezTrain$secchi_depth_m, fitted(DominguezGomezModel))
colnames(DominguezGomezTrainTable) = c('Secchi', 'fitted.values')
DominguezGomezTrainTable

#record training statistics
Table3[6,1] <- "DominguezGomez et al."
Table3[6,2] <- mae(DominguezGomezTrain$secchi_depth_m , fitted(DominguezGomezModel))
Table3[6,3] <- rmse(DominguezGomezTrain$secchi_depth_m , fitted(DominguezGomezModel))
Table3[6,4] <-1 - (sum((DominguezGomezTrainTable$fitted.values - DominguezGomezTrainTable$Secchi) ^ 2)/sum((DominguezGomezTrainTable$Secchi - mean(DominguezGomezTrainTable$Secchi)) ^ 2)) #R2 in NLE isn't valid, but adding in anyway
Table3[6,8] = summary(DominguezGomezModel)$coefficients[2,4] #maybe we should drop p-value too?

#apply model to test data
DominguezGomezTestModel <- predict(DominguezGomezModel, newdata = DominguezGomezTest)
DominguezGomezTestTable <- data.table(DominguezGomezTestModel)
DominguezGomezTestTable$Secchi <- DominguezGomezTest$secchi_depth_m
DominguezGomezTestTable

#record test statistics
Table3[6,5] = mae(DominguezGomezTestTable$Secchi, DominguezGomezTestTable$DominguezGomezTestModel)
Table3[6,6] = rmse(DominguezGomezTestTable$Secchi, DominguezGomezTestTable$DominguezGomezTestModel)
Table3[6,7] = 1 - (sum((DominguezGomezTestTable$DominguezGomezTestModel - DominguezGomezTestTable$Secchi) ^ 2)/sum((DominguezGomezTestTable$Secchi - mean(DominguezGomezTestTable$Secchi)) ^ 2)) #R2 in NLE isn't valid
Table3

#create and save plots of observed-predicted
DominguezGomezFigureTrain <- ggplot(DominguezGomezTrainTable, aes(Secchi,fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dominguez Gomez", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DominguezGomezFigureTest <- ggplot(DominguezGomezTestTable, aes(Secchi,DominguezGomezTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dominguez Gomez", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

DominguezGomezFigureTrain
DominguezGomezFigureTest

#Add columns for predicted/observed summary
Table3[6,9] = summary(lm(DominguezGomezTrainTable$fitted.values~DominguezGomezTrainTable$Secchi))$coefficients[2,1]
Table3[6,10] = summary(lm(DominguezGomezTestTable$DominguezGomezTestModel~DominguezGomezTestTable$Secchi))$coefficients[2,1]
Table3


# #do reality check by doing transformations and treating like a linear equation
# 
# #create new DF and calculate needed variables
# DGcheck <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B2)
# DGcheck$lnSecchi = log(DGcheck$secchi_depth_m)
# DGcheck$lnB2 = log(DGcheck$B2)
# 
# #run DG on training data
# DGcheck_alg = lm(lnSecchi ~ lnB2, data = DGcheck)
# summary(DGcheck_alg)
# #intercept: 3.35
# #slope: -0.31
# #as a linear equation:
# # ln(secchi)= 3.35 - 0.31*ln(B2)
# 
# #transformed to a power equation:
# # secchi = exp(3.34)*B2^-0.31
# # secchi = 28.22*B2^-0.31
# 
# summary(DominguezGomezModel)
# # from nls:
# # secchi = 15.9*B2^-0.19
# 
# DGcheck$lnfitted = predict(DGcheck_alg, DGcheck)
# 
# ggplot(DGcheck, aes(exp(lnSecchi),exp(lnfitted))) + 
#   geom_bin2d(bins = 50 ) +
#   scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
#   theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Dominguez Gomez", size = 8) +
#   geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
#   xlim(0,20) + ylim(0,20)
# 
# lm_b2secchi <- ggplot(DGcheck, aes(x = exp(lnB2), y = exp(lnSecchi))) +
#   geom_point() +
#   geom_line(aes(x = exp(lnB2), y = exp(lnfitted)), size = 1, color = 'red') +
#   labs(x = 'B2', y = 'Secchi depth (m), observed', title = 'transformed linear model')
# 
# nls_b2secchi <- ggplot(DominguezGomezTrain, aes(x = B2, y = secchi_depth_m)) +
#   geom_point() +
#   geom_line(aes(x = B2, y = 15.9*(B2^-0.19)), color = 'blue', size =1)+
#   labs(x = 'B2', y = 'Secchi depth (m), observed', title = 'nls model')
# 
# grid.arrange(nls_b2secchi, lm_b2secchi)
# 
# nls.model.df <- data.frame(x.nls = DominguezGomezTrainTable$Secchi,
#                            y.nls = DominguezGomezTrainTable$fitted.values) %>% 
#   rowid_to_column()
# 
# lm.model.df <- data.frame(x.lm = DGcheck$secchi_depth_m,
#                           y.lm = exp(DGcheck$lnfitted)) %>% 
#   rowid_to_column()
# 
# two.models <- full_join(nls.model.df, lm.model.df)
# 
# ggplot(two.models, aes(x=x.lm, y = x.nls)) +
#   geom_point()
# ggplot(two.models, aes(x=y.lm, y = y.nls)) +
#   geom_point() +
#   labs(x = 'predicted secchi from (log-transformed) linear model',
#        y = 'predicted secchi from nls model') +
#   geom_abline(intercept = 0, slope = 1, size = 1, color = "red") +
#   coord_cartesian(xlim = c(3, 10), ylim = c(3, 10))


#################Giardino et al.#################################################

#create data frame and calculate needed values for algorithm
GiardinoTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B1B2)
GiardinoTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B1B2)

#calculate constants for algorithm based on training data
GiardinoModel <- lm(secchi_depth_m ~ B1B2, data = GiardinoTrain)
GiardinoTrainTable <- data.table(GiardinoModel$model$secchi_depth_m, GiardinoModel$fitted.values)
colnames(GiardinoTrainTable) = c('secchi_depth_m', 'fitted.values')

#record training statistics
Table3[7,1] <- "Giardino et al."
Table3[7,2] <- mae(GiardinoModel$model$secchi_depth_m, GiardinoModel$fitted.values)
Table3[7,3] <- rmse(GiardinoModel$model$secchi_depth_m, GiardinoModel$fitted.values)
Table3[7,4] <-  1 - (sum((GiardinoTrainTable$fitted.values - GiardinoTrainTable$secchi_depth_m) ^ 2)/sum((GiardinoTrainTable$secchi_depth_m - mean(GiardinoTrainTable$secchi_depth_m)) ^ 2))
Table3[7,8] = glance(GiardinoModel)$p.value

#apply model to test data
GiardinoTestModel <- predict(GiardinoModel, newdata = GiardinoTest)
GiardinoTestTable <- data.table(GiardinoTestModel)
GiardinoTestTable$Secchi <- GiardinoTest$secchi_depth_m
GiardinoTestTable

#record test statistics
Table3[7,5] = mae(GiardinoTestTable$Secchi, GiardinoTestTable$GiardinoTestModel)
Table3[7,6] = rmse(GiardinoTestTable$Secchi, GiardinoTestTable$GiardinoTestModel)
Table3[7,7] =  1 - (sum((GiardinoTestTable$GiardinoTestModel - GiardinoTestTable$Secchi) ^ 2)/sum((GiardinoTestTable$Secchi - mean(GiardinoTestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
GiardinoFigureTrain <- ggplot(GiardinoTrainTable, aes(secchi_depth_m, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Giardino et al.", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

GiardinoFigureTest <- ggplot(GiardinoTestTable, aes(Secchi,GiardinoTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Giardino et al.", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)


GiardinoFigureTrain
GiardinoFigureTest

#Add columns for predicted/observed summary
Table3[7,9] = summary(lm(GiardinoTrainTable$fitted.values~GiardinoTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[7,10] = summary(lm(GiardinoTestTable$GiardinoTestModel~GiardinoTestTable$Secchi))$coefficients[2,1]
Table3


#################Kloiber et al.#################################################

#create data frame and calculate needed values for algorithm
KloiberTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B1B3, B1)
KloiberTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B1B3, B1)
KloiberTrain$lnSecchi <- log(KloiberTrain$secchi_depth_m)
KloiberTest$lnSecchi <- log(KloiberTest$secchi_depth_m)

#calculate constants for algorithm based on training data
KloiberModel <- lm(lnSecchi ~ B1B3 + B1, data = KloiberTrain)
KloiberTrainTable <- data.table(KloiberModel$model$lnSecchi, KloiberModel$fitted.values)
colnames(KloiberTrainTable) = c('lnSecchi', 'ln.fitted.values')
KloiberTrainTable

#record training statistics
Table3[8,1] <- "Kloiber et al."
Table3[8,2] <- mae(exp(KloiberModel$model$lnSecchi), exp(KloiberModel$fitted.values))
Table3[8,3] <- rmse(exp(KloiberModel$model$lnSecchi), exp(KloiberModel$fitted.values))
Table3[8,4] <-  1 - (sum((exp(KloiberTrainTable$ln.fitted.values) - exp(KloiberTrainTable$lnSecchi)) ^ 2)/sum((exp(KloiberTrainTable$lnSecchi) - mean(exp(KloiberTrainTable$lnSecchi))) ^ 2))
Table3[8,8] = glance(KloiberModel)$p.value

#apply model to test data
KloiberTestModel <- predict(KloiberModel, newdata = KloiberTest)
KloiberTestTable <- data.table(KloiberTestModel)
colnames(KloiberTestTable) = 'lnKloiberTestModel'
KloiberTestTable$lnSecchi <- log(KloiberTest$secchi_depth_m)
KloiberTestTable

#record test statistics
Table3[8,5] = mae((exp(KloiberTestTable$lnSecchi)), (exp(KloiberTestTable$lnKloiberTestModel)))
Table3[8,6] = rmse((exp(KloiberTestTable$lnSecchi)), (exp(KloiberTestTable$lnKloiberTestModel)))
Table3[8,7] <-  1 - (sum((exp(KloiberTestTable$lnKloiberTestModel) - exp(KloiberTestTable$lnSecchi)) ^ 2)/sum((exp(KloiberTestTable$lnSecchi) - mean(exp(KloiberTestTable$lnSecchi))) ^ 2))
Table3

#create and save plots of observed-predicted
KloiberFigureTrain <- ggplot(KloiberTrainTable, aes((exp(lnSecchi)), exp(ln.fitted.values))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Kloiber et al. **", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

KloiberFigureTest <- ggplot(KloiberTestTable, aes((exp(lnSecchi)), exp(lnKloiberTestModel))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Kloiber et al. **", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

KloiberFigureTrain
KloiberFigureTest

max(exp(KloiberTrainTable$ln.fitted.values))
range(exp(KloiberTestTable$lnKloiberTestModel))

#Kloiber also has astronomically high fitted values (and out of range for the test set, too), hence the terrible stats numbers
KloiberTrainTable$predicted <- exp(KloiberTrainTable$ln.fitted.values)

#Add columns for predicted/observed summary
Table3[8,9] = summary(lm(exp(KloiberTrainTable$ln.fitted.values)~exp(KloiberTrainTable$lnSecchi)))$coefficients[2,1]
Table3[8,10] = summary(lm(exp(KloiberTestTable$lnKloiberTestModel)~exp(KloiberTestTable$lnSecchi)))$coefficients[2,1]
Table3

#################Lathrop and Lillesand#################################################

#create data frame and calculate needed values for algorithm
LathropLillesandTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B2)
LathropLillesandTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B2)
LathropLillesandTrain$lnSecchi <- log(LathropLillesandTrain$secchi_depth_m)
LathropLillesandTest$lnSecchi <- log(LathropLillesandTest$secchi_depth_m)

#calculate constants for algorithm based on training data
LathropLillesandModel <- lm(lnSecchi ~ B2, data = LathropLillesandTrain)
LathropLillesandTrainTable <- data.table(LathropLillesandModel$model$lnSecchi, LathropLillesandModel$fitted.values)
colnames(LathropLillesandTrainTable) = c('lnSecchi', 'ln.fitted.values')
LathropLillesandTrainTable

#record training statistics
Table3[9,1] <- "Lathrop and Lillesand"
Table3[9,2] <- mae((exp(LathropLillesandModel$model$lnSecchi)), (exp(LathropLillesandModel$fitted.values)))
Table3[9,3] <- rmse((exp(LathropLillesandModel$model$lnSecchi)), (exp(LathropLillesandModel$fitted.values)))
Table3[9,4] <-  1 - (sum((exp(LathropLillesandTrainTable$ln.fitted.values) - exp(LathropLillesandTrainTable$lnSecchi)) ^ 2)/sum((exp(LathropLillesandTrainTable$lnSecchi) - mean(exp(LathropLillesandTrainTable$lnSecchi))) ^ 2))
Table3[9,8] = glance(LathropLillesandModel)$p.value

#apply model to test data
LathropLillesandTestModel <- predict(LathropLillesandModel, newdata = LathropLillesandTest)
LathropLillesandTestTable <- data.table(LathropLillesandTestModel)
colnames(LathropLillesandTestTable) = 'lnLathropLillesandTestModel'
LathropLillesandTestTable$lnSecchi <- log(LathropLillesandTest$secchi_depth_m)
LathropLillesandTestTable

#record test statistics
Table3[9,5] = mae((exp(LathropLillesandTestTable$lnSecchi)), (exp(LathropLillesandTestTable$lnLathropLillesandTestModel)))
Table3[9,6] = rmse((exp(LathropLillesandTestTable$lnSecchi)), (exp(LathropLillesandTestTable$lnLathropLillesandTestModel)))
Table3[9,7] <-  1 - (sum((exp(LathropLillesandTestTable$lnLathropLillesandTestModel) - exp(LathropLillesandTestTable$lnSecchi)) ^ 2)/sum((exp(LathropLillesandTestTable$lnSecchi) - mean(exp(LathropLillesandTestTable$lnSecchi))) ^ 2))
Table3

#create and save plots of observed-predicted
LathropLillesandFigureTrain <- ggplot(LathropLillesandTrainTable, aes((exp(lnSecchi)), exp(ln.fitted.values))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Lathrop/Lillesand", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

LathropLillesandFigureTest <- ggplot(LathropLillesandTestTable, aes((exp(lnSecchi)), exp(lnLathropLillesandTestModel))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Lathrop/Lillesand", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

LathropLillesandFigureTrain
LathropLillesandFigureTest

#Add columns for predicted/observed summary
Table3[9,9] = summary(lm(exp(LathropLillesandTrainTable$ln.fitted.values)~exp(LathropLillesandTrainTable$lnSecchi)))$coefficients[2,1]
Table3[9,10] = summary(lm(exp(LathropLillesandTestTable$lnLathropLillesandTestModel)~exp(LathropLillesandTestTable$lnSecchi)))$coefficients[2,1]
Table3

#################Lavery et al.#################################################

#create data frame and calculate needed values for algorithm
LaveryTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3, B1B3)
LaveryTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3, B1B3)

#calculate constants for algorithm based on training data
LaveryModel <- lm(secchi_depth_m ~ B3+B1B3, data = LaveryTrain )
LaveryTrainTable <- data.table(LaveryModel$model$secchi_depth_m, LaveryModel$fitted.values)
colnames(LaveryTrainTable) = c('secchi_depth_m', 'fitted.values')
LaveryTrainTable

#record training statistics
Table3[10,1] <- "Lavery et al."
Table3[10,2] <- mae(LaveryModel$model$secchi_depth_m, LaveryModel$fitted.values)
Table3[10,3] <- rmse(LaveryModel$model$secchi_depth_m, LaveryModel$fitted.values)
Table3[10,4] <- 1 - (sum((LaveryTrainTable$fitted.values - LaveryTrainTable$secchi_depth_m) ^ 2)/sum((LaveryTrainTable$secchi_depth_m - mean(LaveryTrainTable$secchi_depth_m)) ^ 2))
Table3[10,8] = glance(LaveryModel)$p.value

#apply model to test data
LaveryTestModel <- predict(LaveryModel, newdata = LaveryTest)
LaveryTestTable <- data.table(LaveryTestModel)
LaveryTestTable$Secchi <- LaveryTest$secchi_depth_m
LaveryTestTable

#record test statistics
Table3[10,5] = mae(LaveryTestTable$Secchi, LaveryTestTable$LaveryTestModel)
Table3[10,6] = rmse(LaveryTestTable$Secchi, LaveryTestTable$LaveryTestModel)
Table3[10,7] = 1 - (sum((LaveryTestTable$LaveryTestModel - LaveryTestTable$Secchi) ^ 2)/sum((LaveryTestTable$Secchi - mean(LaveryTestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
LaveryFigureTrain <- ggplot(LaveryTrainTable, aes(secchi_depth_m,fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Lavery et al. *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

LaveryFigureTest <- ggplot(LaveryTestTable, aes(Secchi,LaveryTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Lavery et al.", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

LaveryFigureTrain
LaveryFigureTest

#Lavery training creates negative secchi

#Add columns for predicted/observed summary
Table3[10,9] = summary(lm(LaveryTrainTable$fitted.values~LaveryTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[10,10] = summary(lm(LaveryTestTable$LaveryTestModel~LaveryTestTable$Secchi))$coefficients[2,1]
Table3

#################Mancino et al.#################################################

#create data frame and calculate needed values for algorithm
MancinoTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3B2, B1B2, B1)
MancinoTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3B2, B1B2, B1)

#calculate constants for algorithm based on training data
MancinoModel <- lm(secchi_depth_m ~ B3B2+B1B2+B1, data = MancinoTrain )
MancinoTrainTable <- data.table(MancinoModel$model$secchi_depth_m, MancinoModel$fitted.values)
colnames(MancinoTrainTable) = c('secchi_depth_m', 'fitted.values')
MancinoTrainTable

#record training statistics
Table3[11,1] <- "Mancino et al."
Table3[11,2] <- mae(MancinoModel$model$secchi_depth_m, MancinoModel$fitted.values)
Table3[11,3] <- rmse(MancinoModel$model$secchi_depth_m, MancinoModel$fitted.values)
Table3[11,4] <-  1 - (sum((MancinoTrainTable$fitted.values - MancinoTrainTable$secchi_depth_m) ^ 2)/sum((MancinoTrainTable$secchi_depth_m - mean(MancinoTrainTable$secchi_depth_m)) ^ 2))
Table3[11,8] = glance(MancinoModel)$p.value

#apply model to test data
MancinoTestModel <- predict(MancinoModel, newdata = MancinoTest)
MancinoTestTable <- data.table(MancinoTestModel)
MancinoTestTable$Secchi <- MancinoTest$secchi_depth_m
MancinoTestTable

#record test statistics
Table3[11,5] = mae(MancinoTestTable$Secchi, MancinoTestTable$MancinoTestModel)
Table3[11,6] = rmse(MancinoTestTable$Secchi, MancinoTestTable$MancinoTestModel)
Table3[11,7] =  1 - (sum((MancinoTestTable$MancinoTestModel - MancinoTestTable$Secchi) ^ 2)/sum((MancinoTestTable$Secchi - mean(MancinoTestTable$Secchi)) ^ 2))
Table3


#create and save plots of observed-predicted
MancinoFigureTrain <- ggplot(MancinoTrainTable, aes(secchi_depth_m,fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Mancino et al. *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)
MancinoFigureTest <- ggplot(MancinoTestTable, aes(Secchi,MancinoTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Mancino et al. *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

MancinoFigureTrain
MancinoFigureTest

#mancino training creates negaive secchi values
#mancino test creates negative secchi values

#Add columns for predicted/observed summary
Table3[11,9] = summary(lm(MancinoTrainTable$fitted.values~MancinoTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[11,10] = summary(lm(MancinoTestTable$MancinoTestModel~MancinoTestTable$Secchi))$coefficients[2,1]
Table3

#################Wu et al.#################################################

#create data frame and calculate needed values for algorithm
WuTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B3, B1)
WuTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B3, B1)
WuTrain$lnSecchi <- log(WuTrain$secchi_depth_m)
WuTest$lnSecchi <- log(WuTest$secchi_depth_m)

#calculate constants for algorithm based on training data
WuModel <- lm(lnSecchi ~ B3 + B1, data = WuTrain)
WuTrainTable <- data.table(WuModel$model$lnSecchi, WuModel$fitted.values)
colnames(WuTrainTable) = c('lnSecchi', 'ln.fitted.values')
WuTrainTable

#record training statistics
Table3[12,1] <- "Wu et al."
Table3[12,2] <- mae(exp(WuModel$model$lnSecchi), exp(WuModel$fitted.values))
Table3[12,3] <- rmse(exp(WuModel$model$lnSecchi), exp(WuModel$fitted.values))
Table3[12,4] <-  1 - (sum((exp(WuTrainTable$ln.fitted.values) - exp(WuTrainTable$lnSecchi)) ^ 2)/sum((exp(WuTrainTable$lnSecchi) - mean(exp(WuTrainTable$lnSecchi))) ^ 2))
Table3[12,8] = glance(WuModel)$p.value

#calculate constants for algorithm based on training data
WuTestModel <- predict(WuModel, newdata = WuTest)
WuTestTable <- data.table(WuTestModel)
colnames(WuTestTable) = 'lnWuTestModel'
WuTestTable$lnSecchi <- log(WuTest$secchi_depth_m)
WuTestTable

#record test statistics
Table3[12,5] = mae((exp(WuTestTable$lnSecchi)), (exp(WuTestTable$lnWuTestModel)))
Table3[12,6] = rmse((exp(WuTestTable$lnSecchi)), (exp(WuTestTable$lnWuTestModel)))
Table3[12,7] <-  1 - (sum((exp(WuTestTable$lnWuTestModel) - exp(WuTestTable$lnSecchi)) ^ 2)/sum((exp(WuTestTable$lnSecchi) - mean(exp(WuTestTable$lnSecchi))) ^ 2))
Table3

WuFigureTrain <- ggplot(WuTrainTable,aes((exp(lnSecchi)), exp(ln.fitted.values))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Wu et al.", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

WuFigureTest <- ggplot(WuTestTable, aes((exp(lnSecchi)), exp(lnWuTestModel))) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Wu et al.", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

WuFigureTrain
WuFigureTest

#Add columns for predicted/observed summary
Table3[12,9] = summary(lm(exp(WuTrainTable$ln.fitted.values)~exp(WuTrainTable$lnSecchi)))$coefficients[2,1]
Table3[12,10] = summary(lm(exp(WuTestTable$lnWuTestModel)~exp(WuTestTable$lnSecchi)))$coefficients[2,1]
Table3

#################Yip et al.#################################################

#create data frame and calculate needed values for algorithm
YipTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B4, B2, B1)
YipTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B4, B2, B1)

#calculate constants for algorithm based on training data
YipModel <- lm(secchi_depth_m ~ B4+B2+B1, data = YipTrain )
YipTrainTable <- data.table(YipModel$model$secchi_depth_m, YipModel$fitted.values)
colnames(YipTrainTable) = c('secchi_depth_m', 'fitted.values')
YipTrainTable

#record training statistics
Table3[13,1] <- "Yip et al."
Table3[13,2] <- mae(YipModel$model$secchi_depth_m, YipModel$fitted.values)
Table3[13,3] <- rmse(YipModel$model$secchi_depth_m, YipModel$fitted.values)
Table3[13,4] <-  1 - (sum((YipTrainTable$fitted.values - YipTrainTable$secchi_depth_m) ^ 2)/sum((YipTrainTable$secchi_depth_m - mean(YipTrainTable$secchi_depth_m)) ^ 2))
Table3[13,8] = glance(YipModel)$p.value

#apply model to test data
YipTestModel <- predict(YipModel, newdata = YipTest)
YipTestTable <- data.table(YipTestModel)
YipTestTable$Secchi <- YipTest$secchi_depth_m

#record test statistics
Table3[13,5] = mae(YipTestTable$Secchi, YipTestTable$YipTestModel)
Table3[13,6] = rmse(YipTestTable$Secchi, YipTestTable$YipTestModel)
Table3[13,7] = 1 - (sum((YipTestTable$YipTestModel - YipTestTable$Secchi) ^ 2)/sum((YipTestTable$Secchi - mean(YipTestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
YipFigureTrain <- ggplot(YipTrainTable, aes(secchi_depth_m, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
    scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Yip et al. *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

YipFigureTest <- ggplot(YipTestTable, aes(Secchi, YipTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Yip et al. *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

YipFigureTrain
YipFigureTest

#yip training and testing creates negative secchi

#Add columns for predicted/observed summary
Table3[13,9] = summary(lm(YipTrainTable$fitted.values~YipTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[13,10] = summary(lm(YipTestTable$YipTestModel~YipTestTable$Secchi))$coefficients[2,1]
Table3

########################FourBands###########################################

#create data frame and calculate needed values for algorithm
FourBandTrain <- select(training, rowid, date, time, secchi_depth_m, dissPermID, B1, B2, B3, B4)
FourBandTest <- select(testing, rowid, date, time, secchi_depth_m, dissPermID, B1, B2, B3, B4)

#calculate constants for algorithm based on training data
FourBandModel <- lm(secchi_depth_m ~ B1 + B2 + B3 + B4, data = FourBandTrain )
FourBandTrainTable <- data.table(FourBandModel$model$secchi_depth_m, FourBandModel$fitted.values) 
colnames(FourBandTrainTable) = c('secchi_depth_m', 'fitted.values')
FourBandTrainTable

#record training statistics
Table3[14,1] <- "Four Band"
Table3[14,2] <- mae(FourBandModel$model$secchi_depth_m, FourBandModel$fitted.values)
Table3[14,3] <- rmse(FourBandModel$model$secchi_depth_m, FourBandModel$fitted.values)
Table3[14,4] <-  1 - (sum((FourBandTrainTable$fitted.values - FourBandTrainTable$secchi_depth_m) ^ 2)/sum((FourBandTrainTable$secchi_depth_m - mean(FourBandTrainTable$secchi_depth_m)) ^ 2))
Table3[14,8] = glance(FourBandModel)$p.value

#apply model to test data
FourBandTestModel <- predict(FourBandModel, newdata = FourBandTest)
FourBandTestTable <- data.table(FourBandTestModel)
FourBandTestTable$Secchi <- FourBandTest$secchi_depth_m
FourBandTestTable

#record test statistics
Table3[14,5] = mae(FourBandTestTable$Secchi, FourBandTestTable$FourBandTestModel)
Table3[14,6] = rmse(FourBandTestTable$Secchi, FourBandTestTable$FourBandTestModel)
Table3[14,7] = 1 - (sum((FourBandTestTable$FourBandTestModel - FourBandTestTable$Secchi) ^ 2)/sum((FourBandTestTable$Secchi - mean(FourBandTestTable$Secchi)) ^ 2))
Table3

#create and save plots of observed-predicted
FourBandFigureTrain <- ggplot(FourBandTrainTable, aes(secchi_depth_m, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Four Band *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

FourBandFigureTest <- ggplot(FourBandTestTable, aes(Secchi, FourBandTestModel)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Four Band *", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

FourBandFigureTrain
FourBandFigureTest

#four band secchi creates negative secchi

#Add columns for predicted/observed summary
Table3[14,9] = summary(lm(FourBandTrainTable$fitted.values~FourBandTrainTable$secchi_depth_m))$coefficients[2,1]
Table3[14,10] = summary(lm(FourBandTestTable$FourBandTestModel~FourBandTestTable$Secchi))$coefficients[2,1]
Table3

########################Bring in RF Object ######################

RFtraining <- readRDS(paste0('Rubin.MS/RFobjects/HannahFullRF128trees.RDS'))

#Now, predict OOB values using the RF object
TestingPredict <- predict(RFtraining, testing, type="response", predict.all=FALSE)
TestingPredictData <- as.data.frame(TestingPredict)
testing$predicted <- TestingPredictData$TestingPredict

RFTrainingTable <- data.table(RFtraining$y, RFtraining$predicted)
colnames(RFTrainingTable) = c('Secchi', 'fitted.values')
RFTrainingTable

Table3[15,1] <- "Random Forest"
Table3[15,2] <- mae(RFTrainingTable$Secchi, RFTrainingTable$fitted.values)
Table3[15,3] <- rmse(RFTrainingTable$Secchi, RFTrainingTable$fitted.values)
Table3[15,4] <- 1 - (sum((RFTrainingTable$fitted.values - RFTrainingTable$Secchi) ^ 2)/sum((RFTrainingTable$Secchi - mean(RFTrainingTable$Secchi)) ^ 2))
Table3[15,5] = mae(testing$secchi_depth_m, testing$predicted)
Table3[15,6] = rmse(testing$secchi_depth_m, testing$predicted)
Table3[15,7] = 1 - (sum((testing$predicted - testing$secchi_depth_m) ^ 2)/sum((testing$secchi_depth_m - mean(testing$secchi_depth_m)) ^ 2))
Table3[15,8] = NA
Table3

#2d histogram plot build. Here we are building a similar figure as before, but with OOB data
FigureRFTraining <- ggplot(RFTrainingTable, aes(x = Secchi, fitted.values)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,1700)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Random Forest", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

FigureRFTesting <- ggplot(testing, aes(secchi_depth_m, predicted)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,250)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() + theme(text=element_text(size=22)) + annotate("text", hjust = 0, x=1, y = 18, label = "Random Forest", size = 8) +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

FigureRFTraining
FigureRFTesting

#Add columns for predicted/observed summary
Table3[15,9] = summary(lm(RFTrainingTable$fitted.values~RFTrainingTable$Secchi))$coefficients[2,1]
Table3[15,10] = summary(lm(testing$predicted~testing$secchi_depth_m))$coefficients[2,1]
Table3

#########################EXPORT TABLE################################
# all p-values are <0.05, so we can remove that column

Table3 <- Table3 %>% 
  select(-`Training p-value`) %>% 
  mutate_at(vars(2:9),
            ~ round(., digits = 2))
Table3

write_csv(Table3,paste0(dumpdir, "wholedataset_summary/Table3_v19Jan2021.csv"))

#######################MultiPanel Figure###########################

jpeg(paste0(dumpdir,"plots/FullDataScatterTraining_v19Jan2021.jpeg"), units="in", width=30, height=15, res=250)
grid.arrange(AlleeFigureTrain, BabanFigureTrain, ChipmanFigureTrain, DekkerB3FigureTrain, DekkerFigureTrain, DominguezGomezFigureTrain, GiardinoFigureTrain,
             KloiberFigureTrain, LathropLillesandFigureTrain, LaveryFigureTrain, MancinoFigureTrain, 
             WuFigureTrain, YipFigureTrain, FourBandFigureTrain, FigureRFTraining, ncol=5)
dev.off()


jpeg(paste0(dumpdir,"plots/FullDataScatterTesting_v19Jan2021.jpeg"), units="in", width=30, height=15, res=250)
grid.arrange(AlleeFigureTest, BabanFigureTest, ChipmanFigureTest, DekkerB3FigureTest, DekkerFigureTest, DominguezGomezFigureTest, GiardinoFigureTest,
             KloiberFigureTest, LathropLillesandFigureTest, LaveryFigureTest, MancinoFigureTest, 
             WuFigureTest, YipFigureTest, FourBandFigureTest, FigureRFTesting, ncol=5)
dev.off()
