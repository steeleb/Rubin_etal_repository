#Written Nov 2020 by David Lutz for Rubin et al. manuscript
#edited Jan 2021 for final v of MS

library(tidyverse)
# library(dplyr)
# library(tidyr)
# library(lubridate)
# library(ggplot2)
library(data.table)
# library(Metrics)
library(gridExtra)
library(randomForest)
library(readxl)

# create directory paths, commented out for us on Dave's machine
datadir <- '5.a.90.10.datasets/'
RDSdir <- 'Rubin.MS/RFobjects/'
dumpdir <- 'Rubin.MS/summary_data/'

# #local directory for Dave testing
# setwd("D:/NASA IDS/Manuscript 3 Rubin/DaveTestNov")

# Read in data
training = read_csv(paste0(datadir, "90percent_7.23.20.csv"),
                    col_types = cols(.default = col_character())) %>% 
  mutate_at(vars(secchi_depth_m, B1:B4),
            ~ as.numeric(.))
testing = read_csv(paste0(datadir, "OOB_7.23.20.csv"),
                   col_types = cols(.default = col_character())) %>% 
  mutate_at(vars(secchi_depth_m, B1:B4),
            ~ as.numeric(.))

# #rsq function
# rsq <- function(x,y) summary(lm(y~x))$r.squared
# 
# completeFun <- function(data, desiredCols) {
#   completeVec <- complete.cases(data[, desiredCols])
#   return(data[completeVec, ])
# }

######################################################################
#remove single sample from training set with very low B3 value (possibly erroneous but technically within confidence intervals)

training <- training %>% 
  filter(B3 >= 3)

#Random Forest Construction

RFtraining <- randomForest(secchi_depth_m ~ B1 + B2 + B3 + B4,
                                       data = training, ntree=128,
                                       na.action = na.roughfix, importance = TRUE)

saveRDS(RFtraining, file = paste0(RDSdir, 'HannahFullRF128trees_v11Jan2021.RDS')) #save the RF object

#appends RF build data with predicteds from RF model
training$predicted <- RFtraining$predicted

#2d histogram plot build. This is the main way we visualize RF performance
FigureRF <- ggplot(training, aes(secchi_depth_m, predicted)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "midnightblue", high = "cyan", limits = c(0,600)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

#Now, predict OOB values using the RF object
TestingPredict <- predict(RFtraining, testing, type="response", predict.all=FALSE)
TestingPredictData <- as.data.frame(TestingPredict)
testing$predicted <- TestingPredictData$TestingPredict

#2d histogram plot build. Here we are building a similar figure as before, but with OOB data
FigureRFTesting <- ggplot(testing, aes(secchi_depth_m, predicted)) + 
  geom_bin2d(bins = 50 ) +
  scale_fill_gradient(low = "tomato4", high = "darkgoldenrod1", limits = c(0,100)) + xlab("Measured Secchi Depth (m)") + ylab("Predicted Secchi Depth (m)") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, size = 1, color = "red") + 
  xlim(0,20) + ylim(0,20)

#place both plots, 90% and 10% side by side using grid.arrange
RF_fig <- grid.arrange(FigureRF, FigureRFTesting, ncol=2)

# Exports Figure as a high-res tiff. This can be modified to jpeg via replacing 'tiff'. The tiff will
# be of high quality for publication, but for checking subtleties, jpeg is faster and smaller.

jpeg("Rubin.MS/summary_data/plots/Figure2_v11Jan2021.jpeg", units="in", width=12, height=5, res=400)
grid.arrange(FigureRF, FigureRFTesting, ncol=2)
dev.off()

