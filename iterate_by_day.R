#file originally written by Hannah Rubin June 2019
#edited for reproducibility and train/test update Nov 2019-Jan 2021 B. Steele

library(readxl)
library(randomForest)
library(lubridate)
library(tidyverse)
library(Metrics)
library(ggthemes)
library(broom)
# library(dplyr)

#save groms for ggplot
final_theme=theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=16, face='bold', hjust=0.5)) #save as a grom


# create directory paths
datadir <- '5.a.90.10.datasets/'
RFobjectdir <- 'Rubin.MS/RFobjects/'
dumpdir <- 'Rubin.MS/summary_data/daybyday_summary/'
figdir <- 'Rubin.MS/summary_data/plots/'

# Read in secchi and RS data
training <- read_csv(paste0(datadir, '90percent_7.23.20.csv'),
                     col_types = cols(.default = col_character())) %>% 
  select(rowid:date, secchi_depth_m, lat_dd, long_dd, dissPermID, sat_date:Difference, B1B3:B2B1) %>% 
  mutate(dataset = 'train',
         sat_date_doy = yday(as.Date(date)))%>% 
  mutate_at(vars(secchi_depth_m, B1:Difference, B1B3:B2B1),
            ~as.numeric(.))
test <- read_csv(paste0(datadir, 'OOB_7.23.20.csv'),
                 col_types = cols(.default = col_character())) %>% 
  select(rowid:date, secchi_depth_m, lat_dd, long_dd, dissPermID, sat_date:Difference, B1B3:B2B1) %>% 
  mutate(dataset = 'test',
         sat_date_doy = yday(as.Date(date))) %>% 
  mutate_at(vars(secchi_depth_m, B1:Difference, B1B3:B2B1),
            ~as.numeric(.))
data <- full_join(training, test)

#remove the outlier in B3
data <- data %>% 
  filter(B3 >3)

#count number of unique image dates
count_train <- data %>% 
  filter(dataset == 'train') %>% 
  count(., sat_date) %>% 
  rename(n_train = n)

#filter to only include dates with at least 75 measurements
count_train <- count_train[count_train$n_train >= 75, ]

#filter and count test data matched with training dataset with >= 75 obs
count_test <-  data %>% 
  left_join(count_train, .) %>% 
  filter(dataset == 'test') %>% 
  count(., sat_date) %>% 
  rename(n_test = n)


#make sure date is in correct format and create a new dataframe
match_df <- full_join(count_train, count_test) %>% 
  mutate(sat_date = as.Date(sat_date))
match_df

# add stats of the observed dataset #
train_stats <-  left_join(count_train, training) %>% 
  group_by(sat_date, n_train) %>% 
  summarise(min_SD_obs_train = round(min(secchi_depth_m), digits = 1),
            mean_SD_obs_train = round(mean(secchi_depth_m), digits = 1),
            median_SD_obs_train = round(median(secchi_depth_m), digits = 1),
            max_SD_obs_train = round(max(secchi_depth_m), digits = 1),
            sd_SD_obs_train = round(sd(secchi_depth_m), digits = 1))
test_stats <-  left_join(count_test, test) %>% 
  group_by(sat_date, n_test) %>% 
  summarise(min_SD_obs_test = round(min(secchi_depth_m), digits = 1),
            mean_SD_obs_test = round(mean(secchi_depth_m), digits = 1),
            median_SD_obs_test = round(median(secchi_depth_m), digits = 1),
            max_SD_obs_test = round(max(secchi_depth_m), digits = 1),
            sd_SD_obs_test = round(sd(secchi_depth_m), digits = 1))

train_test_dbd_obs_stats <- full_join(train_stats, test_stats)

# write_csv(train_test_dbd_obs_stats, 'Rubin.MS/summary_data/wholedataset_summary/day_by_day_insitu_summarystats.csv')

#make histograms for each sat day
data_filter <- data %>% 
  mutate(sat_date = as.Date(sat_date)) %>% 
  left_join(match_df, .)

for (i in 1:nrow(match_df)) {
  histogram <- data_filter %>% 
    filter(sat_date == match_df$sat_date[i]) %>% 
    ggplot(., aes(x = secchi_depth_m)) +
    geom_histogram(bins = 30) +
    facet_grid(dataset ~ ., scales = 'free_y') +
    labs(x = 'Secchi depth (m)',
         y = 'frequency')  +
    theme_bw() +
    final_theme
  histogram
  ggsave(paste0('Rubin.MS/summary_data/plots/day_by_day/histograms/train_test_histogram_',match_df$sat_date[i],'.jpg'), units = 'in', height = 5, width = 5, dpi = 250)
}

match_template <- match_df %>% 
  mutate(mae_train = as.numeric(''),
         rmse_train = as.numeric(''),
         R2_train = as.numeric(''),
         mae_test = as.numeric(''),
         rmse_test = as.numeric(''),
         R2_test = as.numeric(''),
         pval_train = as.numeric(''),
         po_slope_train = as.numeric(''),
         po_slope_test = as.numeric('')
  )

########################### Random Forest ##############################

#create additional fields for stats reporting
stats_match_df_rf <- match_template

#iterate through unique dates and run random forest for each date then export to a file
for(i in 1:nrow(stats_match_df_rf)){
  date1 = as.Date(stats_match_df_rf$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]

  #run RF on training data
  output1 = randomForest(secchi_depth_m ~ B1 + B2 + B3 + B4,
                         data = dframe_train,
                         ntree = 128,
                         na.action = na.roughfix,
                         importance = TRUE)
  saveRDS(output1, paste0(RFobjectdir, 'day_by_day_RF_', date1, '.RDS')) #save the alg object
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_rf dataframe
  stats_match_df_rf$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p)  #MAE
  stats_match_df_rf$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_rf$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_rf$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  
  #calculate R2
  rss <- sum((dframe_test$dbd_p - dframe_test$secchi_depth_m) ^ 2)
  tss <- sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)
  stats_match_df_rf$R2_test[i] <- 1 - rss/tss  
  #calculate the R2 for train set
  rss <- sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)
  tss <- sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2)
  stats_match_df_rf$R2_train[i] <- 1 - rss/tss
  
  #calculate the slope of predicted and observed
  stats_match_df_rf$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_rf$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Random Forest day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/rf/', 'day_by_day_rf_', date1, '.jpg'))
}


stats_match_df_rf$`Model Name` = 'Random Forest'
write_csv(stats_match_df_rf, paste0(dumpdir, 'day_by_day_rf_summarystats_v', Sys.Date(), '.csv'))


########################### FOUR BAND LINEAR ##############################

#re-initialize the summary table
stats_match_df_fourlinear <- match_template

#create empty object to save models to
fourlinear_dbd_models <- NULL

#iterate through unique dates and run fourlinear alg for each date then export to a file
for(i in 1:nrow(stats_match_df_fourlinear)){
  date1 = as.Date(stats_match_df_fourlinear$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run fourlinear on training data
  output1 = lm(secchi_depth_m ~ B1 + B2 + B3 + B4, data = dframe_train )
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_fourlinear dataframe
  stats_match_df_fourlinear$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p)  #MAE
  stats_match_df_fourlinear$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_fourlinear$R2_test[i] =  1 - (sum((dframe_test$dbd_p - dframe_test$secchi_depth_m) ^ 2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2))#R-squared
  stats_match_df_fourlinear$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_fourlinear$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_fourlinear$R2_train[i] = 1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_fourlinear$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_fourlinear$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_fourlinear$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  fourlinear_dbd_models[[i]] <- summary(output1)
  fourlinear_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Four Bands linear, 1992 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/fourlinear/', 'day_by_day_fourlinear_', date1, '.jpg'))
}


stats_match_df_fourlinear$`Model Name` = 'fourlinear'
write_csv(stats_match_df_fourlinear, paste0(dumpdir, 'day_by_day_fourlinear_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(fourlinear_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_fourlinear_modelsummary_v', Sys.Date(), '.RDS'))


###########################ALLEE AND JOHNSON ##############################
#Allee and Johnson, 1999

#re-initialize the summary table
stats_match_df_allee <- match_template

#create empty object to save models to
allee_dbd_models <- NULL

#iterate through unique dates and run allee alg for each date then export to a file
for(i in 1:nrow(stats_match_df_allee)){
  date1 = as.Date(stats_match_df_allee$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #create necessary values and columns for allee
  B3mean_train <- mean(dframe_train[["B3"]])
  dframe_train$B3dev <- (dframe_train$B3 - B3mean_train)
  dframe_test$B3dev <- (dframe_test$B3 - B3mean_train)
  
  #run allee on training data
  output1 = lm(secchi_depth_m ~ poly(B3dev, 3), data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_allee dataframe
  stats_match_df_allee$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p)  #MAE
  stats_match_df_allee$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_allee$R2_test[i] = 1 - (sum((dframe_test$dbd_p - dframe_test$secchi_depth_m) ^ 2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2))#R-squared
  stats_match_df_allee$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_allee$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_allee$R2_train[i] =   1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_allee$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_allee$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_allee$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  
  allee_dbd_models[[i]] <- summary(output1)
  allee_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Allee and Johnson, 1999 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/allee/', 'day_by_day_allee_', date1, '.jpg'))
  
}

stats_match_df_allee$`Model Name` = 'Allee and Johnson'
write_csv(stats_match_df_allee, paste0(dumpdir, 'day_by_day_allee_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(allee_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_allee_modelsummary_v', Sys.Date(), '.RDS'))

########################### BABAN ##############################

#re-initialize the summary table
stats_match_df_baban <- match_template

#create empty object to save models to
baban_dbd_models <- NULL

#iterate through unique dates and run baban alg for each date then export to a file
for(i in 1:nrow(stats_match_df_baban)){
  date1 = as.Date(stats_match_df_baban$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run baban on training data
  output1 = lm(secchi_depth_m ~ B1, data = dframe_train )
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_baban dataframe
  stats_match_df_baban$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p)  #MAE
  stats_match_df_baban$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_baban$R2_test[i] =  1 - (sum((dframe_test$dbd_p - dframe_test$secchi_depth_m) ^ 2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2))#R-squared
  stats_match_df_baban$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_baban$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_baban$R2_train[i] =  1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_baban$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_baban$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_baban$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  baban_dbd_models[[i]] <- summary(output1)
  baban_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Baban, 1992 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/baban/', 'day_by_day_baban_', date1, '.jpg'))
}


stats_match_df_baban$`Model Name` = 'Baban'
write_csv(stats_match_df_baban, paste0(dumpdir, 'day_by_day_baban_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(baban_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_baban_modelsummary_v', Sys.Date(), '.RDS'))


################# CHIPMAN #################################################

#re-initialize the summary table
stats_match_df_chipman <- match_template

#create empty object to save models to
chipman_dbd_models <- NULL

#iterate through unique dates and run chipman alg for each date then export to a file
for(i in 1:nrow(stats_match_df_chipman)){
  date1 = as.Date(stats_match_df_chipman$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #calculate needed variables
  dframe_train$lnSecchi = log(dframe_train$secchi_depth_m)
  dframe_test$lnSecchi = log(dframe_test$secchi_depth_m)
  
  #run chipman on training data
  output1 = lm(lnSecchi ~ B1B3, data = dframe_train)  
  
  #predict the data based on alg object, save in dataframe
  dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_chipman dataframe
  stats_match_df_chipman$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  stats_match_df_chipman$rmse_test[i] = rmse(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #RMSE
  stats_match_df_chipman$R2_test[i] = 1-(sum((exp(dframe_test$ln_dbd_p) - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_chipman$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  stats_match_df_chipman$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  stats_match_df_chipman$R2_train[i] =  1 - (sum((exp(dframe_train$ln_dbd_p) - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_chipman$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_chipman$po_slope_train[i] <- summary(lm(exp(dframe_train$ln_dbd_p) ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_chipman$po_slope_test[i]<- summary(lm(exp(dframe_test$ln_dbd_p) ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  chipman_dbd_models[[i]] <- summary(output1)
  chipman_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$ln_dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = exp(ln_dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Chipman, et al., 2004 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/chipman/', 'day_by_day_chipman_', date1, '.jpg'))
  
}

stats_match_df_chipman$`Model Name` = 'Chipman, et al.'
write_csv(stats_match_df_chipman, paste0(dumpdir, 'day_by_day_chipman_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(chipman_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_chipman_modelsummary_v', Sys.Date(), '.RDS'))


#################DEKKER and PETERS, ln(B3)#################################################

#re-initialize the summary table
stats_match_df_dekker <- match_template

#create empty object to save models to
dekker_dbd_models <- NULL

#iterate through unique dates and run dekker alg for each date then export to a file
for(i in 1:nrow(stats_match_df_dekker)){
  date1 = as.Date(stats_match_df_dekker$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #create new variables needed for dekker
  dframe_train$lnSecchi <- log(dframe_train$secchi_depth_m)
  dframe_test$lnSecchi <- log(dframe_test$secchi_depth_m)
  dframe_train$lnB3 <- log(dframe_train$B3)
  dframe_test$lnB3 <- log(dframe_test$B3)
  
  #run dekker on training data
  output1 = lm(lnSecchi ~ lnB3, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_dekker dataframe
  stats_match_df_dekker$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  stats_match_df_dekker$rmse_test[i] = rmse(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #RMSE
  stats_match_df_dekker$R2_test[i] = 1-(sum((exp(dframe_test$ln_dbd_p) - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_dekker$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  stats_match_df_dekker$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  stats_match_df_dekker$R2_train[i] =  1 - (sum((exp(dframe_train$ln_dbd_p) - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_dekker$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_dekker$po_slope_train[i] <- summary(lm(exp(dframe_train$ln_dbd_p) ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_dekker$po_slope_test[i]<- summary(lm(exp(dframe_test$ln_dbd_p) ~ dframe_test$secchi_depth_m))$coefficients[2,1]

  
  dekker_dbd_models[[i]] <- summary(output1)
  dekker_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(exp(alldata_dbd$ln_dbd_p)))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = exp(ln_dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Dekker and Peters, 1993 (ln equation) day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/dekker/', 'day_by_day_dekker_', date1, '.jpg'))
}

stats_match_df_dekker$`Model Name` = 'Dekker and Peters 1'
write_csv(stats_match_df_dekker, paste0(dumpdir, 'day_by_day_dekker_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(dekker_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_dekker_modelsummary_v', Sys.Date(), '.RDS'))



#################DEKKER and PETERS, (B3)#################################################

#re-initialize the summary table
stats_match_df_dekker2 <-match_template

#create empty object to save models to
dekker2_dbd_models <- NULL

#iterate through unique dates and run dekker2 alg for each date then export to a file
for(i in 1:nrow(stats_match_df_dekker2)){
  date1 = as.Date(stats_match_df_dekker2$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run dekker2 on training data
  output1 = lm(secchi_depth_m ~ B3, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_dekker2 dataframe
  stats_match_df_dekker2$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p) #MAE
  stats_match_df_dekker2$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_dekker2$R2_test[i] = 1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_dekker2$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_dekker2$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_dekker2$R2_train[i] =  1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_dekker2$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_dekker2$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_dekker2$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  dekker2_dbd_models[[i]] <- summary(output1)
  dekker2_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Dekker and Peters, 1993 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/dekker2/', 'day_by_day_dekker2_', date1, '.jpg'))
}

stats_match_df_dekker2$`Model Name` = 'Dekker and Peters 2'
write_csv(stats_match_df_dekker2, paste0(dumpdir, 'day_by_day_dekker2_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(dekker2_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_dekker2_modelsummary_v', Sys.Date(), '.RDS'))


#################DOMINGUEZ GOMEZ###########################################

#re-initialize the summary table
stats_match_df_DG <- match_template

#create empty object to save models to
DG_dbd_models <- NULL

#iterate through unique dates and run DG alg for each date then export to a file
for(i in 1:nrow(stats_match_df_DG)){
  date1 = as.Date(stats_match_df_DG$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  # transformed to linear
  # #calculate needed variables
  # dframe_train$lnSecchi = log(dframe_train$secchi_depth_m)
  # dframe_test$lnSecchi = log(dframe_test$secchi_depth_m)
  # dframe_train$lnB2 = log(dframe_train$B2)
  # dframe_test$lnB2 = log(dframe_test$B2)
  # 
  # #run DG on training data
  # output1 = lm(lnSecchi ~ lnB2, data = dframe_train)
  # 
  # #predict the data based on alg object, save in dataframe
  # dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  # dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  # 
  # alldata_dbd <- full_join(dframe_train, dframe_test)
  # 
  # #save stats in stats_match_df_DG dataframe
  # stats_match_df_DG$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  # stats_match_df_DG$rmse_test[i] = rmse(dframe_test$secchi_depth_m,exp(dframe_test$ln_dbd_p)) #RMSE
  # stats_match_df_DG$R2_test[i] = 1-(sum((dframe_test$ln_dbd_p - dframe_test$lnSecchi)^2)/sum((dframe_test$lnSecchi - mean(dframe_test$lnSecchi)) ^ 2)) #R-squared
  # stats_match_df_DG$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  # stats_match_df_DG$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  # stats_match_df_DG$R2_train[i] = summary(output1)$r.squared #R-squared
  # stats_match_df_DG$pval_train[i] = glance(output1)$p.value #p-value
  
  # calculated as power function
  #run DG on training data
  output1 = nls(secchi_depth_m ~ b*(B2^z), data = dframe_train, start = list(b = 1, z = 0))

  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))

  alldata_dbd <- full_join(dframe_train, dframe_test)

  #save stats in stats_match_df_DG dataframe
  stats_match_df_DG$mae_test[i] = mae(dframe_test$secchi_depth_m, (dframe_test$dbd_p)) #MAE
  stats_match_df_DG$rmse_test[i] = rmse(dframe_test$secchi_depth_m,(dframe_test$dbd_p)) #RMSE
  stats_match_df_DG$R2_test[i] = 1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_DG$mae_train[i] = mae(dframe_train$secchi_depth_m, (dframe_train$dbd_p))  #MAE
  stats_match_df_DG$rmse_train[i] = rmse(dframe_train$secchi_depth_m, (dframe_train$dbd_p)) #RMSE
  stats_match_df_DG$R2_train[i] = 1-(sum((dframe_train$dbd_p - dframe_train$secchi_depth_m)^2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_DG$pval_train[i] = summary(output1)$coefficients[2,4] #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_DG$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_DG$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  DG_dbd_models[[i]] <- summary(output1)
  DG_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = (dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Dominguez Gomez, et al, 2009 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/DG/', 'day_by_day_DG_', date1, '.jpg'))
  
}

stats_match_df_DG$`Model Name` = 'Dominguez Gomez, et al.'
write_csv(stats_match_df_DG, paste0(dumpdir, 'day_by_day_DG_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(DG_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_DG_modelsummary_v', Sys.Date(), '.RDS'))


#################GIARDINO et al.#################################################

#re-initialize the summary table
stats_match_df_giardino <-match_template

#create empty object to save models to
giardino_dbd_models <- NULL

#iterate through unique dates and run giardino alg for each date then export to a file
for(i in 1:nrow(stats_match_df_giardino)){
  date1 = as.Date(stats_match_df_giardino$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run giardino on training data
  output1 = lm(secchi_depth_m ~ B1B2, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_giardino dataframe
  stats_match_df_giardino$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p) #MAE
  stats_match_df_giardino$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_giardino$R2_test[i] =  1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_giardino$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_giardino$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_giardino$R2_train[i] =  1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_giardino$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_giardino$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_giardino$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  giardino_dbd_models[[i]] <- summary(output1)
  giardino_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Giardino, et al., 2001 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/giardino/', 'day_by_day_giardino_', date1, '.jpg'))
}

stats_match_df_giardino$`Model Name` = 'Giardino, et al.'
write_csv(stats_match_df_giardino, paste0(dumpdir, 'day_by_day_giardino_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(giardino_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_giardino_modelsummary_v', Sys.Date(), '.RDS'))


#################KLOIBER et al.#################################################

#re-initialize the summary table
stats_match_df_kloiber <- match_template

#create empty object to save models to
kloiber_dbd_models <- NULL

#iterate through unique dates and run kloiber alg for each date then export to a file
for(i in 1:nrow(stats_match_df_kloiber)){
  date1 = as.Date(stats_match_df_kloiber$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #create new variables needed for kloiber
  dframe_train$lnSecchi <- log(dframe_train$secchi_depth_m)
  dframe_test$lnSecchi <- log(dframe_test$secchi_depth_m)
  
  #run kloiber on training data
  output1 = lm(lnSecchi ~ B1B2 + B1, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_kloiber dataframe
  stats_match_df_kloiber$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  stats_match_df_kloiber$rmse_test[i] = rmse(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #RMSE
  stats_match_df_kloiber$R2_test[i] = 1-(sum((exp(dframe_test$ln_dbd_p) - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_kloiber$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  stats_match_df_kloiber$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  stats_match_df_kloiber$R2_train[i] =  1 - (sum((exp(dframe_train$ln_dbd_p) - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_kloiber$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_kloiber$po_slope_train[i] <- summary(lm(exp(dframe_train$ln_dbd_p) ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_kloiber$po_slope_test[i]<- summary(lm(exp(dframe_test$ln_dbd_p) ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  
  kloiber_dbd_models[[i]] <- summary(output1)
  kloiber_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              exp(max(alldata_dbd$ln_dbd_p)))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = exp(ln_dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Kloiber, et al., 2002 day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/kloiber/', 'day_by_day_kloiber_', date1, '.jpg'))
}

stats_match_df_kloiber$`Model Name` = 'Kloiber, et al.'
write_csv(stats_match_df_kloiber, paste0(dumpdir, 'day_by_day_kloiber_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(kloiber_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_kloiber_modelsummary_v', Sys.Date(), '.RDS'))




#################Lathrop and Lillesand#################################################

#re-initialize the summary table
stats_match_df_lathrop3 <- match_template


#create empty object to save models to
lathrop3_dbd_models <- NULL

#iterate through unique dates and run lathrop3 alg for each date then export to a file
for(i in 1:nrow(stats_match_df_lathrop3)){
  date1 = as.Date(stats_match_df_lathrop3$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #create new variables needed for lathrop3
  dframe_train$lnSecchi <- log(dframe_train$secchi_depth_m)
  dframe_test$lnSecchi <- log(dframe_test$secchi_depth_m)
  
  #run lathrop3 on training data
  output1 = lm(lnSecchi ~ B2, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_lathrop3 dataframe
  stats_match_df_lathrop3$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  stats_match_df_lathrop3$rmse_test[i] = rmse(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #RMSE
  stats_match_df_lathrop3$R2_test[i] = 1-(sum((exp(dframe_test$ln_dbd_p) - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_lathrop3$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  stats_match_df_lathrop3$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  stats_match_df_lathrop3$R2_train[i] =  1 - (sum((exp(dframe_train$ln_dbd_p) - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_lathrop3$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_lathrop3$po_slope_train[i] <- summary(lm(exp(dframe_train$ln_dbd_p) ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_lathrop3$po_slope_test[i]<- summary(lm(exp(dframe_test$ln_dbd_p) ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  
  lathrop3_dbd_models[[i]] <- summary(output1)
  lathrop3_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              exp(max(alldata_dbd$ln_dbd_p)))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = exp(ln_dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Lathrop and Lilliesand day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/lathrop3/', 'day_by_day_lathrop3_', date1, '.jpg'))
}

stats_match_df_lathrop3$`Model Name` = 'Lathrop and Liliesand'
write_csv(stats_match_df_lathrop3, paste0(dumpdir, 'day_by_day_lathrop3_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(lathrop3_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_lathrop3_modelsummary_v', Sys.Date(), '.RDS'))

#################Lavery et al.#################################################

#re-initialize the summary table
stats_match_df_lavery <- match_template

#create empty object to save models to
lavery_dbd_models <- NULL

#iterate through unique dates and run lavery alg for each date then export to a file
for(i in 1:nrow(stats_match_df_lavery)){
  date1 = as.Date(stats_match_df_lavery$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run lavery on training data
  output1 = lm(secchi_depth_m ~ B3 + B1B3, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_lavery dataframe
  stats_match_df_lavery$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p) #MAE
  stats_match_df_lavery$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_lavery$R2_test[i] =1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_lavery$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_lavery$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_lavery$R2_train[i] = 1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_lavery$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_lavery$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_lavery$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  lavery_dbd_models[[i]] <- summary(output1)
  lavery_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Lavery, et al. day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/lavery/', 'day_by_day_lavery_', date1, '.jpg'))
}

stats_match_df_lavery$`Model Name` = 'Lavery, et al.'
write_csv(stats_match_df_lavery, paste0(dumpdir, 'day_by_day_lavery_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(lavery_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_lavery_modelsummary_v', Sys.Date(), '.RDS'))

#################Mancino et al.#################################################

#re-initialize the summary table
stats_match_df_mancino <- match_template

#create empty object to save models to
mancino_dbd_models <- NULL

#iterate through unique dates and run mancino alg for each date then export to a file
for(i in 1:nrow(stats_match_df_mancino)){
  date1 = as.Date(stats_match_df_mancino$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run mancino on training data
  output1 = lm(secchi_depth_m ~ B3B2 + B1B2 + B1, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_mancino dataframe
  stats_match_df_mancino$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p) #MAE
  stats_match_df_mancino$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_mancino$R2_test[i] = 1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_mancino$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_mancino$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_mancino$R2_train[i] =  1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_mancino$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_mancino$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_mancino$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  mancino_dbd_models[[i]] <- summary(output1)
  mancino_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Mancino, et al. day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/mancino/', 'day_by_day_mancino_', date1, '.jpg'))
}

stats_match_df_mancino$`Model Name` = 'Mancino, et al.'
write_csv(stats_match_df_mancino, paste0(dumpdir, 'day_by_day_mancino_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(mancino_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_mancino_modelsummary_v', Sys.Date(), '.RDS'))



#################wu et al.#################################################

#re-initialize the summary table
stats_match_df_wu <- match_template

#create empty object to save models to
wu_dbd_models <- NULL

#iterate through unique dates and run wu alg for each date then export to a file
for(i in 1:nrow(stats_match_df_wu)){
  date1 = as.Date(stats_match_df_wu$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #create new variables needed for wu
  dframe_train$lnSecchi <- log(dframe_train$secchi_depth_m)
  dframe_test$lnSecchi <- log(dframe_test$secchi_depth_m)
  
  #run wu on training data
  output1 = lm(lnSecchi ~ B3 + B1, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$ln_dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$ln_dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_wu dataframe
  stats_match_df_wu$mae_test[i] = mae(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #MAE
  stats_match_df_wu$rmse_test[i] = rmse(dframe_test$secchi_depth_m, exp(dframe_test$ln_dbd_p)) #RMSE
  stats_match_df_wu$R2_test[i] = 1-(sum((exp(dframe_test$ln_dbd_p) - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared
  stats_match_df_wu$mae_train[i] = mae(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p))  #MAE
  stats_match_df_wu$rmse_train[i] = rmse(dframe_train$secchi_depth_m, exp(dframe_train$ln_dbd_p)) #RMSE
  stats_match_df_wu$R2_train[i] =  1 - (sum((exp(dframe_train$ln_dbd_p) - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_wu$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_wu$po_slope_train[i] <- summary(lm(exp(dframe_train$ln_dbd_p) ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_wu$po_slope_test[i]<- summary(lm(exp(dframe_test$ln_dbd_p) ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  
  wu_dbd_models[[i]] <- summary(output1)
  wu_dbd_models[[i]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              exp(max(alldata_dbd$ln_dbd_p)))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = exp(ln_dbd_p), color = dataset)) +
    geom_point() +
    labs(title = paste0('Wu, et al. day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/wu/', 'day_by_day_wu_', date1, '.jpg'))
}

stats_match_df_wu$`Model Name` = 'Wu, et al.'
write_csv(stats_match_df_wu, paste0(dumpdir, 'day_by_day_wu_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(wu_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_wu_modelsummary_v', Sys.Date(), '.RDS'))


#################Yip et al.#################################################

#re-initialize the summary table
stats_match_df_yip <- match_template

#create a null object for the dbd models
yip_dbd_models <- NULL

#iterate through unique dates and run yip alg for each date then export to a file
for(i in 1:nrow(stats_match_df_yip)){
  date1 = as.Date(stats_match_df_yip$sat_date[i])
  
  dframe_train = training[as.Date(training$sat_date) == as.Date(date1), ]
  dframe_test = test[as.Date(test$sat_date) == as.Date(date1), ]
  
  #run yip on training data
  output1 = lm(secchi_depth_m ~ B4 + B2 + B1, data = dframe_train)
  
  #predict the data based on alg object, save in dataframe
  dframe_train$dbd_p = as.numeric(predict(output1, dframe_train))
  dframe_test$dbd_p = as.numeric(predict(output1, dframe_test))
  
  alldata_dbd <- full_join(dframe_train, dframe_test)
  
  #save stats in stats_match_df_yip dataframe
  stats_match_df_yip$mae_test[i] = mae(dframe_test$secchi_depth_m, dframe_test$dbd_p) #MAE
  stats_match_df_yip$rmse_test[i] = rmse(dframe_test$secchi_depth_m, dframe_test$dbd_p) #RMSE
  stats_match_df_yip$R2_test[i] = 1-(sum((dframe_test$dbd_p - dframe_test$secchi_depth_m)^2)/sum((dframe_test$secchi_depth_m - mean(dframe_test$secchi_depth_m)) ^ 2)) #R-squared 
  stats_match_df_yip$mae_train[i] = mae(dframe_train$secchi_depth_m, dframe_train$dbd_p)  #MAE
  stats_match_df_yip$rmse_train[i] = rmse(dframe_train$secchi_depth_m, dframe_train$dbd_p) #RMSE
  stats_match_df_yip$R2_train[i] = 1 - (sum((dframe_train$dbd_p - dframe_train$secchi_depth_m) ^ 2)/sum((dframe_train$secchi_depth_m - mean(dframe_train$secchi_depth_m)) ^ 2))
  stats_match_df_yip$pval_train[i] = glance(output1)$p.value #p-value
  
  #calculate the slope of predicted and observed
  stats_match_df_yip$po_slope_train[i] <- summary(lm(dframe_train$dbd_p ~ dframe_train$secchi_depth_m))$coefficients[2,1]
  stats_match_df_yip$po_slope_test[i]<- summary(lm(dframe_test$dbd_p ~ dframe_test$secchi_depth_m))$coefficients[2,1]
  
  yip_dbd_models[[date1]] <- summary(output1)
  yip_dbd_models[[date1]]$model_date <- date1
  
  #define max SD for scale
  maxSD = max(max(alldata_dbd$secchi_depth_m), 
              max(alldata_dbd$dbd_p))
  
  ggplot(alldata_dbd, aes(x = secchi_depth_m, y = dbd_p, color = dataset)) +
    geom_point() +
    labs(title = paste0('Yip, et al. day by day\n', date1),
         x = 'observed secchi depth (m)',
         y = 'predicted secchi depth (m)') +
    coord_cartesian(xlim = c(0, maxSD), ylim = c(0, maxSD)) +
    final_theme +
    scale_color_colorblind()
  ggsave(paste0(figdir, 'day_by_day/yip/', 'day_by_day_yip_', date1, '.jpg'))
}

#save stats for tables
stats_match_df_yip$`Model Name` = 'Yip, et al.'
write_csv(stats_match_df_yip, paste0(dumpdir, 'day_by_day_yip_summarystats_v', Sys.Date(), '.csv'))
#save RDS of models
write_rds(yip_dbd_models, paste0(dumpdir, 'model_summary/day_by_day_yip_modelsummary_v', Sys.Date(), '.RDS'))

