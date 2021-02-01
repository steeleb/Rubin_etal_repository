
###file originally writtedn by H. Rubin Aug 2020 'archive/boxplot_algs_1.R'
# rewritten Dec 1 2020 by B. Steele to bring up to FAIR data standards and remove hardcoding, edit figure to be a single two-panel figure instead of 3 single-panel figures.

#libraries
library(plyr)
library(tidyverse)
library(readxl)
library(ggthemes)
library(cowplot)

#set directories
dbd_datadir <- 'C:/Users/steeleb/Dropbox/NASA IDS/insitu_historical/secchi/Rubin.MS/summary_data/daybyday_summary/'
wd_datadir <- 'C:/Users/steeleb/Dropbox/NASA IDS/insitu_historical/secchi/Rubin.MS/summary_data/wholedataset_summary/'
figdir <- 'C:/Users/steeleb/Dropbox/NASA IDS/insitu_historical/secchi/Rubin.MS/summary_data/plots/'


#load mae/rmse data files from day-by-day analysis

#list all of the output files in the datadir folder
filelist <- list.files(dbd_datadir, pattern = '*.csv')
filelist

#create a blank object to save files to
alldata_dbd <- NULL

#loop over all listed files, loading them as lists into the blank object
for (i in 1:length(filelist)) {
  data <- read_csv(paste0(dbd_datadir, filelist[i]))
  alldata_dbd[[i]] <- data
}

#use ldply to convert lists to a single dataframe
alldata_dbd <- ldply(alldata_dbd, data.frame) %>% 
  select(-n_train, -n_test)  %>% 
  mutate(sat_date = as.character(sat_date),
         dbd_whole = 'day by day')

#load mae/rmse from whole dataset analysis
wd_data <- read_csv(paste0(wd_datadir, 'Table3_v19Jan2021.csv'),
                     col_names = c('Model.Name','mae_train', 'rmse_train', 'R2_train',
                                   'mae_test', 'rmse_test', 'R2_test', 'po_slope_train', 'po_slope_test'),
                     skip = 1) %>% 
  mutate(sat_date = 'all data',
         dbd_whole = 'whole dataset') #note this will ellicit the warning message 'na's introduced by coercion' because of the text strin in this for RF.
                     

#format wd_data so that it will merge nicely with the day-by-day
colnames(alldata_dbd)
colnames(wd_data)
str(alldata_dbd)
str(wd_data)

#merge dbd with wd algs
alg_data <- full_join(alldata_dbd, wd_data) %>% 
  arrange(Model.Name)

#get some summary stats for paper
training <- read_csv('5.a.90.10.datasets/90percent_7.23.20.csv')
test <- read_csv('5.a.90.10.datasets/OOB_7.23.20.csv')

tandt <- full_join(training, test)

#number of lakes
length(unique(tandt$dissPermID))
#number of satellite dates
length(unique(tandt$sat_date))
#number of secchi dates
length(unique(tandt$date))
#years
min(tandt$sat_date)
max(tandt$sat_date)
#number of samples
nrow(tandt)

#mean, sd
mean(test$secchi_depth_m)
sd(test$secchi_depth_m)
training <- training %>% 
  filter(B3>3)
mean(training$secchi_depth_m)
sd(training$secchi_depth_m)

#range of MAE, RMSE values for RF dbd
rf_dbd <- alldata_dbd %>% 
  filter(Model.Name == 'Random Forest' & dbd_whole == 'day by day')
min(rf_dbd$mae_train)
max(rf_dbd$mae_train)
min(rf_dbd$rmse_train)
max(rf_dbd$rmse_train)

allother_dbd <- alldata_dbd %>% 
  filter(Model.Name != 'Random Forest' & dbd_whole == 'day by day')
min(allother_dbd$mae_train)
max(allother_dbd$mae_train)
min(allother_dbd$rmse_train)
max(allother_dbd$rmse_train)


#make sure alg names match, add year of pub
unique(alg_data$Model.Name)

alg_data <- alg_data %>% 
  mutate(Model.Name.Year = case_when(grepl('Johnson', Model.Name) ~ 'Allee and Johnson (1999)',
                                grepl('Baban', Model.Name) ~ 'Baban (1992)',
                                grepl('Cox', Model.Name)~ 'Cox, et al. (1998)',
                                grepl('Chipman', Model.Name) ~ 'Chipman, et al. (2004)',
                                Model.Name == 'Dekker and Peters' | Model.Name == 'Dekker and Peters 1' ~ 'Dekker and Peters (1993) (1)',
                                Model.Name == 'Dekker and Peters 2' ~ 'Dekker and Peters (1993) (2)',
                                grepl('Dominguez', Model.Name) ~ 'Dominguez Gomez, et al. (2009)',
                                grepl('Giardino', Model.Name) ~ 'Giardino, et al. (2009)',
                                grepl('Kloiber', Model.Name)~ 'Kloiber, et al. (2002)',
                                Model.Name == 'Lathrop 1' ~ 'Lathrop (1992) (1)',
                                Model.Name == 'Lathrop 2' ~ 'Lathrop (1992) (2)',
                                grepl('Lathrop and', Model.Name) ~ 'Lathrop and Lillesand (1991)',
                                grepl('Lavery', Model.Name) ~ 'Lavery, et al. (1993)',
                                grepl('Mancino', Model.Name)~ 'Mancino, et al. (2009)',
                                grepl('McCullough', Model.Name) ~ 'McCullough, et al. (2012)',
                                grepl('Wu', Model.Name) ~ 'Wu, et al. (2008)',
                                grepl('Yip', Model.Name) ~ 'Yip, et al. (2015)', 
                                grepl('Random', Model.Name) ~ 'Random Forest',
                                grepl('Four', Model.Name, ignore.case = T) ~ 'Four Band Linear',
                                TRUE  ~ Model.Name)) %>% 
  mutate(Model.Name.Year = factor(Model.Name.Year, levels = c('Allee and Johnson (1999)', 'Baban (1992)','Chipman, et al. (2004)','Cox, et al. (1998)',
                                                    'Dekker and Peters (1993) (1)','Dekker and Peters (1993) (2)','Dominguez Gomez, et al. (2009)',
                                                    'Giardino, et al. (2009)','Kloiber, et al. (2002)','Lathrop (1992) (1)',
                                                    'Lathrop (1992) (2)','Lathrop and Lillesand (1991)','Lavery, et al. (1993)', 'Mancino, et al. (2009)',
                                                    'McCullough, et al. (2012)', 'Wu, et al. (2008)','Yip, et al. (2015)',
                                                    'Four Band Linear','Random Forest')))
unique(alg_data$Model.Name.Year)


#plot using ggplot
#make a vertical dataset, reorganize so that train and test can be plotted next to each other
alg_data_vert <- alg_data %>% 
  select(Model.Name.Year, dbd_whole, R2_train, mae_train, rmse_train, R2_test, mae_test, rmse_test, po_slope_train, po_slope_test, sat_date) %>% 
  gather(variable, value, -Model.Name.Year, -dbd_whole, -sat_date) %>% 
  mutate(dataset = case_when(grepl('train', variable) ~ 'train',
                             grepl('test', variable) ~ 'test',
                             TRUE ~ NA_character_),
         variable = case_when(grepl('R2', variable) ~ 'R2',
                              grepl('mae', variable) ~ 'MAE',
                              grepl('rmse', variable) ~ 'RMSE',
                              grepl('po', variable) ~ 'Predicted v Observed slope',
                              TRUE ~ NA_character_)) %>% 
  spread(variable, value)

max_mae <- max(alg_data_vert$MAE)
max_rmse <- max(alg_data_vert$RMSE)
max_R2 <- max(alg_data_vert$R2)
max_slope <- max(alg_data_vert$`Predicted v Observed slope`)

#characterize predicted versus observed slope:
alg_data_vert <- alg_data_vert %>%
  mutate(slope_class = case_when(`Predicted v Observed slope` < 0.1 ~ 'predicted v observed slope < 0.1',
                                 `Predicted v Observed slope` >= 0.1 ~ 'predicted v observed slope >= 0.1'),
         slope_class = factor(slope_class, levels = c('predicted v observed slope < 0.1','predicted v observed slope >= 0.1')))
#add dataset info to the Model Name Year
alg_data_vert <- alg_data_vert %>% 
  mutate(dataset = factor(dataset, levels = c('train', 'test'))) %>% 
  mutate(Model.Name.Year.Dataset = paste(Model.Name.Year, dataset, sep = ' - ')) %>% 
  mutate(Model.Name.Year.Dataset = factor(Model.Name.Year.Dataset, levels = c('Allee and Johnson (1999) - train', 'Allee and Johnson (1999) - test', 
                                                                              'Baban (1992) - train', 'Baban (1992) - test',
                                                                              'Chipman, et al. (2004) - train', 'Chipman, et al. (2004) - test',
                                                                              'Dekker and Peters (1993) (1) - train', 'Dekker and Peters (1993) (1) - test',                                                            
                                                                              'Dekker and Peters (1993) (2) - train', 'Dekker and Peters (1993) (2) - test',
                                                                              'Dominguez Gomez, et al. (2009) - train', 'Dominguez Gomez, et al. (2009) - test',
                                                                              'Giardino, et al. (2009) - train', 'Giardino, et al. (2009) - test',
                                                                              'Kloiber, et al. (2002) - train', 'Kloiber, et al. (2002) - test',
                                                                              'Lathrop and Lillesand (1991) - train',  'Lathrop and Lillesand (1991) - test',
                                                                              'Lavery, et al. (1993) - train',   'Lavery, et al. (1993) - test', 
                                                                              'Mancino, et al. (2009) - train', 'Mancino, et al. (2009) - test',
                                                                              'Wu, et al. (2008) - train', 'Wu, et al. (2008) - test',
                                                                              'Yip, et al. (2015) - train','Yip, et al. (2015) - test',
                                                                              'Four Band Linear - train', 'Four Band Linear - test',
                                                                              'Random Forest - train', 'Random Forest - test')))



# MAE_color <- ggplot(alg_data_vert, aes(x = MAE, y = reorder(Model.Name.Year.Dataset, desc(Model.Name.Year.Dataset)))) +
#   geom_point(aes(fill = slope_class, shape = dbd_whole), size=3) +
#   scale_shape_manual(values = c(21, 24), name = 'dataset') +
#   scale_fill_manual(values = c('#FFFFFF', '#696969'), 
#                     name = 'predicted v. observed slope', 
#                     labels = c('slope < 0.1', 'slope >= 0.1')) +
#   coord_cartesian(xlim = c(0, max_mae)) +
#   labs(y = NULL, x = 'Secchi depth (m)') +
#   theme_bw()+
#   guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))

#MAE_color_facet
ggplot(alg_data_vert, aes(x = MAE, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(fill = slope_class, shape = dbd_whole), size=2.5) +
  facet_grid(.~dataset, scales = 'free_x') +
  expand_limits(x =0) +
  scale_shape_manual(values = c(21, 24), name = 'dataset') +
  scale_fill_manual(values = c('#FFFFFF', '#696969'),
                    name = 'predicted v.\nobserved',
                    labels = c('slope < 0.1', 'slope >= 0.1')) +
  labs(y = NULL, x = 'MAE (m)') +
  theme_bw()+
  theme(legend.position = 'right') +
  guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))
ggsave('Rubin.MS/manuscript/MAE_panels_colorslope_v01Feb2021.jpeg', device = 'jpeg', width = 8, height = 4, dpi = 300)

# MAE_box <- alg_data_vert %>% 
#   filter(dbd_whole == 'day by day')  %>% 
#   ggplot(., aes(x = MAE, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
#   geom_boxplot() +
#   geom_point(data = alg_data_vert[alg_data_vert$dbd_whole == 'whole dataset',], aes(x = MAE), color = 'red', size = 2) +
#   facet_grid(.~dataset) +
#   scale_shape_manual(values = c(21, 24), name = 'dataset') +
#   scale_fill_manual(values = c('#FFFFFF')) +
#   coord_cartesian(xlim = c(0, max_mae)) +
#   labs(y = NULL, x = 'Secchi depth (m)') +
#   theme_bw()+
#   guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))
# 
# RMSE_color <- ggplot(alg_data_vert, aes(x = RMSE, y = reorder(Model.Name.Year.Dataset, desc(Model.Name.Year.Dataset)))) +
#   geom_point(aes(fill = slope_class, shape = dbd_whole), size=3) +
#   scale_shape_manual(values = c(21, 24), name = 'dataset') +
#   scale_fill_manual(values = c('#FFFFFF', '#696969'), 
#                     name = 'predicted v. observed slope', 
#                     labels = c('slope < 0.1', 'slope >= 0.1')) +
#   coord_cartesian(xlim = c(0, max_mae)) +
#   labs(y = NULL, x = 'Secchi depth (m)') +
#   theme_bw()+
#   guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))

#RMSE_color_facet
ggplot(alg_data_vert, aes(x = RMSE, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(fill = slope_class, shape = dbd_whole), size=2.5)+
  facet_grid(.~dataset, scales = 'free_x') +
  expand_limits(x =0) +
  scale_shape_manual(values = c(21, 24), name = 'dataset') +
  scale_fill_manual(values = c('#FFFFFF', '#696969'),
                    name = 'predicted v.\nobserved',
                    labels = c('slope < 0.1', 'slope >= 0.1')) +
  # coord_cartesian(xlim = c(0, max_rmse)) +
  labs(y = NULL, x = 'RMSE (m)') +
  theme_bw()+
  theme(legend.position = 'right') +
  guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))
ggsave('Rubin.MS/manuscript/RMSE_panels_colorslope_v01Feb2021.jpeg', device = 'jpeg', width = 8, height = 4, dpi = 300)



# #MAE
# ggplot(alg_data_vert, aes(x = MAE, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
#   geom_point(aes(shape = dbd_whole, color = dbd_whole), size = 3) +
#   facet_grid(.~dataset, scales = 'free_x') + 
#   labs(y = NULL,
#        x = 'MAE (m)') +
#   scale_shape_discrete(name = NULL)+
#   scale_color_colorblind(name = NULL) +
#   theme_bw() +
#   theme(legend.position = 'bottom')
# ggsave('Rubin.MS/manuscript/MAE_panels_v01Feb2021.jpeg', device = 'jpeg', width = 6.5, height = 5, dpi = 300)
# 
# #RMSE
# ggplot(alg_data_vert, aes(x = RMSE, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
#   geom_point(aes(shape = dbd_whole, color = dbd_whole), size = 3) +
#   facet_grid(.~dataset, scales = 'free_x') + 
#   labs(y = NULL,
#        x = 'RMSE (m)') +
#   scale_shape_discrete(name = NULL)+
#   scale_color_colorblind(name = NULL) +
#   theme_bw() +
#   theme(legend.position = 'bottom')
# ggsave('Rubin.MS/manuscript/RMSE_panels_v01Feb2021.jpeg', device = 'jpeg', width = 6.5, height = 5, dpi = 300)

#R2
ggplot(alg_data_vert, aes(x = R2, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(shape = dbd_whole, color = dbd_whole), size = 3) +
  facet_grid(.~dataset, scales = 'free_x') + 
  labs(y = NULL,
       x = expression('pseudo'-R^2)) +
  scale_shape_discrete(name = NULL)+
  scale_color_colorblind(name = NULL) +
  theme_bw() +
  theme(legend.position = 'bottom')

ggsave('Rubin.MS/manuscript/R2_panels_fullscale_v01Feb2021.jpeg', device = 'jpeg', width = 6.5, height = 5, dpi = 300)

alg_data_vert %>% 
  filter((dataset == 'test' & R2 >-11) | (dataset == 'train' & R2 >-0.25)) %>% 
  ggplot(., aes(x = R2, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(shape = dbd_whole, color = dbd_whole), size = 2.5) +
  facet_grid(.~dataset, scales = 'free_x') + 
  labs(y = NULL,
       x = expression('pseudo'-R^2)) +
  scale_shape_discrete(name = 'dataset')+
  scale_color_colorblind(name = 'dataset') +
  theme_bw() +
  theme(legend.position = 'right')
ggsave('Rubin.MS/manuscript/R2_panels_truncated_v01Feb2021.jpeg', device = 'jpeg', width = 8, height = 4, dpi = 300)


#R2_color_facet
ggplot(alg_data_vert, aes(x = R2, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(fill = slope_class, shape = dbd_whole), size=2.5)+
  facet_grid(.~dataset, scales = 'free_x') +
  expand_limits(x =0) +
  scale_shape_manual(values = c(21, 24), name = 'dataset') +
  scale_fill_manual(values = c('#FFFFFF', '#696969'),
                    name = 'predicted v.\nobserved',
                    labels = c('slope < 0.1', 'slope >= 0.1')) +
  # coord_cartesian(xlim = c(0, max_rmse)) +
  labs(y = NULL, 
       x = expression('pseudo'-R^2)) +
  theme_bw()+
  theme(legend.position = 'right') +
  guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))
ggsave('Rubin.MS/manuscript/R2_panels_colorslope_v01Feb2021.jpeg', device = 'jpeg', width = 8, height = 4, dpi = 300)

#R2_color_facet truncated
alg_data_vert %>% 
  filter((dataset == 'test' & R2 >-11) | (dataset == 'train' & R2 >-0.25)) %>% 
  ggplot(., aes(x = R2, y = reorder(Model.Name.Year, desc(Model.Name.Year)))) +
  geom_point(aes(fill = slope_class, shape = dbd_whole), size=2.5)+
  facet_grid(.~dataset, scales = 'free_x') +
  expand_limits(x =0) +
  scale_shape_manual(values = c(21, 24), name = 'dataset') +
  scale_fill_manual(values = c('#FFFFFF', '#696969'),
                    name = 'predicted v.\nobserved',
                    labels = c('slope < 0.1', 'slope >= 0.1')) +
  # coord_cartesian(xlim = c(0, max_rmse)) +
  labs(y = NULL, 
       x = expression('pseudo'-R^2)) +
  theme_bw()+
  theme(legend.position = 'right') +
  guides(fill=guide_legend(override.aes=list(shape=22, size = 5, color = '#9C99CE')))
ggsave('Rubin.MS/manuscript/R2_panels_colorslope_truncated_v01Feb2021.jpeg', device = 'jpeg', width = 8, height = 4, dpi = 300)



