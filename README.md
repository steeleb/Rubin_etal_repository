# Rubin, et al. Code Repository
 R programs used for the Rubin, et al. manuscript
 Contact: hrubin@vols.utk.edu
 
 This is a repository of the code used to complete analyses in the Rubin, et al. manuscript 'Remote Sensing of Lake Water Clarity: Performance and Transferability of Both Historical Algorithms and Machine Learning'
 
 Raw data are not stored in this repository, but can be made available upon request. 
 
# Descriptions of programs within this repository:
 
band_ratio_architecture update 01.19.R - This program calculates new coefficients for historical algorithms and applies the machine learning algorithm developed in the program 'Random Forest Full Dataset 128 trees 11Jan2021.' It also calculates pseudo-R^2, MAE and RMSE for training and testing sets.
 
Random Forest Full Dataset 128 trees 11Jan2021.R - This program builds a random forest model on the training data for this analysis.

iterate_by_day.R - This program applies each historical and machine learning algorithm over the data on a day-by-day basis. 

boxplot_figure_daybyday_wholedataset_mae_rmse.R - This program creates all non-mapped figures in the manuscript.