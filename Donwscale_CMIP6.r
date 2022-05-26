################################################################################
# STATISTICAL DOWNSCALING OF DAILY CLIMATE DATA USING                          #
# QUANTILE MAPPING TECHNIQUE FOR CMIP6 DATASETS                                #
################################################################################

#Author: Julio Montenegro Gambini, M.Sc.,
#PhD fellow - Technische Universiteit Delft (TU Delft), Netherlands.

#Current version: 1.0

#©Copyright 2013  2021 Julio Montenegro.
#This script is strictly under license GPLv3
#License details: https://www.gnu.org/licenses/gpl-3.0.en.html

# Please, when using this script, cite as: "Montenegro, J. (2021). Statistical
#downscaling of daily climate data using quantile mapping technique for CMIP6
#datasets"

# Installing or loading the required packages ==================================
library(tidyverse)
library(lubridate)
library(qmap)
library(zoo)
library(latticeExtra)
library(readxl)
library(beepr)

#Setting working directory =====================================================
#IMPORTANT!: The main folder has to contain sub-folders of each model
#where 6 files are located.
#The file names are: TASMIN_OBS.xlsx, TASMIN.csv, TASMAX_OBS.xlsx, TASMAX.csv, 
#PR_OBS.xlsx, PR.csv
#Historical observed data: TASMAX_OBS, TASMIN_OBS, PR_OBS
#Future data: TASMIN, TASMAX, PR
#Each file has to contain continous daily data!
#Data cannot contain missing values!

setwd(paste0('C:/Users/monte/Downloads/CMIP6_DOWNSCALING/',"SSP 585"))
dir <- dir()
length(dir)

#Pre-procesing and downscaling (EDIT HERE!) ====================================
#Creating a general loop for processing within each model (folder)
for(ind in 1:1){
  # Name of the folder which contains csv and xlsx excel files
  dir_estation <- dir[ind]

# CHANGE HERE THE FILE NAMES! ==================================================
  #Historical observed data: TASMAX_OBS, TASMIN_OBS, PR_OBS
  #Future data: TASMIN, TASMAX, PR

  # DAILY MINIMUM TEMPERATURE PARAMETERS
  file_tmin_obs <- 'TASMIN_OBS.xlsx'
  file_tmin_mod <- 'TASMIN.csv'
  fecha_in_his_min <- as.Date('1981-01-01')
  fecha_fin_his_min <- as.Date('2019-12-31')
  fecha_in_min <- as.Date('1950-01-01')
  fecha_fin_min <- as.Date('2099-12-31')
  
  # DAILY MAXIMUM TEMPERATURE PARAMETERS
  file_tmax_obs <- 'TASMAX_OBS.xlsx'
  file_tmax_mod <- 'TASMAX.csv'
  fecha_in_his_max <- as.Date('1981-01-01')
  fecha_fin_his_max <- as.Date('2019-12-31')
  fecha_in_max <- as.Date('1950-01-01')
  fecha_fin_max <- as.Date('2099-12-31')
  
  # DAILY PRECIPITATION PARAMETERS
  file_pr_obs <- 'PR_OBS.xlsx'
  file_pr_mod <- 'PR.csv'
  fecha_in_his_pr <- as.Date('1981-01-01')
  fecha_fin_his_pr <- as.Date('2019-12-31')
  fecha_in_pr <- as.Date('1950-01-01')
  fecha_fin_pr <- as.Date('2099-12-31')
  
# GENERATING A QUANTILE MAPPING FUNCTION =======================================
  
  # Downscaling function for each variable
  qp_function <- function(time_ini_his, time_fin_his,
                          time_ini, time_fin,
                          file_his, file_mod,
                          var){
    
    time_ini <- fecha_in_min
    time_fin <- fecha_fin_min
    time_ini_his <- fecha_in_his_min
    time_fin_his <- fecha_fin_his_min
    #file_his <- file_tmin_obs
    #file_mod <- file_tmin_mod
    #var = 'MIN'
    #j <- 2
    
    ## Reading historical baseline and GCM/RCM future time series
    data_historica <- read_excel(paste0(dir_estation,'/',file_his))
    data_modelada <- read.csv(paste0(dir_estation,'/',file_mod))
    
    names_stations <- names(data_historica)[2:length(data_historica)]
    
    df_reg <- c() # empty dataframe for filling with downscaled data
    
    ## Downscaling loop for different time series (station data)
    for (j in 2:ncol(data_historica)) {
      
      ## Historical baseline configuration
      data_his <- data_historica[,c(1,j)]
      colnames(data_his) <- c('FECHA','STATION')
      data_his <- data_his %>%  mutate(FECHA= as.Date(FECHA))
      
      ## GCM/RCM data configuration
      if (var %in% c('MIN','MAX')) {
        correccion <- function(x, na.rm=FALSE) (x-273.15)
      }
      if (var %in% c('PR')) {
        correccion <- function(x, na.rm=FALSE) (x*86400)
      }
      
      data_mod <- data_modelada[,c(1,j)] %>% 
        mutate_if(is.numeric, correccion, na.rm=FALSE)
      colnames(data_mod) <- c('FECHA','STATION')
      
      ####CONDICION IMPORTANTE ADICIONAL COMO CORRECCION
      if(data_modelada[1,1] == "1/1/1950"){
        data_model <- data_mod %>%  mutate(FECHA= as.Date(FECHA, format = "%m/%d/%Y"))
      }else{data_model <- data_mod %>%  mutate(FECHA= as.Date(FECHA))}
      
      ## Date filtering
      
      ### Historical
      var_hist_2 <- data_his %>% 
        filter(FECHA >= time_ini_his & FECHA <= time_fin_his) 
      colnames(var_hist_2) <- c('isodate','hist')
      
      var_hist_2 <- var_hist_2 %>% 
        full_join(
          data.frame(isodate = seq(from=time_ini_his, to=time_fin_his, by ='day')),
          by = 'isodate') 
      
      ### Setting the variable of GCM/RCM data
      var_model_2 <- data_model %>% 
        filter(FECHA >= time_ini & FECHA <= time_fin)
      colnames(var_model_2) <- c('isodate','mode')
      
      var_model_2 <- var_model_2 %>% 
        full_join(
          data.frame(isodate = seq(from=time_ini, to=time_fin, by ='day')),
          by = 'isodate') 
      
      ## Creating variables to be used
      OBS_hist <- var_hist_2 %>% 
        rename(OBS_hist = hist) 
      
      GCM_model <- var_model_2
      
      # Filling missing historical data
      if (var=='PR') {
        OBS_hist <- OBS_hist %>%
          mutate(OBS_hist = ifelse(is.na(OBS_hist),0.1,OBS_hist))
      }
      if (var=='MAX' | var=='MIN') {
        OBS_hist <- OBS_hist %>%
          mutate(OBS_hist = ifelse(is.na(OBS_hist),16,OBS_hist))
      }
      
      # Filling missing GCM/RCM future data
      if (var=='PR') {
        GCM_model <- GCM_model %>%
          mutate(mode = ifelse(is.na(mode),0.1,mode))
      }
      if (var=='MAX' | var=='MIN') {
        GCM_model <- GCM_model %>%
          mutate(mode = ifelse(is.na(mode),15,mode))
      }
      
      GCM_hist <- GCM_model
      
      ## Conversion of datasets to "ts" objects
      
      data_hist <- OBS_hist %>% 
        read.zoo()
      
      data_mod <- GCM_model %>% 
        read.zoo()
      
      data_wt <- GCM_model %>% 
        read.zoo()
      
# SETTING SEASONAL OR MONTHLY ANALYSIS =========================================

      
      seasons_by_year <- list(c("December"),c("January"),c("February"), 
                              c("March"),c("April"),c("May"), 
                              c("June"),c("July"),c("August"),
                              c("September"),c("October"),c("November"))
      #According to system region, month names should be changed (e.g. spanish):
      #seasons_by_year <- list(c("Diciembre"), c("Enero"), c("Febrero"),
      #c("Marzo"), c("Abril"), c("Mayo"), 
      #c("Junio"), c("Julio"), c("Agosto"),
      #c("Septiembre"), c("Octubre"), c("Noviembre"))
      
      for(i in 1:12) {
        
        obs_sl <- data_hist[months(time(data_hist)) %in% seasons_by_year[[i]]]
        mod_sl <- data_mod[months(time(data_mod)) %in% seasons_by_year[[i]]]
        #GCM/RCM data, read!: L. Gudmundsson et al. (2012)
        
        if (sum(mod_sl, na.rm = T)==0) {
          mod_sl[1]<- 0.001
        }
        
        qm_fit <- fitQmapQUANT(obs = coredata(obs_sl),
                               mod = coredata(mod_sl),
                               qstep = 0.001,
                               nboot = 1, 
                               wet.day = 0, # To be changed for temperatures
                               type = "linear")
        
        mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")
        
        data_wt[ months(time(data_wt)) %in% seasons_by_year[[i]]] <- mod_sl_qmapped
        
      }
      
      t <- as.data.frame(data_wt)
      
      if (j == 2) {
        df_reg <- t
      } else{
        df_reg <- cbind(df_reg, t)
      }
      
    }
    
    # Adding column names
    colnames(df_reg) <- names_stations
    
    # Assign the date column (daily)
    df_empty<- data.frame(isodate = seq(from=time_ini, to=time_fin, by ='day'))
    
    df_out <- df_empty %>% 
      cbind(df_reg)
    
    return(df_out)
  }
  
# APPLYING DOWNSCALING =========================================================
  
  # Dowscaling for minimum temperature
  tmin_reg <- qp_function(fecha_in_min, fecha_fin_min,
                          fecha_in_his_min, fecha_fin_his_min,
                          file_tmin_obs, file_tmin_mod, 'MIN')
  
  # Dowscaling for maximum temperature
  tmax_reg <- qp_function(fecha_in_max, fecha_fin_max,
                          fecha_in_his_max, fecha_fin_his_max,
                          file_tmax_obs, file_tmax_mod, 'MAX')
  
  # Dowscaling for precipitation
  pr_reg <- qp_function(fecha_in_pr, fecha_fin_pr,
                        fecha_in_his_pr, fecha_fin_his_pr,
                        file_pr_obs, file_pr_mod, 'PR')
  
  
# EXPORTING DOWNSCALED DATA IN 3 FILES =========================================
  tmin_reg %>%  write.csv(file = paste0(dir_estation,'/tmin_reg.csv'),row.names = F)
  tmax_reg %>%  write.csv(file = paste0(dir_estation,'/tmax_reg.csv'),row.names = F)
  pr_reg %>%  write.csv(file = paste0(dir_estation,'/pr_reg.csv'),row.names = F)
  #Hear the sound when finally the three files were generated
  beep(15)
}
#In case of an error, another sound is played
beep(5)
