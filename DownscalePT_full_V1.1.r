################################################################################
# STATISTICAL DOWNSCALING OF DAILY CLIMATE DATA USING                          #
# QUANTILE MAPPING (QPM)TECHNIQUE                                              #
################################################################################

#Author: Julio Montenegro Gambini, P.E. ASCE, M.Sc.,
#PhD fellow - Technische Universiteit Delft (TU Delft), Netherlands.

#Current version: 1.0
#This script is strictly under license GPLv3
#License details: https://www.gnu.org/licenses/gpl-3.0.en.html
#About quantile mapping technique (Gudmundsson et al., 2012): https://hess.copernicus.org/articles/16/3383/2012/

#' Para este script se necesita:
#' - 3 archivos de la data observada en la carpeta de trabajo
#' - 1 carpeta con nombre de la estacion y dentro de ella, 3 archivos
#'     con la data modelada historica (que coincida con la data observada)
#'     y 3 archivos con la data modelada a futuro

# Librerias
library(tidyverse)
library(lubridate)
library(qmap)
library(zoo)
library(latticeExtra)
library(readxl)

# Indicar la carpeta que contiene a las estaciones
setwd('C:/Users/Julio/Downloads/qmap/qmap')

# llamamos a los archivos historicos
data_pr_his <- read_excel("Histórico_Est.xlsx") 
data_tmax_his <- read_excel("Tmáx_esta.xlsx") 
data_tmin_his <- read_excel("Tmín_esta.xlsx") 

#============================================================================
# Configurar segun la data

## Numero de variables
n_var <- 3

## Nombre del archivo de la data modelada que sigue despues del texto 'RCP4585'
vec_var_2 <- c('TASMAX','TASMIN','PR')
## Nombre resumen de la variable (recomiendo dejar como esta)
vec_var_3 <- c('MAX','MIN','PR')

## Fijar las fechas historicas para cada variable
fecha_in_his_min <- as.Date('1980-01-01')
fecha_fin_his_min <- as.Date('2019-12-31')

fecha_in_his_max <- as.Date('1980-01-01')
fecha_fin_his_max <- as.Date('2019-12-31')

fecha_in_his_pr <- as.Date('1980-01-01')
fecha_fin_his_pr <- as.Date('2016-12-31')

## Fijar las fechas de inicio y fin de la data modelada
fecha_in_min <- as.Date('1950-01-01')
fecha_fin_min <- as.Date('2099-12-31')

fecha_in_max <- as.Date('1950-01-01')
fecha_fin_max <- as.Date('2099-12-31')

fecha_in_pr <- as.Date('1950-01-01')
fecha_fin_pr <- as.Date('2099-12-31')

#============================================================================

# Funci?n para generar la regionalizacion, para cada variable y para un solo bloque de tiempo

ds <- function(time_ini_his, time_fin_his, time_ini, time_fin, data_his, data_model, data_model_his, var){
  
  # time_ini <- as.Date('1980-01-01')
  # time_fin <- as.Date('2099-12-31')
  # # time_ini <- as.Date('2040-01-01')
  # # time_fin <- as.Date('2069-12-31')
  # # time_ini <- as.Date('2070-01-01')
  # # time_fin <- as.Date('2099-12-31')
  # time_ini_his <- as.Date('1980-01-01')
  # time_fin_his <- as.Date('2019-12-31')
  #  data_his <- var_hist
  #  data_model <- var_model
  #  var = 'MAX'
  # data_model_his = var_model_hist

  
  data_model <-data_model %>% 
    dplyr::select(isodate, starts_with('rcp'))
  
  data_model_his <- data_model_his %>% 
    left_join(data_model_his, by = 'isodate') 
  colnames(data_model_his) <- names(data_model)
    
  data_model <- rbind(data_model_his, data_model)
  
  # Creamos vectores con la informacion de los modelos
  name_mod <- sort(names(data_model)[which(substring(names(data_model),1,3)=='rcp')] )
  n_model <- length(name_mod)
  # Creamos df vacios quienes van a contener a los resultados
  df_mod <- c()
  
  # Generamos la regionalizaci?n por tipo de proyeccion
  for (j in 1:n_model) {
    
    # configurar la variable historica
    var_hist_2 <- data_his %>% 
      filter(FECHA >= time_ini_his & FECHA <= time_fin_his) 
    colnames(var_hist_2) <- c('isodate','hist')
    
    var_hist_2 <- var_hist_2 %>% 
      full_join(
        data.frame(isodate = seq(from=time_ini_his, to=time_fin_his, by ='day')),
        by = 'isodate') 
    
    # configurar la variable del modelo
    var_model_2 <- data_model[, c(which(names(data_model) == 'isodate'),
                                  which(names(data_model) == name_mod[j]))] %>% 
      filter(isodate >= time_ini & isodate <= time_fin)
    colnames(var_model_2) <- c('isodate','mode')
    
    var_model_2 <- var_model_2 %>% 
      full_join(
        data.frame(isodate = seq(from=time_ini, to=time_fin, by ='day')),
        by = 'isodate') 

    # creamos las variables a usar (solo para seguir el script de origen)
    OBS_hist <- var_hist_2 %>% 
      rename(OBS_hist = hist) 
    
    GCM_model <- var_model_2
    
    # aqui completo valores faltantes 
    if (var=='PR') {
      OBS_hist <- OBS_hist %>%
        mutate(OBS_hist = ifelse(is.na(OBS_hist),0.1,OBS_hist))
    }
    if (var=='MAX' | var=='MIN') {
      OBS_hist <- OBS_hist %>%
        mutate(OBS_hist = ifelse(is.na(OBS_hist),16,OBS_hist))
    }
    
    # aqui completo valores faltantes 
    if (var=='PR') {
      GCM_model <- GCM_model %>%
        mutate(mode = ifelse(is.na(mode),0.1,mode))
    }
    if (var=='MAX' | var=='MIN') {
      GCM_model <- GCM_model %>%
        mutate(mode = ifelse(is.na(mode),15,mode))
    }
    
    GCM_hist <- GCM_model
  
    data_hist <- OBS_hist %>% 
      read.zoo()
    
    data_mod <- GCM_model %>% 
      read.zoo()

    data_wt <- GCM_model %>% 
      read.zoo()
    
    #============================================================================
    # APLICACIÓN DE LA TÉCNICA DE QUANTILE MAPPING (MÉTODO EMPÍRICO)
    
    seasons_by_year <- list(c("December"), c("January"), c("February"),
                            c("March"), c("April"), c("May"), 
                            c("June"), c("July"), c("August"),
                            c("September"), c("October"), c("November"))
    
    for(i in 1:12) {

      obs_sl <- data_hist[months(time(data_hist)) %in% seasons_by_year[[i]]]
      mod_sl <- data_mod[months(time(data_mod)) %in% seasons_by_year[[i]]]
      #MODEL, read!: L. Gudmundsson et al. (2012)

      if (sum(mod_sl, na.rm = T)==0) {
        mod_sl[1]<- 0.001
      }
      
      qm_fit <- fitQmapQUANT(obs = coredata(obs_sl),
                             mod = coredata(mod_sl),
                             qstep = 0.001,
                             nboot = 1, 
                             wet.day = 0, # Adaptado para temperatuas negativas
                             type = "linear")
      
      mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")

      data_wt[ months(time(data_wt)) %in% seasons_by_year[[i]]] <- mod_sl_qmapped
      
      }
 
    t <- as.data.frame(data_wt)
    
    if (j == 1) {
      df_mod <- t
    } else{
      df_mod <- cbind(df_mod, t)
    }
    
  }
  
  # agregamos losnombres correctas de las columnas
  colnames(df_mod) <- name_mod
  
  # asignamos la columna de tiempo
  df_empty<- data.frame(isodate = seq(from=time_ini, to=time_fin, by ='day'))

  df_out <- df_empty %>% 
    cbind(df_mod)
  
  return(df_out)
  }

#============================================================================
files_eliminar <- c("down.R","down2.R","down3.R", "Histórico_Est.xlsx", "Tmáx_esta.xlsx", "Tmín_esta.xlsx")
estaciones <- setdiff(dir(),files_eliminar)

# Bucle para cada estacion
n_estacion <- length(estaciones)

for (estacion in 1:n_estacion) {
# Bucle para que ejecute la regionalizaci?n para cada bloque y los exporte en excel por cada variable
for (z in 1:3) {
  
  if (vec_var_3[z]=='MAX') {
    
    # Seleccionamos el archivo de una variable
    var_hist <- data_tmax_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    dif.273 <- function(x, na.rm=FALSE) (x-273.15)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585TASMAX.csv')) %>%
      mutate_if(is.numeric, dif.273, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
    
    var_model_hist <- read.csv(paste0(estaciones[estacion],'/HISTASMAX.csv')) %>% 
      dplyr::select(isodate, starts_with('hist')) %>% 
      mutate(isodate = as.Date(isodate)) %>% 
      mutate_if(is.numeric, dif.273, na.rm=FALSE)
    
    list_time <- ds(time_ini_his = fecha_in_his_max,
                    time_fin_his = fecha_fin_his_max,
                    time_ini = fecha_in_max,
                    time_fin = fecha_fin_max,
                    data_his = var_hist,
                    data_model = var_model,
                    data_model_his = var_model_hist,
                    var = 'MAX')
  }
  
  if (vec_var_3[z]=='MIN') {
    # Seleccionamos el archivo de una variable
    var_hist <- data_tmin_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    dif.273 <- function(x, na.rm=FALSE) (x-273.15)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585TASMIN.csv')) %>%
      mutate_if(is.numeric, dif.273, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
    
    var_model_hist <- read.csv(paste0(estaciones[estacion],'/HISTASMIN.csv')) %>% 
      dplyr::select(isodate, starts_with('hist')) %>% 
      mutate(isodate = as.Date(isodate)) %>% 
      mutate_if(is.numeric, dif.273, na.rm=FALSE)
    
    list_time <- ds(time_ini_his = fecha_in_his_min,
                    time_fin_his = fecha_fin_his_min,
                    time_ini = fecha_in_min,
                    time_fin = fecha_fin_min,
                    data_his = var_hist,
                    data_model = var_model,
                    data_model_his = var_model_hist,
                    var = 'MIN')
    
  }
  
  if (vec_var_3[z]=='PR') {
    # Seleccionamos el archivo de una variable
    
    var_hist <- data_pr_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    mult.86400 <- function(x, na.rm=FALSE) (x*86400)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585PR.csv')) %>%
      mutate_if(is.numeric, mult.86400, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
    
    var_model_hist <- read.csv(paste0(estaciones[estacion],'/HISTPR.csv')) %>% 
      dplyr::select(isodate, starts_with('hist')) %>% 
      mutate(isodate = as.Date(isodate)) %>% 
      mutate_if(is.numeric, mult.86400, na.rm=FALSE)
      
    list_time <- ds(time_ini_his = fecha_in_his_pr,
                    time_fin_his = fecha_fin_his_pr,
                    time_ini = fecha_in_pr,
                    time_fin = fecha_fin_pr,
                    data_his = var_hist,
                    data_model = var_model,
                    data_model_his = var_model_hist,
                    var = 'PR')
  }
  
  # Exportamos
  name_file_out <- paste0(estaciones[estacion],'/',estaciones[estacion],'_',vec_var_3[z],'.csv')

  write.csv(list_time, file = name_file_out,row.names = F)

  rm(list_time)
  
  }

}

