# Packages to be used
library(dplyr)
library(lubridate)
library(qmap)
library(zoo)
library(latticeExtra)
library(xlsx)
library(readxl)

# Indicar la carpeta que contiene a las estaciones
setwd('D:/downscaling/')

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

## Fijar las fechas historicas
fecha_in_his <- as.Date('1980-01-01')
fecha_fin_his <- as.Date('2009-12-31')

## Fijar las fechas de los bloques
fecha_in_1 <- as.Date('2010-01-01')
fecha_fin_1 <- as.Date('2039-12-31')

fecha_in_2 <- as.Date('2040-01-01')
fecha_fin_2 <- as.Date('2069-12-31')

fecha_in_3 <- as.Date('2070-01-01')
fecha_fin_3 <- as.Date('2099-12-31')

## Agrupar en vector (dejar como esta)
fecha_in <- c(fecha_in_1, fecha_in_2, fecha_in_3)
fecha_fin <- c(fecha_fin_1, fecha_fin_2, fecha_fin_3)
#============================================================================

# Funci?n para generar la regionalizacion, para cada variable y para un solo bloque de tiempo

ds <- function(time_ini_his, time_fin_his, time_ini, time_fin, data_his, data_model, var){
  
  #time_ini <- as.Date('2010-01-01')
  #time_fin <- as.Date('2039-12-31')
  # time_ini <- as.Date('2040-01-01')
  # time_fin <- as.Date('2069-12-31')
  # time_ini <- as.Date('2070-01-01')
  # time_fin <- as.Date('2099-12-31')
  # time_ini_his <- as.Date('1980-01-01')
  # time_fin_his <- as.Date('2009-12-31')
  # data_his <- var_hist
  # data_model <- var_model
  # var = 'PR'
  
  # Creamos vectores con la informacion de los modelos
  #name_his <- sort(names(data_his)[which(substring(names(data_his),1,3)=='his')] )
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
      filter(isodate <= time_fin & isodate >= time_ini) %>% 
      full_join(
        data.frame(isodate = seq(from=time_ini, to=time_fin, by ='day')),
        by = 'isodate') 

    # aqui completo valores faltantes con 0.1
    if (var=='PR') {
      OBS_hist <- OBS_hist %>%
        mutate(OBS_hist = ifelse(is.na(OBS_hist),0.1,OBS_hist))
    }
    if (var=='MAX' | var=='MIN') {
      OBS_hist <- OBS_hist %>%
        mutate(OBS_hist = ifelse(is.na(OBS_hist),15,OBS_hist))
    }
    
    # forzamos a que coincida el numero de dias para ambas series
    if (nrow(var_hist_2) > nrow(var_model_2)) {
      var_model_2 <- rbind(var_model_2, data.frame(isodate=Sys.Date(),mode=var_model_2[nrow(var_model_2),2]))
    } 
    if (nrow(var_model_2) > nrow(var_hist_2)) {
      var_model_2 <- var_model_2 %>% 
        slice(1:nrow(var_hist_2))
    } 
    
    # creamos las variables a usar (solo para seguir el script de origen)
    OBS_hist <- var_hist_2 %>% 
      rename(OBS_hist = hist) 
    
    GCM_model <- var_model_2
    
    # aqui completo valores faltantes con 0.1
    if (var=='PR') {
      GCM_model <- GCM_model %>%
        mutate(mode = ifelse(is.na(mode),0.1,mode))
    }
    if (var=='MAX' | var=='MIN') {
      GCM_model <- GCM_model %>%
        mutate(mode = ifelse(is.na(mode),15,mode))
    }
    
    GCM_hist <- GCM_model

    data_at <- cbind(OBS_hist, GCM_model = GCM_model[,2]) %>% 
      read.zoo()
    data_wt <- cbind(OBS_hist, GCM_hist = GCM_hist[,2]) %>% 
      read.zoo()
    
    #============================================================================
    # APLICACIÓN DE LA TÉCNICA DE QUANTILE MAPPING (MÉTODO EMPÍRICO)
    
    data_wt$gcm_downscaled <- data_wt$GCM_hist
    
    seasons_by_year <- list(c("Diciembre","Enero","Febrero"), 
                            c("Marzo","Abril","Mayo"), 
                            c("Junio","Julio","Agosto"),
                            c("Setiembre","Octubre","Noviembre"))
    
    seasonal_qm_fit_model <- list()
    
    ## funcion de quantil maping    0.1 es a correcion, es linea
    for(i in 1:4) {
      obs_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$OBS_hist
      mod_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$GCM_hist
      #MODEL, read!: L. Gudmundsson et al. (2012)
      
      if (sum(mod_sl, na.rm = T)==0) {
        mod_sl[1]<- 0.001
      }
      
      qm_fit <- fitQmapQUANT(obs = coredata(obs_sl),
                             coredata(mod_sl),
                             qstep = 0.01,
                             nboot = 1, 
                             wet.day = 0, # Adaptado para temperatuas negativas
                             type = "linear")
      
      mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")
      data_wt$gcm_downscaled[ months(time(data_wt)) %in%
                                seasons_by_year[[i]]] <- mod_sl_qmapped
      
      seasonal_qm_fit_model[[i]] <- qm_fit 
    }
    
    #============================================================================
    #INTERPOLANDO LA INFORMACION EN CASO DE DATOS MENSUALES
    
    data_at$GCM_downscaled <- data_wt$gcm_downscaled
    
    t <- as.data.frame(data_at[,3])
    
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
   if (nrow(df_mod)> nrow(df_empty)) {
     df_mod <- df_mod %>% 
       slice(1:nrow(df_empty))
   }
  
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
  
  if (vec_var_3[z]=='PR') {
    # Seleccionamos el archivo de una variable
    var_hist <- data_pr_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    mult.86400 <- function(x, na.rm=FALSE) (x*86400)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585',vec_var_2[z],'.csv')) %>%
      mutate_if(is.numeric, mult.86400, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
  }
  
  if (vec_var_3[z]=='MAX') {
    # Seleccionamos el archivo de una variable
    var_hist <- data_tmax_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    dif.273 <- function(x, na.rm=FALSE) (x-273.15)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585',vec_var_2[z],'.csv')) %>%
      mutate_if(is.numeric, dif.273, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
  }

  if (vec_var_3[z]=='MIN') {
    # Seleccionamos el archivo de una variable
    var_hist <- data_tmax_his[,c('FECHA',estaciones[estacion])] %>% 
      mutate(FECHA= as.Date(FECHA))
    
    dif.273 <- function(x, na.rm=FALSE) (x-273.15)
    
    var_model <- read.csv(paste0(estaciones[estacion],'/RCP4585',vec_var_2[z],'.csv')) %>%
      mutate_if(is.numeric, dif.273, na.rm=FALSE) %>% 
      mutate(isodate = as.Date(isodate))
  }
  
  # Hacemos un bucle para generar la tabla en cada bloque de tiempo
  list_time <- list()
  for (i in 1:3) {
    list_time[[i]] <- ds(fecha_in_his, fecha_fin_his, fecha_in[i], fecha_fin[i], var_hist, var_model, var = vec_var_3[z])
  }

  # Exportamos
  name_file_out <- paste0(estaciones[estacion],'/',estaciones[estacion],'_',vec_var_3[z],'.xls')

  write.xlsx(list_time[[1]], file = name_file_out,row.names = F,
             sheetName = paste0(fecha_in_1,'-',fecha_fin_1), append = FALSE)
  
  write.xlsx(list_time[[2]], file = name_file_out, row.names = F,
             sheetName= paste0(fecha_in_2,'-',fecha_fin_2), append=TRUE)
  
  write.xlsx(list_time[[3]], file = name_file_out,row.names = F,
             sheetName= paste0(fecha_in_3,'-',fecha_fin_3), append=TRUE)

  rm(list_time)
  
  }

}
