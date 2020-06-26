###############################################################################
DOWNSCALLING ESTADÍSTICO DE DATOS CLIMÁTICOS USANDO MAPEO DE QUANTILES (QMAP)##
###############################################################################
Referencias:
Technical Note: Downscaling RCM precipitation to the station scale using statistical transformations – a comparison of methods 
https://www.hydrol-earth-syst-sci.net/16/3383/2012/
Qmap package information - CRAN:
https://cran.r-project.org/web/packages/qmap/qmap.pdf
#==============================================================================
# INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
setwd("C:/Users/Julio/Documents/INVESTIGACION/R SCRIPTS/Downscaling/Diario")
getwd()
ls()
rm(list=ls())
por ejemplo:
#install.packages(c("qmap", "zoo", "latticeExtra"))

# CARGAR LAS LIBRERIAS
library(qmap)
library(zoo)
library(latticeExtra)

#=============================================================================
#CARGANDO LOS ARCHIVOS NECESARIOS

# Observacion historica PISCO
OBS_hist <- read.zoo("Data_historica.csv", header = TRUE, sep = ",",
                     format = "%Y-%m-%d")
# Observacion del modelo GCM
GCM_model <- read.zoo("Data_Descargada.csv", header = TRUE, sep = ",", 
                      format = "%Y-%m-%d")

# Ajustando la serie hasta diciembre 2005
GCM_hist <- window(GCM_model, end = "2010-12-31") 

data_at <- cbind(OBS_hist, GCM_model)
data_wt <- cbind(OBS_hist, GCM_hist)

#============================================================================
#GRAFICANDO LAS SERIES DE TIEMPO
plot(data_at, plot.type = "single", col = c(1, 2), lwd = 1, 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)", xlab = "years")

plot(data_wt, plot.type = "single", col = c(1, 2),  lwd = 1, 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)", xlab = "years")

# months(time(data_wt))   correr y luego cambiar de acuerdo al nombre de los meses 
plot(data_wt[months(time(data_wt)) %in% c("December","January","February")], 
     plot.type = "single", col = c(1, 2),  lwd = 1, type='p', 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)", xlab = "years")

# GRAFICANDO SCATTERPLOT
plot(GCM_hist~OBS_hist, coredata(data_wt), col = c(1,2))

plot(GCM_hist~OBS_hist,
     coredata(data_wt[months(time(data_wt)) %in% 
                        c("December","January","February")]), col = c(1,2))

# GRAFICANDO ECDF
ecdfplot(~ OBS_hist +  GCM_hist, data = data.frame(data_wt), 
         lwd = 2, col = c(1, 2))

ecdfplot(~ OBS_hist +  GCM_hist,
         data = data.frame(data_wt[months(time(data_wt)) %in% 
                                     c("December","January","February")]), 
         lwd = 2, col = c(1, 2))

#============================================================================
# APLICACIÓN DE LA TÉCNICA DE QUANTILE MAPPING (EMPÍRICO)

data_wt$gcm_dowscaled <- data_wt$GCM_hist

seasons_by_year <- list(c("December","January","February"), 
                        c("March","April","May"), 
                        c("June","July","August"),
                        c("September","October","November"))

seasonal_qm_fit_model <- list()

## funcion de quantil maping    0.1 es a correcion, es linea
for(i in 1:4) {
  obs_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$OBS_hist
  mod_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$GCM_hist
  #MODEL, read!: L. Gudmundsson et al. (2012)
  qm_fit <- fitQmapQUANT(obs = coredata(obs_sl),
                         coredata(mod_sl),
                         qstep = 0.01,
                         nboot = 1, 
                         wet.day = TRUE, # Adaptado para temperatuas negativas
                         type = "linear")
  
  mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")
  data_wt$gcm_dowscaled[ months(time(data_wt)) %in%
                           seasons_by_year[[i]]] <- mod_sl_qmapped
  
  seasonal_qm_fit_model[[i]] <- qm_fit 
}

#============================================================================
# GRÁFICOS DE SERIES TEMPORALES
plot(data_wt, plot.type = "single", col = c(1, 2, 4),  lwd = 2, 
     main = c("OBS vs GCM vs GCM downscaled"),
     xlab = "Years", ylab = "pp (mm/mes)")

plot(data_wt[months(time(data_wt)) %in% c("December",
                                          "January",
                                          "February")], 
     plot.type = "single", col = c(1, 2, 4),  lwd = 2, type = "p",
     main = c("OBS vs GCM vs GCM downscaled"), 
     xlab = "Years", ylab = "pp (mm/mes)")

plot(gcm_dowscaled~OBS_hist, coredata(data_wt), col = c(4,1))

plot(gcm_dowscaled~OBS_hist, coredata(data_wt[months(time(data_wt)) %in%
                                                c("December",
                                                  "January",
                                                  "February")]), 
     col = c(4,1))

ecdfplot(~ OBS_hist +  GCM_hist + gcm_dowscaled, data = data.frame(data_wt), 
         lwd = 3, col = c(1, 2, 4))

ecdfplot(~ OBS_hist +  GCM_hist + gcm_dowscaled,
         data = data.frame(data_wt[months(time(data_wt)) %in%
                                     c("December","January","February")]), 
         lwd = 3, col = c(1, 2, 4))

ecdfplot(~ OBS_hist + gcm_dowscaled, data = data.frame(data_wt), 
         lwd = 3, col = c(1, 2, 4))

#============================================================================
#INTERPOLANDO LA INFORMACION

data_at$GCM_downscaled <- data_at$GCM_model

for(i in 1:4) {
  mod_sl <- data_at[months(time(data_at)) %in% seasons_by_year[[i]]]$GCM_model
  mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), seasonal_qm_fit_model[[i]], type = "linear")
  data_at$GCM_downscaled[ months(time(data_at)) %in% seasons_by_year[[i]]] <- mod_sl_qmapped
}
# verificar la nueva tabla y el grafico
View(data_at)

plot(data_at, plot.type = "single", col = c(1, 2, 4), lwd = 1, 
     main = c("OBS vs GCM VS GCM DOWNSCALED"),
     ylab = "pp (mm/mes)", xlab = "years")

# ..............................................................................
# GUARDAR DATOS ESCALADOS EN UN ARCHIVO .CSV

write.zoo(data_at[,3], file = "GCM_DOWNSCALED_1981_2100_RCP85.csv", sep = ",")
print("El proceso se ha completado satisfactoriamente!")
