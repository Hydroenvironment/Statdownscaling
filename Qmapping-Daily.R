###############################################################################
#DOWNSCALLING ESTADÍSTICO DE DATOS CLIMÁTICOS USANDO MAPEO DE QUANTILES (QMAP)#
###############################################################################
#Referencias:
#Technical Note: Downscaling RCM precipitation to the station scale using statistical transformations - a comparison of methods 
#https://www.hydrol-earth-syst-sci.net/16/3383/2012/
#Qmap package information - CRAN:
#https://cran.r-project.org/web/packages/qmap/qmap.pdf
#==============================================================================
# INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
setwd("C:/Users/Julio/Documents/INVESTIGACION/R SCRIPTS/Downscaling/Diario")
getwd()
ls()
rm(list=ls())
#por ejemplo:
#install.packages(c("qmap", "zoo","latticeExtra"))
  
# CARGAR LAS LIBRERIAS
library(qmap)
library(zoo)
library(latticeExtra)

#=============================================================================
#CARGANDO LOS ARCHIVOS NECESARIOS

# Observacion historica PISCO
OBS_hist <- read.zoo("0Estacion_C2.csv", header = TRUE, sep = ",",
                     format = "%Y-%m-%d")
# Observacion del modelo GCM
GCM_model <- read.zoo("ACCESS1-0-RCP45.csv", header = TRUE, sep = ",", 
                      format = "%Y-%m-%d")

# Ajustando la serie hasta el final de la data hist?rica
GCM_hist <- window(GCM_model, end = "2039-12-31") 

data_at <- cbind(OBS_hist, GCM_model)
data_wt <- cbind(OBS_hist, GCM_hist)

#============================================================================
#GRAFICANDO LAS SERIES DE TIEMPO
plot(data_at, plot.type = "single", col = c(1, 2), lwd = 0.01
     , ylab = "pp (mm/day)", xlab = "years")
title(expression(bold("OBSERVED DAILY HISTORICAL DATA vs"*phantom("GCM DAILY HISTORICAL DATA"))), col.main = "black")
title(expression(bold(phantom("OBSERVED DAILY HISTORICAL DATA vs ")*"GCM DAILY HISTORICAL DATA")), col.main = "red")

plot(data_wt, plot.type = "single", col = c(1, 2),  lwd = 0.01
     , ylab = "pp (mm/day)", xlab = "years")
title(expression(bold("OBSERVED DAILY DATA vs"*phantom("GCM DAILY DATA"))), col.main = "black")
title(expression(bold(phantom("OBSERVED DAILY DATA vs ")*"GCM DAILY DATA")), col.main = "red")


# months(time(data_wt))   correr y luego cambiar de acuerdo al nombre de los meses 
plot(data_wt[months(time(data_wt)) %in% c("December","January","February","March")], 
     plot.type = "single", col = c(1, 2),  lwd = 1, type='p'
     , ylab = "pp (mm/day)", xlab = "years")
title(expression(bold("OBSERVED DAILY DATA VS "*phantom("GCM DAILY DATA")*"DURING RAINY SEASON (DEC-JAN-FEB-MAR)")), col.main = "black", cex.main=1)
title(expression(bold(phantom("OBSERVED DAILY DATA vs ")*"GCM DAILY DATA"*phantom("DURING RAINY SEASON (DEC-JAN-FEB-MAR)"))), col.main = "red", cex.main=1)


# GRAFICANDO SCATTERPLOT
plot(GCM_hist~OBS_hist, coredata(data_wt), col = c(1,2), ylab = "PP GCM (mm)", xlab = "PP OBSERVED (mm)")
title(expression(bold("SCATTER PLOT OF OBSERVED DAILY DATA VS "*phantom("GCM DAILY DATA"))), col.main = "black", cex.main=1)
title(expression(bold(phantom("SCATTER PLOT OF OBSERVED DAILY DATA vs ")*"GCM DAILY DATA")), col.main = "red", cex.main=1)

plot(GCM_hist~OBS_hist,
     coredata(data_wt[months(time(data_wt)) %in% 
                        c("December","January","February","March")]), col = c(1,2),ylab = "PP GCM (mm)", xlab = "PP OBSERVED (mm)")
title(expression(bold("SCATTER PLOT OF OBSERVED DAILY DATA VS"*phantom(" GCM DAILY DATA ")* "DURING RAINY SEASON (DEC-JAN-FEB-MAR)")), col.main = "black", cex.main=0.8)
title(expression(bold(phantom("SCATTER PLOT OF OBSERVED DAILY DATA vs ")*" GCM DAILY DATA "* phantom("DURING RAINY SEASON (DEC-JAN-FEB-MAR)"))), col.main = "red", cex.main=0.8)

# GRAFICANDO ECDF
ecdfplot(~ OBS_hist +  GCM_hist, data = data.frame(data_wt), lwd = 2, col = c(1, 2),main="CDF PLOT OF OBSERVED DAILY DATA VS GCM DAILY DATA",ylab = "Empirical CDF", xlab = "Daily GCM (red) and Observed (black) data")

ecdfplot(~ OBS_hist +  GCM_hist,
         data = data.frame(data_wt[months(time(data_wt)) %in% 
                                     c("December","January","February")]), 
         lwd = 2, col = c(1, 2),main="CDF PLOT OF DAILY OBS. VS GCM DATA DURING RAINFALL SEASON",ylab = "Empirical CDF", xlab = "Daily GCM (red) and Observed (black) data")

#============================================================================
# APLICACIÓN DE LA TÉCNICA DE QUANTILE MAPPING (MÉTODO EMPÍRICO)

data_wt$gcm_downscaled <- data_wt$GCM_hist

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
                         wet.day = 0, # Adaptado para temperatuas negativas
                         type = "linear")
  
  mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")
  data_wt$gcm_downscaled[ months(time(data_wt)) %in%
                           seasons_by_year[[i]]] <- mod_sl_qmapped
  
  seasonal_qm_fit_model[[i]] <- qm_fit 
}

#============================================================================
# GR?FICOS DE SERIES TEMPORALES
plot(data_wt, plot.type = "single", col = c(1, 2, 4),  lwd = 1,xlab = "Years", ylab = "pp (mm/day)")
title(expression(bold("OBSERVED DAILY HISTORICAL DATA vs" * phantom(" GCM DAILY DATA ") * phantom("vs DOWNSCALED GCM DAILY DATA"))), col.main = "black",cex.main=0.9)
title(expression(bold(phantom("OBSERVED DAILY HISTORICAL DATA vs ") * " GCM DAILY DATA " * phantom("vs DOWNSCALED GCM DAILY DATA"))), col.main = "red",cex.main=0.9)
title(expression(bold(phantom("OBSERVED DAILY HISTORICAL DATA vs ") * phantom(" GCM DAILY DATA ") *"vs DOWNSCALED GCM DAILY DATA")), col.main = "dodgerblue3",cex.main=0.9)

plot(data_wt[months(time(data_wt)) %in% c("December","January","February")], plot.type = "single", col = c(1, 2, 4),  
     lwd = 2, type = "p",xlab = "Years", ylab = "pp (mm/mes)")
title(expression(bold("OBSERVED DAILY HISTORICAL DATA vs" * phantom(" GCM DAILY DATA ") * phantom(" vs DOWNSCALED GCM DAILY DATA ")*"DURING RAINY SEASON")), col.main = "black",cex.main=0.73)
title(expression(bold(phantom("OBSERVED DAILY HISTORICAL DATA vs ") * " GCM DAILY DATA " * phantom(" vs DOWNSCALED GCM DAILY DATA ")*phantom("DURING RAINY SEASON"))), col.main = "red",cex.main=0.73)
title(expression(bold(phantom("OBSERVED DAILY HISTORICAL DATA vs ") * phantom(" GCM DAILY DATA ")*" vs DOWNSCALED GCM DAILY DATA "*phantom("DURING RAINY SEASON"))), col.main = "dodgerblue3",cex.main=0.73)

plot(gcm_downscaled~OBS_hist, coredata(data_wt), col = c(4,1))
title(expression(bold("SCATTER PLOT OF OBSERVED DAILY DATA VS "*phantom("GCM DAILY DATA"))), col.main = "black", cex.main=1)
title(expression(bold(phantom("SCATTER PLOT OF OBSERVED DAILY DATA vs ")*"GCM DAILY DATA")), col.main = "red", cex.main=1)

plot(gcm_dowscaled~OBS_hist, coredata(data_wt[months(time(data_wt))%in%c("December","January","February")]), col = c(4,1),ylab = "PP GCM DOWNSCALED (mm)", xlab = "PP OBSERVED (mm)")
title(expression(bold("SCATTER PLOT OF OBSERVED DAILY DATA VS "*phantom("GCM DOWNSCALED DAILY DATA")*"DURING RAINY SEASON")), col.main = "black", cex.main=0.8)
title(expression(bold(phantom("SCATTER PLOT OF OBSERVED DAILY DATA vs ")*"GCM DOWNSCALED DAILY DATA"*phantom("DURING RAINY SEASON"))), col.main = "dodgerblue3", cex.main=0.8)


ecdfplot(~ OBS_hist +  GCM_hist + gcm_downscaled, data = data.frame(data_wt), lwd = 3, col = c(1, 2, 4),
         main="CDF PLOT OF OBS. DAILY DATA VS GCM AND GCM DOWNSCALED DAILY DATA",
         ylab = "Empirical CDF", xlab = "Daily GCM (red), GCM downscaled (blue) and Observed (black) data")

ecdfplot(~ OBS_hist +  GCM_hist + gcm_downscaled,data = data.frame(data_wt[months(time(data_wt))
                                                                          %in%c("December","January","February")]),
         lwd = 3, col = c(1, 2, 4), main="CDF PLOT OF OBS.vsGCMvsGCM DOWNSC. DAILY DATA - RAINY SEASON",
ylab = "Empirical CDF", xlab = "Daily GCM (red), GCM downscaled (blue) and Observed (black) data")

ecdfplot(~ OBS_hist + gcm_downscaled, data = data.frame(data_wt), lwd = 3, col = c(1, 2, 4), 
         main="CDF PLOT OF OBS. DAILY DATA VS GCM DOWNSCALED DAILY DATA", ylab = "Empirical CDF", 
         xlab = "Daily GCM downscaled (red) and Observed (black) data")

#============================================================================
#INTERPOLANDO LA INFORMACION EN CASO DE DATOS MENSUALES

data_at$GCM_downscaled <- data_wt$gcm_downscaled

#for(i in 1:4) {
#mod_sl <- data_at[months(time(data_at)) %in% seasons_by_year[[i]]]$GCM_model
#mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), seasonal_qm_fit_model[[i]], type = "linear")
#data_at$GCM_downscaled[ months(time(data_at)) %in% seasons_by_year[[i]]] <- mod_sl_qmapped
  #}
# verificar la nueva tabla y el grafico
View(data_at)
View(data_wt)

plot(data_at, plot.type = "single", col = c(1, 2, 4), lwd = 1, 
     main = c("OBS vs GCM VS GCM DOWNSCALED"),
     ylab = "pp (mm/mes)", xlab = "years")

# ..............................................................................
# GUARDAR DATOS ESCALADOS EN UN ARCHIVO .CSV

write.zoo(data_at[,3],file = "OUT.csv", sep = ",")
print("El proceso se ha completado satisfactoriamente!")
plot(data_at, plot.type = "single", col = c(1, 2, 4), lwd = 1, 
     main = c("OBS vs GCM VS GCM DOWNSCALED"),
     ylab = "pp (mm/mes)", xlab = "years")

# ..............................................................................
# GUARDAR DATOS ESCALADOS EN UN ARCHIVO .CSV

write.zoo(data_at[,3],file = "Out2.csv", sep = ",")
print("El proceso se ha completado satisfactoriamente!")
