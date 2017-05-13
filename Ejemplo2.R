rm(list=ls()) # Eliminamos lo que haya en el espacio de trabajo
library("tseries") # para el test ADF

# Cargamos serie con tendencia y estacionalidad
serie<-scan("ModeloAR2_TendyEst_simulado.dat")
plot.ts(serie) # la mostramos
print("Mostrando la serie con tendencia y estacionalidad")
acf(serie) # Mostramos ACF
print("Mostrando el ACF")
pacf(serie) # Mostramos PACF
print("Mostrando el PACF")
resul <- adf.test(serie) # Pasamos el test de ADF
cat("El resultado del test de Dickey-Fuller aumentado es un p-value=", resul$p.value, "\n")


#----------------------------------------------------
# Eliminar tendencia con aproximac. funcional
#----------------------------------------------------
# Hipótesis 1: Modelo Lineal x= a + b*t (x=serie; t=tiempo; a,b=parámetros a estimar)
tiempo<- 1:length(serie) # Creamos la variable que modela al tiempo
parametros <- lm (serie ~ tiempo) # Ajustamos modelo lineal
# Calculamos la estimación de la tendencia
TendEstimada<-parametros$coefficients[1]+tiempo*parametros$coefficients[2] 
# Mostramos en la misma figura la serie y la tendencia estimada
series<-matrix(c(t(serie), t(TendEstimada)), ncol=2); # Mostramos resultado
matplot(series, pch=1, type= "l")
print("Metodo de eliminacion de tendencia por aproximacion funcional.")
print("La hipotesis de modelado lineal parece que no sirve aqui.")

# Hipótesis 2: Modelo polinomico x= a + b*t +c*t^2 (x=serie; t=tiempo; a,b,c=parámetros a estimar)
tiempo<- 1:length(serie) # Creamos la variable que modela al tiempo
parametros <- lm (serie ~ tiempo + I(tiempo^2)) # Ajustamos modelo lineal
# Calculamos la estimación de la tendencia
TendEstimada<-parametros$coefficients[1]+tiempo*parametros$coefficients[2] + (tiempo^2)*parametros$coefficients[3] 
# Mostramos en la misma figura la serie y la tendencia estimada
series<-matrix(c(t(serie), t(TendEstimada)), ncol=2); # Mostramos resultado
matplot(series, pch=1, type= "l")
print("Metodo de eliminacion de tendencia por aproximacion funcional.")
print("La hipotesis de modelado polinomico de 2o. grado parece que va mejor.")

# Eliminamos la tendencia
SerSinTend<-serie-TendEstimada;
plot.ts(SerSinTend)
print("Serie sin la tendencia.")

#----------------------------------------------------
# Eliminar tendencia con filtrado
#----------------------------------------------------
# Estimación de tendencia por filtrado de medias moviles de orden k
for (k in 3:5) {
  filtro<-rep(1/k, k); # Creamos el filtro
  # Filtramos señal
  SerFiltrada<-filter(serie,filter=filtro,sides=2,method="convolution")
  # Mostramos en la misma figura la serie y la tendencia estimada
  series<-matrix(c(t(serie), t(SerFiltrada)), ncol=2);
  matplot(series, pch=1, type= "l")
  cat("Calculo tendencia con filtro de orden k=", k, "\n")
  print("Pulse una tecla para continuar...")
  #pause<-readline(); # para pausar la ejecución
}
print("Vemos que, a mayor k, más se suaviza la serie. Vamos a continuar hasta que tengamos una tendencia clara.")
for (k in 6:15) {
  filtro<-rep(1/k, k); # Creamos el filtro
  # Filtramos señal
  SerFiltrada<-filter(serie,filter=filtro,sides=2,method="convolution")
  # Mostramos en la misma figura la serie y la tendencia estimada
  series<-matrix(c(t(serie), t(SerFiltrada)), ncol=2);
  matplot(series, pch=1, type= "l")
  cat("Calculo tendencia con filtro de orden k=", k, "\n")
  print("Pulse una tecla para continuar...")
  #pause<-readline(); # para pausar la ejecución
}
print("Paramos, por ejemplo, con k=15.")
# Eliminamos la tendencia
SerSinTend2<-serie-SerFiltrada;
plot.ts(SerSinTend2)
print("Serie sin la tendencia (calculada mediante filtros).")
# Mostramos en la misma figura la serie sin tendencia mediante los 2 metodos
series<-matrix(c(t(SerSinTend), t(SerSinTend2)), ncol=2);
matplot(series, pch=1, type= "l")
print("Aqui vemos como queda la serie tras quitar la tendencia con los 2 metodos.")

# Quitamos los NA de la serie sin tendencia con filtro
aux<- SerSinTend2[!is.na(SerSinTend2)];

acf(SerSinTend) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia con el método de aprox. funcional")

acf(aux) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia con el método de filtrado")

pacf(SerSinTend) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia con el método de aprox. funcional")

pacf(aux) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia con el método de filtrado")
