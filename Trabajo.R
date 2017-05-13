rm(list=ls()) # Eliminamos lo que haya en el espacio de trabajo
library("tseries") # para el test ADF

#Leemos la serie
serie<-scan("SerieTrabajoPractico.dat")
# Echamos un primer vistazo a la serie.
# Como se trata de series anuales, de haber estacionalidad, posiblemente aparezca cada 12 meses.
# Por tanto, inicialmente creamos la serie temporal con periodo de estacionalidad 12, para visualizar
serie.ts<- ts(serie, frequency = 12)
# Visualizamos la descomposición
plot(decompose(serie.ts))
cat("Visualizamos primero la descomposición de la serie, buscando patrones visuales que nos den idea de por dónde empezar (tendencia, estacionalidad)\n");
cat("Vemos 3 cosas: \n");
cat("  1. Que la variabilidad de la estacionalidad no aumenta ni decrece. ");
cat("  2. Que puede que la tendencia no juegue un papel importante.\n")
cat("  3. Que posiblemente hay una estacionalidad.\n")
print("Como con esto no estamos seguros de si la serie es estacionaria o no, vamos a ver el ACF, PACF y pasaremos el test ADF:")
print("Mostrando el ACF")
acf(serie) # Mostramos ACF
print("Mostrando el PACF")
pacf(serie) # Mostramos PACF
resul <- adf.test(serie) # Pasamos el test de ADF
cat("El resultado del test de Dickey-Fuller aumentado es un p-value=", resul$p.value, "\n")
print("Podemos decir debido a este resultado, que la serie es Estacionaria con un nivel de confianza del 99%. Además en el gráfico ACF vemos como con un Lag = 6, que la serie va a ser Estacionaria")

# Dividimos la serie en training y test (nos quedamos con los NTest últimos para el test)
NPred= 12; # Valores a predecir
NTest= 12; # Valores que vamos a dejar para test
serieTr<- serie[1:(length(serie)-NTest)];
tiempoTr<- 1:length(serieTr)
serieTs<- serie[(length(serie)-NTest+1):length(serie)];
tiempoTs<- (tiempoTr[length(tiempoTr)]+1):(tiempoTr[length(tiempoTr)]+NTest);
plot.ts(serieTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, serieTs, col="red")
cat("Dividimos el conjunto de datos en entrenamiento (para ajuste, negro) y test (para comprobar los modelos, rojo).\n");

#######################################
#Tendencia
#######################################
# Modelado y eliminación de tendencia
# Hipótesis: Modelo Lineal x= a + b*t (x=serie; t=tiempo; a,b=parámetros a estimar)
parametros.H1 <- lm (serieTr ~ tiempoTr) # Ajustamos modelo
# Calculamos la estimación de la tendencia
TendEstimadaTr.H1<-parametros.H1$coefficients[1]+tiempoTr*parametros.H1$coefficients[2] 
TendEstimadaTs.H1<-parametros.H1$coefficients[1]+tiempoTs*parametros.H1$coefficients[2] 
# Mostramos en la misma figura la serie y la tendencia estimada
plot.ts(serieTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTr, TendEstimadaTr.H1, col="blue")
lines(tiempoTs, serieTs, col="red")
lines(tiempoTs, TendEstimadaTs.H1, col="green")
RSS.tendencia.H1= sum( (parametros.H1$residuals)^2); # Calculamos suma de errores al cuadrado
print("Metodo de eliminacion de tendencia por aproximacion lineal.")
cat(c("El modelo presenta un RSS (Residual Sum of Squares)=", RSS.tendencia.H1, "\n"));
JB <- jarque.bera.test(parametros.H1$residuals)
JB$p.value
JB <- jarque.bera.test((TendEstimadaTs.H1-serieTs))
JB$p.value
TT <- t.test(c(parametros.H1$residuals, TendEstimadaTs.H1-serieTs))
TT$p.value
print("Viendo estos resultados, podemos pensar, que al tener todos un p-value > 0.05, no existen diferencias significativas en los datos y la hipótesis lineal es factible")

print("Si observamos los resultados, no parece que sea el modelo correcto para este caso, por lo que voy a probar con un flitrado")
# Estimación de tendencia por filtrado de medias moviles de orden k
for (k in 1:4) {
  filtro<-rep(1/k, k); # Creamos el filtro
  # Filtramos señal
  SerFiltradaTr<-filter(serieTr,filter=filtro,sides=2,method="convolution")
  # Mostramos en la misma figura la serie y la tendencia estimada
  series<-matrix(c(t(serieTr), t(SerFiltradaTr)), ncol=2);
  matplot(series, pch=1, type= "l")
  cat("Calculo tendencia con filtro de orden k=", k, "\n")
  print("Pulse una tecla para continuar...")
  pause<-readline(); # para pausar la ejecución
}
for (k in 1:3) {
  filtro<-rep(1/k, k); # Creamos el filtro
  # Filtramos señal
  SerFiltradaTs<-filter(serieTs,filter=filtro,sides=2,method="convolution")
  # Mostramos en la misma figura la serie y la tendencia estimada
  series<-matrix(c(t(serieTs), t(SerFiltradaTs)), ncol=2);
  matplot(series, pch=1, type= "l")
  cat("Calculo tendencia con filtro de orden k=", k, "\n")
  print("Pulse una tecla para continuar...")
  pause<-readline(); # para pausar la ejecución
}
print("Vemos que, a mayor k, más se suaviza la serie. Vamos a quedarnos con el ultimo valor.")
# Eliminamos la tendencia
SerSinTend2Tr<-serieTr-SerFiltradaTr
SerSinTend2Ts<-serieTs-SerFiltradaTs
plot.ts(SerSinTend2Tr, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, SerSinTend2Ts, col="red")
print("Serie sin la tendencia (calculada mediante filtros).")

print("Vemos que la eliminación de tendencia es una posibilidad, pero vamos a pasar al siguiente paso sin eliminarla, ya que puedo asumir que no es elemental")


#######################################
#Tendencia
#######################################
# Calculamos y eliminamos la estacionalidad
k<- 12; # Asumimos periodo de estacionalidad k= 12
estacionalidad.H1<- decompose(serie.ts)$seasonal[1:k];
#Eliminamos estacionalidad para el modelo
aux<-rep(estacionalidad.H1, length(serieTr)/length(estacionalidad.H1));
serieTr.SinEst.H1<- serieTr-aux;
serieTs.SinEst.H1<- serieTs-estacionalidad.H1;
plot.ts(serieTr.SinEst.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, serieTs.SinEst.H1, col="red")
print("Serie sin la estacionalidad\n")


#Vamos a ver si es estacionaria con el Test de Dickey-Fuller aumentado
adftest.H1<- adf.test(serieTr.SinEst.H1);
cat(c("Resultados del test ADF (aproximación lineal sin tendencia ni estacionalidad): ", adftest.H1$p.value, "\n"));
cat("No es necesario diferenciar.\n")

# Vemos los ACF y PACF
acf(serieTr.SinEst.H1) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia lineal y sin estacionalidad\n")
pacf(serieTr.SinEst.H1) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia lineal y sin estacionalidad")


# Mostramos las series sin tendencia ni estacionalidad
plot.ts(serieTr.SinEst.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs[1:(length(tiempoTs))], serieTs.SinEst.H1, col="red")
print("Serie sin la estacionalidad.")
cat("Con esto, podemos comenzar probando un modelo ARIMA(1, 0, 1):\n");








# Ajustamos el modelo
modelo.H1<- arima(serieTr.SinEst.H1, order=c(1, 0, 1))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados.H1<- serieTr.SinEst.H1+modelo.H1$residuals; 
# Calculamos las predicciones 
Predicciones.H1<- predict(modelo.H1, n.ahead = NPred); 
valoresPredichos.H1<- Predicciones.H1$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr.H1<- sum((modelo.H1$residuals)^2);
errorTs.H1<- sum((valoresPredichos.H1-serieTs.SinEst.H1)^2);
cat("Error en ajuste con ARIMA(1, 0, 1): ", errorTr.H1, "\n")
cat("Error en la predicción de test con ARIMA(1, 0, 1): ", errorTs.H1, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(valoresAjustados.H1, col="blue")
lines(tiempoTs, serieTs.SinEst.H1, col="red")
lines(tiempoTs, valoresPredichos.H1, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtestM1<- Box.test(modelo.H1$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtestM1$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB.H1<- jarque.bera.test(modelo.H1$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB.H1$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW.H1<- shapiro.test(modelo.H1$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW.H1$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo.H1$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo.H1$residuals))

################################
# Si tuviésemos más modelos, deberíamos compararlos de la siguiente forma:
#  1. Calculando el criterio de AIC
#  2. Calculando la predicción de la serie en test, y comprobando cuál modelo
# produce mejores resultados al predecir los valores del test
#  3. Seleccionando el modelo más apropiado (el más simple, el que mejor predicción proporcione...)
################################
