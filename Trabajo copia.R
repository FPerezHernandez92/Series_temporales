#######################################
# 2.2 Pre-procesameinto
#######################################
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
print("Vemos como la frecuencia no es de 12 meses sinó cada 6 meses, así que vamos a modificar esta serie")
serie.ts <- ts(serie, frequency = 6)
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
title("Serie original")
lines(tiempoTs, serieTs, col="red")
cat("Dividimos el conjunto de datos en entrenamiento (para ajuste, negro) y test (para comprobar los modelos, rojo).\n");

#######################################
#2.3 Tendencia
#######################################
# Modelado y eliminación de tendencia
# Hipótesis: Modelo Lineal x= a + b*t (x=serie; t=tiempo; a,b=parámetros a estimar)
parametros.MLineal <- lm (serieTr ~ tiempoTr) # Ajustamos modelo
# Calculamos la estimación de la tendencia
TendEstimadaTr.MLineal<-parametros.MLineal$coefficients[1]+tiempoTr*parametros.MLineal$coefficients[2] 
TendEstimadaTs.MLineal<-parametros.MLineal$coefficients[1]+tiempoTs*parametros.MLineal$coefficients[2] 
# Mostramos en la misma figura la serie y la tendencia estimada
plot.ts(serieTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Modelo lineal")
lines(tiempoTr, TendEstimadaTr.MLineal, col="blue")
lines(tiempoTs, serieTs, col="red")
lines(tiempoTs, TendEstimadaTs.MLineal, col="green")
RSS.tendencia.MLineal= sum( (parametros.MLineal$residuals)^2); # Calculamos suma de errores al cuadrado
print("Metodo de eliminacion de tendencia por aproximacion lineal.")
cat(c("El modelo presenta un RSS (Residual Sum of Squares)=", RSS.tendencia.MLineal, "\n"));
JB <- jarque.bera.test(parametros.MLineal$residuals)
JB$p.value
JB <- jarque.bera.test((TendEstimadaTs.MLineal-serieTs))
JB$p.value
TT <- t.test(c(parametros.MLineal$residuals, TendEstimadaTs.MLineal-serieTs))
TT$p.value
print("Viendo estos resultados, podemos pensar, que al tener todos un p-value > 0.05, no existen diferencias significativas en los datos y la hipótesis lineal es factible")

print("Si observamos los resultados, no parece que sea el modelo más correcto para este caso, por lo que voy a probar con un flitrado")
# Estimación de tendencia por filtrado de medias moviles de orden k
for (k in 1:4) {
  filtro<-rep(1/k, k); # Creamos el filtro
  # Filtramos señal
  SerFiltradaTr<-filter(serieTr,filter=filtro,sides=2,method="convolution")
  # Mostramos en la misma figura la serie y la tendencia estimada
  series<-matrix(c(t(serieTr), t(SerFiltradaTr)), ncol=2);
  matplot(series, pch=1, type= "l")
  title("Método Filtrado train")
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
  title("Método Filtrado test")
  cat("Calculo tendencia con filtro de orden k=", k, "\n")
  print("Pulse una tecla para continuar...")
  pause<-readline(); # para pausar la ejecución
}
print("Vemos que, a mayor k, más se suaviza la serie. Vamos a quedarnos con el ultimo valor.")

# Eliminamos la tendencia
SerSinTendMLinTr <- serieTr - TendEstimadaTr.MLineal
SerSinTendMlinTs <- serieTs - TendEstimadaTs.MLineal
SerSinTendFilTr<-serieTr-SerFiltradaTr
SerSinTendFilTs<-serieTs-SerFiltradaTs
par(mfrow=c(2,2))  
plot.ts(SerSinTendFilTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Filtrado")
lines(tiempoTs, SerSinTendFilTs, col="red")
plot.ts(serieTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Original")
lines(tiempoTs, serieTs, col="red")
plot.ts(SerSinTendMLinTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Método Lineal")
lines(tiempoTs, SerSinTendMlinTs, col="red")
par(mfrow=c(1,1))
print("Vemos que la eliminación de tendencia es una posibilidad, pero vamos a pasar al siguiente paso sin eliminarla, ya que puedo asumir que no es elemental")


#######################################
#2.4 Estacionalidad
#######################################
# Calculamos y eliminamos la estacionalidad
k<- 6; # Asumimos periodo de estacionalidad k= 6
estacionalidad<- decompose(serie.ts)$seasonal[1:k];
#Eliminamos estacionalidad para el modelo ya que vemos que hay
aux<-rep(estacionalidad, length(serieTr)/length(estacionalidad));
serieTr.SinEst<- serieTr-aux;
serieTs.SinEst<- serieTs-estacionalidad;
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Serie sin estacionalidad")
lines(tiempoTs, serieTs.SinEst, col="red")
print("Serie sin la estacionalidad\n")

#######################################
#2.5 Estacionaridad
#######################################
#Vamos a ver si es estacionaria con el Test de Dickey-Fuller aumentado
adftest<- adf.test(serieTr.SinEst);
cat(c("Resultados del test ADF (aproximación lineal sin tendencia ni estacionalidad): ", adftest$p.value, "\n"));
cat("No es necesario diferenciar.\n")


#######################################
#2.6 Modelos ARIMA
#######################################
# Vemos los ACF y PACF
acf(serieTr.SinEst) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia lineal y sin estacionalidad\n")
pacf(serieTr.SinEst) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia lineal y sin estacionalidad")

# Mostramos las series sin estacionalidad
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Serie sin estacionalidad")
lines(tiempoTs[1:(length(tiempoTs))], serieTs.SinEst, col="red")
print("Serie sin la estacionalidad.")
cat("Con esto, podemos comenzar probando un modelo ARIMA(2, 0, 2):\n");

#######################################
#2.6.1 Modelo ARIMA
#######################################
# Ajustamos el modelo ARIMA(2,0,2)
modelo1<- arima(serieTr.SinEst, order=c(2, 0, 2))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados1<- serieTr.SinEst+modelo1$residuals; 
# Calculamos las predicciones 
Predicciones1<- predict(modelo1, n.ahead = NPred); 
valoresPredichos1<- Predicciones1$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr1<- sum((modelo1$residuals)^2);
errorTs1<- sum((valoresPredichos1-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(2, 0, 2): ", errorTr1, "\n")
cat("Error en la predicción de test con ARIMA(2, 0, 2): ", errorTs1, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados1, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos1, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest1<- Box.test(modelo1$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest1$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB1<- jarque.bera.test(modelo1$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB1$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW1<- shapiro.test(modelo1$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW1$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo1$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo1$residuals))

#######################################
#2.6.2 Modelos ARIMA
#######################################
# Ajustamos el modelo ARIMA(0,0,1)
modelo2<- arima(serieTr.SinEst, order=c(0, 0, 1))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados2<- serieTr.SinEst+modelo2$residuals; 
# Calculamos las predicciones 
Predicciones2<- predict(modelo2, n.ahead = NPred); 
valoresPredichos2<- Predicciones2$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr2<- sum((modelo2$residuals)^2);
errorTs2<- sum((valoresPredichos2-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(0, 0, 1): ", errorTr2, "\n")
cat("Error en la predicción de test con ARIMA(0, 0, 1): ", errorTs2, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados2, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos2, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest2<- Box.test(modelo2$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest2$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB2<- jarque.bera.test(modelo2$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB2$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW2<- shapiro.test(modelo2$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW2$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo2$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo2$residuals))

#######################################
#2.6.3 Modelos ARIMA
#######################################
# Ajustamos el modelo ARIMA(1,0,0)
modelo3<- arima(serieTr.SinEst, order=c(1, 0, 0))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados3<- serieTr.SinEst+modelo3$residuals; 
# Calculamos las predicciones 
Predicciones3<- predict(modelo3, n.ahead = NPred); 
valoresPredichos3<- Predicciones3$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr3<- sum((modelo3$residuals)^2);
errorTs3<- sum((valoresPredichos3-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(1, 0, 0): ", errorTr3, "\n")
cat("Error en la predicción de test con ARIMA(1, 0, 0): ", errorTs3, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados3, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos3, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest3<- Box.test(modelo3$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest3$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB3<- jarque.bera.test(modelo3$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB3$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW3<- shapiro.test(modelo3$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW3$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo3$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo3$residuals))

#######################################
#2.6.4 Modelos ARIMA
#######################################
# Ajustamos el modelo ARIMA(2,0,1)
modelo4<- arima(serieTr.SinEst, order=c(2, 0, 1))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados4<- serieTr.SinEst+modelo4$residuals; 
# Calculamos las predicciones 
Predicciones4<- predict(modelo4, n.ahead = NPred); 
valoresPredichos4<- Predicciones4$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr4<- sum((modelo4$residuals)^2);
errorTs4<- sum((valoresPredichos4-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(2, 0, 1): ", errorTr4, "\n")
cat("Error en la predicción de test con ARIMA(2, 0, 1): ", errorTs4, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados4, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos4, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest4<- Box.test(modelo4$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest4$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB4<- jarque.bera.test(modelo4$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB4$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW4<- shapiro.test(modelo4$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW4$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo4$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo4$residuals))

#######################################
#2.6.5 Modelos ARIMA
#######################################
# Ajustamos el modelo ARIMA(1,0,2)
modelo5<- arima(serieTr.SinEst, order=c(1, 0, 2))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados5<- serieTr.SinEst+modelo5$residuals; 
# Calculamos las predicciones 
Predicciones5<- predict(modelo5, n.ahead = NPred); 
valoresPredichos5<- Predicciones5$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr5<- sum((modelo5$residuals)^2);
errorTs5<- sum((valoresPredichos5-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(1, 0, 2): ", errorTr5, "\n")
cat("Error en la predicción de test con ARIMA(1, 0, 2): ", errorTs5, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados5, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos5, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest5<- Box.test(modelo5$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest5$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB5<- jarque.bera.test(modelo5$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB5$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW5<- shapiro.test(modelo5$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW5$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo5$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo5$residuals))

#######################################
#2.6.6 Modelos ARIMA
#######################################
# Ajustamos el modelo ARIMA(1,0,1)
modelo6<- arima(serieTr.SinEst, order=c(1, 0, 1))
# Cogemos los valores del ajuste y las predicciones
# Cogemos los valores que se han ajustado de la serie 
valoresAjustados6<- serieTr.SinEst+modelo6$residuals; 
# Calculamos las predicciones 
Predicciones6<- predict(modelo6, n.ahead = NPred); 
valoresPredichos6<- Predicciones6$pred; # Cogemos las predicciones
# Calculamos el error cuadrático acumulado del ajuste, en ajuste y en test
errorTr6<- sum((modelo6$residuals)^2);
errorTs6<- sum((valoresPredichos6-serieTs.SinEst)^2);
cat("Error en ajuste con ARIMA(1, 0, 1): ", errorTr6, "\n")
cat("Error en la predicción de test con ARIMA(1, 0, 1): ", errorTs6, "\n")

# Mostramos las gráficas del ajuste y predicción en test
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]))
title("Ajuste y predicción")
lines(valoresAjustados6, col="blue")
lines(tiempoTs, serieTs.SinEst, col="red")
lines(tiempoTs, valoresPredichos6, col="blue")
cat("Predicción con el modelo\n");

# Tests para la selección del modelo y su validación
boxtest6<- Box.test(modelo6$residuals) # Test de aleatoriedad de Box-Pierce
cat(c("El test de Box-Pierce da un p-value=", boxtest6$p.value, " para el modelo. Lo pasa (los errores son aleatorios)\n"))
JB6<- jarque.bera.test(modelo6$residuals); # Test de normalidad de Jarque Bera
cat(c("El test de Jarque Bera da un p-value=", JB6$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))
SW6<- shapiro.test(modelo6$residuals); # Test de normalidad de Shapiro-Wilk
cat(c("El test de Shapiro-Wilk da un p-value=", SW6$p.value, " para el modelo. Lo pasa (p-value>0.05)\n"))

# Mostramos histograma de residuos 
cat("Mostramos el histograma de los residuos del modelo.\n")
hist(modelo6$residuals, col="blue", prob=T,ylim=c(0,20),xlim=c(-0.2,0.2))
lines(density(modelo6$residuals))


#######################################
#2.7 Selección del modelo
#######################################
print("Probamos con 6 modelos: ARIMA(2,0,2), ARIMA(0,0,1), ARIMA(1,0,0), ARIMA(2,0,1), ARIMA(1,0,2) y ARIMA(1,0,1).");
print("Calculamos las predicciones y el modelo en si");

#Ajuste de los modelos
par(mfrow=c(2,3)) 
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 1")
lines(valoresAjustados1, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos1, col="blue")
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 2")
lines(valoresAjustados2, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos2, col="blue")
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 3")
lines(valoresAjustados3, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos3, col="blue")
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 4")
lines(valoresAjustados4, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos4, col="blue")
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 5")
lines(valoresAjustados5, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos5, col="blue")
plot.ts(serieTr.SinEst, xlim=c(1, tiempoTs[length(tiempoTs)]));
title("Ajuste Modelo 6")
lines(valoresAjustados6, col="red");
lines(tiempoTs, serieTs.SinEst, col="green")
lines(valoresPredichos6, col="blue")
print("Serie original (negro y verde), ajustada y predicha con los diferentes modelos")
par(mfrow=c(1,1)) 

#Histograma de los residuos
par(mfrow=c(2,3)) 
hist(modelo1$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo1$residuals))
hist(modelo2$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo2$residuals))
hist(modelo3$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo3$residuals))
hist(modelo4$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo4$residuals))
hist(modelo5$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo5$residuals))
hist(modelo6$residuals, col="blue", prob=T,ylim=c(0,4),xlim=c(-0.2,0.2))
lines(density(modelo6$residuals))
par(mfrow=c(1,1))
# Comparamos ambos modelos por el criterio de AIC
resultsAIC<-AIC(modelo1, modelo2, modelo3, modelo4, modelo5, modelo6)
print(resultsAIC)

#######################################
#2.8 Predicción de la serie
#######################################
# Probamos ahora a calibrar el mejor modelo con la serie completa
tiempo<- 1:length(serie)
# Calculamos estacionalidad
k<-6
estacionalidad<-rep(0, k);
for (i in 1:k) {
  secuencia<-seq(i, length(serie), by=k);
  for (j in secuencia) {
    estacionalidad[i]<- estacionalidad[i] + serie[j];
  }
  
  estacionalidad[i]<-estacionalidad[i]/length(secuencia);
}
aux<-rep(estacionalidad, length(serie)/length(estacionalidad));
#Eliminamos estacionalidad
serieSinEst<- serie-aux;
# Ajustamos el modelo que hemos seleccionado
modelo<- arima(serieSinEst, order=c(0, 0, 1))
# Obtenemos ajuste y predicción
valoresAjustados<- serieSinEst+modelo$residuals;
Predicciones<- predict(modelo, n.ahead = NPred); 
valoresPredichos<- Predicciones$pred; # Cogemos las predicciones
# Por último, deshacemos cambios
valoresAjustados<- valoresAjustados+aux; # Estacionalidad
valoresPredichos<- valoresPredichos+estacionalidad;
tiempoPred<- (tiempo[length(tiempo)]+(1:NPred));
plot.ts(serie, xlim=c(1, max(tiempoPred)), ylim=c(-0.2, 1.5))
title("Predicción Final")
lines(valoresAjustados, col="blue")
lines(valoresPredichos, col="red")
