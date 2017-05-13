# Serie: AirPassengers (fuente: data(AirPassengers))
# Monthly totals of international airline passengers, 1949 to 1960
# La serie son medidas anuales, muestreadas mensualmente, del número de pasajeros de avión (en miles) desde 1949 hasta 1960
# Se proporciona: Datos entre 1949 y 1959
# Se pide: Proporcionar la predicción del número de pasajeros para todos los meses de 1960
rm(list=ls()) # Eliminamos lo que haya en el espacio de trabajo
library("tseries") # para el test ADF

NPred= 12; # Valores a predecir
NTest= 12; # Valores que vamos a dejar para test
# Cargamos serie y datos reales a predecir 
serie<-scan("pasajeros_1949_1959.dat")
# Echamos un primer vistazo a la serie.
# Como se trata de series anuales, de haber estacionalidad, posiblemente aparezca cada 12 meses.
# Por tanto, inicialmente creamos la serie temporal con periodo de estacionalidad 12, para visualizar
serie.ts<- ts(serie, frequency = 12)
# Visualizamos la descomposición
plot(decompose(serie.ts))
cat("Visualizamos primero la descomposición de la serie, buscando patrones visuales ");
cat("que nos den idea de por dónde empezar (tendencia, estacionalidad)\n");
cat("Vemos 3 cosas: \n");
cat("  1. Que la variabilidad de la estacionalidad aumenta, provocando que la serie sea no estacionaria por varianza. ");
cat("Para solucionarlo, aplicamos log() a la serie.\n")
cat("  2. Que hay una tendencia creciente.\n")
cat("  3. Que hay una estacionalidad.\n")
# Aplicamos logaritmo para suavizar el problema de no estacionaridad por varianza
serie.ts<- log(serie.ts);
serie.log<- log(serie);
# Visualizamos de nuevo
plot(decompose(serie.ts))
cat("Solucionado el problema de la no estacionaridad en varianza.\n");
# Dividimos la serie en training y test (nos quedamos con los NTest últimos para el test)
serieTr<- serie.log[1:(length(serie.log)-NTest)];
tiempoTr<- 1:length(serieTr)
serieTs<- serie.log[(length(serie.log)-NTest+1):length(serie)];
tiempoTs<- (tiempoTr[length(tiempoTr)]+1):(tiempoTr[length(tiempoTr)]+NTest);
plot.ts(serieTr, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, serieTs, col="red")
cat("Dividimos el conjunto de datos en entrenamiento (para ajuste, negro) y test (para comprobar los modelos, rojo).\n");

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
# Eliminamos la tendencia
serieTr.SinTend.H1<- serieTr-TendEstimadaTr.H1;
serieTs.SinTend.H1<- serieTs-TendEstimadaTs.H1;
plot.ts(serieTr.SinTend.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, serieTs.SinTend.H1, col="red")
print("Serie sin la tendencia.")




 

# Calculamos y eliminamos la estacionalidad
k<- 12; # Asumimos periodo de estacionalidad k= 12
estacionalidad.H1<- decompose(serie.ts)$seasonal[1:k];

#Eliminamos estacionalidad para el modelo
aux<-rep(estacionalidad.H1, length(serieTr)/length(estacionalidad.H1));
serieTr.SinTendEst.H1<- serieTr.SinTend.H1-aux;
serieTs.SinTendEst.H1<- serieTs.SinTend.H1-estacionalidad.H1;
plot.ts(serieTr.SinTendEst.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs, serieTs.SinTendEst.H1, col="red")
print("Serie sin la estacionalidad\n")



# Vemos los ACF y PACF
acf(serieTr.SinTendEst.H1) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia lineal y sin estacionalidad\n")

pacf(serieTr.SinTendEst.H1) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia lineal y sin estacionalidad")



adftest.H1<- adf.test(serieTr.SinTendEst.H1);
cat(c("Resultados del test ADF (aproximación lineal sin tendencia ni estacionalidad): ", adftest.H1$p.value, "\n"));
cat("Hay que diferenciar.\n")

# Diferenciamos la serie
serieTr.SinTendEstDiff.H1<- diff(serieTr.SinTendEst.H1);
serieTs.SinTendEstDiff.H1<- diff(serieTs.SinTendEst.H1);

# Volvemos a aplicar el test
adftest.H1<- adf.test(serieTr.SinTendEstDiff.H1);
cat(c("Resultados del test ADF (aproximación lineal sin tendencia ni estacionalidad, y diferenciada): ", adftest.H1$p.value, ". Ya sí es estacionaria\n"));

# Volvemos a mostrar los ACF y PACF de la serie ya diferenciada
acf(serieTr.SinTendEstDiff.H1) # Mostramos ACF
print("Mostrando el ACF de la serie sin tendencia lineal, sin estacionalidad y diferenciada\n")

pacf(serieTr.SinTendEstDiff.H1) # Mostramos PACF
print("Mostrando el PACF de la serie sin tendencia lineal, sin estacionalidad y diferenciada")


# Mostramos las series sin tendencia ni estacionalidad
plot.ts(serieTr.SinTendEstDiff.H1, xlim=c(1, tiempoTs[length(tiempoTs)]))
lines(tiempoTs[1:(length(tiempoTs)-1)], serieTs.SinTendEstDiff.H1, col="red")
print("Serie sin la tendencia, ni la estacionalidad, y diferenciada (modelo lineal).")
print("Pulse una tecla para continuar...")
pause<-readline(); # para pausar la ejecución


cat("Con esto, podemos comenzar probando un modelo ARIMA(4, 1, 0):\n");
print("Pulse una tecla para continuar...")
pause<-readline(); # para pausar la ejecución
