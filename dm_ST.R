install.packages("fUnitRoots")
library(fUnitRoots)
library(dplyr)
library(lubridate)
library(zoo)
require(tseries)


serie <- read.csv("valeurs_mensuelles.csv",sep=";")
serie <- serie %>% select(-Codes)
serie <- serie[-c(1:3),]
colnames(serie) <- c("Date", "Valeur")

serie$Valeur <- as.numeric(serie$Valeur)
dates <- as.yearmon(serie$Date, "%Y-%m")
value_serie <- zoo(serie$Valeur, order.by = dates)

plot(value_serie)

#La série semble exhiber une tendance déterministe linéaire, en effet la pente est assez linéaire
trend <- 1:length(value_serie)
linear_trend <- lm(value_serie ~ trend))
resid <- linear_trend$residuals

pacf(resid)
acf(resid,lag=60)
#la série est extrêmement persistante, il y a de l'autocorrélation jusqu'à plus de 50 lags
#la PACF s'arrête à 5, donc q < 6  

## Dans le TD5, ils commentent la significativité de la régression : peut-être faire un commentaire
#"Le coefficient associ´e `a la tendance lin´eaire (dates) est bien n´egative, et peut-ˆetre significative (on ne peut pas
#vraiment le confirmer car le test n’est pas valide en pr´esence de r´esidus possiblement autocorr´el´es"

# VERIFICATION DE LA STATIONNARITE (test de racine unitaire), on sélectionne le test sans constante ni trend pusqu'il s'agit des résidus
# or les résidus d'une régression sont déjà centrés 

# Un test Dickey-Fuller simple avec lag = 0 ne fonctionnerait pas puisqu'on soupçonne un processus AR(p) avec p >1
plot(resid)
adfTest(resid, lag = 1, type = "nc") #j'ai pas compris il fallait mettre lag cb ici 

#
pp.test(resid)
kpss.test(resid, null = "Level")

# On peut faire un test de racine unitaire sur la série de base 
adf <-adfTest(value_serie, lag = 1, type = "ct")


# pas trop compris ce qu'il fallait faire
Qtests(adf@test$lm$residuals, 24, length(adf@test$lm$coefficients))


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
