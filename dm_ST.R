
library(zoo)
library(tseries)
library(fUnitRoots)

# IMPORTATION ET NETTOYAGE
serie <- read.csv("valeurs_mensuelles.csv", sep=";", skip=3, header=TRUE)
colnames(serie) <- c("Date", "IPI", "Codes")

# Nettoyage et tri
serie <- serie[!is.na(serie$IPI), ]
serie <- serie[order(serie$Date), ]
serie$IPI <- as.numeric(serie$IPI)

# Indexation temporelle et création de l'objet zoo
dates <- as.yearmon(serie$Date, "%Y-%m")
ipi_complet <- zoo(serie$IPI, order.by = dates)

# Troncature avant le COVID-19 (fin 2019)
ipi <- window(ipi_complet, end = as.yearmon("2019-12"))

plot(ipi, main="IPI - Industrie Pharmaceutique (Pré-COVID)", ylab="Indice", xlab="Temps")

# FONCTIONS DE VALIDATION DU TD5 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}

# --- TESTS DE STATIONNARITÉ ---
# La série a une tendance visible, on teste avec constante et tendance (type="ct")
cat("\n--- Test ADF sur la série en niveau ---\n")
adf_niveau <- adfTest_valid(ipi, 24, "ct")
print(adf_niveau)


# La boucle a automatiquement sélectionné 11 retards (Lag Order: 11) pour 
# garantir l'absence d'autocorrélation des résidus (validé par Ljung-Box).
#
# Statistique Dickey-Fuller : -1.8054
# P-VALUE : 0.659
#
# INTERPRÉTATION :
# La p-value (0.659) est très largement supérieure au seuil de 5% (0.05).
# Dans le cadre du test ADF, l'hypothèse nulle (H0) postule la présence 
# d'une racine unitaire.
# -> Nous ne pouvons donc PAS rejeter H0.
# -> CONCLUSION : La série de l'IPI en niveau est NON-STATIONNAIRE.

# Différenciation première pour stationnariser
d_ipi <- diff(ipi, 1)
plot(d_ipi, main="IPI Pharmaceutique - Différence première", ylab="Δ IPI")

# Test ADF sur la série différenciée (sans tendance, type="c")
cat("\n--- Test ADF sur la série différenciée ---\n")
adf_diff <- adfTest_valid(d_ipi, 24, "c")
print(adf_diff)

# INTERPRÉTATION :
# La p-value (0.01) est strictement inférieure au seuil de 5% (0.05).
# -> Nous REJETONS l'hypothèse nulle (H0) de présence d'une racine unitaire.
# -> CONCLUSION : La série différenciée en différence première est STATIONNAIRE.
#
# BILAN DE L'ÉTUDE DE STATIONNARITÉ (PARTIE I) :
# La série en niveau étant non-stationnaire et la série en différence 
# première étant stationnaire, on conclut que la série de l'IPI 
# de l'industrie pharmaceutique est intégrée d'ordre 1, noté I(1).
#