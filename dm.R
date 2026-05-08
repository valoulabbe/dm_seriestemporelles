library(zoo)
library(tseries)
library(fUnitRoots)

# --- 1. IMPORTATION ET NETTOYAGE ---
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

# Nouveau titre pour le graphique
plot(ipi, main="IPI - Production de l'industrie alimentaire (Pré-COVID)", ylab="Indice", xlab="Temps")

# --- 2. FONCTIONS DE VALIDATION ------

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
  # AJOUT : on force la boucle à s'arrêter si k dépasse kmax
  while (noautocorr==0 & k <= kmax){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    
    # On met un tryCatch pour éviter que R ne plante si le lag est trop grand
    adf <- try(adfTest(series, lags=k, type=adftype), silent=TRUE)
    
    if (inherits(adf, "try-error")) {
      cat("erreur de calcul\n")
      k <- k + 1
      next
    }
    
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1
      cat("OK \n")
    } else {
      cat("nope \n")
      k <- k+1
    }
  }
  
  if (noautocorr == 0) {
    cat("\nATTENTION : Aucun lag <= kmax n'a permis de blanchir totalement les résidus.\nLe test renvoyé est celui avec kmax retards.\n")
  }
  
  return(adf)
}

# --- 3. TESTS DE STATIONNARITÉ ---
# On teste avec constante et tendance (type="ct") pour la série en niveau
cat("\n--- Test ADF sur la série en niveau ---\n")
adf_niveau <- adfTest_valid(ipi, 24, "ct")
print(adf_niveau)

# -------------------------------------------------------------------------
# INTERPRÉTATION DU TEST ADF SUR LA SÉRIE EN NIVEAU :
# La boucle a automatiquement sélectionné 8 retards (Lag Order: 8) pour 
# garantir l'absence d'autocorrélation des résidus (validé par Ljung-Box).
#
# Statistique Dickey-Fuller : -2.6276
# P-VALUE : 0.312
#
# INTERPRÉTATION :
# La p-value (0.312) est supérieure au seuil de 5% (0.05).
# Dans le cadre du test ADF, l'hypothèse nulle (H0) postule la présence 
# d'une racine unitaire.
# -> Nous ne pouvons donc PAS rejeter H0.
# -> CONCLUSION : La série de l'IPI en niveau est NON-STATIONNAIRE.
# -------------------------------------------------------------------------

# Différenciation première pour stationnariser
d_ipi <- diff(ipi, 1)

# Nouveau titre pour la série différenciée
plot(d_ipi, main="IPI Production industrie alimentaire - Différence première", ylab="Δ IPI")

# Test ADF sur la série différenciée (sans tendance, type="c")
cat("\n--- Test ADF sur la série différenciée ---\n")
adf_diff <- adfTest_valid(d_ipi, 24, "c")
print(adf_diff)

# -------------------------------------------------------------------------
# INTERPRÉTATION DU TEST ADF SUR LA SÉRIE DIFFÉRENCIÉE :
# La boucle a automatiquement sélectionné 7 retards (Lag Order: 10).
# Statistique Dickey-Fuller : -10.1529
# P-VALUE : < 0.01
#
# La p-value est strictement inférieure au seuil de 5% (0.05).
# -> Nous REJETONS l'hypothèse nulle (H0) de présence d'une racine unitaire.
# -> CONCLUSION : La série en différence première est STATIONNAIRE.
#
# BILAN DE L'ÉTUDE DE STATIONNARITÉ (PARTIE I) :
# La série en niveau étant non-stationnaire et la série en différence 
# première étant stationnaire, on conclut que la série d'Extraction 
# est intégrée d'ordre 1, noté I(1).
# -------------------------------------------------------------------------

# --- 4. ETUDE ACF ET PACF en vue du CHOIX DE l'ARMA adapté ---
# Etude des fonctions d'autocorrélations et d'autocorrélations partielles
acf(d_ipi)
pacf(d_ipi)

# Les fonctions d'autocorrélations sont signicatives jusqu'à q*=2
# Les fonctions d'autocorrélations partielles elles, jusqu'à p*=3 
# Si la série suit un processus ARIMA, elle suit "au plus",
# un processus ARIMA(p*=3,d=1,q*=2)


#On vérifie la validité d'un ARIMA(3,1,2) en regardant l'autocorrélation des résidus 
arima312 <- arima(ipi,c(3,1,2)) 
Box.test(arima312$residuals, lag=6, type="Ljung-Box", fitdf=5) 


#L'hypothèse nulle n'est pas rejetée au seuil de 5%, on peut supposer 
# que les résidus jusqu'à l'horizon 6 ne sont pas autocorrélés
# On peut vérifier pour des horizons plus long en traçant les ACF/PACF des résidus 
pacf(arima312$residuals)
acf(arima312$residuals)

#On observe de l'autocorrélations partielles faible (-0.15) 
#pour un horizon 24, ce qui est sûrement un résidu des 
# traitements corrigeants la saisonnalité, 

#On peut vérifier cela en regardant la nullité jointe des coefficients d'autocorrélations des résidus 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima312$residuals, 24, 5)

#Si on teste à l'horizon 7, la nullité conjointe des coefficients d'autocorrélation est rejetée
#on peut cependant tout de même considérer que le modèle demeure valide 


# --- 5. SIMPLIFICATION DU MODELE ARMA ---

#Vérifions si les coefficients de l'ARIMA(312) sont significatifs afin
#d'éventuellement prendre un modèle plus simple

signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t))) * 2
  return(rbind(coef,se,t,pval))
}
signif(arima312) #tests de significativité  des coefficients de l’ARIMA(3,1,2)

# INTERPRÉTATION DU TEST de nullité des COEFFICIENTS de l'ARIMA(3,1,2) :
# Les coefficients de ar1, ma1 et ma2 sont significatifs au seuil de 1%
# P-VALUE : < 0.01

# En revanche les coefficients de ar2 et ar3 ne le sont pas au seuil de 5%
# On ne peut donc pas rejeter l'hypothèse (H0) que ces coefficients soient nuls 



# Testons l'ajustement et la validité de tous les sous modèles possibles

test_model <- function(p, q, series = ipi) {
  model_name <- paste0("arima", p, "1", q)
  estim <- arima(series, c(p, 1, q))
  assign(model_name, estim, envir = .GlobalEnv)
  fitdf <- p + q
  lb_lag <- fitdf + 1
  print(Box.test(estim$residuals, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf))
  print(Qtests(estim$residuals, 24, fitdf))
  print(signif(estim))
  invisible(estim)
}

test_model(p=1 ,q=0, series=ipi)
test_model(p=2 ,q=0, series=ipi)
test_model(p=3 ,q=0, series=ipi)
# Tous ces modèles s'ils sont très bien ajustés ne sont pas valides, on les exclut


test_model(p=1 ,q=1, series=ipi)
test_model(p=2 ,q=1, series=ipi)
test_model(p=3 ,q=1, series=ipi)
#Tous ces modèles sont valides, mais seul ARIMA(111) est relativement bien ajusté, 
# même si son coefficient AR(1) n'est pas significatif 
# au seuil de 5% (p-value : 0.067), donc on ne peut pas rejeter qu'il soit nul


test_model(p=1 ,q=2, series=ipi)
test_model(p=2 ,q=2, series=ipi)
#Ces modèles sont valides mais très mal ajustés, on les écarte

#CONCLUSION : Les deux modèles à peu près valides et ajustés sont ARIMA312 et ARIMA101
AIC(arima312)
BIC(arima312)

AIC(arima111)
BIC(arima111)


