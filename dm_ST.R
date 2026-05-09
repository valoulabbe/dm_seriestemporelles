### ===========================================================
### DM Séries temporelles — IPI Industries alimentaires (CVS-CJO)
### 
### ===========================================================

# --- Q1. IMPORTATION DES DONNÉES ---
require(zoo)
require(tseries)
require(fUnitRoots)

datafile <- "valeurs_mensuelles.csv"
data <- read.csv(datafile, sep = ";", skip = 3, header = TRUE)
colnames(data) <- c("Date", "IPI", "Codes")

# Nettoyage et tri
data <- data[!is.na(data$IPI), ]
data <- data[order(data$Date), ]
data$IPI <- as.numeric(data$IPI)

# Construction de l'objet zoo, troncature pré-COVID
dates <- as.yearmon(data$Date, "%Y-%m")
ipi_complet <- zoo(data$IPI, order.by = dates)
ipi <- window(ipi_complet, end = as.yearmon("2019-12"))


# --- Q2. REPRÉSENTATION GRAPHIQUE DE LA SÉRIE ET DE SA DIFFÉRENCE ---
d_ipi <- diff(ipi, 1)
plot(cbind(ipi, d_ipi), main = "IPI Industries alimentaires : niveau et différence première")
# La série en niveau présente une tendance croissante apparente.
# La différence première semble stationnaire autour d'une constante proche de 0.
# La série est probablement I(1).


# --- FONCTIONS ---
Qtests <- function(series, k, fitdf = 0) {
  pvals <- apply(matrix(1:k), 1, FUN = function(l) {
    pval <- if (l <= fitdf) NA else Box.test(series, lag = l, type = "Ljung-Box", fitdf = fitdf)$p.value
    return(c("lag" = l, "pval" = pval))
  })
  return(t(pvals))
}

adfTest_valid <- function(series, kmax, type) {
  k <- 0
  noautocorr <- 0
  while (noautocorr == 0) {
    cat(paste0("ADF with ", k, " lags: residuals OK? "))
    adf <- adfTest(series, lags = k, type = type)
    pvals <- Qtests(adf@test$lm$residuals, 24,
                    fitdf = length(adf@test$lm$coefficients))[, 2]
    if (sum(pvals < 0.05, na.rm = TRUE) == 0) {
      noautocorr <- 1; cat("OK \n")
    } else {
      cat("nope \n")
    }
    k <- k + 1
  }
  return(adf)
}

signif <- function(estim) {
  coef <- estim$coef
  se   <- sqrt(diag(estim$var.coef))
  t    <- coef / se
  pval <- (1 - pnorm(abs(t))) * 2
  return(rbind(coef, se, pval))
}


# --- Q3. TESTS DE RACINE UNITAIRE ---

# (a) Justification du type d'ADF sur la série EN NIVEAU
# Avant de choisir le type d'ADF, on vérifie via une régression linéaire
# si la série exhibe une constante et/ou une tendance significative.

summary(lm(coredata(ipi) ~ as.numeric(index(ipi))))



# Le coefficient associé aux dates est positif et significatif (0.297)
# => présence d'une tendance + d'une constante => type = "ct"

cat("\n--- ADF sur la série en niveau (type = ct) ---\n")
adf_niveau <- adfTest_valid(ipi, 24, type = "ct")
print(adf_niveau)
# Lag Order = 8 ; statistique = -2.6276 ; p-value = 0.312 > 5%
# => On NE rejette PAS H0 : la série en niveau a une racine unitaire.

# (b) Justification du type d'ADF sur la série DIFFÉRENCIÉE
summary(lm(coredata(d_ipi) ~ as.numeric(index(d_ipi))))
# Coefficients constante et tendance non significatifs (p-value 0.545)
# => type = "nc" (sans constante ni tendance)

cat("\n--- ADF sur la série différenciée (type = nc) ---\n")
adf_diff <- adfTest_valid(d_ipi, 24, type = "nc")
print(adf_diff)
# p-value = 0.01 < 0.05 => on rejette H0
# => La série différenciée est stationnaire.

# CONCLUSION Q3 : la série de l'IPI alimentaire est I(1).
# On modélisera ARIMA(p, 1, q) avec d* = 1.


# --- Q4. ACF / PACF DE LA SÉRIE STATIONNAIRE ---
x <- d_ipi
par(mfrow = c(1, 2))
acf(x,  main = "ACF de Δ IPI")
pacf(x, main = "PACF de Δ IPI")
par(mfrow = c(1, 1))
# ACF significative jusqu'à l'ordre 2  -> q* = 2
# PACF significative jusqu'à l'ordre 3 -> p* = 3
# Ordres maximaux : p* = 3, q* = 2.

# Fonction d'affichage type TD 4 : combine signif + Qtests
arimafit <- function(estim) {
  adjust <- round(signif(estim), 3)
  pvals  <- Qtests(estim$residuals, 24, fitdf = length(estim$coef))
  pvals  <- matrix(apply(matrix(1:24, nrow = 6), 2,
                         function(c) round(pvals[c, ], 3)), nrow = 6)
  colnames(pvals) <- rep(c("lag", "pval"), 4)
  cat("Tests de nullité des coefficients :\n")
  print(adjust)
  cat("\nTests d'absence d'autocorrélation des résidus :\n")
  print(pvals)
}


# --- Estimation et examen de chaque modèle ---

cat("\n========== ARIMA(1,1,0) ==========\n")
estim <- arima(ipi, c(1, 1, 0), include.mean = FALSE); arima110 <- estim
arimafit(estim)
# AR(1) significatif mais Ljung-Box rejeté à tous les lags => NON VALIDE.

cat("\n========== ARIMA(2,1,0) ==========\n")
estim <- arima(ipi, c(2, 1, 0), include.mean = FALSE); arima210 <- estim
arimafit(estim)
# AR(1) et AR(2) significatifs mais Ljung-Box rejeté => NON VALIDE.

cat("\n========== ARIMA(3,1,0) ==========\n")
estim <- arima(ipi, c(3, 1, 0), include.mean = FALSE); arima310 <- estim
arimafit(estim)
# Tous coefs AR significatifs mais Ljung-Box rejeté à plusieurs lags => NON VALIDE.

cat("\n========== ARIMA(0,1,1) ==========\n")
estim <- arima(ipi, c(0, 1, 1), include.mean = FALSE); arima011 <- estim
arimafit(estim)
# MA(1) significatif au seuil 1% (p<0.001), 
# Ljung-Box non rejeté aux lags canoniques (6, 12, 18, 24) ; rejets borderline aux lags 9-11.

cat("\n========== ARIMA(0,1,2) ==========\n")
estim <- arima(ipi, c(0, 1, 2), include.mean = FALSE); arima012 <- estim
arimafit(estim)
# MA(1) sig mais MA(2) p=0.104, non significatif au seuil 5% => MAL AJUSTÉ.

cat("\n========== ARIMA(1,1,1) ==========\n")
estim <- arima(ipi, c(1, 1, 1), include.mean = FALSE); arima111 <- estim
arimafit(estim)
# MA(1) sig mais AR(1) p=0.067, non significatif au seuil 5% => MAL AJUSTÉ.

cat("\n========== ARIMA(2,1,1) ==========\n")
estim <- arima(ipi, c(2, 1, 1), include.mean = FALSE); arima211 <- estim
arimafit(estim)
# AR(1) et MA(1) sig mais AR(2) p=0.426, non significatif => MAL AJUSTÉ.

cat("\n========== ARIMA(3,1,1) ==========\n")
estim <- arima(ipi, c(3, 1, 1), include.mean = FALSE); arima311 <- estim
arimafit(estim)
# MA(1) sig, mais AR(1), AR(2), AR(3) tous non significatifs => MAL AJUSTÉ.

cat("\n========== ARIMA(1,1,2) ==========\n")
estim <- arima(ipi, c(1, 1, 2), include.mean = FALSE); arima112 <- estim
arimafit(estim)
# MA(1) sig, mais AR(1) p=0.213 et MA(2) p=0.457 non significatifs => MAL AJUSTÉ.

cat("\n========== ARIMA(2,1,2) ==========\n")
estim <- arima(ipi, c(2, 1, 2), include.mean = FALSE); arima212 <- estim
arimafit(estim)
# AUCUN coefficient significatif => MAL AJUSTÉ.

cat("\n========== ARIMA(3,1,2) ==========\n")
estim <- arima(ipi, c(3, 1, 2), include.mean = FALSE); arima312 <- estim
arimafit(estim)
# AR(1), MA(1), MA(2) sig mais AR(2) p=0.085 et AR(3) p=0.108 non sig => MAL AJUSTÉ.


# --- Bilan final ---
cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("BILAN : un seul modèle est valide ET bien ajusté : ARIMA(0,1,1)\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")


# --- Confirmation par AIC/BIC sur tous les modèles ---
modeles <- list(
  "ARIMA(1,1,0)" = arima110, "ARIMA(2,1,0)" = arima210, "ARIMA(3,1,0)" = arima310,
  "ARIMA(0,1,1)" = arima011, "ARIMA(0,1,2)" = arima012,
  "ARIMA(1,1,1)" = arima111, "ARIMA(2,1,1)" = arima211, "ARIMA(3,1,1)" = arima311,
  "ARIMA(1,1,2)" = arima112, "ARIMA(2,1,2)" = arima212, "ARIMA(3,1,2)" = arima312
)
tab <- sapply(modeles, function(m) c(AIC = AIC(m), BIC = BIC(m)))
cat("\n--- Tableau récapitulatif AIC / BIC ---\n")
print(round(t(tab), 2))

# Le BIC est minimisé par ARIMA(0,1,1), confirmant la sélection.


# --- Modèle retenu ---
cat("\n========== MODÈLE RETENU : ARIMA(0,1,1) ==========\n")
print(arima011)
cat("\nSignificativité des coefficients :\n")
print(signif(arima011))
cat("\nLjung-Box sur les résidus (24 lags) :\n")
print(Qtests(arima011$residuals, 24, fitdf = 1))

# Équation du modèle :
#   (1 - L) X_t = (1 + theta_1 L) eps_t
#   <=> Δ X_t = eps_t + theta_1 * eps_{t-1}
#   avec theta_1 ≈ -0.61, significatif au seuil 1%.

# --- Chargement des librairies supplémentaires pour la Partie III ---
require(ggplot2)
require(ellipse)

# ================== PARTIE III : PRÉVISION ===================================

# --- Calcul des prévisions et de la matrice Sigma ---
# Récupération des paramètres du modèle ARIMA(0,1,1)

theta <- coef(arima011)["ma1"]
sigma2 <- arima011$sigma2

# Définition du paramètre psi1 
psi1 <- 1 + theta

# Dernier point observé (X_T) et dernier résidu estimé (eps_T)
X_T <- as.numeric(tail(ipi, 1))
eps_T <- as.numeric(tail(arima011$residuals, 1))

# Prévisions ponctuelles
X_T1_hat <- X_T + theta * eps_T
X_T2_hat <- X_T1_hat # Constant car ARIMA sans dérive

# --- Représentation graphique de l'ellipse  ---
# Matrice de variance-covariance 
Sigma <- sigma2 * matrix(c(1,    psi1,
                           psi1, 1 + psi1^2),
                         nrow = 2)

# Centre de l'ellipse (prévisions ponctuelles)
centre_df <- data.frame(x = X_T1_hat, y = X_T2_hat)

# Ellipse de confiance à 95%
ell    <- ellipse(Sigma, centre = c(X_T1_hat, X_T2_hat), level = 0.95)
ell_df <- as.data.frame(ell)

# Graphique
p_ellipse <- ggplot() +
  
  geom_polygon(
    data = ell_df,
    aes(x = x, y = y),
    fill = "grey70", alpha = 0.4, color = "black", linewidth = 0.5
  ) +
  
  geom_point(
    data = centre_df,
    aes(x = x, y = y),
    color = "red", fill = "red", shape = 21, size = 3
  ) +
  
  geom_vline(xintercept = X_T1_hat,
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  
  geom_hline(yintercept = X_T2_hat,
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  
  labs(
    title = "Région de confiance à 95% pour (X_{T+1}, X_{T+2})",
    x = expression(X[T+1]),
    y = expression(X[T+2])
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

print(p_ellipse)
