# Projekt ma na celu zbadanie wspó³zale¿noœci i podobieñstw zmiennoœci 
# indeksów gie³dowych grupy BRICS, która za cel swojego powstwania 
# uzna³a detronizacjê dolara jako g³ównej waluty transakcyjnej na rynkach œwiatowych

library(quantmod)
library(rugarch)
library(e1071)
library(tseries)
library(lmtest)
library(FinTS)
library(zoo)
library(TTR)
library(xts)
library(vars)
# A quantmod to bardzo potezna biblioteka

getSymbols("^BVSP", src = "yahoo", auto.assign = TRUE) # Brazylia
getSymbols("IMOEX.ME", src = "yahoo", auto.assign = TRUE) # Rosja
getSymbols("^BSESN", src = "yahoo", auto.assign = TRUE) # Indie
getSymbols("399001.SZ", src = "yahoo", auto.assign = TRUE) # Chiny
getSymbols("^JN0U.JO", src = "yahoo", auto.assign = TRUE) # Republika Po³udniowej Afryki
Szangh=`399001.SZ`
BVSP_clean <- na.omit(BVSP)
IMOEX.ME_clean <- na.omit(IMOEX.ME)
BSESN_clean <- na.omit(BSESN)
SSEC_clean <- na.omit(Szangh)
JN0U.JO_clean <- na.omit(JN0U.JO)

# Przesuñ datê rozpoczêcia dla indeksu RPA
start_date <- as.Date("2017-09-26")

# Wycinamy dane tylko od tej daty
BVSP_clean <- window(BVSP_clean, start = start_date)
IMOEX.ME_clean <- window(IMOEX.ME_clean, start = start_date)
BSESN_clean <- window(BSESN_clean, start = start_date)
SSEC_clean <- window(Szangh, start = start_date)
JN0U.JO_clean <- window(JN0U.JO_clean, start = start_date)

kraje=c("Brazylia","Rosja","Indie","Chiny","RPA")

# £¹czymy ramki danych w oparciu o wspóln¹ datê
ceny <- merge(BVSP_clean$BVSP.High, IMOEX.ME_clean$IMOEX.ME.High, BSESN_clean$BSESN.High,
              SSEC_clean$`399001.SZ.High`, JN0U.JO_clean$JN0U.JO.High, all = FALSE)
colnames(ceny) <- kraje
layout_matrix <- matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE)
layout(layout_matrix)

# Wykresy dla samych cen
for (i in 1:ncol(ceny)) {
  print(plot(ceny[, i], main = paste0("Cena g³ównego indeksu gie³dowego: ", colnames(ceny)[i]),col=i))
}

##################################### Pora na logarytmy#######################
log_ret <- function(x) {
  return(diff(log(x)))
}

log_returns <- apply(ceny, 2, log_ret)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
for (i in 1:ncol(log_returns)) {
  plot(log_returns[, i], type = "l", main = paste0("Stopy zwrotu logarytmiczne dla: ", colnames(log_returns)[i]), ylab = "Logarytmiczne stopy zwrotu", xlab = "Czas",col=i)
}



# Obliczanie skoœnoœci i kurtozy dla ka¿dego kraju
skosnosc <- apply(log_returns, 2, skewness)
kurtosis <- apply(log_returns, 2, kurtosis)
wyniki_rozklady <- data.frame(Kraj = colnames(log_returns), Skosnosc = skosnosc, Kurtosis = kurtosis)
print(wyniki_rozklady)

# Wykresy QQ-plot dla ka¿dego kraju
par(mfrow = c(3, 2), mar = c(1, 1, 1, 1))
for (i in 1:ncol(log_returns)) {
  qqnorm(log_returns[, i], main = paste0("Wykres QQ dla: ", colnames(log_returns)[i]))
  qqline(log_returns[, i], col = "red")
}
#### Wszystkie kraje charakteryzuj¹ siê asymetri¹ lewostronn¹, ponadto 
# wszystkie szeregi czasowe s¹ leptokurtyczne (grube ogony)

X11(width = 8, height = 12) 
par(mfrow = c(5, 2), mar = c(5, 5, 3, 1))
for (i in 1:ncol(log_returns)) {
  acf(log_returns[, i], main = paste0("ACF dla: ", colnames(log_returns)[i]), lag.max = 20)
  pacf(log_returns[, i], main = paste0("PACF dla: ", colnames(log_returns)[i]), lag.max = 20)
}

log_returns_squared <- log_returns^2
X11(width = 8, height = 12) 
par(mfrow = c(5, 2), mar = c(5, 5, 3, 1))
for (i in 1:ncol(log_returns_squared)) {
  pacf(log_returns_squared[, i], main = paste0("PACF dla kwadratów log. stóp zwrotu: ", colnames(log_returns_squared)[i]), lag.max = 20)
}
### WyraŸny efekt ARCH dla ka¿dego szeregu 

#### Sprawdzmy efekt ARCH, stacjonarnoœc i autokorelacje

autocorrelation_tests <- numeric(ncol(log_returns))
stationarity_tests <- numeric(ncol(log_returns))
arch_tests <- numeric(ncol(log_returns))

for (i in 1:ncol(log_returns)) {
  autocorrelation_tests[i] <- Box.test(log_returns[, i], lag = 20, type = "Ljung-Box")$p.value
  stationarity_tests[i] <- kpss.test(log_returns[, i], null = "Level")$p.value
  arch_tests[i] <- ArchTest(log_returns[, i], lags = 12)$p.value
}

results.assump <- data.frame(Autocorrelation = autocorrelation_tests,
                      Stationarity = stationarity_tests,
                      ARCH = arch_tests,
                      row.names = colnames(log_returns))
results.assump
#### wyniki sugeruj¹, ¿e dla wszystkich krajów wystêpuje silna autokorelacja
#(ma³e wartoœci p w teœcie na autokorelacjê) oraz efekt ARCH (ma³e wartoœci p
#w teœcie ARCH). To oznacza, ¿e zmiennoœæ szeregów czasowych mo¿e byæ zwi¹zana
#z przesz³ymi wartoœciami, a efekt ARCH wskazuje na istnienie zjawiska tzw.
#skupienia siê zmiennoœci.


### Wyniki testu KPSS sugeruj¹ równie¿, ¿e dla wszystkich krajów szeregi czasowe
#s¹ stacjonarne (du¿e wartoœci p w teœcie KPSS), co oznacza, ¿e maj¹ sta³¹
#wartoœæ oczekiwan¹ i wariancjê w czasie.


##### Zobaczmy ruchome korelacje dla wszystkich kombinacji indeksów
# Ustalamy okno czasowe na oko³o miesi¹c
window_size <- 31

# Obliczamy ruchome korelacje pomiêdzy wszystkimi kombinacjami krajów
rolling_correlations <- list()
for (i in 1:(ncol(log_returns) - 1)) {
  for (j in (i + 1):ncol(log_returns)) {
    rolling_correlations[[paste0(colnames(log_returns)[i], "_", colnames(log_returns)[j])]] <- 
      rollapply(log_returns[, i], window_size, function(x) cor(x, log_returns[(which(index(log_returns) %in% index(x)))[1]:((which(index(log_returns) %in% index(x)))[1] + window_size - 1), j]), by.column = FALSE, align = "right")
  }
}

if (capabilities("X11")) {
  X11()
} else {
  dev.new()
}
par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))
for (i in 1:length(rolling_correlations)) {
  plot(index(rolling_correlations[[i]]), coredata(rolling_correlations[[i]]), type = "l", main = names(rolling_correlations)[i], xlab = "Data", ylab = "Korelacja", ylim = c(-1, 1))
  abline(h = 0, col = "red", lty = 2)
}
#### Widzimy ¿e dla ka¿dej kombinacji korelacje ruchome s¹ w pewnych okresach 
# wiêksze a w pewnych mniejsze lub praktyczne zerowe


### Pora na VAR-MGARCH

# zmienne
y <- log_returns
n <- nrow(y)
k <- ncol(y)

# Wyznaczanie optymalnej liczby opóŸnieñ w modelu VAR
var_select <- VARselect(y, lag.max = 12, type = "const", season = NULL, exogen = NULL)

p <- var_sel$selection[1] # liczba opóŸnieñ

# Model VAR
var_fit <- VAR(y, p = p, type = "const", lag.max = NULL, ic = "AIC", season = NULL)

# Reszty modelu VAR
resid_var <- serial.test(var_fit)
resid_var
roots(var_fit)
X11()
par(mfrow = c(5, 2), mar = c(5, 5, 2, 1))
pacf(residuals(var_fit))

# Autokorelacja niestety dalej wystepuje w modelu VAR(31), model ma bardzo
# d³ug¹ pamiêæ, mo¿na pomyœleæ o modelach VARFIMA


# Model MGARCH, ustawi³em rozk³ad jako skoœny t student,
# jako ¿e skoœnoœæ dla ka¿dego szeregu by³a ujemna, m¹drzejsza coœ takiego zrobiæ

spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                   distribution.model = "std", 
                   start.pars = list(mu = c(0, 0, 0, 0, 0), omega = c(0.1, 0.1, 0.1, 0.1, 0.1), alpha = c(0.1, 0.1, 0.1, 0.1, 0.1), beta = c(0.8, 0.8, 0.8, 0.8, 0.8), shape = -1))


fit <- ugarchfit(spec = spec, data = residuals(var_fit), solver = "hybrid")
show(fit)
X11()
par(mfrow = c(5, 2), mar = c(5, 5, 2, 1))
plot(fit)

# jest lepiej ale problem autokorelacji nie znikn¹³, model VARFIMA jako 
# model œredniej warunkowej mo¿e zrobiæ robotê

