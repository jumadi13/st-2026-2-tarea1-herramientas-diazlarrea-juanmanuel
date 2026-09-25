suppressPackageStartupMessages({
  library(tibble); library(ggplot2); library(patchwork); library(dplyr)
})

source("R/00-lectura.R")
source("R/01-graficos.R")
source("R/02-metodos.R")
source("R/03-evaluacion.R")
source("ejemplos/driver.R")

dir.create("figs", showWarnings = FALSE)
resultados <- list()

# 1) Media simple 
resultados$media <- ejecutar_ejemplo(
  nombre = "1. Media simple - discoveries",
  serie_ts = discoveries,
  fuente = "Number of \"great\" inventions and scientific discoveries, 1860-1959",
  unidad = "conteo de descubrimientos",
  construir_ajuste = function(y_est) {
    list(ajuste = ajustar_media(y_est), parametros_elegidos = NULL, opt = NULL)
  },
  p_metodo = 1, k_dw = 1, prefijo_fig = "01_media"
)

# 2) Media movil 
resultados$mm <- ejecutar_ejemplo(
  nombre = "2. Media movil - Nile",
  serie_ts = Nile,
  fuente = "Measurements of the annual flow of the river Nile at Aswan, 1871-1970",
  unidad = "10^8 m^3",
  construir_ajuste = function(y_est) {
    opt <- optimizar(y_est, "mm", rejilla = 2:12)
    list(ajuste = ajustar_mm(y_est, opt$optimo$k),
         parametros_elegidos = list(k = opt$optimo$k), opt = opt)
  },
  p_metodo = 0, k_dw = 1, prefijo_fig = "02_mm"
)

# 3) SES 
resultados$ses <- ejecutar_ejemplo(
  nombre = "3. SES - discoveries",
  serie_ts = discoveries,
  fuente = "Number of \"great\" inventions and scientific discoveries, 1860-1959",
  unidad = "conteo de descubrimientos",
  construir_ajuste = function(y_est) {
    opt <- optimizar(y_est, "ses", rejilla = seq(0.02, 0.98, by = 0.02))
    list(ajuste = ajustar_ses(y_est, opt$optimo$alpha),
         parametros_elegidos = list(alpha = opt$optimo$alpha), opt = opt)
  },
  p_metodo = 1, k_dw = 1, prefijo_fig = "03_ses"
)

# 3-bis) Contraejemplo: SES sobre nottem  
resultados$ses_contraejemplo <- ejecutar_ejemplo(
  nombre = "3b. CONTRAEJEMPLO: SES - nottem (estacional)",
  serie_ts = nottem,
  fuente = "Average monthly temperatures at Nottingham Castle, 1920-1939",
  unidad = "grados Fahrenheit",
  construir_ajuste = function(y_est) {
    opt <- optimizar(y_est, "ses", rejilla = seq(0.02, 0.98, by = 0.02))
    list(ajuste = ajustar_ses(y_est, opt$optimo$alpha),
         parametros_elegidos = list(alpha = opt$optimo$alpha), opt = opt)
  },
  p_metodo = 1, k_dw = 1, estacional = TRUE, s = 12, prefijo_fig = "03b_ses_contraejemplo"
)

# 4) Doble media movil 
resultados$dmm <- ejecutar_ejemplo(
  nombre = "4. Doble media movil - austres",
  serie_ts = austres,
  fuente = "Quarterly Time Series of the Number of Australian Residents",
  unidad = "miles de personas",
  construir_ajuste = function(y_est) {
    opt <- optimizar(y_est, "dmm", rejilla = 2:12)
    list(ajuste = ajustar_dmm(y_est, opt$optimo$k),
         parametros_elegidos = list(k = opt$optimo$k), opt = opt)
  },
  p_metodo = 0, k_dw = 1, prefijo_fig = "04_dmm"
)

# 5) Tendencia lineal 
resultados$tend_lineal <- ejecutar_ejemplo(
  nombre = "5. Tendencia lineal - LakeHuron",
  serie_ts = LakeHuron,
  fuente = "Level of Lake Huron 1875-1972",
  unidad = "pies (nivel del lago)",
  construir_ajuste = function(y_est) {
    a <- ajustar_tendencia(y_est, "lineal")
    nombres_beta <- a$parametros$tabla_coeficientes$parametro
    list(ajuste = a, parametros_elegidos = as.list(setNames(a$parametros$beta, nombres_beta)), opt = NULL)
  },
  p_metodo = 2, k_dw = 1, prefijo_fig = "05_tend_lineal"
)

# 6) Tendencia cuadratica 
resultados$tend_cuadratica <- ejecutar_ejemplo(
  nombre = "6. Tendencia cuadratica - co2",
  serie_ts = co2,
  fuente = "Mauna Loa Atmospheric CO2 Concentration, 1959-1997",
  unidad = "ppm",
  construir_ajuste = function(y_est) {
    a <- ajustar_tendencia(y_est, "cuadratica")
    nombres_beta <- a$parametros$tabla_coeficientes$parametro
    list(ajuste = a, parametros_elegidos = as.list(setNames(a$parametros$beta, nombres_beta)), opt = NULL)
  },
  p_metodo = 3, k_dw = 2, prefijo_fig = "06_tend_cuadratica"
)

# 7) Tendencia exponencial 
resultados$tend_exponencial <- ejecutar_ejemplo(
  nombre = "7. Tendencia exponencial - airmiles",
  serie_ts = airmiles,
  fuente = "US Airline Passenger Miles, 1937-1960",
  unidad = "millones de millas",
  construir_ajuste = function(y_est) {
    a <- ajustar_tendencia(y_est, "exponencial", corregir_sesgo = TRUE)
    nombres_beta <- a$parametros$tabla_coeficientes$parametro
    list(ajuste = a, parametros_elegidos = as.list(setNames(a$parametros$beta, nombres_beta)), opt = NULL)
  },
  p_metodo = 2, k_dw = 1, prefijo_fig = "07_tend_exponencial"
)

# 8) Holt 
resultados$holt <- ejecutar_ejemplo(
  nombre = "8. Holt lineal - WWWusage",
  serie_ts = WWWusage,
  fuente = "Internet Usage per Minute (server connections)",
  unidad = "usuarios conectados",
  construir_ajuste = function(y_est) {
    opt <- optimizar(y_est, "holt", rejilla = list(alpha = seq(0.05, 0.95, by = 0.05),
                                                     beta  = seq(0.05, 0.95, by = 0.05)))
    list(ajuste = ajustar_holt(y_est, opt$optimo$alpha, opt$optimo$beta),
         parametros_elegidos = list(alpha = opt$optimo$alpha, beta = opt$optimo$beta), opt = opt)
  },
  p_metodo = 2, k_dw = 1, prefijo_fig = "08_holt"
)

# resumen 
resumen <- purrr::map_dfr(resultados, function(r) {
  pe <- unlist(r$parametros_elegidos)
  texto_parametros <- if (is.null(pe)) {
    "(sin constante que ajustar)"
  } else {
    paste(names(pe), signif(pe, 4), sep = "=", collapse = ", ")
  }
  tibble(
    ejemplo = r$nombre, n = r$Tn, h = r$h,
    parametros = texto_parametros,
    MASE_metodo = round(r$med_validacion$MASE, 3),
    MASE_ingenuo = round(r$med_ingenuo_val$MASE, 3),
    gana_al_ingenuo = r$med_validacion$MASE < r$med_ingenuo_val$MASE
  )
})

cat("\n================= RESUMEN =================\n")
print(as.data.frame(resumen), row.names = FALSE)

saveRDS(resultados, "figs/resultados_completos.rds")
cat("\nfigs/ contiene", length(list.files("figs", pattern = "png$")), "imagenes.\n")
cat("resultados completos guardados en figs/resultados_completos.rds (para informe.qmd)\n")
