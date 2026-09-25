# 03-evaluacion.R
# ljung_box(), jarque_bera(), durbin_watson(): las tres pruebas sobre una serie o
# sobre los errores de un metodo.
# medidas(): MSE, MAD, MAPE, MASE.
# validar_errores(): empaqueta todo el diagnostico de los errores de un paso.

ljung_box <- function(r, T, m, p = 0) {
  stopifnot(is.numeric(r), length(r) >= m, m >= 1, is.numeric(T), T > m,
            is.numeric(p), p >= 0, m - p > 0)
  r_m <- r[seq_len(m)]
  h <- seq_len(m)
  Q <- T * (T + 2) * sum(r_m^2 / (T - h))
  gl <- m - p
  list(
    estadistico = Q,
    gl = gl,
    valor_critico = stats::qchisq(0.95, df = gl),
    valor_p = 1 - stats::pchisq(Q, df = gl)
  )
}

jarque_bera <- function(e) {
  stopifnot(is.numeric(e), !anyNA(e), length(e) >= 2)
  N <- length(e)
  m <- mean(e)
  s2 <- mean((e - m)^2)
  A <- mean((e - m)^3) / s2^1.5   # asimetria
  K <- mean((e - m)^4) / s2^2     # curtosis (normal = 3)
  JB <- (N / 6) * (A^2 + (K - 3)^2 / 4)
  list(
    estadistico = JB,
    gl = 2,
    valor_critico = stats::qchisq(0.95, df = 2),
    valor_p = 1 - stats::pchisq(JB, df = 2),
    asimetria = A,
    curtosis = K,
    N = N,
    aviso_muestra_pequena = N < 20
  )
}

# Tabla de valores criticos de Durbin-Watson al 5%, k = numero de regresores
# SIN el intercepto.
# Cubre k=1 y k=2.
.tabla_dw <- list(
  "1" = list(
    n  = c(15,16,17,18,19,20,21,22,23,24,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,200),
    dL = c(1.08,1.10,1.13,1.16,1.18,1.20,1.22,1.24,1.26,1.27,1.29,1.35,1.40,1.44,1.48,1.50,1.53,1.55,1.57,1.58,1.60,1.61,1.62,1.63,1.64,1.65,1.72,1.76),
    dU = c(1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.45,1.49,1.52,1.54,1.57,1.59,1.60,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.69,1.75,1.78)
  ),
  "2" = list(
    n  = c(15,16,17,18,19,20,21,22,23,24,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,200),
    dL = c(0.95,0.98,1.02,1.05,1.08,1.10,1.13,1.15,1.17,1.19,1.21,1.28,1.34,1.39,1.43,1.46,1.49,1.51,1.54,1.55,1.57,1.59,1.60,1.61,1.62,1.63,1.71,1.75),
    dU = c(1.54,1.54,1.54,1.53,1.53,1.54,1.54,1.54,1.54,1.55,1.55,1.57,1.58,1.60,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.70,1.71,1.72,1.76,1.79)
  )
)

durbin_watson <- function(e, k = 1) {
  stopifnot(is.numeric(e), !anyNA(e), length(e) >= 2, k %in% c(1, 2))
  Tn <- length(e)
  d <- sum(diff(e)^2) / sum(e^2)
  tabla <- .tabla_dw[[as.character(k)]]

  if (Tn <= max(tabla$n)) {
    # dentro del rango tabulado: cotas dL/dU reales
    dL <- stats::approx(tabla$n, tabla$dL, xout = Tn, rule = 2)$y
    dU <- stats::approx(tabla$n, tabla$dU, xout = Tn, rule = 2)$y
    decision <- if (d < dL) {
      "rechaza H0: autocorrelacion positiva"
    } else if (d > 4 - dL) {
      "rechaza H0: autocorrelacion negativa"
    } else if (d > dU && d < 4 - dU) {
      "no rechaza H0"
    } else {
      "zona inconclusa"
    }
    list(estadistico = d, metodo = "tabla", n = Tn, k = k,
         dL = dL, dU = dU, decision = decision)
  } else {
    # Tn excede el maximo tabulado (no existe tabla publicada mas alla de
    # n=200): se usa la aproximacion asintotica estandar d ~= 2(1-r1),
    # r1 ~ N(0,1/T) bajo H0 (este es el procedimiento documentado para
    # muestras grandes)
    r1 <- 1 - d / 2
    z <- r1 / (1 / sqrt(Tn))
    valor_p <- 2 * (1 - stats::pnorm(abs(z)))
    decision <- if (valor_p < 0.05) "rechaza H0: hay autocorrelacion" else "no rechaza H0"
    list(estadistico = d, metodo = "asintotico (T fuera de tabla)", n = Tn, k = k,
         r1_equivalente = r1, z = z, valor_critico_z = stats::qnorm(0.975),
         valor_p = valor_p, decision = decision)
  }
}

mae_ingenuo <- function(y, estacional = FALSE, s = 1) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) > s)
  if (estacional) {
    e <- y[(s + 1):length(y)] - y[1:(length(y) - s)]
  } else {
    e <- diff(y)
  }
  mean(abs(e))
}

medidas <- function(y, yhat, mae_ref) {
  stopifnot(length(y) == length(yhat), is.numeric(mae_ref), mae_ref > 0)
  e <- y - yhat
  MSE <- mean(e^2, na.rm = TRUE)
  MAD <- mean(abs(e), na.rm = TRUE)
  MAPE <- mean(abs(e / y), na.rm = TRUE) * 100
  MASE <- MAD / mae_ref
  list(MSE = MSE, MAD = MAD, MAPE = MAPE, MASE = MASE,
       n = sum(!is.na(e)))
}

validar_errores <- function(e, p = 0, m = NULL, k = 1) {
  e <- e[!is.na(e)]
  n <- length(e)
  if (is.null(m)) m <- min(floor(n / 4), 24)
  stopifnot(m - p > 0)

  grafico_errores <- ggplot2::ggplot(
    tibble::tibble(t = seq_len(n), e = e), ggplot2::aes(x = t, y = e)
  ) +
    ggplot2::geom_line(color = "#1f3b57") +
    ggplot2::geom_point(size = 0.8, color = "#1f3b57") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
    ggplot2::labs(title = "Errores de un paso", x = "t", y = "e_t") +
    ggplot2::theme_minimal(base_size = 12)

  corr <- correlograma(e, m = m)

  media_e <- mean(e)
  s_e <- stats::sd(e)
  t_stat <- media_e / (s_e / sqrt(n))
  t_gl <- n - 1

  list(
    n = n, m = m, p = p,
    grafico_errores = grafico_errores,
    correlograma = corr,
    media_cero = list(
      estadistico = t_stat,
      gl = t_gl,
      valor_critico = stats::qt(0.975, df = t_gl),
      valor_p = 2 * (1 - stats::pt(abs(t_stat), df = t_gl)),
      media_observada = media_e
    ),
    ljung_box = ljung_box(corr$acf, T = n, m = m, p = p),
    jarque_bera = jarque_bera(e),
    durbin_watson = durbin_watson(e, k = k)
  )
}
