# 02-metodos.R: metodos de nivel
# ajustar_media(), ajustar_mm(), ajustar_ses(), ajustar_dmm()
# cada uno devuelve list(yhat, pronosticar, parametros).

ajustar_media <- function(y) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) >= 2)
  Tn <- length(y)
  yhat <- rep(NA_real_, Tn)

  suma <- y[1]
  n_acum <- 1
  for (t in 2:Tn) {
    yhat[t] <- suma / n_acum          # Ŷ_t = media recursiva de y[1..t-1]
    suma <- suma + y[t]
    n_acum <- n_acum + 1
  }
  ybar_T <- suma / n_acum             # media de toda la muestra

  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    rep(ybar_T, h)
  }

  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(beta0_hat = ybar_T))
}

ajustar_mm <- function(y, k) {
  stopifnot(is.numeric(y), !anyNA(y),
            length(k) == 1, k == round(k), k >= 2, k <= length(y))
  Tn <- length(y)
  yhat <- rep(NA_real_, Tn)

  mm <- sum(y[1:k]) / k                       # MM_k(k)
  if (k + 1 <= Tn) yhat[k + 1] <- mm

  if (Tn > k) {
    for (t in (k + 1):Tn) {
      mm <- mm + (y[t] - y[t - k]) / k        # MM_t(k), un solo recorrido
      if (t + 1 <= Tn) yhat[t + 1] <- mm
    }
  }
  # al salir del bucle, mm = MM_T(k) (usa TODA la muestra)

  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    rep(mm, h)
  }

  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(k = k, MM_T = mm))
}

ajustar_ses <- function(y, alpha) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) >= 2,
            length(alpha) == 1, is.numeric(alpha), alpha > 0, alpha < 1)
  Tn <- length(y)
  yhat <- rep(NA_real_, Tn)
  yhat[2] <- y[1]                              # calentamiento: Ŷ2 = Y1

  if (Tn >= 3) {
    for (t in 2:(Tn - 1)) {
      e_t <- y[t] - yhat[t]
      yhat[t + 1] <- yhat[t] + alpha * e_t     # forma de correccion de error
    }
  }
  e_T <- y[Tn] - yhat[Tn]
  nivel_final <- yhat[Tn] + alpha * e_T        # = alpha*Y_T + (1-alpha)*Yhat_T

  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    rep(nivel_final, h)
  }

  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(alpha = alpha, nivel_final = nivel_final))
}

ajustar_dmm <- function(y, k) {
  stopifnot(is.numeric(y), !anyNA(y),
            length(k) == 1, k == round(k), k >= 2, 2 * k - 1 <= length(y))
  Tn <- length(y)

  # trayectoria completa de MM_t(k), t = k, ..., T
  mm <- rep(NA_real_, Tn)
  s <- sum(y[1:k])
  mm[k] <- s / k
  if (Tn > k) {
    for (t in (k + 1):Tn) {
      s <- s + y[t] - y[t - k]
      mm[t] <- s / k
    }
  }

  # trayectoria completa de DMM_t(k) = MM de la serie mm, mismo k, t = 2k-1..T
  dmm <- rep(NA_real_, Tn)
  s2 <- sum(mm[k:(2 * k - 1)])
  dmm[2 * k - 1] <- s2 / k
  if (Tn > 2 * k - 1) {
    for (t in (2 * k):Tn) {
      s2 <- s2 + mm[t] - mm[t - k]
      dmm[t] <- s2 / k
    }
  }

  Et <- 2 * mm - dmm
  beta1t <- (2 / (k - 1)) * (mm - dmm)

  yhat <- rep(NA_real_, Tn)
  idx <- (2 * k - 1):(Tn - 1)
  if (length(idx) > 0 && idx[1] <= idx[length(idx)]) yhat[idx + 1] <- Et[idx] + beta1t[idx]

  E_T <- Et[Tn]
  beta1_T <- beta1t[Tn]

  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    E_T + beta1_T * seq_len(h)
  }

  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(k = k, E_T = E_T, beta1_T = beta1_T,
                          E_t = Et, beta1_t = beta1t, MM_t = mm, DMM_t = dmm))
}
# 02-metodos.R: tendencias por minimos cuadrados
# ajustar_tendencia(y, tipo, corregir_sesgo)

ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                               corregir_sesgo = FALSE) {
  tipo <- match.arg(tipo)
  stopifnot(is.numeric(y), !anyNA(y))
  Tn <- length(y)
  t <- seq_len(Tn)

  if (tipo == "exponencial") {
    if (any(y <= 0)) {
      stop("ajustar_tendencia('exponencial') requiere y > 0 en todas las observaciones (se estima sobre log(y))")
    }
    z <- log(y)
    X <- cbind(1, t)
    nombres <- c("a", "theta")
  } else if (tipo == "lineal") {
    z <- y
    X <- cbind(1, t)
    nombres <- c("beta0", "beta1")
  } else {
    z <- y
    X <- cbind(1, t, t^2)
    nombres <- c("beta0", "beta1", "beta2")
  }
  p <- ncol(X)
  stopifnot(Tn > p + 1)  # deja grados de libertad para sigma2 y HAC

  # ecuaciones normales
  XtX <- crossprod(X)
  beta <- as.numeric(solve(XtX, crossprod(X, z)))
  ajustado_z <- as.numeric(X %*% beta)
  residuos <- z - ajustado_z              # residuos en la escala de z (log si es exponencial)
  gl <- Tn - p
  sigma2 <- sum(residuos^2) / gl

  XtX_inv <- solve(XtX)
  se_ord <- sqrt(diag(XtX_inv) * sigma2)

  # errores estandar HAC (Newey-West, nucleo de Bartlett)
  L <- floor(4 * (Tn / 100)^(2 / 9))
  XE <- X * residuos                       # T x p, fila t = e_t * x_t'
  Omega <- crossprod(XE)                   # suma_t e_t^2 x_t x_t'
  if (L >= 1) {
    for (l in 1:L) {
      w <- 1 - l / (L + 1)
      Gl <- crossprod(XE[(l + 1):Tn, , drop = FALSE], XE[1:(Tn - l), , drop = FALSE])
      Omega <- Omega + w * (Gl + t(Gl))
    }
  }
  V_hac <- XtX_inv %*% Omega %*% XtX_inv
  se_hac <- sqrt(diag(V_hac))

  t_ord <- beta / se_ord
  p_ord <- 2 * (1 - stats::pt(abs(t_ord), df = gl))
  t_hac <- beta / se_hac
  p_hac <- 2 * (1 - stats::pt(abs(t_hac), df = gl))

  R2 <- 1 - sum(residuos^2) / sum((z - mean(z))^2)
  dw_estadistico <- sum(diff(residuos)^2) / sum(residuos^2)

  tabla_coeficientes <- data.frame(
    parametro     = nombres,
    estimacion    = beta,
    se_ordinario  = se_ord,
    se_hac        = se_hac,
    t_ordinario   = t_ord,
    p_ordinario   = p_ord,
    t_hac         = t_hac,
    p_hac         = p_hac
  )

  # valores ajustados en la escala ORIGINAL de y
  if (tipo == "exponencial") {
    correccion <- if (corregir_sesgo) exp(sigma2 / 2) else 1
    yhat <- exp(ajustado_z) * correccion
  } else {
    yhat <- ajustado_z
    correccion <- 1
  }

  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    t_nuevo <- Tn + seq_len(h)
    Xnuevo <- if (tipo == "cuadratica") cbind(1, t_nuevo, t_nuevo^2) else cbind(1, t_nuevo)
    z_pred <- as.numeric(Xnuevo %*% beta)
    if (tipo == "exponencial") exp(z_pred) * correccion else z_pred
  }

  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = list(
      tipo = tipo,
      beta = beta,
      tabla_coeficientes = tabla_coeficientes,
      R2 = R2,
      sigma2 = sigma2,
      durbin_watson = dw_estadistico,
      L_hac = L,
      k_regresores = p - 1,       # para durbin_watson(residuos, k = ...) de 03-evaluacion.R
      corregir_sesgo = corregir_sesgo,
      residuos = residuos,        
      nota_exponencial = if (tipo == "exponencial") {
        "e^(a+theta*t) estima la MEDIANA condicional, no la media; corregir_sesgo=TRUE multiplica por e^(sigma2/2) para aproximar la media"
      } else {
        NA_character_
      }
    )
  )
}

# 02-metodos.R: ajustar_holt() y optimizar()


ajustar_holt <- function(y, alpha, beta) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) >= 3,
            length(alpha) == 1, alpha > 0, alpha < 1,
            length(beta) == 1, beta > 0, beta < 1)
  Tn <- length(y)
  L <- rep(NA_real_, Tn)
  Th <- rep(NA_real_, Tn)
  yhat <- rep(NA_real_, Tn)

  L[1] <- y[1]
  Th[1] <- 0
  yhat[2] <- L[1] + Th[1]                 # calentamiento: L1=Y1, That1=0

  for (t in 2:Tn) {
    Yhat_t <- L[t - 1] + Th[t - 1]        # = yhat[t]
    L[t]  <- alpha * y[t] + (1 - alpha) * Yhat_t
    Th[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Th[t - 1]
    if (t + 1 <= Tn) yhat[t + 1] <- L[t] + Th[t]
  }

  L_T <- L[Tn]; Th_T <- Th[Tn]
  pronosticar <- function(h) {
    stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == round(h))
    L_T + Th_T * seq_len(h)
  }

  list(yhat = yhat, pronosticar = pronosticar,
       parametros = list(alpha = alpha, beta = beta, L_T = L_T, Th_T = Th_T,
                          L_t = L, Th_t = Th))
}

optimizar <- function(y, metodo = c("mm", "dmm", "ses", "holt"), rejilla) {
  metodo <- match.arg(metodo)
  stopifnot(is.numeric(y), !anyNA(y))

  mse_de <- function(ajuste_fn) {
    res <- tryCatch(ajuste_fn(), error = function(e) NULL)
    if (is.null(res)) return(NA_real_)
    mean((y - res$yhat)^2, na.rm = TRUE)
  }

  if (metodo %in% c("mm", "dmm")) {
    stopifnot(is.numeric(rejilla))
    ajustar_fn <- if (metodo == "mm") ajustar_mm else ajustar_dmm
    mse <- vapply(rejilla, function(k) mse_de(function() ajustar_fn(y, k)), numeric(1))
    tabla <- data.frame(k = rejilla, MSE = mse)

  } else if (metodo == "ses") {
    stopifnot(is.numeric(rejilla))
    mse <- vapply(rejilla, function(a) mse_de(function() ajustar_ses(y, a)), numeric(1))
    tabla <- data.frame(alpha = rejilla, MSE = mse)

  } else {  # holt
    stopifnot(is.list(rejilla), all(c("alpha", "beta") %in% names(rejilla)))
    combos <- expand.grid(alpha = rejilla$alpha, beta = rejilla$beta)
    combos$MSE <- mapply(function(a, b) mse_de(function() ajustar_holt(y, a, b)),
                          combos$alpha, combos$beta)
    tabla <- combos
  }

  if (all(is.na(tabla$MSE))) stop("ningun punto de la rejilla produjo un ajuste valido")
  optimo <- tabla[which.min(tabla$MSE), ]
  list(tabla = tabla, optimo = optimo)
}
