# ejecutar_ejemplo(): corre el protocolo de 7 pasos sobre una serie

# construir_ajuste: function(y_estimacion) -> list(ajuste = <lista con
#   yhat/pronosticar/parametros>, parametros_elegidos = <lista para el
#   resumen>, opt = <resultado de optimizar(), o NULL si el metodo no
#   necesita rejilla>). Aqui es donde cada ejemplo corre su propio
#   optimizar() usando solo el tramo de estimacion.

ejecutar_ejemplo <- function(nombre, serie_ts, fuente, unidad, construir_ajuste,
                              p_metodo, k_dw = 1, estacional = FALSE, s = 1,
                              prefijo_fig = NULL) {

  datos <- leer_serie(serie_ts, fuente = fuente, unidad = unidad)
  Tn <- nrow(datos)
  h <- min(12, floor(0.2 * Tn))
  if (estacional) h <- max(h, s)

  idx_est <- 1:(Tn - h)
  idx_val <- (Tn - h + 1):Tn
  y_est <- datos$y[idx_est]
  y_val <- datos$y[idx_val]

  corr_serie <- correlograma(datos)
  lb_serie <- ljung_box(corr_serie$acf, T = Tn, m = corr_serie$m, p = 0)

  ai <- construir_ajuste(y_est)
  ajuste <- ai$ajuste

  mae_ref <- mae_ingenuo(y_est, estacional = estacional, s = s)

  med_estimacion <- medidas(y_est, ajuste$yhat, mae_ref)

  pron_val <- ajuste$pronosticar(h)
  med_validacion <- medidas(y_val, pron_val, mae_ref)

  pron_ingenuo <- if (estacional) {
    rep(tail(y_est, s), length.out = h)
  } else {
    rep(tail(y_est, 1), h)
  }
  med_ingenuo_val <- medidas(y_val, pron_ingenuo, mae_ref)

  if (!is.null(ajuste$parametros$residuos)) {
    errores_1paso <- ajuste$parametros$residuos          # tendencias (escala de z)
  } else {
    errores_1paso <- ajuste$yhat[!is.na(ajuste$yhat)]
    errores_1paso <- y_est[!is.na(ajuste$yhat)] - errores_1paso
  }
  val_err <- validar_errores(errores_1paso, p = p_metodo, k = k_dw)

  df_est <- tibble::tibble(fecha = datos$fecha[idx_est], y = y_est, ajuste = ajuste$yhat)
  df_val <- tibble::tibble(fecha = datos$fecha[idx_val], y = y_val, pronostico = pron_val)

  g_final <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = fecha, y = y), data = df_est, color = "grey40") +
    ggplot2::geom_line(ggplot2::aes(x = fecha, y = y), data = df_val, color = "grey40") +
    ggplot2::geom_line(ggplot2::aes(x = fecha, y = ajuste), data = df_est,
                        color = "#1f3b57", na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(x = fecha, y = pronostico), data = df_val,
                        color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = nombre, x = "Fecha", y = unidad,
                  caption = sprintf("gris = observado | azul = ajuste (estimacion) | rojo = pronostico (validacion) | n=%d, h=%d",
                                     Tn, h)) +
    ggplot2::theme_minimal(base_size = 11)

  if (!is.null(prefijo_fig)) {
    ggplot2::ggsave(sprintf("figs/%s_serie.png", prefijo_fig), graficar_serie(datos, nombre),
                     width = 7, height = 4, dpi = 110)
    ggplot2::ggsave(sprintf("figs/%s_correlograma.png", prefijo_fig), corr_serie$grafico,
                     width = 6, height = 6, dpi = 110)
    ggplot2::ggsave(sprintf("figs/%s_ajuste.png", prefijo_fig), g_final,
                     width = 7, height = 4, dpi = 110)
    ggplot2::ggsave(sprintf("figs/%s_errores.png", prefijo_fig), val_err$grafico_errores,
                     width = 7, height = 3, dpi = 110)
    ggplot2::ggsave(sprintf("figs/%s_correlograma_errores.png", prefijo_fig), val_err$correlograma$grafico,
                     width = 6, height = 6, dpi = 110)
  }

  list(
    nombre = nombre, datos = datos, Tn = Tn, h = h,
    corr_serie = corr_serie, lb_serie = lb_serie,
    ajuste = ajuste, parametros_elegidos = ai$parametros_elegidos, opt = ai$opt,
    mae_ref = mae_ref,
    med_estimacion = med_estimacion, med_validacion = med_validacion,
    med_ingenuo_val = med_ingenuo_val,
    val_err = val_err, grafico_final = g_final
  )
}
