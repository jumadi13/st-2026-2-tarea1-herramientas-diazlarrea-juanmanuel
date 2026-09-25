# 01-graficos.R
# graficar_serie(): grafico de la serie en el tiempo con ggplot2
# correlograma(): ACF + PACF vía stats::pacf()

graficar_serie <- function(datos, titulo) {
  stopifnot(is.data.frame(datos), all(c("fecha", "y") %in% names(datos)))
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  n <- nrow(datos)

  ggplot2::ggplot(datos, ggplot2::aes(x = fecha, y = y)) +
    ggplot2::geom_line(color = "#1f3b57", linewidth = 0.6) +
    ggplot2::scale_x_date(date_labels = "%Y") +
    ggplot2::labs(
      title = titulo,
      x = "Fecha",
      y = if (is.null(unidad)) "y" else unidad,
      caption = sprintf(
        "Fuente: %s. n = %d observaciones.",
        if (is.null(fuente)) "no declarada" else fuente, n
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, face = "italic")
    )
}

correlograma <- function(datos, m = NULL) {
  # acepta un tibble con columna y, o un vector numerico plano
  if (is.list(datos) && "y" %in% names(datos)) {
    y <- datos$y
  } else {
    y <- as.numeric(datos)
  }
  stopifnot(is.numeric(y), !anyNA(y))

  Tn <- length(y)
  if (is.null(m)) m <- min(floor(Tn / 4), 24)
  stopifnot(m >= 1, m < Tn)

  ybar <- mean(y)
  denom <- sum((y - ybar)^2)  

  r <- vapply(seq_len(m), function(h) {
    num <- sum((y[(h + 1):Tn] - ybar) * (y[1:(Tn - h)] - ybar))
    num / denom
  }, numeric(1))

  pacf_vals <- as.numeric(stats::pacf(y, lag.max = m, plot = FALSE)$acf)

  banda <- stats::qnorm(0.975) / sqrt(Tn)  

  df_acf  <- tibble::tibble(rezago = seq_len(m), valor = r)
  df_pacf <- tibble::tibble(rezago = seq_len(m), valor = pacf_vals)

  tema <- ggplot2::theme_minimal(base_size = 12)

  p_acf <- ggplot2::ggplot(df_acf, ggplot2::aes(x = rezago, y = valor)) +
    ggplot2::geom_col(width = 0.15, fill = "#1f3b57") +
    ggplot2::geom_hline(yintercept = c(-banda, banda), linetype = "dashed", color = "firebrick") +
    ggplot2::geom_hline(yintercept = 0, color = "grey40") +
    ggplot2::labs(title = "ACF", x = NULL, y = "r_h") +
    tema

  p_pacf <- ggplot2::ggplot(df_pacf, ggplot2::aes(x = rezago, y = valor)) +
    ggplot2::geom_col(width = 0.15, fill = "#1f3b57") +
    ggplot2::geom_hline(yintercept = c(-banda, banda), linetype = "dashed", color = "firebrick") +
    ggplot2::geom_hline(yintercept = 0, color = "grey40") +
    ggplot2::labs(title = "PACF", x = "Rezago", y = "phi_hh") +
    tema

  list(
    grafico = patchwork::wrap_plots(p_acf, p_pacf, ncol = 1),
    acf = r,
    pacf = pacf_vals,
    banda = banda,
    m = m,
    n = Tn
  )
}
