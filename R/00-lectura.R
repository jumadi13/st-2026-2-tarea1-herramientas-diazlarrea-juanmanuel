# 00-lectura.R
# leer_serie(): acepta un objeto ts o la ruta a un csv (columnas fecha, valor)
# y devuelve una tibble (t, fecha, y) con los atributos frecuencia, fuente, unidad.

.fecha_desde_ts <- function(inicio, frecuencia, n) {
  # inicio = c(año, periodo_1indexado)
  anio0 <- inicio[1]
  periodo0 <- inicio[2] - 1L  
  idx0 <- anio0 * frecuencia + periodo0
  idx <- idx0 + (0:(n - 1))
  anio <- idx %/% frecuencia
  periodo <- idx %% frecuencia            
  mes <- floor(periodo * 12 / frecuencia) + 1L
  as.Date(sprintf("%04d-%02d-01", anio, mes))
}

.inferir_frecuencia <- function(fecha) {
  dias <- as.numeric(diff(fecha))
  paso_medio <- stats::median(dias)
  candidatos <- c(1, 7, 30.44, 91.31, 365.25)
  frecuencias <- c(365, 52, 12, 4, 1)
  frecuencias[which.min(abs(candidatos - paso_medio))]
}

.verificar_fechas <- function(fecha, frecuencia) {
  if (is.unsorted(fecha, strictly = TRUE)) {
    stop("las fechas deben ser estrictamente crecientes")
  }
  n <- length(fecha)
  if (n < 2) return(invisible(TRUE))
  dias <- as.numeric(diff(fecha))
  paso_esperado <- 365.25 / frecuencia
  tol <- max(paso_esperado * 0.15, 3)  # al menos 3 dias de tolerancia
  if (any(abs(dias - paso_esperado) > tol)) {
    stop(sprintf(
      "las fechas no son equiespaciadas segun la frecuencia declarada (%s obs/anio): paso esperado ~%.1f dias, maxima desviacion observada %.1f dias",
      frecuencia, paso_esperado, max(abs(dias - paso_esperado))
    ))
  }
  invisible(TRUE)
}

leer_serie <- function(x, fuente, unidad) {
  stopifnot(is.character(fuente), length(fuente) == 1, nzchar(fuente))
  stopifnot(is.character(unidad), length(unidad) == 1, nzchar(unidad))

  if (stats::is.ts(x)) {

    y <- as.numeric(x)
    n <- length(y)
    frecuencia <- stats::frequency(x)
    fecha <- .fecha_desde_ts(stats::start(x), frecuencia, n)

  } else if (is.character(x) && length(x) == 1) {

    stopifnot(file.exists(x))
    crudo <- utils::read.csv(x, stringsAsFactors = FALSE)
    stopifnot(all(c("fecha", "valor") %in% names(crudo)))
    fecha <- as.Date(crudo$fecha)
    y <- as.numeric(crudo$valor)
    n <- length(y)
    frecuencia <- .inferir_frecuencia(fecha)

  } else {
    stop("x debe ser un objeto ts o la ruta a un archivo csv con columnas fecha y valor")
  }

  stopifnot(!anyNA(y))
  .verificar_fechas(fecha, frecuencia)

  datos <- tibble::tibble(t = seq_len(n), fecha = fecha, y = y)
  attr(datos, "frecuencia") <- frecuencia
  attr(datos, "fuente") <- fuente
  attr(datos, "unidad") <- unidad
  datos
}
