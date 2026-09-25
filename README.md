<https://github.com/jumadi13/st-2026-2-tarea1-herramientas-diazlarrea-juanmanuel>

# Tarea 1 — Caja de herramientas de pronóstico

## Qué contiene el repositorio

| Archivo | Contenido | Dependencias |
|------------------------|------------------------|------------------------|
| `R/00-lectura.R` | `leer_serie()` | R base |
| `R/01-graficos.R` | `graficar_serie()`, `correlograma()` | ggplot2, patchwork, tibble |
| `R/02-metodos.R` | Los ocho métodos de pronóstico y `optimizar()` | R base |
| `R/03-evaluacion.R` | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `validar_errores()` | ggplot2, tibble |
| `ejemplos/driver.R` | `ejecutar_ejemplo()`: corre el protocolo de 7 pasos sobre una serie (ayudante interno, no exigido por nombre) | tibble, ggplot2 |
| `ejemplos/ejemplos.R` | Corre los 8 ejemplos + el contraejemplo de punta a punta | dplyr, purrr, tibble, ggplot2, patchwork |
| `informe/informe.qmd` | Análisis narrado de los 9 casos | \+ knitr (vía Quarto) |
| `figs/` | Imágenes generadas por `ejemplos.R` | — |

## Cómo se corre

``` bash
git clone https://github.com/jumadi13/st-2026-2-tarea1-herramientas-diazlarrea-juanmanuel.git
cd st-2026-2-tarea1-herramientas-diazlarrea-juanmanuel
Rscript -e 'install.packages(c("dplyr","tidyr","purrr","tibble","ggplot2","patchwork"))'
Rscript ejemplos/ejemplos.R
```

Para el informe: abrir `informe/informe.qmd` en RStudio y usar el botón **Render**, o desde terminal:

``` bash
quarto render informe/informe.qmd
```

## Cómo se usan las funciones

Ejemplo mínimo con `ajustar_ses()`:

``` r
source("R/00-lectura.R")
datos <- leer_serie(discoveries, fuente = "The World Almanac and Book of Facts, 1975", unidad = "conteo")
source("R/02-metodos.R")
ajuste <- ajustar_ses(datos$y, alpha = 0.22)
ajuste$pronosticar(h = 5)   # 5 pronosticos extramuestrales
```

## Convenciones que fijan los números

-   **ACF (`correlograma()`)**: divisor único $T$ (Definición 2.5), no $T-h$. Verificado contra `acf(plot=FALSE)$acf` con diferencia máxima $<10^{-12}$.
-   **Calentamiento**: media simple y SES, $\hat Y_2=Y_1$; media móvil y DMM, $k$ y $2k-1$ períodos respectivamente; tendencias, sin calentamiento (ajuste global); Holt, $L_1=Y_1$, $\hat T_1=0$.
-   **Grados de libertad de Ljung-Box** ($p$ en $\chi^2_{m-p}$): $p=1$ para media simple y SES; $p=0$ para media móvil y DMM (el ancho $k$ se trata como banda de suavizamiento fija, no como parámetro ajustado por mínimos cuadrados); $p=2$ para tendencia lineal/exponencial y Holt; $p=3$ para tendencia cuadrática.
-   **Durbin-Watson**: se usa la tabla clásica de cotas $d_L,d_U$ al 5% (transcrita en `R/03-evaluacion.R`) cuando el tamaño de muestra está cubierto por ella ($T\le200$); para T=456 (ejemplo de `co2`) ninguna tabla publicada cubre ese tamaño, así que se usa la aproximación asintótica estándar $d\approx2(1-\hat\rho_1)$, $\hat\rho_1\sim N(0,1/T)$ bajo $H_0$ — no es una sustitución arbitraria, es la práctica documentada para T fuera de rango.
-   **Tendencia exponencial**: se estima sobre $\ln(Y_t)$; el pronóstico $e^{\hat a+\hat\theta t}$ estima la mediana condicional, no la media; `corregir_sesgo=TRUE` multiplica por $e^{\hat\sigma^2/2}$.
-   **Partición**: $h=\min(12,\lfloor 0.2T\rfloor)$; en la serie estacional (contraejemplo, `nottem`) se exige además al menos un ciclo completo ($s=12$), condición que aquí ya se cumple con ese mismo $h$.

## Resumen de resultados

| Ejemplo | Parámetros | MASE método | MASE ingenuo |
|------------------|------------------|------------------|------------------|
| 1\. Media simple — discoveries | (sin constante) | 0.919 | 1.136 |
| 2\. Media móvil — Nile | k=9 | 0.825 | 0.835 |
| 3\. SES — discoveries | α=0.22 | 0.615 | 1.136 |
| 3b. Contraejemplo: SES — nottem | α=0.98 | 3.741 | 0.619 |
| 4\. Doble media móvil — austres | k=2 | 1.863 | 6.091 |
| 5\. Tendencia lineal — LakeHuron | β₀=580.3, β₁=−0.0285 | 2.067 | 2.164 |
| 6\. Tendencia cuadrática — co2 | β₀=314.9, β₁=0.0657, β₂=9.33e-05 | 1.904 | 1.916 |
| 7\. Tendencia exponencial — airmiles | a=6.032, θ=0.212 | 20.807 | 4.496 |
| 8\. Holt — WWWusage | α=0.95, β=0.95 | 1.890 | 7.625 |

## Declaración de uso de IA

Usé Claude (Anthropic) como asistente durante el desarrollo de esta tarea. En términos generales: - **Qué pedí**: ayuda para planear el orden de trabajo de la tarea, resolver dudas puntuales sobre las fórmulas de las notas (ACF, HAC, Durbin-Watson), y para escribir/depurar el código de `R/*.R`, `ejemplos/*.R` e `informe.qmd` principalmente. - **Qué recibí**: explicaciones de conceptos (diferencia ACF/PACF, por qué Durbin-Watson no tiene tablas para T grande), código para las funciones de `R/`, y una plantilla de redacción para el informe. - **Qué verifiqué por cuenta propia**: corrí personalmente `Rscript ejemplos/ejemplos.R` en mi propio computador y confirmé que los números coinciden con los reportados; revisé cada fórmula contra las notas de clase.
