#' Network meta-analysis in diabetes
#'
#' Network meta-analysis in diabetes comparing effects of a
#' number of drugs on the HbA1c value (%).
#' These data are used as an example in Senn et al. (2013) and have
#' been sourced from the R package netmeta, converted into long format.
#'
#' @docType data
#'
#' @usage data(Senn2013)
#'
#' @format A data frame with 53 rows and 4 variables:
#' \describe{
#'   \item{y}{Mean difference versus placebo in HbA1c (\%)}
#'   \item{se}{Standard error of mean difference}
#'   \item{trt}{Treatment}
#'   \item{study}{Study name}
#' }
#'
#' @keywords datasets
#'
#' @references Senn S, Gavini F, Magrez D, Scheen A (2013):
#' Issues in performing a network meta-analysis. Statistical Methods
#' in Medical Research, 22, 169–89
#'
#' Balduzzi S, Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G,
#' Efthimiou O, Schwarzer G (2023). “netmeta: An R Package for Network
#' Meta-Analysis Using Frequentist Methods.” Journal of Statistical Software,
#' 106(2), 1–40. doi:10.18637/jss.v106.i02.
#'
#' @examples
#' data(Senn2013)
#' head(Senn2013)
"Senn2013"
