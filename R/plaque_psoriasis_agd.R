#' Plaque psoriasis aggregate-level data
#'
#' Aggregate data from five studies
#' (Phillippo, 2019).
#'
#' @docType data
#'
#' @format A data frame with 15 rows and 26 variables:
#' \describe{
#'   \item{studyc}{study name}
#'   \item{trtc_long}{treatment name (long format)}
#'   \item{trtc}{treatment name}
#'   \item{trtn}{numeric treatment code}
#'   \item{pasi75_r}{PASI 75 outcome count}
#'   \item{pasi75_n}{PASI 75 outcome denominator}
#'   \item{pasi90_r}{PASI 90 outcome count}
#'   \item{pasi90_n}{PASI 90 outcome denominator}
#'   \item{pasi100_r}{PASI 100 outcome count}
#'   \item{pasi100_n}{PASI 100 outcome denominator}
#'   \item{sample_size_w0}{sample size at week zero}
#'   \item{age_mean}{mean age (years)}
#'   \item{age_sd}{standard deviation of age (years)}
#'   \item{bmi_mean}{mean BMI}
#'   \item{bmi_sd}{standard deviation of BMI}
#'   \item{pasi_w0_mean}{mean PASI score at week 0}
#'   \item{pasi_w0_sd}{standard deviation of PASI score at week 0}
#'   \item{male}{percentage of males}
#'   \item{bsa_mean}{mean body surface area (percent)}
#'   \item{bsa_sd}{standard deviation of body surface area (percent)}
#'   \item{weight_mean}{mean weight (kilograms)}
#'   \item{weight_sd}{standard deviation of weight (kilograms)}
#'   \item{durnpso_mean}{mean duration of psoriasis (years)}
#'   \item{durnpso_sd}{standard deviation of duration of psoriasis (years)}
#'   \item{prevsys}{percentage of individuals with previous systemic treatment}
#'   \item{psa}{percentage of individuals with psoriatic arthritis}
#' }
#'
#' @references
#' Phillippo DM (2019). Calibration of Treatment Effects in Network
#' Meta-Analysis using Individual Patient Data. PhD thesis,
#' University of Bristol.
"plaque_psoriasis_agd"
