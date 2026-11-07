#' Plaque psoriasis aggregate-level data
#'
#' Two data frames, plaque_psoriasis_ipd and plaque_psoriasis_agd, containing (simulated) individual patient data from four studies and aggregate data from five studies (Phillippo 2019). Outcomes are binary success/failure to achieve 75%, 90%, or 100% reduction in symptoms on the Psoriasis Area and Severity Index (PASI) scale.
#'
#' @docType data
#'
#' @usage
#'
#' data(plaque_psoriasis_ipd)
#'
#' data(plaque_psoriasis_agd))
#'
#' @format The individual patient data are contained in a data frame plaque_psoriasis_ipd with 4118 rows, one per individual, and 16 variables:
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
#' @source ...
"plaque_psoriasis_agd"
