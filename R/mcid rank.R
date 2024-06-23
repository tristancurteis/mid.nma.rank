rm(list = ls())


#' Treatment rankings and rank probabilities
#'
#' Produce posterior treatment rankings and rank probabilities from a fitted NMA
#' model. Specify a minimally important difference (MID) specific to the outcome
#' and its scale to require ranking to be based on an MID.
#' When a meta-regression is fitted with effect modifier interactions
#' with treatment, posterior ranks will differ by study population.
#'
#' @param x A `stan_nma` object created by [nma()]
#' @param newdata Only used if a regression model is fitted. A data frame of
#'   study details, one row per study, giving the covariate values at which to
#'   produce relative effects. Column names must match variables in the
#'   regression model. If `NULL`, relative effects are produced for all studies
#'   in the network.
#' @param study Column of `newdata` which specifies study names, otherwise
#'   studies will be labelled by row number.
#' @param lower_better Logical, are lower treatment effects better (`TRUE`;
#'   default) or higher better (`FALSE`)? See details.
#' @param probs Numeric vector of quantiles of interest to present in computed
#'   summary, default `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#' @param summary Logical, calculate posterior summaries? Default `TRUE`.
#' @param sucra Logical, calculate the surface under the cumulative ranking
#'   curve (SUCRA) for each treatment? Default `FALSE`.
#' @param mid Numeric, a positive value, account for minimally important difference?
#' Default of 0, i.e. no minimally important difference.
#'
#' @return A [nma_summary] object if `summary = TRUE`, otherwise a list
#'   containing a 3D MCMC array of samples and (for regression models) a data
#'   frame of study information.
#' @export
#'
#' @seealso [plot.nma_summary()] for plotting the ranks and rank probabilities.
#'
#' @details The function `posterior_ranks_mid()` produces posterior rankings, which
#'   have a distribution (e.g. mean/median rank and 95% Credible Interval). The
#'   function `posterior_rank_probs_mid()` produces rank probabilities, which give
#'   the posterior probabilities of being ranked first, second, etc. out of all
#'   treatments.
#'
#'   When a minimally important difference (MID) is specified, treatments will
#'   only be ranked as superior to one another (i.e. ranked higher or lower than
#'   one another) only if an MID exists between them.
#'
#'   The argument `lower_better` specifies whether lower treatment
#'   effects or higher treatment effects are preferred. For example, with a
#'   negative binary outcome lower (more negative) log odds ratios are
#'   preferred, so `lower_better = TRUE`. Conversely, for example, if treatments
#'   aim to increase the rate of a positive outcome then `lower_better = FALSE`.
#'
#'   @examples
#'   # example code
#'
#'   library(multinma)
#'   library(mid.rank.nma)
#'
#'   #load data
#'   data(Senn2013)
#'   head(Senn2013)
#'
#'   Senn2013$se[Senn2013$study == "Willms1999"][2] <- 0.4 #Impute standard error for reference arm
#'
#'   #prepare network
#'   diabetes_net <- set_agd_contrast(Senn2013,
#'                                 study = study,
#'                                 trt = trt,
#'                                 y = y,
#'                                 se = se,
#'                                 trt_ref = "Placebo")
#'
#'   #fit model
#'   diabetes_fit_FE <- nma(diabetes_net,
#'                       trt_effects = "fixed",
#'                       prior_intercept = normal(scale = 100),
#'                       prior_trt = normal(scale = 100))
#'
#'
#'
#'   # With minimally important difference
#'   mid <- 0.5
#'
#'   rank_mid <- posterior_ranks_mid(diabetes_fit_FE, summary = TRUE, sucra = TRUE, mid = mid)
#'   rank_probs_mid <- posterior_rank_probs_mid(diabetes_fit_FE, mid = mid)
#'   rank_probs_cum_mid <- posterior_rank_probs_mid(diabetes_fit_FE, cumulative = T, mid = mid)
#'
#'   plot(rank_mid)
#'   plot(rank_probs_mid)
#'   plot(rank_probs_cum_mid)
#'
#'

posterior_ranks_mid <- function(x,
                            newdata = NULL,
                            study = NULL,
                            lower_better = TRUE,
                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                            sucra = FALSE,
                            summary = TRUE,
                            mid = 0) {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")

  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")

  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")

  if (!is.numeric(mid))
    abort("`mid` should be numeric.")

  check_probs(probs)

  # Cannot produce relative effects for inconsistency models
  if (x$consistency != "consistency")
    abort(glue::glue("Cannot produce ranks under inconsistency '{x$consistency}' model."))

  # Get reference treatment, number of treatments
  trt_ref <- levels(x$network$treatments)[1]
  ntrt <- nlevels(x$network$treatments)

  # All other checks handled by relative_effects()
  rel_eff <- relative_effects(x = x, newdata = newdata, study = {{ study }},
                              all_contrasts = FALSE, summary = FALSE)

  studies <- rel_eff$studies

  if (is.null(studies)) { # No study-specific treatment effects
    # Add zeros for d[1]
    dim_d <- dim(rel_eff$sim)
    dim_d[3] <- dim_d[3] + 1
    dimnames_d <- dimnames(rel_eff$sim)
    dimnames_d[[3]] <- c(paste0("d[", trt_ref, "]"), dimnames_d[[3]])
    d <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
    d[ , , 1] <- 0
    d[ , , 2:ntrt] <- rel_eff$sim

    # Calculate ranks
    if (mid > 0) {    #If ranking based on mid

      rk <- mid_based_ranks(d, mid = mid, dim_d = dim_d,
                             dimnames_d = dimnames_d, lower_better = lower_better)

    }else if(mid==0){                 #If ranking not based on mid

      # Get ranks at each iteration
      rk <- aperm(apply(d, 1:2, rank, ties.method = "average"),
                  c("iterations", "chains", "parameters"))

      # If higher treatment effects are better
      if (!lower_better) rk <- ntrt + 1 - rk

    }

    # Rename parameters
    dimnames(rk)[[3]] <- paste0("rank[", levels(x$network$treatments), "]")

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs) %>%
        tibble::add_column(.trt = x$network$treatments, .before = 1)

      if (sucra) {
        # Calculate SUCRA using scaled mean rank relation of Rucker and Schwarzer (2015)
        sucras <- unname((ntrt - rk_summary$mean) / (ntrt - 1))
        rk_summary <- tibble::add_column(rk_summary, sucra = sucras, .after = "sd")
      }

      out <- list(summary = rk_summary, sims = rk)
    } else {
      out <- list(sims = rk)
    }

  } else { # Study-specific treatment effects   ###CMi - to be updated
    nstudy <- nrow(studies)

    d <- rel_eff$sim
    d_names <- dimnames(d)[[3]]

    # Calculate ranks within each study population
    dim_rk <- dim(d)
    dim_rk[3] <- dim_rk[3] + nstudy
    rk_names <- vector("character", dim_rk[[3]])
    dimnames_rk <- dimnames(d)
    dimnames_rk[[3]] <- rk_names
    rk <- array(NA_real_, dim = dim_rk, dimnames = dimnames_rk)

    dim_temp_d <- dim(d)
    dim_temp_d[[3]] <- ntrt
    dimnames_temp_d <- dimnames(d)
    dimnames_temp_d[[3]] <- paste0("d[", levels(x$network$treatments), "]")
    temp_d <- array(NA_real_, dim = dim_temp_d, dimnames = dimnames_temp_d)
    temp_d[ , , 1] <- 0

    for (i in seq_len(nstudy)) {
      rk_names[(i - 1)*ntrt + 1:ntrt] <-
        c(paste0("d[", studies$.study[i], ": ", trt_ref, "]"),
          d_names[(i - 1)*(ntrt - 1) + 1:(ntrt - 1)])

      temp_d[ , , 2:ntrt] <- d[ , , (i - 1)*(ntrt - 1) + 1:(ntrt - 1)]

      # Calculate ranks
      if (mid > 0) {  #If ranking based on mid

        rk[ , , (i - 1)*ntrt + 1:ntrt] <- mid_based_ranks(temp_d, mid = mid, dim_d = dim_temp_d, dimnames_d = dimnames_temp_d, lower_better = lower_better)

      }else if(mid == 0){                 #If ranking not based on mid

        # Get ranks at each iteration
        rk[ , , (i - 1)*ntrt + 1:ntrt] <-
          aperm(apply(temp_d, 1:2, rank, ties.method = "average"),
                c("iterations", "chains", "parameters"))
      }


      # If higher treatment effects are better   #TCu: Moved inside else statement (2023-11-12)
      if (!lower_better) rk <- ntrt + 1 - rk
    }


    # Rename parameters
    dimnames(rk)[[3]] <- gsub("^d\\[", "rank\\[", rk_names)

    # Get summaries
    if (summary) {
      rk_summary <- summary_mcmc_array(rk, probs = probs) %>%
        # Add in study info
        tibble::add_column(.study = rep(studies$.study, each = ntrt),
                           .trt = rep(x$network$treatments, times = nstudy),
                           .before = 1)

      if (sucra) {
        # Calculate SUCRA using scaled mean rank relation of Rucker and Schwarzer (2015)
        sucras <- unname((ntrt - rk_summary$mean) / (ntrt - 1))
        rk_summary <- tibble::add_column(rk_summary, sucra = sucras, .after = "sd")
      }

      out <- list(summary = rk_summary, sims = rk, studies = studies)
    } else {
      out <- list(sims = rk, studies = studies)
    }
  }

  if (summary) {
    class(out) <- c("nma_ranks", "nma_summary")
    attr(out, "xlab") <- "Treatment"
    attr(out, "ylab") <- "Posterior Rank"
  }
  return(out)
}


####posterior rank probs function ####

#' @param cumulative Logical, return cumulative rank probabilities? Default is
#'   `FALSE`, return posterior probabilities of each treatment having a given
#'   rank. If `TRUE`, cumulative posterior rank probabilities are returned for
#'   each treatment having a given rank or better.
#' @export
#' @rdname posterior_ranks_mid

posterior_rank_probs_mid <- function(x,
                                 newdata = NULL,
                                 study = NULL,
                                 lower_better = TRUE,
                                 cumulative = FALSE,
                                 sucra = FALSE,
                                 mid = NULL) {
  # Checks
  if (!rlang::is_bool(cumulative))
    abort("`cumulative` should be TRUE or FALSE.")

  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")

  if (!is.numeric(mid))
    abort("`mid` should be numeric.")

  # All other checks handled by posterior_ranks_mid()
  rk <- posterior_ranks_mid(x = x, newdata = newdata, study = {{ study }},
                        lower_better = lower_better, summary = FALSE, mid = mid)

  ntrt <- nlevels(x$network$treatments)
  studies <- rk$studies

  if (is.null(studies)) { # No study-specific treatment effects

    p_rank <- apply(apply(rk$sims, 1:3, `==`, seq(1,ntrt, 0.5)), c(4, 1), mean)  #TCu (2024-03-23: this requires updating to allow for .5s) #CMi (24Apr2024 - adjusted accordingly)
    if (cumulative) p_rank <- t(apply(p_rank, 1, cumsum))

    if (sucra) {
      if (cumulative) sucras <- rowMeans(p_rank[, -ntrt, drop = FALSE])
      else sucras <- rowMeans(t(apply(p_rank, 1, cumsum))[, -ntrt, drop = FALSE])
    }

    rownames(p_rank) <- stringr::str_replace(rownames(p_rank), "^rank\\[", "d\\[")
    colnames(p_rank) <- paste0("p_rank[", seq(1,ntrt, 0.5), "]") #TCu (2024-03-23: this requires updating to allow for .5s) #CMi (24Apr2024 - adjusted accordingly)

    p_rank <- tibble::as_tibble(p_rank, rownames = "parameter") %>%
      tibble::add_column(.trt = x$network$treatments, .before = 1)

    #Calculate p best or equal, i.e. lowest rank (not necessarily 1st)
    min_rank <- apply(rk$sims, 1:2, min)  #lowest rank for each iteration
    min_cond <- array(NA_real_, dim = dim(rk$sims), dimnames = dimnames(rk$sims))

    for(i in 1:dim(rk$sims)[1]){ #for each iteration

      for(j in 1:dim(rk$sims)[2]){ # for each chain

        min_cond[i,j,] <- 1*(rk$sims[i,j,] == min_rank[i,j])

      }
    }

    p_rank$p_mid_best <- apply(min_cond, 3, mean) #probability within an MID of lowest rank

    if (sucra) p_rank$sucra <- unname(sucras)

    out <- list(summary = p_rank)

  } else { # Study-specific treatment effects

    nstudy <- nrow(studies)

    p_rank <- matrix(NA_real_, nrow = ntrt * nstudy, ncol = length(seq(1, ntrt, 0.5)))   #TCu (2024-03-23: this requires updating to allow for .5s)
    rk_sims <- rk$sims
    for (i in seq_len(nstudy)) {
      p_rank[(i - 1)*ntrt + 1:ntrt, ] <-
        apply(apply(rk_sims[ , , (i - 1)*ntrt + 1:ntrt], 1:3, `==`, seq(1, ntrt, 0.5)), c(4, 1), mean)
    }

    if (cumulative) p_rank <- t(apply(p_rank, 1, cumsum))

    if (sucra) {
      if (cumulative) sucras <- rowMeans(p_rank[, -ntrt, drop = FALSE])
      else sucras <- rowMeans(t(apply(p_rank, 1, cumsum))[, -ntrt, drop = FALSE])
    }

    rownames(p_rank) <- stringr::str_replace(dimnames(rk_sims)[[3]], "^rank\\[", "d\\[")
    colnames(p_rank) <- paste0("p_rank[", seq(1, ntrt, 0.5), "]")  #TCu (2024-03-23: this requires updating to allow for .5s)

    p_rank <- tibble::as_tibble(p_rank, rownames = "parameter") %>%
      # Add in study info
      tibble::add_column(.study = rep(studies$.study, each = ntrt),
                         .trt = rep(x$network$treatments, times = nstudy),
                         .before = 1)

    if (sucra) p_rank$sucra <- unname(sucras)

    out <- list(summary = p_rank, studies = studies)

  }

  class(out) <- c("nma_rank_probs", "nma_summary")
  attr(out, "xlab") <- "Rank"
  attr(out, "ylab") <- if (cumulative) "Cumulative Rank Probability" else "Rank Probability"
  return(out)
}

####mid based ranks ranks function ####


#' @param d Array of relative treatment effects as defined above, based on relative_effects

#' Calculate mcid adjusted ranks
#' @noRd
mid_based_ranks <- function(x, mid, dim_d, dimnames_d, lower_better){

  if (!is.array(x) || length(dim(x)) != 3) abort("Error in mid_based_ranks: First parameter is not a 3D MCMC array, [Iterations, Chains, Parameters]")

  #Intermediary storage for mid ranking calculations
  test <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
  nsup <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)

  #Sequential calculations to reduce memory burden

  for(i in 1:dim_d[3]){    #For trt i

    for(j in 1:dim_d[3]){  #For trt j


      if (!lower_better){
        test[,,j] <- 1*(x[,,i] - x[,,j] - mid > 0)     # If positive orientation, test whether trt i is mid superior to trt j
        # Test difference from mid in mean difference, assuming higher mean diff is better; Events are good => +ve xi-xj and xi-xj-delta>0, means i more effective than j
      }

      if (lower_better){
        test[,,j] <- 1*(x[,,j] - x[,,i] - mid > 0)     #If negative orientation, test whether trt i is mid superior to trt j
        # Test difference from mid in mean difference, assuming lower mean diff is better; Events are bad =>  -ve xi-xj and -(xi-xj)-delta>=0 means j more effective than i
      }


    }

    nsup[,,i]  <-  apply(test[,,], 1:2, sum)    #Calculate number of mid superiorities for trt i

    #Note if mid = 0, and decision rule includes = sign, then require next line (only for case mid==0)
    #nsup[,,i]  <-  apply(test[,,], 1:2, sum) - test[,,i]

  }

  #Get ranks at each iteration
  rk <- aperm(apply(nsup, 1:2, rank, ties.method = "average"),
              c("iterations", "chains", "parameters"))

  rk <- dim_d[3] + 1 - rk    #since based on number of superiorities

  return(rk)

}

####summary mcmc array function ####

summary_mcmc_array <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  if (!is.array(x) || length(dim(x)) != 3) abort("Not a 3D MCMC array, [Iterations, Chains, Parameters]")
  check_probs(probs)

  pars <- dimnames(x)[[3]]
  p_mean <- apply(x, 3, mean)
  p_sd <- apply(x, 3, sd)
  p_ess_bulk <- apply(x, 3, rstan::ess_bulk)
  p_ess_tail <- apply(x, 3, rstan::ess_tail)
  p_rhat <- apply(x, 3, rstan::Rhat)
  # p_se_mean <- p_sd / sqrt(apply(x, 3, rstan:::ess_mean))

  p_quan <- apply(x, 3, quantile, probs = probs)
  if (length(probs) == 1) {
    p_quan <- tibble::tibble(!! paste0(probs*100, "%") := p_quan)
  } else {
    p_quan <- as.data.frame(t(p_quan))
  }

  ss <- tibble::tibble(
    parameter = pars,
    mean = p_mean,
    # se_mean = p_se_mean,
    sd = p_sd,
    !!! p_quan,
    Bulk_ESS = p_ess_bulk, Tail_ESS = p_ess_tail, Rhat = p_rhat)

  return(ss)
}

#' Validate probs argument
#' @noRd
check_probs <- function(probs) {
  if (!rlang::is_double(probs, finite = TRUE) || any(probs < 0) || any(probs > 1))
    rlang::abort("`probs` must be a numeric vector of probabilities between 0 and 1.")
}


