rm(list = ls())


#' Treatment rankings and rank probabilities
#'
#' Produce posterior treatment rankings and rank probabilities from a fitted NMA or ML-NMR
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
#'   1. example code for application in NMA
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
#'   2. example code for application in ML-NMR
#'
#'   library(multinma)
#'   library(mid.rank.nma)
#'
#' data(plaque_psoriasis_ipd)
#' data(plaque_psoriasis_agd)
#' 
#' pso_ipd <- plaque_psoriasis_ipd
#' pso_agd <- plaque_psoriasis_agd
#' 
#' 
#' # IPD studies
#' pso_ipd <- plaque_psoriasis_ipd %>% 
#'   mutate(
#'     # Variable transformations
#'     bsa = bsa / 100,
#'     weight = weight / 10,
#'     durnpso = durnpso / 10,
#'     prevsys = as.numeric(prevsys),
#'     psa = as.numeric(psa),
#'     # Treatment classes
#'     trtclass = case_when(trtn == 1 ~ "Placebo",
#'                          trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
#'                          trtn == 4 ~ "TNFa blocker",
#'                          trtn == 7 ~ "IL-12/23 blocker"),
#'     # Check complete cases for covariates of interest
#'     is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
#'   ) %>% 
#'   arrange(studyc, trtn)
#' 
#' # AgD studies
#' pso_agd <- plaque_psoriasis_agd %>% 
#'   mutate(
#'     # Variable transformations
#'     bsa_mean = bsa_mean / 100, 
#'     bsa_sd = bsa_sd / 100,
#'     weight_mean = weight_mean / 10,
#'     weight_sd = weight_sd / 10,
#'     durnpso_mean = durnpso_mean / 10,
#'     durnpso_sd = durnpso_sd / 10,
#'     prevsys = prevsys / 100,
#'     psa = psa / 100,
#'     # Treatment classes
#'     trtclass = case_when(trtn == 1 ~ "Placebo",
#'                          trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
#'                          trtn == 4 ~ "TNFa blocker",
#'                          trtn == 7 ~ "IL-12/23 blocker")
#'   ) %>% 
#'   arrange(studyc, trtn)
#' 
#' 
#' pso_ipd %>% 
#'   group_by(studyc) %>% 
#'   summarise(n_total = n(),
#'             n_missing = sum(!is_complete), 
#'             pct_missing = mean(!is_complete) * 100)
#' 
#' #Remove patients with missing data (very few patients)
#' pso_ipd <- filter(pso_ipd, is_complete)
#' 
#' #===============================#
#' #####  Set up Network       #####
#' #===============================#
#' 
#' pso_net <- combine_network(
#'   set_ipd(pso_ipd,
#'           study = studyc,
#'           trt = trtc,
#'           r = pasi90,
#'           trt_class = trtclass),
#'   set_agd_arm(pso_agd,
#'               study = studyc,
#'               trt = trtc,
#'               r = pasi90_r,
#'               n = pasi90_n,
#'               trt_class = trtclass)
#' )
#' 
#' pso_net
#' 
#' #===============================#
#' #####     Integration    #####
#' #===============================#
#' 
#' pso_net <- add_integration(pso_net,
#'                            durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
#'                            prevsys = distr(qbern, prob = prevsys),
#'                            bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
#'                            weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
#'                            psa = distr(qbern, prob = psa),
#'                            n_int = 64)
#' 
#' 
#' #========================================#
#' #####          Fit Model             #####
#' #========================================#
#' 
#' #binomial logit
#' pso_fit_FE <- nma(pso_net,
#'                   trt_effects = "fixed",
#'                   link = "logit",
#'                   likelihood = "bernoulli2",
#'                   regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
#'                   class_interactions = "common",
#'                   prior_intercept = normal(scale = 10),
#'                   prior_trt = normal(scale = 10),
#'                   prior_reg = normal(scale = 10),
#'                   init_r = 0.1,
#'                   QR = TRUE)
#' 
#' #==========================================#
#' #####             Plot Output          #####
#' #==========================================#
#' 
#' beta_plot <- plot(pso_fit_FE,
#'      pars = "beta",
#'      stat = "halfeye",
#'      ref_line = 0)
#' 
#' #Relative effects on logit scale
#' pso_releff_FE <- relative_effects(pso_fit_FE)
#' forest_logit <- plot(pso_releff_FE, ref_line = 0)
#' 
#' #Relative effects on probability scale
#' pso_pred_FE <- predict(pso_fit_FE, type = "response")
#' forest_probs <- plot(pso_pred_FE, ref_line = c(0, 1))
#' 
#' #Relative effects on odds ratio
#' pso_or_FE <- predict(pso_fit_FE, scale = "oddsratio")
#' forest_OR <- plot(pso_or_FE, ref_line = c(0, 1))
#' 
#' #=====================#
#' #####     MID     #####
#' #=====================#
#' 
#' pso_agd$avg_pasi75perc <- pso_agd$pasi75_r/pso_agd$pasi75_n
#' mean(pso_agd$avg_pasi75perc)
#' 
#' #For aggregate studies, the average pasi75 rate for placebo arms is 3%
#' avg_pasi75_trn1 <- mean(pso_agd$pasi75_r[pso_agd$trtn == 1]/pso_agd$pasi75_n[pso_agd$trtn == 1], na.rm = TRUE)
#' avg_pasi75_trn1 
#' #probit scale
#' 
#' # midRD <- 0.1                #%10 MID on risk difference scale
#' # p1 <- 0.5 #mean(pso_agd$avg_pasi75perc) #avg_pasi75_trn1       # baseline probability
#' # p2 <- p1 + midRD            # 10 percentage-point increase from baseline probability
#' # delta_eta <- qnorm(p2) - qnorm(p1)
#' # mid <- delta_eta            #Minimally importance difference versus placebo on probit scale
#' # mid
#' 
#' #Log odds ratio scale
#' #mid <- log(1.5) #assuming at least 1.5 multiplicative increase in odds (i.e. OR of 1.5)
#' 
#' #=====================================#
#' #####     Motivated by IXORA-S     #####
#' #=====================================#
#' 
#' #Unique study and treatment names, IPD
#' uniq_combinations <- unique(pso_ipd[c("studyc", "trtc_long")])
#' uniq_combinations
#' 
#' #Research question:
#' #https://academic.oup.com/bjd/article/177/4/1014/6674120?login=false
#' #IXORA-S
#' #Assumed PASI 90% response rate of 43% for ustekinumab
#' #12.6 risk difference
#' 
#' #How does IXE_Q2W rank if we require a 12.6% risk difference to establish superior ranking?
#' 
#' #with baseline probability as well
#' p1 <- 0.43      # baseline probability
#' dp <- 0.126     # 10% absolute risk difference
#' p2 <- p1 + dp
#' logit <- function(p) log(p / (1 - p))
#' delta_logit <- logit(p2) - logit(p1)
#' mid <- delta_logit
#' exp(mid) # Corresponds to exceeding all treatments by an OR of:
#' 
#' #===============================#
#' #####     Ranking           #####
#' #===============================#
#' 
#' # Without minimally important difference, # With ties.method = minimum
#' rank <- posterior_ranks(pso_fit_FE, summary = TRUE, sucra = TRUE, lower_better = FALSE)
#' rank_probs <- posterior_rank_probs(pso_fit_FE, lower_better = FALSE)
#' rank_probs_cum <- posterior_rank_probs(pso_fit_FE, cumulative = T, lower_better = FALSE)
#' 
#' rank_0 <- plot(rank)
#' probs_0 <- plot(rank_probs)
#' cum_0 <- plot(rank_probs_cum)
#' 
#' rank_0
#' probs_0
#' cum_0
#' 
#' # With minimally important difference
#' 
#' rank_mid <- posterior_ranks_mid(pso_fit_FE, summary = TRUE, sucra = TRUE, mid = mid, lower_better = FALSE)
#' rank_probs_mid <- posterior_rank_probs_mid(pso_fit_FE, mid = mid, lower_better = FALSE)
#' rank_probs_cum_mid <- posterior_rank_probs_mid(pso_fit_FE, cumulative = T, mid = mid, lower_better = FALSE)
#' 
#' # rank_mid$summary %>% 
#' #   filter(.study == "IXORA-S") 
#' # 
#' # #select parameter column of rank_probs_mid$summary and filter for IXORA-S study
#' # rank_probs_mid$summary %>% 
#' #   filter(.study == "IXORA-S") %>% select(parameter) 
#' 
#' 
#' rank_mid_plot <- plot(rank_mid)
#' probs_mid <- plot(rank_probs_mid)
#' cum_mid <- plot(rank_probs_cum_mid)
#' 
#' rank_mid_plot
#' probs_mid
#' cum_mid 



posterior_ranks_mid <- function(x,
                                newdata = NULL,
                                study = NULL,
                                lower_better = TRUE,
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                sucra = FALSE,
                                summary = TRUE,
                                mid = 0,
                                method = "midpoint") {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")
  
  if (!is.numeric(mid))
    abort("`mid` should be numeric.")
  
  if (!method %in% c("midpoint", "min", "max"))
    abort("parameter method must be one of: \"midpoint\", \"min\", or \"max\".")
  
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
      
      rk <- mid_based_ranks(d, mid = mid, 
                            dim_d = dim_d,
                            dimnames_d = dimnames_d, 
                            lower_better = lower_better, 
                            method = method)
      
    }else if(mid==0){                 #If ranking not based on mid (no ties possible)
      
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
    
  } else { # Study-specific treatment effects   ###TCu: updates 08/07/2025
    
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
        
        rk[ , , (i - 1)*ntrt + 1:ntrt] <- mid_based_ranks(temp_d, 
                                                          mid = mid, 
                                                          dim_d = dim_temp_d, 
                                                          dimnames_d = dimnames_temp_d, 
                                                          lower_better = lower_better, 
                                                          method = method)
        
      }else if(mid == 0){                 #If ranking not based on mid
        
        # Get ranks at each iteration
        rk[ , , (i - 1)*ntrt + 1:ntrt] <-
          aperm(apply(temp_d, 1:2, rank, ties.method = "average"),
                c("iterations", "chains", "parameters"))
      
        # If higher treatment effects are better   
        if (!lower_better) rk <- ntrt + 1 - rk         
      }
      
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

posterior_ranks_mid <- function(x,
                                newdata = NULL,
                                study = NULL,
                                lower_better = TRUE,
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                sucra = FALSE,
                                summary = TRUE,
                                mid = 0,
                                method = "midpoint") {
  # Checks
  if (!rlang::is_bool(lower_better))
    abort("`lower_better` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(summary))
    abort("`summary` should be TRUE or FALSE.")
  
  if (!rlang::is_bool(sucra))
    abort("`sucra` should be TRUE or FALSE.")
  
  if (!is.numeric(mid))
    abort("`mid` should be numeric.")
  
  if (!method %in% c("midpoint", "min", "max"))
    abort("parameter method must be one of: \"midpoint\", \"min\", or \"max\".")
  
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
      
      rk <- mid_based_ranks(d, mid = mid, 
                            dim_d = dim_d,
                            dimnames_d = dimnames_d, 
                            lower_better = lower_better, 
                            method = method)
      
    }else if(mid==0){                 #If ranking not based on mid (no ties possible)
      
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
    
  } else { # Study-specific treatment effects   ###TCu: updates 08/07/2025
    
      #TCu (2026-05-13): Carla to QC from here to end of function (line 191)
      #Carla, the function mid_based_ranks (line 150) has been QCed as part of the standard NMA output (i.e. where it features on line 80)
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
        
        rk[ , , (i - 1)*ntrt + 1:ntrt] <- mid_based_ranks(temp_d, 
                                                          mid = mid, 
                                                          dim_d = dim_temp_d, 
                                                          dimnames_d = dimnames_temp_d, 
                                                          lower_better = lower_better, 
                                                          method = method)
        
      }else if(mid == 0){                 #If ranking not based on mid
        
        # Get ranks at each iteration
        rk[ , , (i - 1)*ntrt + 1:ntrt] <-
          aperm(apply(temp_d, 1:2, rank, ties.method = "average"),
                c("iterations", "chains", "parameters"))
      
        # If higher treatment effects are better   
        if (!lower_better) rk <- ntrt + 1 - rk         
      }
      
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
}####mid based ranks ranks function ####


#' @param d Array of relative treatment effects as defined above, based on relative_effects

#' Calculate mid adjusted ranks
#' @noRd
mid_based_ranks <- function(x, mid, dim_d, dimnames_d, lower_better, method){
  
  if (!is.array(x) || length(dim(x)) != 3) abort("Error in mid_based_ranks: First parameter is not a 3D MCMC array, [Iterations, Chains, Parameters]")
  
  #Intermediary storage for mid ranking calculations
  test <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
  nsup <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
  ninf <- array(NA_real_, dim = dim_d, dimnames = dimnames_d)
  
  #Sequential calculations to reduce memory burden
  
  for(i in 1:dim_d[3]){    #For trt i
    
    for(j in 1:dim_d[3]){  #For trt j
      
      
      if (!lower_better){
        test[,,j] <- 1*(x[,,i] - x[,,j] - mid >= 0)     # If positive orientation, test whether trt i is mid superior to trt j
        # Test difference from mid in mean difference, assuming higher mean diff is better; Events are good => +ve xi-xj and xi-xj-delta>0, means i more effective than j
      }
      
      if (lower_better){
        test[,,j] <- 1*(x[,,j] - x[,,i] - mid >= 0)     #If negative orientation, test whether trt i is mid superior to trt j
        # Test difference from mid in mean difference, assuming lower mean diff is better; Events are bad =>  -ve xi-xj and -(xi-xj)-delta>=0 means j more effective than i
      }
      
      
    }
    
    #Calculate number of mid superiorities for trt i
    
    nsup[,,i]  <- apply(test[,,], 1:2, sum) - test[,,i]  #Subtract comparison vs itself (this final term test[,,i] is always zero unless mid=0)
    
  }
  
  #rank based on the number of (non-self) equal or inferior treatments
  rk_min <- aperm(apply((dim_d[3] - 1) - nsup, 1:2, rank, ties.method = "min"),
                  c("iterations", "chains", "parameters"))
  rk_max <- aperm(apply((dim_d[3] - 1) - nsup, 1:2, rank, ties.method = "max"),
                  c("iterations", "chains", "parameters"))
  rk_ave <- aperm(apply((dim_d[3] - 1) - nsup, 1:2, rank, ties.method = "average"),
                  c("iterations", "chains", "parameters"))
  
  if(method == "midpoint"){
    rk <- rk_ave  #take average approach
  }
  if(method == "min"){
    rk <- rk_min  #take minimum approach
  }
  if(method == "max"){
    rk <- rk_max  #take maximum approach
  }
  
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

