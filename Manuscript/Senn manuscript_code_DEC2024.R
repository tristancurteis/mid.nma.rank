#manuscript results

#Clear workspace
rm(list = ls())

library(devtools)
#install_github("tristancurteis/mid.nma.rank")
#devtools::install_github("dmphillippo/multinma")
#install.packages("multinma", repos = c("https://dmphillippo.r-universe.dev", getOption("repos")))

library(ggplot2)
library(multinma)
library(mid.rank.nma)
library(rprojroot)
library(dplyr)
library(openxlsx)

# Plotting function

#' Treatment rankings and rank probabilities
#'
#' The below code allows for the correct production of plots of MID-adjusted probability jth and cumulative probability jth
#' Refer to manuscript code for example usage.
#' plot.nma_rank_probs below is to replace the equivalent version in multinma, if MID-adjusted plots are required.

#' @rdname plot.nma_summary
#' @export
plot.nma_rank_probs <- function(x, ...) {
  # Get axis labels from attributes, if available
  p_xlab <- attr(x, "xlab", TRUE)
  if (is.null(p_xlab)) p_xlab <- ""

  p_ylab <- attr(x, "ylab", TRUE)
  if (is.null(p_ylab)) p_ylab <- ""

  dat <- as.data.frame(x)

  ntrt <- nrow(dat)

  if (has_studies <- rlang::has_name(dat, ".study")) {
    dat$Study <- forcats::fct_inorder(factor(dat$.study))
    dat$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(dat$parameter, "(?<=\\: ).+(?=\\])")))
  } else {
    dat$Treatment <- forcats::fct_inorder(factor(
      stringr::str_extract(dat$parameter, "(?<=\\[).+(?=\\])")))
  }

  dat <- tidyr::pivot_longer(dat, cols = dplyr::starts_with("p_rank"),
                             names_to = "rank", values_to = "probability")
  #names_pattern = "^p_rank\\[?\\d+([.]5)?\\]$", #Update from multinma is here
  #names_transform = list(rank = as.numeric))    #Update from multinma is here

  dat$rank <- gsub("p_rank\\[","", dat$rank)   #Update from multinma is here
  dat$rank <- gsub("\\]","", dat$rank)         #Update from multinma is here
  dat$rank <- as.numeric(dat$rank)             #Update from multinma is here

  p <- ggplot2::ggplot(dat,
                       ggplot2::aes(x = .data$rank, y = .data$probability)) +
    ggplot2::geom_line(...) +
    ggplot2::ylab(p_ylab) +
    ggplot2::scale_x_continuous(p_xlab, breaks = 1:ntrt, minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_multinma()

  if (has_studies) {
    p <- p + ggplot2::facet_grid(Study~Treatment)
  } else {
    p <- p + ggplot2::facet_wrap(~Treatment)
  }

  return(p)
}

#load data
data(Senn2013, package = "mid.rank.nma")
head(Senn2013)

Senn2013$se[Senn2013$study == "Willms1999"][2] <- 0.4 #Impute standard error for reference arm

options(mc.cores = parallel::detectCores())


#======================#
#####  Directories #####
#======================#

# Directory where the .RProj file is
w_dir <- find_root(has_file("Manuscript.Rproj"))

res_dir <- paste0(w_dir, "/Results/Senn 2013/")
res_dir_corrections <- paste0(w_dir, "/Results/Corrections/Senn 2013/")


#===============================#
#####  Fit Model            #####
#===============================#

diabetes_net <- set_agd_contrast(Senn2013,
                                 study = study,
                                 trt = trt,
                                 y = y,
                                 se = se,
                                 trt_ref = "Placebo")

diabetes_fit_FE <- nma(diabetes_net,
                       trt_effects = "fixed",
                       prior_intercept = normal(scale = 1000),
                       prior_trt = normal(scale = 1000),
                       iter = 50000,   #Update to 50000
                       chains = 4)

diabetes_fit_RE <- nma(diabetes_net,
                       trt_effects = "random",
                       prior_intercept = normal(scale = 1000),
                       prior_trt = normal(scale = 1000),
                       iter = 50000,  #Update to 50000
                       chains = 4)

#===================================================#
#####  Model Selection and Densities            #####
#===================================================#

dic(diabetes_fit_RE)
dic(diabetes_fit_FE)

plot_prior_posterior(diabetes_fit_RE, prior = c("trt", "het"))


#===============================#
#####        Network         #####
#===============================#

net <- plot(diabetes_net)
net

#========================================#
#####     Treatment Effects          #####
#========================================#

rel_eff <- relative_effects(diabetes_fit_RE, trt_ref = "Placebo")
# rel_eff$summary <- rel_eff$summary[order(rel_eff$summary$mean),]
# rel_eff$sims <- rel_eff$sims[,,order(rel_eff$summary$mean)]

forest <- plot(rel_eff, ref_line = 0) + labs(x = "Mean Difference Vs Placebo, HbA1c (%)")
forest

#===============================#
#####     Ranking           #####
#===============================#

# Without minimally important difference, # With ties.method = minimum
rank <- posterior_ranks(diabetes_fit_RE, summary = TRUE, sucra = TRUE)
rank_probs <- posterior_rank_probs(diabetes_fit_RE)
rank_probs_cum <- posterior_rank_probs(diabetes_fit_RE, cumulative = T)

rank_0 <- plot(rank)
probs_0 <- plot(rank_probs)
cum_0 <- plot(rank_probs_cum)

rank_0

# With minimally important difference
mid <- 0.3 # CMi updated to new value

rank_mid <- posterior_ranks_mid(diabetes_fit_RE, summary = TRUE, sucra = TRUE, mid = mid)
rank_probs_mid <- posterior_rank_probs_mid(diabetes_fit_RE, mid = mid)
rank_probs_cum_mid <- posterior_rank_probs_mid(diabetes_fit_RE, cumulative = T, mid = mid)

rank_mid_plot <- plot(rank_mid)
probs_mid <- plot(rank_probs_mid)
cum_mid <- plot(rank_probs_cum_mid)

rank_mid_plot

p_mid_best <- rank_probs_mid$summary$p_best_or_equal #CMi updated from  rank_probs_mid$summary$p_mid_best as column didn't exist

# Plots
rank_mid_plot <- plot(rank_mid)
probs_mid <- plot(rank_probs_mid)
cum_mid <- plot(rank_probs_cum_mid)

rank_mid_plot
probs_mid
cum_mid

#===============================#
#####         Table         #####
#===============================#

# round up numbers to specified number of digits
digits <- function(x, dig = 1, adj = 0) {
  # x: numeric value (or converted to numeric if string)
  # dig: number of digits for output, default 1 but can adjust
  # adj: add small number (e.g. adj = 0.000001) if you want to round up (like excel) - default 0 (round to nearest even)
  sprintf(paste0("%." , dig, "f"), as.numeric(x)+adj)
}

rank_table <- cbind(rank$summary$`50%`,
                    rank_mid$summary$`50%`,
                    digits(rank$summary$sucra, 2),
                    digits(rank_mid$summary$sucra, 2),
                    digits(rank_probs$summary$`p_rank[1]`,2),
                    digits(rank_probs_mid$summary$`p_rank[1]`,2),
                    digits(p_mid_best,2))
rank_table <- as.data.frame(rank_table)
colnames(rank_table) <- c("Median Rank", "MID-Median Rank",
                          "SUCRA", "MID-SUCRA",
                          "P Best", "P MID Best", "P MID Equal Best")
rownames(rank_table) <- rank$summary$.trt
rank_table

#===============================#
#####        Export         #####
#===============================#

#Export forest plot
forest

#Export forest plot
loc <- paste0(res_dir_corrections, "Senn 2013 Forest.png")
png(loc)
print(forest)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 Forest.pdf")
pdf(loc)
print(forest)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 Forest.jpg")
jpeg(loc)
print(forest)
dev.off()

#Export network plot
loc <- paste0(res_dir_corrections, "Senn 2013 Network.png")
png(loc)
print(net)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 Network.pdf")
pdf(loc)
print(net)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 Network.jpg")
jpeg(loc)
print(net)
dev.off()

#export ranking plots
loc <- paste0(res_dir_corrections, "Senn 2013 Median Rank MID = 0.pdf")
pdf(loc)
print(rank_0)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 Median Rank MID.pdf")
pdf(loc)
print(rank_mid_plot)
dev.off()

#export probs plots
loc <- paste0(res_dir_corrections, "Senn 2013 Probs MID.pdf")
pdf(loc)
print(probs_mid)
dev.off()

loc <- paste0(res_dir_corrections, "Senn 2013 CUM MID.pdf")
pdf(loc)
print(cum_mid)
dev.off()

# loc <- paste0(res_dir, "Senn 2013 Median Prob MID = 0.pdf")
# pdf(loc)
# print(probs_0)
# dev.off()
#
# loc <- paste0(res_dir, "Senn 2013 Median Prob MID.pdf")
# pdf(loc)
# print(probs_mid)
# dev.off()
#
# loc <- paste0(res_dir, "Senn 2013 Median Cum MID = 0.pdf")
# pdf(loc)
# print(cum_0)
# dev.off()
#
# loc <- paste0(res_dir, "Senn 2013 Median Cum MID.pdf")
# pdf(loc)
# print(cum_mid)
# dev.off()

#export ranking table

filename <- "Ranking Table"                                                                 #File name
exp_file <- paste0(res_dir_corrections, filename, ".xlsx")    #File path
write.xlsx(rank_table, exp_file, overwrite = T, rowNames=T)


