
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

