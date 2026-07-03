  # Charité - Universitätsmedizin Berlin
  # Institute of Public Health
  # Hans-Aloys Wischmann
  # July 2nd, 2026

  # ensure consistency across systems
  Sys.setlocale("LC_ALL", 'en_US.UTF-8')
  Sys.setenv(LANG = "en_US.UTF-8")
  knitr::opts_chunk$set(echo = FALSE, fig.width = 6, dpi = 1200, comment = NA)
  knitr::opts_knit$set(root.dir = getwd())

  # load required libraries
  library(tidyverse)
  library(flextable)
  library(data.table)
  library(marginaleffects)
  library(readxl)
  library(table1)
  library(broom)
  library(lavaan)
  library(MoEClust)
  library(ggh4x)
  library(pROC)
  library(mgcv)
  library(car)

  # utility to format a p.value
  format.p <- function(p) {
    case_when(p < 0.001 ~ "<0.001",
              p < 0.10  ~ sprintf("%.03f", p),
              TRUE      ~ sprintf("%.02f", p))
  }

  # utilities to compute Weighted mean, weighed sd, and weighted p.value from weighted paired t.test
  wmean  <- function(x, w) { return(sum(w * x) / sum(w)) }
  wvar   <- function(x, w) { wm <- wmean(x, w); return(sum(w * (x - wm) * (x - wm)) / (sum(w) - 1)) }
  wsd    <- function(x, w) { return(sqrt(wvar(x, w))) }
  wttest <- function(x, y, w) {
    # Weighted means and variances
    mx <- wmean(x, w)
    my <- wmean(y, w)
    vx <- wvar(x, w)
    vy <- wvar(y, w)
  
    # Weighted t-statistic, degrees of freedom, p.value
    t_stat <- (mx - my) / sqrt(vx / length(x) + vy / length(y))
    dof <- (vx / length(x) + vy / length(y))^2 / ((vx / length(x))^2 / (length(x) - 1) + (vy / length(y))^2 / (length(y) - 1))
    return(2.0 * pt(-abs(t_stat), dof))
  }

  # utility function to plot to *.pdf file and *.png file, set default size to 7.25 in x 5.25 in
  plotPngPdf <- function(file, object, width = 7.25, height = 5.25) {
    themed <- object + theme(text = element_text(size = 10), plot.title = element_text(size = 10))
    png(paste(file, ".png", sep = ""), width, height, units = "in", res = 1200)
    print(themed)
    invisible(capture.output(dev.off()))
    pdf(paste(file, ".pdf", sep = ""), width, height, paper = "a4")
    print(themed)
    invisible(capture.output(dev.off()))
    knitr::include_graphics(paste(file, ".png", sep = ""), dpi = 1200)
  }

  # default theme for flextables
  theme_flextable <- function(x, ..., j = NULL, align = NULL) {
    if (!inherits(x, "flextable")) {
      stop("Function theme_netkoh() supports only flextable objects.")
    }

    defaults <- get_flextable_defaults()

    big_border <- officer::fp_border(
      width = defaults$border.width * 2,
      color = defaults$border.color
    )
    std_border <- officer::fp_border(
      width = defaults$border.width,
      color = defaults$border.color
    )
  
    h_nrow <- nrow_part(x, "header")
    f_nrow <- nrow_part(x, "footer")
    b_nrow <- nrow_part(x, "body")
  
    x <- border_remove(x)
    x <- hline_top(x,    border = big_border, part = "header")
    x <- hline_bottom(x, border = big_border, part = "header")
    x <- hline_bottom(x, border = big_border, part = "body")
    x <- hline_bottom(x, border = big_border, part = "footer")
    fix_border_issues(x)

    x <- bold(x, bold = TRUE, part = "header")

    x <- flextable::font(x, fontname = "Calibri", part = "all")
    x <- fontsize(x, size = 10, part = "header")
    x <- fontsize(x, size = 10, part = "body")
    x <- fontsize(x, size =  9, part = "footer")

    x <- padding(x, padding.top = 1, padding.bottom = 2, part = "header")
    x <- padding(x, padding.top = 2, padding.bottom = 2, part = "body")
    x <- padding(x, padding.top = 0, padding.bottom = 1, part = "footer")

    x <- align(x, align = "center", part = "all")
    x <- align(x, j = 1, align = "left", part = "all")
    x <- align(x, align = "left", part = "footer")
    if(!is.null(j) && !is.null(align)) {
      x <- align(x, j = j, align = align, part = "header")
      x <- align(x, j = j, align = align, part = "body")
    }

    x <- italic(x, i = grepl("Fehlende Angaben", x$body$dataset[,1], fixed = TRUE), j = 1)
    x <- autofit(x)
  }
  set_flextable_defaults(theme_fun = theme_flextable, font.family = "Calibri", font.size = 10)
