
#' Actual fitting function
#'
#' \code{logisticFit} returns fit parameters for an IC or GR curve
#'
#' returns fit parameters
#'
#' @param concs concentration values
#' @param normValues normalized response values
#' @param x_0 upper limit; (=1 by default)
#' @param curve_type response curve: either IC (for values in \code{[0,1]}) or GR (for values in \code{[-1,1]})
#' @param force force a sigmoidal fit even if the fit is not significantly better than a flat fit
#' @param cap cap values at (x_0 + cap)
#' @return vector of parameters
#' @importFrom drc drm drmc LL.3u
#' @export
logisticFit <-
  function(concs,
           normValues,
           x_0 = 1,
           curve_type = c("IC", "GR"),
           force = FALSE,
           cap = 0.1) {
    
    # define variables and prepare data
    log10concs <- log10(concs)
    df_ <- data.frame(log10conc = log10concs,
                     normValues = pmin(normValues, x_0 + cap))

    fit_para <- c("h", "x_inf", "c50")

    out <- array(NA, length(get_header("response_metrics")))
    names(out) <- get_header("response_metrics")
    out["maxlog10Concentration"] <- max(log10concs)
    out["N_conc"] <- length(unique(log10concs))
    out["x_0"] <- x_0

    # fit parameters and boundaries
    if (curve_type == "IC") {
      priors <- c(2, 0.4, median(concs))
      lower <- c(.1, 0, min(concs) / 10)
    } else if (curve_type == "GR") {
      priors <- c(2, 0.1, median(concs))
      lower <- c(.1, -1, min(concs) / 10)
    }
    upper <- c(5, min(x_0 + .1, 1), max(concs) * 10)

    controls <- drc::drmc()
    controls$relTol <- 1e-06
    controls$errorm <- FALSE
    controls$noMessage <- TRUE
    controls$rmNA <- TRUE

    ######################################
    # IC curve fitting
    output_model_new <- try(drc::drm(
      normValues ~ log10conc,
      data = df_,
      logDose = 10,
      fct = drc::LL.3u(upper = x_0, names = fit_para),
      start = priors,
      lowerl = lower,
      upperl = upper,
      control = controls,
      na.action = na.omit
    ))

    # assuming proper fit result
    if (class(output_model_new) != "try-error") {
      for (p in fit_para) {
        out[p] <- stats::coef(output_model_new)[paste0(p, ":(Intercept)")]
      }
      # F-test for the significance of the sigmoidal fit
      Npara <- 3 # N of parameters in the growth curve
      Npara_flat <- 1 # F-test for the models
      RSS2 <- sum(stats::residuals(output_model_new) ^ 2, na.rm = TRUE)
      RSS1 <- sum((df_$normValues - mean(df_$normValues,
                                        na.rm = TRUE)) ^ 2, na.rm = TRUE)
      df1 <- (Npara - Npara_flat)
      df2 <- (length(na.omit(df_$normValues)) - Npara + 1)
      f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
      f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
      out["r2"] <- 1 - RSS2 / RSS1
    }


    # non-fitted metrics
    xAvg <- aggregate(
      df_$normValues,
      by = list(log10conc = df_$log10conc),
      FUN = function(x)
        mean(x, na.rm = T)
    )
    colnames(xAvg)[2] <- "normValues"
    l <- dim(xAvg)[1]

    out["x_max"] <- min(xAvg$normValues[c(l, l - 1)], na.rm = TRUE)

    out["x_mean"] <- mean(xAvg$normValues)
    out["x_AOC"] <- 1 - mean(xAvg$normValues)

    # analytical solution for ic50
    out["xc50"] <- out["c50"] * ((x_0 - out["x_inf"]) / (0.5 - out["x_inf"]) - 1) ^
      (1 / out["h"])

    # testing the significance of the fit and replacing with flat function if required
    pcutoff <- ifelse(force, 1, .05)
    if (!is.na(f_pval)) {
      out["flat_fit"] <- ifelse(f_pval >= pcutoff |
                                 is.na(out["c50"]), 1, 0)
    } else {
      out["flat_fit"] <- ifelse(is.na(out["c50"]), 1, 0)
    }

    # Replace values for flat fits: c50 = 0, h = 0.01 and xc50 = +/- Inf
    if (out["flat_fit"] == 1) {
      out["c50"] <- 0
      out["h"] <- 0.0001
      out["xc50"] <- ifelse(mean(xAvg$normValues) > .5, Inf, -Inf)
      out["x_inf"] <- mean(xAvg$normValues)
    }

    # Add xc50 = +/-Inf for any curves that don"t reach RelativeViability = 0.5
    if (is.na(out["xc50"])) {
      out["xc50"] <- ifelse(out["x_inf"] > .5, Inf, -Inf)
    }
    return(out)
  }



#' logistic function for fitting drug-dose response curve
#'
#' \code{logistic_4parameters} returns values based on concentration and fit parameters
#'
#' returns values based on concentration and fit parameters
#' 
#' @param c concentration (can be an array)
#' @param Vinf asymptotic value at high concentration
#' @param V0 asymptotic value at low concentration
#' @param EC50 mid-point of the curve
#' @param h Hill coefficient
#' @return array of response values
#' @export
logistic_4parameters <- function(c, Vinf, V0, EC50, h) {
  Vinf + (V0 - Vinf) / (1 + (c / EC50) ^ h)
}
