
#' PK_to_conc_profile
#' Calculate the serum concentration profile based on the drug PK, dose, and schedule
#'
#' Calculate the serum concentration profile based on the drug PK, dose, and schedule.
#' The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e 
#' provided in the input variable PK_para. 
#'
#' @param PK_para numeric array for the PK properties of the drug (required fields: MW, VF, k_a, k_e)
#' \itemize{
#'   \item{'MW' molecular weight in g/mol}
#'   \item{'VF' distribution volume in L/kg}
#'   \item{'k_a' absoption rate in 1/h}
#'   \item{'k_e' elimination rate in 1/h}
#'   }
#' @param Dose  numeric values for the dose of the drug given in mg/kg
#' @param Schedule  string for the frequency of the treatment given (QC, BID or EOD)
#' @param Duration  numeric value for the duration of the treatment (in days)
#' @param Dose_uM_0 numeric value for the initial serum concentration (in µM)
#'  Defaults to \code{0}.
#'
#' @return data.frame with \code{Time} (in days) and \code{Conc} serum concentration (in µM)
#' 
#' @export
#' 

PK_to_conc_profile = function(PK_para, Dose, Schedule, Duration, Dose_uM_0 = 0) {
  # Calculate the PK serum profile based on the drug, dose, and schedule
  #PK_para are the properties of the drug kept in a file
  #Dose is the drug conc. given for each in vivo trial
  #Schedule is the frequency of injections supplied
  #Duration is the last day of the trial
  Dose_uM = Dose*1000/PK_para$MW/PK_para$VF
  if (Schedule == "QD"){
    T <- 1
  } else if (Schedule == "BID"){
    T <- 0.5
  } else if (Schedule == "EOD"){
    T <- 2
  } else if (Schedule == "Q3D"){
    T <- 3
  }
  # do the integration based on a 1 compartment PK model
  k_a <- PK_para$k_a #check if given in hr^-1
  k_e <- PK_para$k_e #check if given in hr^-1
  
  #solve explicitly by creating matrix, finding e-val, e-vec, and
  #solving linear system of ODEs
  A <- 24*matrix(c(-k_a, 0, k_a, -k_e), 2, 2, byrow=TRUE)
  ev <- eigen(A)
  eigval <- ev$values
  eigvec <- ev$vectors
  initvec <- matrix(c(Dose_uM, Dose_uM_0), 2, 1, byrow=TRUE)
  consts <- solve(eigvec,initvec)
  
  #create pulsetrain of multiple doses given ever T days
  
  conc <- function(t) (consts[1]*eigvec[2,1]*exp(eigval[1]*(t)) + consts[2]*eigvec[2,2]*exp(eigval[2]*(t)))
  tp <- seq(0,Duration,.integration_step())
  pp <- conc(tp)
  d <- seq(0,Duration,T)
  
  serum_conc <- function(t) pulstran(t,d,pp,1/.integration_step())
  return(list(Time = tp, Conc = serum_conc))
}


#' expanding GR function for rate based on concentration
relk_fct = function(c,GR_para) {
  GR_c = logistic_4parameters(pmax(0,c), GR_para$GR_inf, 1, GR_para$GEC50, GR_para$h_GR)
  return ( log2(GR_c + 1) )
}


#' relk_over_time
#' Calculate the serum concentration profile based on the drug PK, dose, and schedule
#'
#' Calculate the serum concentration profile based on the drug PK, dose, and schedule.
#' The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e 
#' provided in the input variable PK_para. 
#'
#' @param GR_para  numeric array for the in vitro GR parameters of the drug (required fields: GR_inf, GEC50, h_GR)
#' Can calculate the average of multiple parameters at the same time given multiple rows in \code{GR_para}
#' @param conc_profile  data.frame with \code{Time} (in days) and \code{Conc} (serum concentration in µM)
#' can be the output of PK_to_conc_profile
#'
#' @return data.frame with\code{Time} (in days) and relative growth rate \code{k} over time
#' 
#' @export
#' 

relk_over_time = function(GR_para, conc_profile) {
  k_inhibition = data.frame(Time = conc_profile$Time, relk = 0)
  
  if (nrow(GR_para)>1) {
    temp_relk = matrix(0, nrow(conc_profile), nrow(GR_para))
    for (i in 1:nrow(GR_para)) {
      temp_relk[,i] = relk_fct(conc_profile$Conc(conc_profile$Time), GR_para[i,])
    }
    k_inhibition$relk = rowMeans(temp_relk)
  } else {
    k_inhibition$relk = relk_fct(conc_profile$Conc(conc_profile$Time), GR_para)
  }
  
  return(k_inhibition)
}


#' GR_inVitro_integration
#' Calculate relative tumor growth inhibition (TGR) based on conc_profile
#'
#' Calculate the serum concentration profile based on the drug PK, dose, and schedule.
#' The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e 
#' provided in the input variable PK_para. 
#'
#' @param conc_profile  data.frame with \code{Time} (in days) and relative growth rate \code{k} over time
#' can be the output of \code{PK_to_conc_profile}
#' @param GR_para  numeric array for the in vitro GR parameters of the drug (required fields: GR_inf, GEC50, h_GR)
#' Can calculate the average of multiple parameters at the same time given multiple rows in \code{GR_para}
#' @param int_method string to define the type of integration
#'
#' @return numeric value for relative tumor growth rate (TGR)
#' 
#' @export
#' 
#' @examples
#' 
#' ## regular QD treatment
#' 
#' PK_para = data.frame(k_a = 2, k_e = 0.25, MW = 450, VF = 1)
#' GR_para = data.frame(GR_inf = -0.5, GEC50 = 0.3, h_GR = 1.5)
#' Schedule = "QD"
#' Duration = 21
#' Dose = 5
#' conc_profile = PK_to_conc_profile(PK_para, Dose, Schedule, Duration)
#' relk_1 = relk_over_time(GR_para, conc_profile)
#' GR_predicted = GR_inVitro_integration(conc_profile, GR_para)

GR_inVitro_integration = function(conc_profile, GR_para, int_method = 'mean') {
  
  Duration = max(conc_profile$Time)
  # integration based on different methods (parameter  int_method  )
  if (int_method == 'mean' ) {
    k_inhibition = relk_over_time(GR_para, conc_profile)
    k_mean = mean( k_inhibition$relk[ k_inhibition$Time >= (Duration-2)*T &
                                        k_inhibition$Time < (Duration-1)*T ])
    
  } else if (int_method == 'cont_int' ) {
    integrand = function(x) {
      # integration should be over time --> convert in c then in GR then in relk
      if (nrow(GR_para)>1) {
        GR_c = matrix(0, length(x), nrow(GR_para))
        for (i in 1:nrow(GR_para)) {
          GR_c[,i] = relk_fct(conc_profile$Conc(x), GR_para[i,])
        }
        return(log2( rowMeans((2**GR_c) - 1) + 1))
      } else {
        GR_c = relk_fct(conc_profile$Conc(x), GR_para)
        return ( log2(GR_c + 1) )
      }
    }
    int_val = stats::integrate(integrand,
                        lower = (Duration-2)*T,
                        upper = (Duration-1)*T,
                        subdivisions=20*T/.integration_step())
    k_mean = int_val$value/T
  }
  
  # returns a single GR value for the integration over the Duration of the experiment
  GR_mean = 2 ** k_mean - 1
  return( GR_mean )
}

#' step for the numerical integration
.integration_step = function() {
  return( 0.0005 )
}
