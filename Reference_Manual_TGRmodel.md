<!-- toc -->

May 18, 2022

# DESCRIPTION

```
Package: TGRmodel
Title: Functions to predict Tumor Growth Rate inhibition using in vitro data
Version: 0.0.1.0
Authors@R: c(
    person("Rocky", "Diegmiller", , "therock5493@gmail.com", role = c("aut", "cre"))
    person("Marc", "Hafner", , "hafner.marc@gene.com", role = c("aut", "cre"), comment = c(ORCID = "0000-0003-1337-7598")))
Description: Functions to predict Tumor Growth Rate inhibition using in vitro data. 
    Model described in the manuscript "Growth-rate model predicts in vivo tumor response from in vitro data" (2022, CPT: Pharmacometrics & Systems Pharmacology).
    Please cite "Diegmiller et al., in revision at CPT Pharmacometrics Syst Pharmacol, 2022"
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.0
```


# `.integration_step`

step for the numerical integration


## Description

step for the numerical integration


## Usage

```r
.integration_step()
```


# `GR_inVitro_integration`

GR_inVitro_integration
 Calculate relative tumor growth inhibition (TGR) based on conc_profile


## Description

Calculate the serum concentration profile based on the drug PK, dose, and schedule.
 The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e
 provided in the input variable PK_para.


## Usage

```r
GR_inVitro_integration(conc_profile, GR_para, int_method = "mean")
```


## Arguments

Argument      |Description
------------- |----------------
`conc_profile`     |     data.frame with `Time` (in days) and relative growth rate `k` over time can be the output of `PK_to_conc_profile`
`GR_para`     |     numeric array for the in vitro GR parameters of the drug (required fields: GR_inf, GEC50, h_GR) Can calculate the average of multiple parameters at the same time given multiple rows in `GR_para`
`int_method`     |     string to define the type of integration


## Value

numeric value for relative tumor growth rate (TGR)


## Examples

```r
## regular QD treatment

PK_para = data.frame(k_a = 2, k_e = 0.25, MW = 450, VF = 1)
GR_para = data.frame(GR_inf = -0.5, GEC50 = 0.3, h_GR = 1.5)
Schedule = "QD"
Duration = 21
Dose = 5
conc_profile = PK_to_conc_profile(PK_para, Dose, Schedule, Duration)
relk_1 = relk_over_time(GR_para, conc_profile)
GR_predicted = GR_inVitro_integration(conc_profile, GR_para)
```


# `logistic_4parameters`

logistic function for fitting drug-dose response curve


## Description

`logistic_4parameters` returns values based on concentration and fit parameters


## Usage

```r
logistic_4parameters(c, Vinf, V0, EC50, h)
```


## Arguments

Argument      |Description
------------- |----------------
`c`     |     concentration (can be an array)
`Vinf`     |     asymptotic value at high concentration
`V0`     |     asymptotic value at low concentration
`EC50`     |     mid-point of the curve
`h`     |     Hill coefficient


## Details

returns values based on concentration and fit parameters


## Value

array of response values


# `logisticFit`

Actual fitting function


## Description

`logisticFit` returns fit parameters for an IC or GR curve


## Usage

```r
logisticFit(
  concs,
  normValues,
  x_0 = 1,
  curve_type = c("IC", "GR"),
  force = FALSE,
  cap = 0.1
)
```


## Arguments

Argument      |Description
------------- |----------------
`concs`     |     concentration values
`normValues`     |     normalized response values
`x_0`     |     upper limit; (=1 by default)
`curve_type`     |     response curve: either IC ( [0,1](#0,1) ) or GR( [-1,1](#-1,1) )
`force`     |     force a sigmoidal fit even if the fit is not significantly better than a flat fit
`cap`     |     cap values at (x_0 + cap)


## Details

returns fit parameters


## Value

vector of parameters


# `PK_to_conc_profile`

PK_to_conc_profile
 Calculate the serum concentration profile based on the drug PK, dose, and schedule


## Description

Calculate the serum concentration profile based on the drug PK, dose, and schedule.
 The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e
 provided in the input variable PK_para.


## Usage

```r
PK_to_conc_profile(PK_para, Dose, Schedule, Duration, Dose_uM_0 = 0)
```


## Arguments

Argument      |Description
------------- |----------------
`PK_para`     |     numeric array for the PK properties of the drug (required fields: MW, VF, k_a, k_e)  

*  'MW' molecular weight in g/mol  

*  'VF' distribution volume in L/kg  

*  'k_a' absoption rate in 1/h  

*  'k_e' elimination rate in 1/h 
`Dose`     |     numeric values for the dose of the drug given in mg/kg
`Schedule`     |     string for the frequency of the treatment given (QC, BID or EOD)
`Duration`     |     numeric value for the duration of the treatment (in days)
`Dose_uM_0`     |     numeric value for the initial serum concentration (in µM) Defaults to `0` .


## Value

list with `Time` (in days) and `Conc` serum concentration (in µM)
 and `dose_period` for the dosing schedule


# `relk_fct`

expanding GR function for rate based on concentration


## Description

expanding GR function for rate based on concentration


## Usage

```r
relk_fct(c, GR_para)
```


# `relk_over_time`

relk_over_time
 Calculate the serum concentration profile based on the drug PK, dose, and schedule


## Description

Calculate the serum concentration profile based on the drug PK, dose, and schedule.
 The concentration (in µM) is calculated based on one-compartment PK model with parameters k_a and k_e
 provided in the input variable PK_para.


## Usage

```r
relk_over_time(GR_para, conc_profile)
```


## Arguments

Argument      |Description
------------- |----------------
`GR_para`     |     numeric array for the in vitro GR parameters of the drug (required fields: `GR_inf` , `GEC50` , `h_GR` )
`conc_profile`     |     data.frame with `Time` (in days) and `Conc` (serum concentration in µM) can be the output of PK_to_conc_profile


## Value

data.frame with `Time` (in days) and relative growth rate `k` over time


