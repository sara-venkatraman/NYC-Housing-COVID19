Poisson regression models (Table 2 of manuscript)
================
Sara Venkatraman
1/12/2021

#### Model setup

First we run a script which reads the syndromic surveillance data, and
we also load a few libraries.

``` r
# Load script which reads the syndromic surveillance data and sets up the design matrices used for modeling
source("ModelingSetup.R")

# Packages for obtaining robust standard errors and VIFs
library(lmtest);  library(sandwich);  library(car);  library(MASS)

# Packages for spatiotemporal modeling
library(maptools);  library(spdep);  library(INLA)

# Packages for plotting and printing tables
library(ggplot2);  library(gridExtra);  library(knitr)
```

Now we read in a zip code-level NYC shapefile that will later enable us
to construct spatiotemporal models of case counts over time and over 173
zip codes.

    ## Reading layer `tl_2010_36_zcta510NYC' from data source `/Users/saravenkatraman/Documents/Cornell University (PhD)/Research/COVID-19 Demographics/Paper 1 Surveillance Analysis/Spatiotemporal Model Files/tl_2010_36_zcta510NYC.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 226 features and 12 fields
    ## geometry type:  MULTIPOLYGON
    ## dimension:      XY
    ## bbox:           xmin: 913090.7 ymin: 120053.5 xmax: 1080968 ymax: 283594.7
    ## projected CRS:  NAD83 / New York Long Island (ftUS)

The next few lines of code produce a dataframe (“design matrix”) of the
following form. Below, “Case count” refers to the suspected cases,
i.e. the total number of ILI + pneumonia emergency department
presentations observed on that day. In this dataframe, both
overcrowdedness and multigenerational housing are binned into quartiles.

<img src="https://user-images.githubusercontent.com/6864298/104348006-3c34f880-54cf-11eb-8775-310a7b2f9e7a.png" style="width:70.0%" />

``` r
# Design matrix construction
designResponse.ili <- Concatenate.Zipcode.Data(zctaOrder, "influenza", "2020-03-16", variablesToDiscretize=c("PctOvercrowded", "PctMultigen"), quartile=T)
designResponse.pneu <- Concatenate.Zipcode.Data(zctaOrder, "pneumonia", "2020-03-16", variablesToDiscretize=c("PctOvercrowded", "PctMultigen"), quartile=T)
designResponse <- cbind(designResponse.ili$Count + designResponse.pneu$Count, designResponse.ili[,-1])
colnames(designResponse)[1] <- "Count";  remove(designResponse.ili);  remove(designResponse.pneu)

# Get sum of essential employment percentages
designResponse$PctEssEmpl <- rowSums(designResponse[,18:23])
```

Define functions for neatly printing the model coefficients and
confidence intervals. The first function applies to GLMs and the second
applies to INLA models.

``` r
Print.Model.Results.GLM <- function(modelGLM, numDecimal) {
  # Get coefficient estimates and confidence intervals. Combine them (along with the)
  # p-value) into one table, called 'modelResults' - results rounded to 'numDecimal'
  modelCoef <- coeftest(modelGLM)
  modelCI <- coefci(modelGLM, vcov=vcovHC(modelGLM, type="HC3"))
  modelResults <- cbind(round(exp(modelCoef[,1]), numDecimal), 
                        modelCoef[,4], 
                        round(exp(modelCI[,1:2]), numDecimal))
  colnames(modelResults)[1:2] <- c("exp(Estimate)", "p-value")
  
  # Create a 1-column table called 'resultsSummary' that stores model results in
  # the following format: "exp(Estimate),  (ciLower, ciUpper)"
  resultsSummary <- matrix("", nrow=nrow(modelResults), ncol=1)
  rownames(resultsSummary) <- rownames(modelResults)
  for(i in 1:nrow(modelResults)) {
    string.i <- paste(modelResults[i,1], "  (", modelResults[i,3], ", ", modelResults[i,4], ")", sep="")
    resultsSummary[i,1] <- string.i
  }
  resultsSummary <- resultsSummary[c(3:nrow(resultsSummary), 2, 1), ]
  return(modelResults)
} 

Print.Model.Results.INLA <- function(modelINLA, numDecimal) {
  # Get coefficient estimates and confidence intervals. Combine them (along with the)
  # p-value) into one table, called 'modelResults' - results rounded to 'numDecimal'
  modelResults <- cbind(round(exp(model5.INLA$summary.fixed[,1]), numDecimal), 
                        model5.INLA$summary.fixed[,2], 
                        round(exp(model5.INLA$summary.fixed[,-c(1:2,4,7)]), numDecimal))
  colnames(modelResults)[1:2] <- c("exp(Estimate)", "SD")
  
  # Create a 1-column table called 'resultsSummary' that stores model results in
  # the following format: "exp(Estimate),  (ciLower, ciUpper)"
  resultsSummary <- matrix("", nrow=nrow(modelResults), ncol=1)
  rownames(resultsSummary) <- rownames(modelResults)
  for(i in 1:nrow(modelResults)) {
    string.i <- paste(modelResults[i,1], "  (", modelResults[i,3], ", ", modelResults[i,4], ")", sep="")
    resultsSummary[i,1] <- string.i
  }
  resultsSummary <- resultsSummary[c(3:nrow(resultsSummary), 2, 1), ]
  return(modelResults)
}
```

#### Model 1: Housing-related exposure covariates only (quasi-Poisson GLM)

``` r
model1.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. None need to be removed (based on VIF > 10 criterion)
kable(vif(model1.ILIpneu))
```

|                |     GVIF |  Df | GVIF^(1/(2\*Df)) |
|:---------------|---------:|----:|-----------------:|
| Time           | 1.000000 |   1 |         1.000000 |
| PctOvercrowded | 1.432021 |   3 |         1.061675 |
| PctMultigen    | 1.432021 |   3 |         1.061675 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model1.ILIpneu, 2))
```

|                  | exp(Estimate) | p-value | 2.5 % | 97.5 % |
|:-----------------|--------------:|--------:|------:|-------:|
| (Intercept)      |          0.37 |   0e+00 |  0.35 |   0.39 |
| Time             |          1.04 |   0e+00 |  1.04 |   1.05 |
| PctOvercrowdedQ2 |          1.56 |   0e+00 |  1.47 |   1.65 |
| PctOvercrowdedQ3 |          1.75 |   0e+00 |  1.66 |   1.85 |
| PctOvercrowdedQ4 |          2.05 |   0e+00 |  1.92 |   2.18 |
| PctMultigenQ2    |          1.19 |   1e-07 |  1.12 |   1.26 |
| PctMultigenQ3    |          1.30 |   0e+00 |  1.22 |   1.37 |
| PctMultigenQ4    |          1.59 |   0e+00 |  1.49 |   1.69 |

#### Model 2: Add clinical risk factors for COVID-19 to model 1

``` r
model2.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + COPD_CrudePrev + CSMOKING_CrudePrev + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. COPD has the largest VIF.
kable(vif(model2.ILIpneu))
```

|                     |      GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|----------:|----:|-----------------:|
| Time                |  1.000000 |   1 |         1.000000 |
| PctOvercrowded      |  2.858385 |   3 |         1.191297 |
| PctMultigen         |  3.513347 |   3 |         1.232973 |
| BPHIGH\_CrudePrev   | 20.428846 |   1 |         4.519828 |
| DIABETES\_CrudePrev | 25.608136 |   1 |         5.060448 |
| CHD\_CrudePrev      | 24.034236 |   1 |         4.902473 |
| OBESITY\_CrudePrev  |  8.539429 |   1 |         2.922230 |
| COPD\_CrudePrev     | 44.229821 |   1 |         6.650550 |
| CSMOKING\_CrudePrev | 25.429741 |   1 |         5.042791 |
| PctOver65           |  2.435163 |   1 |         1.560501 |

``` r
# Re-fit model without COPD
model2.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model2.ILIpneu))
```

|                     |      GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|----------:|----:|-----------------:|
| Time                |  1.000000 |   1 |         1.000000 |
| PctOvercrowded      |  2.726129 |   3 |         1.181928 |
| PctMultigen         |  3.345774 |   3 |         1.222971 |
| BPHIGH\_CrudePrev   | 11.225443 |   1 |         3.350439 |
| DIABETES\_CrudePrev |  9.730038 |   1 |         3.119301 |
| CHD\_CrudePrev      |  4.818064 |   1 |         2.195009 |
| OBESITY\_CrudePrev  |  8.356499 |   1 |         2.890761 |
| CSMOKING\_CrudePrev |  6.755660 |   1 |         2.599165 |
| PctOver65           |  2.308810 |   1 |         1.519477 |

``` r
# Re-fit model without COPD and hypertension
model2.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model2.ILIpneu))
```

|                     |     GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|---------:|----:|-----------------:|
| Time                | 1.000000 |   1 |         1.000000 |
| PctOvercrowded      | 2.315100 |   3 |         1.150169 |
| PctMultigen         | 3.501300 |   3 |         1.232267 |
| DIABETES\_CrudePrev | 6.986645 |   1 |         2.643226 |
| CHD\_CrudePrev      | 3.673233 |   1 |         1.916568 |
| OBESITY\_CrudePrev  | 3.684709 |   1 |         1.919560 |
| CSMOKING\_CrudePrev | 5.326226 |   1 |         2.307862 |
| PctOver65           | 2.335673 |   1 |         1.528291 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model2.ILIpneu, 2))
```

|                     | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
|:--------------------|--------------:|----------:|------:|-------:|
| (Intercept)         |          0.40 | 0.0000000 |  0.32 |   0.50 |
| Time                |          1.04 | 0.0000000 |  1.04 |   1.05 |
| PctOvercrowdedQ2    |          1.22 | 0.0000000 |  1.15 |   1.30 |
| PctOvercrowdedQ3    |          1.47 | 0.0000000 |  1.39 |   1.56 |
| PctOvercrowdedQ4    |          1.84 | 0.0000000 |  1.72 |   1.97 |
| PctMultigenQ2       |          1.07 | 0.0487494 |  1.00 |   1.15 |
| PctMultigenQ3       |          1.09 | 0.0321660 |  1.01 |   1.18 |
| PctMultigenQ4       |          1.02 | 0.6089763 |  0.93 |   1.12 |
| DIABETES\_CrudePrev |          1.16 | 0.0000000 |  1.14 |   1.18 |
| CHD\_CrudePrev      |          0.70 | 0.0000000 |  0.68 |   0.73 |
| OBESITY\_CrudePrev  |          1.01 | 0.0000000 |  1.01 |   1.02 |
| CSMOKING\_CrudePrev |          1.00 | 0.6588549 |  0.99 |   1.01 |
| PctOver65           |          1.01 | 0.0004611 |  1.00 |   1.02 |

#### Model 3: Add socioeconomic covariates to model 1

``` r
model3.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No variables need to be removed.
kable(vif(model3.ILIpneu))
```

|                   |     GVIF |  Df | GVIF^(1/(2\*Df)) |
|:------------------|---------:|----:|-----------------:|
| Time              | 1.000000 |   1 |         1.000000 |
| PctOvercrowded    | 2.596283 |   3 |         1.172354 |
| PctMultigen       | 4.301477 |   3 |         1.275272 |
| PctWhite          | 2.963210 |   1 |         1.721398 |
| PctBelowPovThresh | 6.727794 |   1 |         2.593799 |
| MedianIncome      | 6.218066 |   1 |         2.493605 |
| PctEssEmpl        | 2.274603 |   1 |         1.508179 |
| PctOver65         | 1.835834 |   1 |         1.354930 |
| PopDensity        | 1.620499 |   1 |         1.272988 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model3.ILIpneu, 2))
```

|                   | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
|:------------------|--------------:|----------:|------:|-------:|
| (Intercept)       |          6.79 | 0.0000000 |  4.14 |  11.16 |
| Time              |          1.04 | 0.0000000 |  1.04 |   1.05 |
| PctOvercrowdedQ2  |          1.35 | 0.0000000 |  1.27 |   1.44 |
| PctOvercrowdedQ3  |          1.43 | 0.0000000 |  1.35 |   1.52 |
| PctOvercrowdedQ4  |          1.53 | 0.0000000 |  1.42 |   1.65 |
| PctMultigenQ2     |          1.02 | 0.6406188 |  0.95 |   1.09 |
| PctMultigenQ3     |          1.08 | 0.0855993 |  1.00 |   1.16 |
| PctMultigenQ4     |          1.25 | 0.0000008 |  1.14 |   1.37 |
| PctWhite          |          1.00 | 0.0000000 |  0.99 |   1.00 |
| PctBelowPovThresh |          0.98 | 0.0000000 |  0.97 |   0.98 |
| MedianIncome      |          1.00 | 0.0000000 |  1.00 |   1.00 |
| PctEssEmpl        |          0.97 | 0.0000000 |  0.96 |   0.98 |
| PctOver65         |          0.97 | 0.0000000 |  0.96 |   0.98 |
| PopDensity        |          1.00 | 0.2833195 |  1.00 |   1.00 |

#### Model 4: Add both clinical and socioeconomic covariates to model 1

``` r
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + COPD_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. COPD has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|----------:|----:|-----------------:|
| Time                |  1.000000 |   1 |         1.000000 |
| PctOvercrowded      |  4.130818 |   3 |         1.266697 |
| PctMultigen         |  6.414250 |   3 |         1.363089 |
| BPHIGH\_CrudePrev   | 34.844689 |   1 |         5.902939 |
| DIABETES\_CrudePrev | 31.851449 |   1 |         5.643709 |
| CHD\_CrudePrev      | 28.450327 |   1 |         5.333885 |
| OBESITY\_CrudePrev  | 13.729488 |   1 |         3.705332 |
| COPD\_CrudePrev     | 45.443651 |   1 |         6.741191 |
| CSMOKING\_CrudePrev | 25.124885 |   1 |         5.012473 |
| PctWhite            | 11.137346 |   1 |         3.337266 |
| PctBelowPovThresh   | 12.061307 |   1 |         3.472939 |
| MedianIncome        |  9.125432 |   1 |         3.020833 |
| PctEssEmpl          |  2.950281 |   1 |         1.717638 |
| PctOver65           |  2.720163 |   1 |         1.649292 |
| PopDensity          |  2.280230 |   1 |         1.510043 |

``` r
# Re-fit model without COPD
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|----------:|----:|-----------------:|
| Time                |  1.000000 |   1 |         1.000000 |
| PctOvercrowded      |  4.009476 |   3 |         1.260418 |
| PctMultigen         |  6.091621 |   3 |         1.351415 |
| BPHIGH\_CrudePrev   | 22.338237 |   1 |         4.726334 |
| DIABETES\_CrudePrev | 19.193303 |   1 |         4.381016 |
| CHD\_CrudePrev      | 12.476584 |   1 |         3.532221 |
| OBESITY\_CrudePrev  | 12.915354 |   1 |         3.593794 |
| CSMOKING\_CrudePrev |  7.882084 |   1 |         2.807505 |
| PctWhite            | 10.960946 |   1 |         3.310732 |
| PctBelowPovThresh   | 11.199253 |   1 |         3.346528 |
| MedianIncome        |  9.052417 |   1 |         3.008723 |
| PctEssEmpl          |  2.884778 |   1 |         1.698463 |
| PctOver65           |  2.675371 |   1 |         1.635656 |
| PopDensity          |  2.248778 |   1 |         1.499593 |

``` r
# Re-fit model without COPD and hypertension
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Diabetes now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|----------:|----:|-----------------:|
| Time                |  1.000000 |   1 |         1.000000 |
| PctOvercrowded      |  3.170583 |   3 |         1.212057 |
| PctMultigen         |  5.625667 |   3 |         1.333611 |
| DIABETES\_CrudePrev | 19.208616 |   1 |         4.382764 |
| CHD\_CrudePrev      |  5.818530 |   1 |         2.412163 |
| OBESITY\_CrudePrev  |  4.679825 |   1 |         2.163290 |
| CSMOKING\_CrudePrev |  7.868768 |   1 |         2.805133 |
| PctWhite            |  7.514403 |   1 |         2.741241 |
| PctBelowPovThresh   |  9.907965 |   1 |         3.147692 |
| MedianIncome        |  8.934041 |   1 |         2.988987 |
| PctEssEmpl          |  2.916710 |   1 |         1.707838 |
| PctOver65           |  2.711916 |   1 |         1.646790 |
| PopDensity          |  1.958756 |   1 |         1.399555 |

``` r
# Re-fit model without COPD, hypertension, and diabetes
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model4.ILIpneu))
```

|                     |     GVIF |  Df | GVIF^(1/(2\*Df)) |
|:--------------------|---------:|----:|-----------------:|
| Time                | 1.000000 |   1 |         1.000000 |
| PctOvercrowded      | 3.202031 |   3 |         1.214053 |
| PctMultigen         | 5.501666 |   3 |         1.328666 |
| CHD\_CrudePrev      | 3.240267 |   1 |         1.800074 |
| OBESITY\_CrudePrev  | 4.492285 |   1 |         2.119501 |
| CSMOKING\_CrudePrev | 7.677644 |   1 |         2.770856 |
| PctWhite            | 3.432937 |   1 |         1.852819 |
| PctBelowPovThresh   | 9.455367 |   1 |         3.074958 |
| MedianIncome        | 8.333537 |   1 |         2.886787 |
| PctEssEmpl          | 2.634561 |   1 |         1.623133 |
| PctOver65           | 2.680648 |   1 |         1.637268 |
| PopDensity          | 2.009401 |   1 |         1.417534 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model4.ILIpneu, 2))
```

|                     | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
|:--------------------|--------------:|----------:|------:|-------:|
| (Intercept)         |         16.33 | 0.0000000 | 10.05 |  26.52 |
| Time                |          1.04 | 0.0000000 |  1.04 |   1.05 |
| PctOvercrowdedQ2    |          1.20 | 0.0000001 |  1.13 |   1.28 |
| PctOvercrowdedQ3    |          1.43 | 0.0000000 |  1.34 |   1.53 |
| PctOvercrowdedQ4    |          1.60 | 0.0000000 |  1.48 |   1.73 |
| PctMultigenQ2       |          1.18 | 0.0000191 |  1.09 |   1.27 |
| PctMultigenQ3       |          1.39 | 0.0000000 |  1.27 |   1.51 |
| PctMultigenQ4       |          1.55 | 0.0000000 |  1.41 |   1.71 |
| CHD\_CrudePrev      |          0.73 | 0.0000000 |  0.71 |   0.75 |
| OBESITY\_CrudePrev  |          1.03 | 0.0000000 |  1.03 |   1.04 |
| CSMOKING\_CrudePrev |          0.98 | 0.0022827 |  0.96 |   0.99 |
| PctWhite            |          1.00 | 0.7518419 |  1.00 |   1.00 |
| PctBelowPovThresh   |          0.98 | 0.0000000 |  0.97 |   0.98 |
| MedianIncome        |          1.00 | 0.0000000 |  1.00 |   1.00 |
| PctEssEmpl          |          0.96 | 0.0000000 |  0.96 |   0.97 |
| PctOver65           |          1.01 | 0.0001067 |  1.00 |   1.02 |
| PopDensity          |          1.00 | 0.4569164 |  1.00 |   1.00 |

#### Model 5: Bayesian spatiotemporal model, using model 4 covariates

``` r
# Add zip code ID number to design matrix (needed for spatial and temporal random effects)
zipcodeID <- sort(rep(1:length(allZipcodes), 30))
designResponse$ZipID <- zipcodeID
designResponse$ZipID2 <- zipcodeID

# Construct spatiotemporal model using same set of covariates in (reduced) model 4
model5.INLAformula <- Count ~ 1 + f(ZipID, model="bym", offset(Population/10000), graph=NYCadj) + f(ZipID2, Time, model="rw1") + Time + PctOvercrowded + PctMultigen + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + PopDensity
model5.INLA <- inla(model5.INLAformula, family="poisson", data=designResponse, control.compute=list(dic=TRUE,cpo=TRUE))

# Print model results
kable(Print.Model.Results.INLA(model5.INLA, 2))
```

|                     | exp(Estimate) |        SD | 0.025quant | 0.975quant | mode |
|:--------------------|--------------:|----------:|-----------:|-----------:|-----:|
| (Intercept)         |          8.38 | 0.8969562 |       1.42 |      48.27 | 8.47 |
| Time                |          1.04 | 0.0008182 |       1.04 |       1.04 | 1.04 |
| PctOvercrowdedQ2    |          1.19 | 0.1276568 |       0.93 |       1.54 | 1.19 |
| PctOvercrowdedQ3    |          1.17 | 0.1370162 |       0.89 |       1.53 | 1.17 |
| PctOvercrowdedQ4    |          1.28 | 0.2039914 |       0.86 |       1.91 | 1.27 |
| PctMultigenQ2       |          1.48 | 0.1633655 |       1.07 |       2.04 | 1.48 |
| PctMultigenQ3       |          1.96 | 0.1882489 |       1.35 |       2.84 | 1.97 |
| PctMultigenQ4       |          1.93 | 0.2242021 |       1.24 |       2.99 | 1.94 |
| CHD\_CrudePrev      |          0.90 | 0.0559462 |       0.80 |       1.00 | 0.90 |
| OBESITY\_CrudePrev  |          1.00 | 0.0152268 |       0.97 |       1.03 | 1.00 |
| CSMOKING\_CrudePrev |          1.08 | 0.0401657 |       1.00 |       1.17 | 1.08 |
| PctWhite            |          1.00 | 0.0035215 |       0.99 |       1.00 | 1.00 |
| PctBelowPovThresh   |          0.96 | 0.0134162 |       0.93 |       0.98 | 0.96 |
| MedianIncome        |          1.00 | 0.0000028 |       1.00 |       1.00 | 1.00 |
| PctEssEmpl          |          0.98 | 0.0145039 |       0.95 |       1.01 | 0.98 |
| PctOver65           |          1.02 | 0.0194660 |       0.98 |       1.06 | 1.02 |
| PopDensity          |          1.00 | 0.0000018 |       1.00 |       1.00 | 1.00 |
