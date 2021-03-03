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

![](https://user-images.githubusercontent.com/6864298/104348006-3c34f880-54cf-11eb-8775-310a7b2f9e7a.png)

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

|                |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------- | -------: | -: | ---------------: |
| Time           | 1.000000 |  1 |         1.000000 |
| PctOvercrowded | 1.432021 |  3 |         1.061675 |
| PctMultigen    | 1.432021 |  3 |         1.061675 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model1.ILIpneu, 2))
```

|                  | exp(Estimate) | p-value | 2.5 % | 97.5 % |
| :--------------- | ------------: | ------: | ----: | -----: |
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

|                     |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | --------: | -: | ---------------: |
| Time                |  1.000000 |  1 |         1.000000 |
| PctOvercrowded      |  2.858385 |  3 |         1.191297 |
| PctMultigen         |  3.513347 |  3 |         1.232973 |
| BPHIGH\_CrudePrev   | 20.428846 |  1 |         4.519828 |
| DIABETES\_CrudePrev | 25.608136 |  1 |         5.060448 |
| CHD\_CrudePrev      | 24.034236 |  1 |         4.902473 |
| OBESITY\_CrudePrev  |  8.539429 |  1 |         2.922230 |
| COPD\_CrudePrev     | 44.229821 |  1 |         6.650550 |
| CSMOKING\_CrudePrev | 25.429741 |  1 |         5.042791 |
| PctOver65           |  2.435163 |  1 |         1.560501 |

``` r
# Re-fit model without COPD
model2.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model2.ILIpneu))
```

|                     |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | --------: | -: | ---------------: |
| Time                |  1.000000 |  1 |         1.000000 |
| PctOvercrowded      |  2.726129 |  3 |         1.181928 |
| PctMultigen         |  3.345774 |  3 |         1.222971 |
| BPHIGH\_CrudePrev   | 11.225443 |  1 |         3.350439 |
| DIABETES\_CrudePrev |  9.730038 |  1 |         3.119301 |
| CHD\_CrudePrev      |  4.818064 |  1 |         2.195009 |
| OBESITY\_CrudePrev  |  8.356499 |  1 |         2.890761 |
| CSMOKING\_CrudePrev |  6.755660 |  1 |         2.599165 |
| PctOver65           |  2.308810 |  1 |         1.519477 |

``` r
# Re-fit model without COPD and hypertension
model2.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model2.ILIpneu))
```

|                     |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | -------: | -: | ---------------: |
| Time                | 1.000000 |  1 |         1.000000 |
| PctOvercrowded      | 2.315100 |  3 |         1.150169 |
| PctMultigen         | 3.501300 |  3 |         1.232267 |
| DIABETES\_CrudePrev | 6.986645 |  1 |         2.643226 |
| CHD\_CrudePrev      | 3.673233 |  1 |         1.916568 |
| OBESITY\_CrudePrev  | 3.684709 |  1 |         1.919560 |
| CSMOKING\_CrudePrev | 5.326226 |  1 |         2.307862 |
| PctOver65           | 2.335673 |  1 |         1.528291 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model2.ILIpneu, 2))
```

|                     | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
| :------------------ | ------------: | --------: | ----: | -----: |
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
model3.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No variables need to be removed.
kable(vif(model3.ILIpneu))
```

|                   |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------- | -------: | -: | ---------------: |
| Time              | 1.000000 |  1 |         1.000000 |
| PctOvercrowded    | 2.421367 |  3 |         1.158804 |
| PctMultigen       | 3.934837 |  3 |         1.256477 |
| PctWhite          | 2.824732 |  1 |         1.680694 |
| PctBelowPovThresh | 6.357139 |  1 |         2.521337 |
| MedianIncome      | 6.077919 |  1 |         2.465343 |
| PctEssEmpl        | 2.232690 |  1 |         1.494219 |
| PctOver65         | 1.836148 |  1 |         1.355045 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model3.ILIpneu, 2))
```

|                   | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
| :---------------- | ------------: | --------: | ----: | -----: |
| (Intercept)       |          7.00 | 0.0000000 |  4.28 |  11.43 |
| Time              |          1.04 | 0.0000000 |  1.04 |   1.05 |
| PctOvercrowdedQ2  |          1.34 | 0.0000000 |  1.27 |   1.43 |
| PctOvercrowdedQ3  |          1.43 | 0.0000000 |  1.34 |   1.52 |
| PctOvercrowdedQ4  |          1.53 | 0.0000000 |  1.43 |   1.65 |
| PctMultigenQ2     |          1.01 | 0.7945176 |  0.94 |   1.08 |
| PctMultigenQ3     |          1.06 | 0.1362862 |  0.99 |   1.14 |
| PctMultigenQ4     |          1.23 | 0.0000014 |  1.13 |   1.35 |
| PctWhite          |          1.00 | 0.0000000 |  0.99 |   1.00 |
| PctBelowPovThresh |          0.98 | 0.0000000 |  0.97 |   0.98 |
| MedianIncome      |          1.00 | 0.0000000 |  1.00 |   1.00 |
| PctEssEmpl        |          0.97 | 0.0000000 |  0.96 |   0.98 |
| PctOver65         |          0.97 | 0.0000000 |  0.96 |   0.98 |

#### Model 4: Add both clinical and socioeconomic covariates to model 1

``` r
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + COPD_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. COPD has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | --------: | -: | ---------------: |
| Time                |  1.000000 |  1 |         1.000000 |
| PctOvercrowded      |  3.789372 |  3 |         1.248613 |
| PctMultigen         |  6.241148 |  3 |         1.356888 |
| BPHIGH\_CrudePrev   | 31.020219 |  1 |         5.569580 |
| DIABETES\_CrudePrev | 31.506003 |  1 |         5.613021 |
| CHD\_CrudePrev      | 27.923966 |  1 |         5.284313 |
| OBESITY\_CrudePrev  | 11.426929 |  1 |         3.380374 |
| COPD\_CrudePrev     | 44.873278 |  1 |         6.698752 |
| CSMOKING\_CrudePrev | 24.971914 |  1 |         4.997191 |
| PctWhite            | 10.953904 |  1 |         3.309668 |
| PctBelowPovThresh   |  9.398259 |  1 |         3.065658 |
| MedianIncome        |  9.023125 |  1 |         3.003852 |
| PctEssEmpl          |  2.970259 |  1 |         1.723444 |
| PctOver65           |  2.710832 |  1 |         1.646461 |

``` r
# Re-fit model without COPD
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | --------: | -: | ---------------: |
| Time                |  1.000000 |  1 |         1.000000 |
| PctOvercrowded      |  3.639142 |  3 |         1.240223 |
| PctMultigen         |  5.916089 |  3 |         1.344846 |
| BPHIGH\_CrudePrev   | 19.786771 |  1 |         4.448232 |
| DIABETES\_CrudePrev | 19.267942 |  1 |         4.389526 |
| CHD\_CrudePrev      | 11.197568 |  1 |         3.346277 |
| OBESITY\_CrudePrev  | 10.860046 |  1 |         3.295458 |
| CSMOKING\_CrudePrev |  7.596813 |  1 |         2.756232 |
| PctWhite            | 10.793511 |  1 |         3.285348 |
| PctBelowPovThresh   |  8.839773 |  1 |         2.973175 |
| MedianIncome        |  8.926680 |  1 |         2.987755 |
| PctEssEmpl          |  2.896746 |  1 |         1.701983 |
| PctOver65           |  2.665769 |  1 |         1.632718 |

``` r
# Re-fit model without COPD and hypertension
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Diabetes now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                     |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | --------: | -: | ---------------: |
| Time                |  1.000000 |  1 |         1.000000 |
| PctOvercrowded      |  3.013370 |  3 |         1.201827 |
| PctMultigen         |  5.570868 |  3 |         1.331437 |
| DIABETES\_CrudePrev | 19.186090 |  1 |         4.380193 |
| CHD\_CrudePrev      |  5.753161 |  1 |         2.398575 |
| OBESITY\_CrudePrev  |  4.492364 |  1 |         2.119520 |
| CSMOKING\_CrudePrev |  7.455840 |  1 |         2.730538 |
| PctWhite            |  7.308061 |  1 |         2.703343 |
| PctBelowPovThresh   |  8.490240 |  1 |         2.913802 |
| MedianIncome        |  8.905563 |  1 |         2.984219 |
| PctEssEmpl          |  2.916488 |  1 |         1.707773 |
| PctOver65           |  2.712402 |  1 |         1.646937 |

``` r
# Re-fit model without COPD, hypertension, and diabetes
model4.ILIpneu <- glm(Count ~ Time + PctOvercrowded + PctMultigen + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65 + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model4.ILIpneu))
```

|                     |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :------------------ | -------: | -: | ---------------: |
| Time                | 1.000000 |  1 |         1.000000 |
| PctOvercrowded      | 3.048808 |  3 |         1.204171 |
| PctMultigen         | 5.450791 |  3 |         1.326610 |
| CHD\_CrudePrev      | 3.155221 |  1 |         1.776294 |
| OBESITY\_CrudePrev  | 4.255068 |  1 |         2.062782 |
| CSMOKING\_CrudePrev | 7.238084 |  1 |         2.690369 |
| PctWhite            | 3.276167 |  1 |         1.810018 |
| PctBelowPovThresh   | 7.875766 |  1 |         2.806380 |
| MedianIncome        | 8.316568 |  1 |         2.883846 |
| PctEssEmpl          | 2.632747 |  1 |         1.622574 |
| PctOver65           | 2.682575 |  1 |         1.637857 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model4.ILIpneu, 2))
```

|                     | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
| :------------------ | ------------: | --------: | ----: | -----: |
| (Intercept)         |         15.51 | 0.0000000 |  9.70 |  24.82 |
| Time                |          1.04 | 0.0000000 |  1.04 |   1.05 |
| PctOvercrowdedQ2    |          1.21 | 0.0000000 |  1.13 |   1.29 |
| PctOvercrowdedQ3    |          1.43 | 0.0000000 |  1.34 |   1.53 |
| PctOvercrowdedQ4    |          1.60 | 0.0000000 |  1.49 |   1.73 |
| PctMultigenQ2       |          1.18 | 0.0000127 |  1.10 |   1.27 |
| PctMultigenQ3       |          1.39 | 0.0000000 |  1.28 |   1.52 |
| PctMultigenQ4       |          1.56 | 0.0000000 |  1.42 |   1.71 |
| CHD\_CrudePrev      |          0.73 | 0.0000000 |  0.71 |   0.75 |
| OBESITY\_CrudePrev  |          1.03 | 0.0000000 |  1.03 |   1.04 |
| CSMOKING\_CrudePrev |          0.98 | 0.0030895 |  0.97 |   0.99 |
| PctWhite            |          1.00 | 0.6247960 |  1.00 |   1.00 |
| PctBelowPovThresh   |          0.98 | 0.0000000 |  0.97 |   0.98 |
| MedianIncome        |          1.00 | 0.0000000 |  1.00 |   1.00 |
| PctEssEmpl          |          0.96 | 0.0000000 |  0.96 |   0.97 |
| PctOver65           |          1.01 | 0.0001021 |  1.00 |   1.02 |

#### Model 5: Bayesian spatiotemporal model, using model 4 covariates

``` r
# Add zip code ID number to design matrix (needed for spatial and temporal random effects)
zipcodeID <- sort(rep(1:length(allZipcodes), 30))
designResponse$ZipID <- zipcodeID
designResponse$ZipID2 <- zipcodeID

# Construct spatiotemporal model using same set of covariates in (reduced) model 4
model5.INLAformula <- Count ~ 1 + f(ZipID, model="bym", offset(Population/10000), graph=NYCadj) + f(ZipID2, Time, model="rw1") + Time + PctOvercrowded + PctMultigen + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + PctOver65
model5.INLA <- inla(model5.INLAformula, family="poisson", data=designResponse, control.compute=list(dic=TRUE,cpo=TRUE))

# Print model results
kable(Print.Model.Results.INLA(model5.INLA, 2))
```

|                     | exp(Estimate) |        SD | 0.025quant | 0.975quant | mode |
| :------------------ | ------------: | --------: | ---------: | ---------: | ---: |
| (Intercept)         |          3.92 | 0.9585759 |       0.58 |      25.15 | 4.03 |
| Time                |          1.04 | 0.0008282 |       1.04 |       1.04 | 1.04 |
| PctOvercrowdedQ2    |          1.42 | 0.1307378 |       1.11 |       1.85 | 1.42 |
| PctOvercrowdedQ3    |          1.41 | 0.1396338 |       1.07 |       1.85 | 1.41 |
| PctOvercrowdedQ4    |          1.73 | 0.2100929 |       1.15 |       2.62 | 1.71 |
| PctMultigenQ2       |          1.51 | 0.1734788 |       1.07 |       2.12 | 1.51 |
| PctMultigenQ3       |          1.94 | 0.1995872 |       1.31 |       2.86 | 1.94 |
| PctMultigenQ4       |          1.84 | 0.2372049 |       1.15 |       2.92 | 1.85 |
| CHD\_CrudePrev      |          0.92 | 0.0581923 |       0.83 |       1.04 | 0.92 |
| OBESITY\_CrudePrev  |          0.97 | 0.0159467 |       0.94 |       1.00 | 0.97 |
| CSMOKING\_CrudePrev |          1.06 | 0.0421120 |       0.98 |       1.16 | 1.06 |
| PctWhite            |          1.00 | 0.0036870 |       0.99 |       1.00 | 1.00 |
| PctBelowPovThresh   |          0.98 | 0.0135750 |       0.95 |       1.01 | 0.98 |
| MedianIncome        |          1.00 | 0.0000028 |       1.00 |       1.00 | 1.00 |
| PctEssEmpl          |          1.01 | 0.0141053 |       0.98 |       1.04 | 1.01 |
| PctOver65           |          1.03 | 0.0208511 |       0.99 |       1.07 | 1.03 |
