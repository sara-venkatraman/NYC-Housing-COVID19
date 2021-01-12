Models for Table 2
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
i.e. the total number of ILI + pneumonia cases observed on that day. In
this dataframe, both overcrowdedness and multigenerational housing are
binned into quartiles.

![](DesignMatrixSetup.png)

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
model1.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. None need to be removed (based on VIF > 10 criterion)
kable(vif(model1.ILIpneu))
```

|                         |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | -------: | -: | ---------------: |
| Time                    | 1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) | 1.432021 |  3 |         1.061675 |
| ordered(PctMultigen)    | 1.432021 |  3 |         1.061675 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model1.ILIpneu, 2))
```

|                           | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
| :------------------------ | ------------: | --------: | ----: | -----: |
| (Intercept)               |          0.71 | 0.0000000 |  0.68 |   0.74 |
| Time                      |          1.04 | 0.0000000 |  1.04 |   1.05 |
| ordered(PctOvercrowded).L |          1.66 | 0.0000000 |  1.59 |   1.73 |
| ordered(PctOvercrowded).Q |          0.87 | 0.0000000 |  0.83 |   0.90 |
| ordered(PctOvercrowded).C |          1.09 | 0.0000038 |  1.05 |   1.12 |
| ordered(PctMultigen).L    |          1.39 | 0.0000000 |  1.33 |   1.45 |
| ordered(PctMultigen).Q    |          1.01 | 0.4810760 |  0.97 |   1.05 |
| ordered(PctMultigen).C    |          1.05 | 0.0056912 |  1.01 |   1.09 |

#### Model 2: Add clinical risk factors for COVID-19 to model 1

``` r
model2.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + COPD_CrudePrev + CSMOKING_CrudePrev + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. COPD has the largest VIF.
kable(vif(model2.ILIpneu))
```

|                         |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | --------: | -: | ---------------: |
| Time                    |  1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) |  2.584312 |  3 |         1.171451 |
| ordered(PctMultigen)    |  3.535624 |  3 |         1.234272 |
| BPHIGH\_CrudePrev       | 19.596931 |  1 |         4.426842 |
| DIABETES\_CrudePrev     | 25.024530 |  1 |         5.002452 |
| CHD\_CrudePrev          | 21.534671 |  1 |         4.640546 |
| OBESITY\_CrudePrev      |  8.007129 |  1 |         2.829687 |
| COPD\_CrudePrev         | 42.218345 |  1 |         6.497565 |
| CSMOKING\_CrudePrev     | 25.581867 |  1 |         5.057852 |

``` r
# Re-fit model without COPD
model2.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model2.ILIpneu))
```

|                         |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | --------: | -: | ---------------: |
| Time                    |  1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) |  2.493952 |  3 |         1.164523 |
| ordered(PctMultigen)    |  3.385975 |  3 |         1.225408 |
| BPHIGH\_CrudePrev       | 11.208735 |  1 |         3.347945 |
| DIABETES\_CrudePrev     |  9.748967 |  1 |         3.122333 |
| CHD\_CrudePrev          |  4.412565 |  1 |         2.100611 |
| OBESITY\_CrudePrev      |  8.007108 |  1 |         2.829683 |
| CSMOKING\_CrudePrev     |  6.169262 |  1 |         2.483800 |

``` r
# Re-fit model without COPD and hypertension
model2.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model2.ILIpneu))
```

|                         |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | -------: | -: | ---------------: |
| Time                    | 1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) | 1.973983 |  3 |         1.120015 |
| ordered(PctMultigen)    | 3.519230 |  3 |         1.233317 |
| DIABETES\_CrudePrev     | 6.924940 |  1 |         2.631528 |
| CHD\_CrudePrev          | 3.155821 |  1 |         1.776463 |
| OBESITY\_CrudePrev      | 3.406973 |  1 |         1.845799 |
| CSMOKING\_CrudePrev     | 4.520681 |  1 |         2.126189 |

``` r
# Print model results
kable(Print.Model.Results.GLM(model2.ILIpneu, 2))
```

|                           | exp(Estimate) |   p-value | 2.5 % | 97.5 % |
| :------------------------ | ------------: | --------: | ----: | -----: |
| (Intercept)               |          0.70 | 0.0000009 |  0.59 |   0.83 |
| Time                      |          1.04 | 0.0000000 |  1.04 |   1.05 |
| ordered(PctOvercrowded).L |          1.54 | 0.0000000 |  1.47 |   1.61 |
| ordered(PctOvercrowded).Q |          1.00 | 0.9600630 |  0.96 |   1.04 |
| ordered(PctOvercrowded).C |          1.01 | 0.6401624 |  0.98 |   1.04 |
| ordered(PctMultigen).L    |          1.02 | 0.5161690 |  0.95 |   1.09 |
| ordered(PctMultigen).Q    |          0.93 | 0.0007162 |  0.90 |   0.97 |
| ordered(PctMultigen).C    |          1.00 | 0.8268683 |  0.96 |   1.03 |
| DIABETES\_CrudePrev       |          1.16 | 0.0000000 |  1.14 |   1.18 |
| CHD\_CrudePrev            |          0.72 | 0.0000000 |  0.69 |   0.74 |
| OBESITY\_CrudePrev        |          1.01 | 0.0000004 |  1.01 |   1.02 |
| CSMOKING\_CrudePrev       |          0.99 | 0.3003368 |  0.98 |   1.01 |

#### Model 3: Add socioeconomic covariates to model 1

``` r
model3.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No variables need to be removed.
kable(vif(model3.ILIpneu))
```

|                         |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | -------: | -: | ---------------: |
| Time                    | 1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) | 2.276243 |  3 |         1.146929 |
| ordered(PctMultigen)    | 3.865827 |  3 |         1.252777 |
| PctWhite                | 2.798164 |  1 |         1.672771 |
| PctBelowPovThresh       | 4.802450 |  1 |         2.191449 |
| MedianIncome            | 5.671430 |  1 |         2.381476 |
| PctEssEmpl              | 2.085771 |  1 |         1.444220 |

#### Model 4: Add both clinical and socioeconomic covariates to model 1

``` r
model4.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + COPD_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. COPD has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                         |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | --------: | -: | ---------------: |
| Time                    |  1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) |  3.504577 |  3 |         1.232459 |
| ordered(PctMultigen)    |  6.205648 |  3 |         1.355599 |
| BPHIGH\_CrudePrev       | 31.117292 |  1 |         5.578288 |
| DIABETES\_CrudePrev     | 31.073302 |  1 |         5.574343 |
| CHD\_CrudePrev          | 25.264831 |  1 |         5.026413 |
| OBESITY\_CrudePrev      | 11.400263 |  1 |         3.376428 |
| COPD\_CrudePrev         | 44.143526 |  1 |         6.644059 |
| CSMOKING\_CrudePrev     | 25.081367 |  1 |         5.008130 |
| PctWhite                | 10.769059 |  1 |         3.281624 |
| PctBelowPovThresh       |  8.739480 |  1 |         2.956261 |
| MedianIncome            |  9.005577 |  1 |         3.000929 |
| PctEssEmpl              |  2.924882 |  1 |         1.710229 |

``` r
# Re-fit model without COPD
model4.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + BPHIGH_CrudePrev + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Hypertension now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                         |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | --------: | -: | ---------------: |
| Time                    |  1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) |  3.380065 |  3 |         1.225051 |
| ordered(PctMultigen)    |  5.902429 |  3 |         1.344328 |
| BPHIGH\_CrudePrev       | 19.787088 |  1 |         4.448268 |
| DIABETES\_CrudePrev     | 19.211002 |  1 |         4.383036 |
| CHD\_CrudePrev          |  9.945668 |  1 |         3.153675 |
| OBESITY\_CrudePrev      | 10.882241 |  1 |         3.298824 |
| CSMOKING\_CrudePrev     |  7.436230 |  1 |         2.726945 |
| PctWhite                | 10.541991 |  1 |         3.246843 |
| PctBelowPovThresh       |  7.954881 |  1 |         2.820440 |
| MedianIncome            |  8.905729 |  1 |         2.984247 |
| PctEssEmpl              |  2.835891 |  1 |         1.684010 |

``` r
# Re-fit model without COPD and hypertension
model4.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + DIABETES_CrudePrev + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. Diabetes now has the largest VIF.
kable(vif(model4.ILIpneu))
```

|                         |      GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | --------: | -: | ---------------: |
| Time                    |  1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) |  2.793404 |  3 |         1.186740 |
| ordered(PctMultigen)    |  5.584043 |  3 |         1.331961 |
| DIABETES\_CrudePrev     | 19.117271 |  1 |         4.372330 |
| CHD\_CrudePrev          |  4.726752 |  1 |         2.174110 |
| OBESITY\_CrudePrev      |  4.428040 |  1 |         2.104291 |
| CSMOKING\_CrudePrev     |  7.302073 |  1 |         2.702235 |
| PctWhite                |  7.132344 |  1 |         2.670645 |
| PctBelowPovThresh       |  7.693449 |  1 |         2.773707 |
| MedianIncome            |  8.875607 |  1 |         2.979196 |
| PctEssEmpl              |  2.845218 |  1 |         1.686778 |

``` r
# Re-fit model without COPD, hypertension, and diabetes
model4.ILIpneu <- glm(Count ~ Time + ordered(PctOvercrowded) + ordered(PctMultigen) + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite  + PctBelowPovThresh + MedianIncome + PctEssEmpl + offset(log(Population/10000)), family=quasipoisson, data=designResponse)

# Check variance inflation factors. No more variables need to be removed.
kable(vif(model4.ILIpneu))
```

|                         |     GVIF | Df | GVIF^(1/(2\*Df)) |
| :---------------------- | -------: | -: | ---------------: |
| Time                    | 1.000000 |  1 |         1.000000 |
| ordered(PctOvercrowded) | 2.790353 |  3 |         1.186524 |
| ordered(PctMultigen)    | 5.446228 |  3 |         1.326425 |
| CHD\_CrudePrev          | 2.263282 |  1 |         1.504421 |
| OBESITY\_CrudePrev      | 4.172275 |  1 |         2.042615 |
| CSMOKING\_CrudePrev     | 7.062900 |  1 |         2.657612 |
| PctWhite                | 3.174338 |  1 |         1.781667 |
| PctBelowPovThresh       | 7.222658 |  1 |         2.687500 |
| MedianIncome            | 8.308585 |  1 |         2.882462 |
| PctEssEmpl              | 2.576829 |  1 |         1.605250 |

#### Model 5: Bayesian spatiotemporal model, using model 4 covariates

``` r
# Add zip code ID number to design matrix (needed for spatial and temporal random effects)
zipcodeID <- sort(rep(1:length(allZipcodes), 30))
designResponse$ZipID <- zipcodeID
designResponse$ZipID2 <- zipcodeID

# Construct spatiotemporal model using same set of covariates in (reduced) model 4
model5.INLAformula <- Count ~ 1 + f(ZipID, model="bym", offset(Population/10000), graph=NYCadj) + f(ZipID2, Time, model="rw1") + Time + ordered(PctOvercrowded) + ordered(PctMultigen) + CHD_CrudePrev + OBESITY_CrudePrev + CSMOKING_CrudePrev + PctWhite + PctBelowPovThresh + MedianIncome + PctEssEmpl
model5.INLA <- inla(model5.INLAformula, family="poisson", data=designResponse, control.compute=list(dic=TRUE,cpo=TRUE))
kable(Print.Model.Results.INLA(model5.INLA, 6))
```

|                           | exp(Estimate) |        SD | 0.025quant | 0.975quant |      mode |
| :------------------------ | ------------: | --------: | ---------: | ---------: | --------: |
| (Intercept)               |     15.613213 | 0.8217829 |   3.073092 |  77.607986 | 15.795646 |
| Time                      |      1.040350 | 0.0008313 |   1.038642 |   1.042037 |  1.040362 |
| ordered(PctOvercrowded).L |      1.401975 | 0.1445271 |   1.058248 |   1.867164 |  1.397868 |
| ordered(PctOvercrowded).Q |      0.902654 | 0.0943710 |   0.749950 |   1.086741 |  0.902435 |
| ordered(PctOvercrowded).C |      1.140009 | 0.0878159 |   0.961780 |   1.358134 |  1.136941 |
| ordered(PctMultigen).L    |      1.549097 | 0.1667621 |   1.111431 |   2.141446 |  1.555782 |
| ordered(PctMultigen).Q    |      0.797066 | 0.1076032 |   0.644002 |   0.983245 |  0.798414 |
| ordered(PctMultigen).C    |      0.963156 | 0.0930177 |   0.801567 |   1.155610 |  0.963817 |
| CHD\_CrudePrev            |      0.922837 | 0.0576141 |   0.824977 |   1.034373 |  0.921792 |
| OBESITY\_CrudePrev        |      0.967290 | 0.0154396 |   0.937756 |   0.996402 |  0.968032 |
| CSMOKING\_CrudePrev       |      1.076362 | 0.0409342 |   0.995350 |   1.168822 |  1.073878 |
| PctWhite                  |      0.995736 | 0.0036233 |   0.988593 |   1.002780 |  0.995813 |
| PctBelowPovThresh         |      0.973244 | 0.0127931 |   0.948814 |   0.997758 |  0.973521 |
| MedianIncome              |      0.999980 | 0.0000027 |   0.999974 |   0.999985 |  0.999980 |
| PctEssEmpl                |      1.008430 | 0.0138817 |   0.981267 |   1.036257 |  1.008460 |
