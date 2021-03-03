# --- Read data ---

# Read data, set rownames to zipcodes, remove "X" from dates in column names
respiratory <- read.csv("../Paper 1 Surveillance Analysis/Datasets Produced in R/RespiratoryInterpolated.csv")
influenza <- read.csv("../Paper 1 Surveillance Analysis/Datasets Produced in R/InfluenzaInterpolated.csv")
pneumonia <- read.csv("../Paper 1 Surveillance Analysis/Datasets Produced in R/PneumoniaInterpolated.csv")

respiratory$zipcode <- as.character(respiratory$zipcode)
influenza$zipcode <- as.character(influenza$zipcode)
pneumonia$zipcode <- as.character(pneumonia$zipcode)
rownames(respiratory) <- respiratory$zipcode
rownames(influenza) <- influenza$zipcode
rownames(pneumonia) <- pneumonia$zipcode
colnames(respiratory) <- gsub("X", "", colnames(respiratory))
colnames(influenza) <- gsub("X", "", colnames(influenza))
colnames(pneumonia) <- gsub("X", "", colnames(pneumonia))

# --- Define model covariates, observation dates, and relevant zipcode sets ---

# Define covariates to include in model
modelCovariates <- c("PctOvercrowded", "PctMultigen", "ACCESS2_CrudePrev", 
                     "BPHIGH_CrudePrev", "CASTHMA_CrudePrev", "CHD_CrudePrev", 
                     "COPD_CrudePrev", "CSMOKING_CrudePrev", "DIABETES_CrudePrev",
                     "OBESITY_CrudePrev", "PctBelowPovThresh", "MedianIncome", 
                     "PctWhite", "Construction", "Manufacturing", "WholesaleTrade", 
                     "RetailTrade", "TransportUtil", "EducHealthSoc", "PctOver65")

# Define dates in the observational window (March 2020)
observationDates <- as.Date(gsub("\\.", "-", colnames(respiratory[2:31])))

# Get list of all zipcodes
allZipcodes <- respiratory$zipcode

# Get list of zipcodes in 3rd quartile of overcrowding percentage 
overcrowdedZips <- respiratory$zipcode[respiratory$PctOvercrowded > quantile(respiratory$PctOvercrowded, na.rm=T)[4]]
nonovercrowdedZips <- respiratory$zipcode[!(respiratory$zipcode %in% overcrowdedZips)]

# Get a list of zipcodes in each quartile of overcrowding percentage
zipCrowdQ1 <- respiratory$zipcode[respiratory$PctOvercrowded <= quantile(respiratory$PctOvercrowded, na.rm=T)[2]]
zipCrowdQ2 <- respiratory$zipcode[respiratory$PctOvercrowded > quantile(respiratory$PctOvercrowded, na.rm=T)[2] & respiratory$PctOvercrowded <= quantile(respiratory$PctOvercrowded, na.rm=T)[3]]
zipCrowdQ3 <- respiratory$zipcode[respiratory$PctOvercrowded > quantile(respiratory$PctOvercrowded, na.rm=T)[3] & respiratory$PctOvercrowded <= quantile(respiratory$PctOvercrowded, na.rm=T)[4]]
zipCrowdQ4 <- respiratory$zipcode[respiratory$PctOvercrowded > quantile(respiratory$PctOvercrowded, na.rm=T)[4]]

# Get a list of zipcodes in each quartile of multigen. housing percentage
zipMultigenQ1 <- respiratory$zipcode[respiratory$PctMultigen <= quantile(respiratory$PctMultigen, na.rm=T)[2]]
zipMultigenQ2 <- respiratory$zipcode[respiratory$PctMultigen > quantile(respiratory$PctMultigen, na.rm=T)[2] & respiratory$PctMultigen <= quantile(respiratory$PctMultigen, na.rm=T)[3]]
zipMultigenQ3 <- respiratory$zipcode[respiratory$PctMultigen > quantile(respiratory$PctMultigen, na.rm=T)[3] & respiratory$PctMultigen <= quantile(respiratory$PctMultigen, na.rm=T)[4]]
zipMultigenQ4 <- respiratory$zipcode[respiratory$PctMultigen > quantile(respiratory$PctMultigen, na.rm=T)[4]]

# --- Function definitions ---

# Construct segmented regression design matrix + response vector for one zipcode 
# for one illness. illness is one of "respiratory", "influenza", or "pneumonia".
# interventionDate is a string ("YYYY-MM-DD") indicating when intervention began, if 
# interrupted time series models are being run. variablesToDiscretize is a vector of
# variable names that will be dichotomized if quartile=FALSE or binned into quartiles
# if quartile=TRUE.
Zipcode.Design.Matrix.Response <- function(zipcode, illness, interventionDate, variablesToDiscretize=c(), quartile=F) {
  if(illness == "respiratory") { 
    zipcodeData <- respiratory[zipcode, ] 
  } else if(illness == "influenza") { 
    zipcodeData <- influenza[zipcode, ] 
  } else if(illness == "pneumonia") { 
    zipcodeData <- pneumonia[zipcode, ] 
  }  
  response <- t(zipcodeData[1,2:31])
  design <- as.data.frame(matrix(0, nrow=nrow(response), ncol=5+length(modelCovariates)))
  
  # 1st column: time since start of observational window
  design[,1] <- as.numeric(observationDates - observationDates[1] + 1)
  
  # 2nd column: {0,1} for pre-/post-intervention
  interventionDate <- as.Date(interventionDate)
  design[,2] <- as.numeric(observationDates >= interventionDate)
  
  # 3rd column: time since start of intervention
  design[,3] <- as.numeric(pmax(observationDates - interventionDate + 1, rep(0, nrow(response))))
  
  # 4th through third to last columns: zipcode-level covariates
  design[,4:(ncol(design)-2)] <- unname(as.matrix(zipcodeData[,modelCovariates])[rep(1,nrow(response)),])
  
  # Second to last column: population
  design[,ncol(design)-1] <- zipcodeData$Population
  
  # Last column: zipcode
  design[,ncol(design)] <- zipcode
  
  # Concatenate response and design matrix into one data frame
  designAndResponse <- cbind(as.data.frame(response), as.data.frame(design))
  
  # Set column names
  colnames(designAndResponse) <- c("Count", "Time", "D", "P", modelCovariates, "Population", "Zipcode")
  rownames(designAndResponse) <- paste(rownames(designAndResponse), rep(paste("_", zipcode, sep=""), 30), sep="")
  
  # Discretize variables
  if(length(variablesToDiscretize) > 0) {
    if(quartile == FALSE) {
      for(var in variablesToDiscretize) {
        designAndResponse[, var] <- cut(designAndResponse[, var], quantile(respiratory[, var])[c(1,3,5)], include.lowest=T, labels=c("BelowMed", "AboveMed"))
      }
    } else {
      for(var in variablesToDiscretize) {
        designAndResponse[, var] <- cut(designAndResponse[, var], quantile(respiratory[, var]), include.lowest=T, labels=c("Q1", "Q2", "Q3", "Q4"))
      }
    }
  }
  
  return(designAndResponse)
}

# Stack the response vectors and design matrices for several zipcodes into one
# dataframe, for a specified illness and intervention date. The variablesToDiscretize
# and quartile arguments are the same as above.
Concatenate.Zipcode.Data <- function(zipcodes, illness, interventionDate, variablesToDiscretize=c(), quartile=F) {
  designAndResponse <- c()
  for(i in 1:length(zipcodes)) {
    designAndResponse_thisZipcode <- Zipcode.Design.Matrix.Response(zipcodes[i], illness, interventionDate, variablesToDiscretize, quartile)
    designAndResponse <- rbind(designAndResponse, designAndResponse_thisZipcode)
  }
  return(designAndResponse)
}

