# Preprocess the tropical cyclone dataset for model training.
# 
# @param dta = the raw data frame from the envData dataset.
#
# @returns a list of 2 elements: the preprocessed data excluding the response, and
# the response

preprocess = function(dta){
  
  # encode missing data from 'basin' as an actual category
  #dta$basin = sapply(dta$basin, function(cat){if(is.na(cat)){"NA"}else{cat}})
  dta$basin[is.na(dta$basin)] = "NA" # better readability and efficiency
  
  # might be a good idea to keep basin as a factor for modelling purposes!
  dta$basin = as.factor(dta$basin)
  
  # encode categories in 'basin' using one hot encoding
  dta.encoded = cbind(dta, model.matrix(~ 0 + basin, data = dta))
  
  # get data frame of predictors only: we can remove the original 'basin' column
  dta.predictors = subset(dta.encoded, select = -c(basin, TC_genesis))
  
  # get data frame of the response only
  dta.response = dta.encoded$TC_genesis
  
  return(list(predictors = dta.predictors,
              response = dta.response))
}