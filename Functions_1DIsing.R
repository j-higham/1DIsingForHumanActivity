getZ <- function(mu1, mu2, isingStringLength, isingLookupLocation){
  
  load(file = paste(isingLookupLocation, isingStringLength, ".Robject", sep=""))
  Z <- lookupTable_Z[which(rownames(lookupTable_Z)==round(mu1,2)),  which(colnames(lookupTable_Z)==round(mu2,2))]
  
  return(Z)
}

getMu1Mu2 <- function(isingStringSet, isingLookupLocation, mu1=0, mu2=0){
  
  numStrings <- length(isingStringSet)
  if(numStrings == 0){
    return(c(NA,NA))
  }
  
  stepSize <- 0.1
  
  i <- 0
  notMaxMu1 <- TRUE
  notMaxMu2 <- TRUE
  while(notMaxMu1 || notMaxMu2){
    i <- i + 1
    currentProb <- calculateStringSetProbability(mu1, mu2, isingLookupLocation, isingStringSet)
    if(notMaxMu1){
      newProbN <- calculateStringSetProbability((mu1+stepSize), mu2, isingLookupLocation, isingStringSet)
      newProbS <- calculateStringSetProbability((mu1-stepSize), mu2, isingLookupLocation, isingStringSet)
    }
    else{
      newProbN <= currentProb
      newProbS <= currentProb
    }
    if(notMaxMu2){
      newProbE <- calculateStringSetProbability(mu1, (mu2+stepSize), isingLookupLocation, isingStringSet)
      newProbW <- calculateStringSetProbability(mu1, (mu2-stepSize), isingLookupLocation, isingStringSet)
    }
    else{
      newProbE <= currentProb
      newProbW <= currentProb
    }
    
#    print(paste(mu1, " ", mu2))
    if(max(currentProb, newProbN, newProbE, newProbS, newProbW) == currentProb){
      notMaxMu1 = FALSE
      notMaxMu2 = FALSE
    }
    if(max(currentProb, newProbN, newProbE, newProbS, newProbW) == newProbN){
      mu1 <- mu1+stepSize
    }
    if(max(currentProb, newProbN, newProbE, newProbS, newProbW) == newProbE){
      mu2 <- mu2+stepSize
    }
    if(max(currentProb, newProbN, newProbE, newProbS, newProbW) == newProbS){
      mu1 <- mu1-stepSize
    }
    if(max(currentProb, newProbN, newProbE, newProbS, newProbW) == newProbW){
      mu2 <- mu2-stepSize
    }
    if(mu1 > 2.0 || mu1 < -2.0){
      notMaxMu1 = FALSE
    } 
    if(mu2 > 2.0 || mu2 < -2.0){
      notMaxMu2 = FALSE
    } 
  }
  
  return(c(mu1, mu2))
}

calculateStringSetProbability <- function(mu1, mu2, isingLookupLocation, isingStringSet){

  numStrings <- length(isingStringSet)
  if(numStrings == 0){
    return(NA)
  }
  
  isingStringSetEnergy <- rep(0, numStrings)
  isingStringSetProbability <- rep(0, numStrings)
  
  for(string in 1:numStrings){
  
    stringLength <- length(isingStringSet[[string]])
    Z <- getZ(mu1,mu2,stringLength,isingLookupLocation)
    
    isingString <- isingStringSet[[string]]
    isingStringSetProbability[string] <- calculateStringProbability(mu1, mu2, Z, isingString)
  }
  totalProbability <- sum(log(isingStringSetProbability))
  
  return(totalProbability)
}

calculateStringProbability <- function(mu1, mu2, Z, isingString){
  
  isingStringEnergy <- 0
  isingStringProbability <- 0
  
  stringLength <- length(isingString)

  for(i in 1:stringLength){
    isingStringEnergy <- isingStringEnergy - mu1*isingString[i]
  }
  for(j in 1:(stringLength-1)){
    isingStringEnergy <- isingStringEnergy - mu2*isingString[j]*isingString[j+1]
  }
  isingStringProbability <- exp(-isingStringEnergy)/Z
    
  return(isingStringProbability)
}
