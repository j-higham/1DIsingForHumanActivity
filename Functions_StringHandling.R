getHourlyStringSets <- function(isingData, hour){
  
  startTime <- min(isingData$t)
  endTime <- max(isingData$t)
  startDay <- floor((startTime)/seconds_per_day)
  endDay <- floor((endTime)/seconds_per_day)
  numDays <- endDay-startDay
  
  isingStringSets <- list(list())
  
  for(day in 1:numDays){
    
    isingStringSet <- list()
    
    dayStart <- (startTime+((day-1)*seconds_per_day))
    hourStart <- dayStart + ((hour-1)*seconds_per_hour)
    hourEnd <- hourStart + seconds_per_hour
    
    if(!is.na(which(isingData$t > hourEnd)[1]) && !is.na(which(isingData$t > hourStart)[1])){
      hourString <- isingData$isingString[(which(isingData$t > hourStart)[1]):((which(isingData$t > hourEnd)[1]))]
    }
    if(is.na(which(isingData$t > hourEnd)[1]) && !is.na(which(isingData$t > hourStart)[1])){
      hourString <- isingData$isingString[(which(isingData$t > hourStart)[1]):(length(isingData$t))]
    }
    if(is.na(which(isingData$t > hourEnd)[1]) && is.na(which(isingData$t > hourStart)[1])){
      hourString <- NULL
    }
    
    while(length(hourString) > 2*ising_string_length){
      isingStringSet[[length(isingStringSet)+1]] <- hourString[1:ising_string_length]
      hourString <- hourString[(ising_string_length+1):length(hourString)]
    }
    if(length(hourString) > ising_string_length){
      isingStringSet[[length(isingStringSet)+1]] <- hourString[1:(floor(length(hourString)/2))]
      isingStringSet[[length(isingStringSet)+1]] <- hourString[(ceiling(length(hourString)/2)):(length(hourString))]
    }
    else if (length(hourString) > 4){
      isingStringSet[[length(isingStringSet)+1]] <- hourString
    }
    isingStringSets[[day]] <- isingStringSet
  }
  
  return(isingStringSets)
}

getLikelihoodOfSegmentGivenMu <- function(isingStringSet, isingLookupLocation, mu1=0, mu2=0){
  
  totalProbability <- calculateStringSetProbability(mu1, mu2, isingLookupLocation, isingStringSet)
  
  return(totalProbability)
}

getLikelihoodOfSegment <- function(isingStringSet, isingLookupLocation){
  
  mu <- getMu1Mu2(isingStringSet, isingLookupLocation)
  totalProbability <- calculateStringSetProbability(mu[1], mu[2], isingLookupLocation, isingStringSet)
  
  return(totalProbability)
}

getAllStringSetProbabilities <- function(isingStringSets, isingLookupLocation){
  
  allStringSetProbabilities <- rep(0, length(isingStringSets))
  for(i in 1:length(isingStringSets)){
    mu <- getMu1Mu2(unlist(isingStringSets[i], recursive = FALSE), isingLookupLocation)
    allStringSetProbabilities[i] <- calculateStringSetProbability(mu[1], mu[2], isingLookupLocation, unlist(isingStringSets[i], recursive = FALSE))
  }
  
  return(allStringSetProbabilities)
}