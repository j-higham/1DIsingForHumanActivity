getChangePoint <- function(isingStringSets, isingLookupLocation, minimumSegmentLength){

  n <- length(isingStringSets)
  segment1Likelihood=rep(NA,n)
  segment2Likelihood=rep(NA,n)
  changepointLikelihood=rep(NA,n)

  seg1mu <- c(0,0)
  seg2mu <- c(0,0)
  
  for(i in minimumSegmentLength:(n-minimumSegmentLength)){
    seg1mu <- getMu1Mu2(unlist(isingStringSets[1:i], recursive = FALSE), isingLookupLocation, seg1mu[1], seg1mu[2])
    segment1Likelihood[i] <- (1/i)*getLikelihoodOfSegmentGivenMu(unlist(isingStringSets[1:i], recursive = FALSE), isingLookupLocation, seg1mu[1], seg1mu[2])
    seg2mu <- getMu1Mu2(unlist(isingStringSets[1:i], recursive = FALSE), isingLookupLocation, seg2mu[1], seg2mu[2])
    segment2Likelihood[i] <- (1/(n-i))*getLikelihoodOfSegmentGivenMu(unlist(isingStringSets[(i+1):n], recursive = FALSE), isingLookupLocation, seg2mu[1], seg2mu[2])
  }
  changepointLikelihood <- segment1Likelihood+segment2Likelihood
  
  return(which(changepointLikelihood==max(changepointLikelihood, na.rm = TRUE)))
}

getPELTChangePoint <- function(isingStringSets, isingLookupLocation, minimumSegmentLength){
  
  n <- length(isingStringSets)
  penalty <- 3*log(n)
  
  lastChangeCpts=rep(0,n+1)
  lastChangeLikelihood=-penalty
  checklist=NULL
  
  for(i in minimumSegmentLength:(2*minimumSegmentLength-1)){
    lastChangeLikelihood[i+1]=-2*getLikelihoodOfSegment(unlist(isingStringSets[1:i], recursive = FALSE), isingLookupLocation)
    lastChangeCpts[i+1]=0
  }
    
  checklist=c(0,minimumSegmentLength)
  for(tstar in (2*minimumSegmentLength):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      if(tmpt==0){ # needs to be separate as not adding 1
        return(lastChangeLikelihood[tmpt+1]-2*getLikelihoodOfSegment(unlist(isingStringSets[1:tstar], recursive = FALSE), isingLookupLocation)+penalty)
      }
      return(lastChangeLikelihood[tmpt+1]-2*getLikelihoodOfSegment(unlist(isingStringSets[(tmpt+1):tstar], recursive = FALSE), isingLookupLocation)+penalty)
      # no + 1 on the tmpt as otherwise we miss a day of observations
    }))
    lastChangeLikelihood[tstar+1]=min(tmplike,na.rm=TRUE)
    lastChangeCpts[tstar+1]=checklist[which.min(tmplike)[1]]
    checklist=checklist[tmplike<=(lastChangeLikelihood[tstar+1]+penalty)]
    checklist=c(checklist,tstar-minimumSegmentLength+1)
    #checklist=checklist[!is.na(checklist)]
  }
  
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastChangeCpts[last+1])
    last=lastChangeCpts[last+1]
  }
  
  return(sort(fcpt)[-1])
}
