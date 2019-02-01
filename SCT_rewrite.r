#!/usr/bin/env Rscript
#-----------------------------------------------------------------------------
options(warn = 2, scipen = 999)
#------------------------------------------------------------------------------
# FUNCTIONS

#+ SCT - spatial consistency test
spatial_consistency_test<-function(t2,
                                   x,
                                   y,
                                   z,
                                   t,
                                   vp) {
#------------------------------------------------------------------------------
  n<-length(x)
  # distance matrices
  disth<-(outer(x,x,FUN="-")**2.+ outer(y,y,FUN="-")**2.)**0.5
  distz<-abs(outer(z,z,FUN="-"))
  print(disth)
  print(distz)
  # station and all the others
  Dh<-max(10000,
          mean(apply(cbind(1:nrow(disth),disth),
                     MARGIN=1,
                     FUN=function(x){
       as.numeric(quantile(x[which(!((1:length(x))%in%c(1,(x[1]+1))))],
                          probs=0.1))})))
  Dz<-200
  # background error correlation matrix
  S<-exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
  Dh = 11600
  print("Dh")
  print(Dh)
  # S+eps2I
  diag(S)<-diag(S)+0.5
  print(round(S,2))
  # innvoation
  d<-t-vp
  # init station flags
  index_ok<-1:n
  stationFlags<-rep(0,n)
  first<-T
  # loop over SCT iterations
  while (length(which(stationFlags==0))>1) {
    # update selection
    # first iteration, inver the matrix
    indx<-which(stationFlags==1)
#    if (length(indx)>0) print(stationFlags)
    if (first) {
      SRinv<-chol2inv(chol(S))
      print(round(SRinv,2))
      # from S+R go back to S
      diag(S)<-diag(S)-0.5
      first<-F
    } else if (length(indx)>1) {
      S<-S[-indx,-indx]
      diag(S)<-diag(S)+0.5
      SRinv<-chol2inv(chol(S))
      # from S+R go back to S
      diag(S)<-diag(S)-0.5
      d<-d[-indx]
      stationFlags<-stationFlags[-indx]
      index_ok<-index_ok[-indx]
    } else {
      # Update inverse matrix (Uboldi et al 2008, Appendix AND erratum!)
      aux<-SRinv
      SRinv<-aux[-indx,-indx]-
             (tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
      S<-S[-indx,-indx]
      rm(aux)
      d<-d[-indx]
      stationFlags<-stationFlags[-indx]
      index_ok<-index_ok[-indx]
    }
    # next tree lines: compute cvres=(obs - CrossValidation_prediction)
    Zinv<-1/diag(SRinv)
    print(round(Zinv,2))
    SRinv.d<-crossprod(SRinv,d)
    print(round(SRinv.d,2))
    ares<-crossprod(S,SRinv.d)-d      #   a-Obs
    print("d")
    print(round(d,2))
    print("ares")
    print(round(-ares,2))
    cvres<--Zinv*SRinv.d              # CVa-Obs, Lussana et al 2010 Eq(A13)
    sig2o<-mean(d*(-ares))            # Lussana et al 2010, Eq(32)
    print("sig2o")
    print(sig2o)
    if (sig2o<0.01) sig2o<-0.01       # safe threshold
    # pog=cvres/(sig2obs+sig2CVpred), Lussana et al 2010 Eq(20)
    pog<-(ares*cvres)/sig2o
    # check if any obs fails the test
    if (any(pog>t2)) {
      stationFlags[which(pog>t2)]<-1
    } else {
      break
    }
  } # end cycle SCT model
  return(index_ok)
}

data<-read.table(file="testdata.txt",header=T,sep=";",stringsAsFactors=F,strip.white=T)

index_ok<-spatial_consistency_test(t2=3,
                                   x=data$x,
                                   y=data$y,
                                   z=data$z,
                                   t=data$t,
                                   vp=data$vp)
print("station id, obs, background")
print(cbind(1:length(data$x),data$t,data$vp))
print(paste("good observations:",toString(index_ok)))
q()
