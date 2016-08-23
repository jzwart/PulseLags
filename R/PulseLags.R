# JAZ, code for analyses described in Zwart et al. doi: 


# the following code (lines 8-409) was used to detect extreme precipitation events and duration and magnitude of response 
#    for t-DOC load and metabolic metrics for 3 different lakes (MO, EL, and CR) 
# The .Rdata file on GitHub has compiled all the analyses using this code. lake specific data are organized in data frames 
# as the lakeID followed by Data (e.g. 'elData') 

events<-data$dateTime[data$precip>(mean(data$precip)+2*sd(data$precip))]
# if events are too close together (i.e. 5 days) then just take earliest one 
eventsNew<-events[1]
for(i in 2:length(events)){
  if(abs((DOY(events[i])-DOY(events[i-1])))<5){
    next
  }else{
    eventsNew<-c(eventsNew,events[i])
  }
}

# using window size of 5 days to find a threshold limit on the response 
#need to remove seasonal trend in data to find appropriate responses 
# using quadratic eqation to remove this 
n=5
temp<-data.frame(dateTime=data$dateTime,R=data$R,GPP=data$GPP,NEP=data$NEP,RInt=data$RInt,
                 GPPInt=data$GPPInt,NEPInt=data$NEPInt,GPPweight=data$GPPweight,Rweight=data$Rweight,NEPweight=data$NEPweight)
temp2<-data.frame()
for(i in 1:length(years)){
  cur<-temp[strftime(temp$dateTime,'%Y')==years[i],]  
  rQuad<-lm(cur$Rweight~DOY(cur$dateTime)+I(DOY(cur$dateTime)^2))
  gppQuad<-lm(cur$GPPweight~DOY(cur$dateTime)+I(DOY(cur$dateTime)^2))
  nepQuad<-lm(cur$NEPweight~DOY(cur$dateTime)+I(DOY(cur$dateTime)^2))
  
  cur$rQuad<-(rQuad$coefficients[2]*DOY(cur$dateTime)+rQuad$coefficients[3]*DOY(cur$dateTime)^2+rQuad$coefficients[1])
  cur$gppQuad<-(gppQuad$coefficients[2]*DOY(cur$dateTime)+gppQuad$coefficients[3]*DOY(cur$dateTime)^2+gppQuad$coefficients[1])
  cur$nepQuad<-(nepQuad$coefficients[2]*DOY(cur$dateTime)+nepQuad$coefficients[3]*DOY(cur$dateTime)^2+nepQuad$coefficients[1])
  temp2<-rbind(temp2,cur)
}
temp<-temp2


#added on the min or max value of filtered data to constrain metabolism positive or negative
temp$rFilter<-temp$Rweight-temp$rQuad-min(temp$Rweight-temp$rQuad,na.rm=T) 
temp$gppFilter<-temp$GPPweight-temp$gppQuad-min(temp$GPPweight-temp$gppQuad,na.rm=T)
temp$nepFilter<-temp$NEPweight-temp$nepQuad-max(temp$NEPweight-temp$nepQuad,na.rm=T)

# rho 
## what is the maximum displacement of NEP after an event? 2014-11-17
eventResponseR<-data.frame(dateTime=eventsNew,threshold=rep(NA,length(eventsNew)),
                           Rtotal=rep(NA,length(eventsNew)),backgroundR=rep(NA,length(eventsNew)),
                           excessR=rep(NA,length(eventsNew)),maxDisplace=rep(NA,length(eventsNew)),
                           duration=rep(NA,length(eventsNew)),lag=rep(NA,length(eventsNew)))
window<-5
for(i in 1:length(eventsNew)){
  cur<-temp[strftime(temp$dateTime,'%Y')==strftime(eventsNew[i],'%Y'),] #have to split up data by years 
  eventResponseR$threshold[i]<-mean(cur$rFilter[DOY(cur$dateTime)%in%
                                                  (DOY(eventResponseR$dateTime[i])-window-1):(DOY(eventResponseR$dateTime[i])-1)],na.rm=T)+
    sd(cur$rFilter[DOY(cur$dateTime)%in%
                     (DOY(eventResponseR$dateTime[i])-window-1):(DOY(eventResponseR$dateTime[i])-1)],na.rm=T)*2
  curMean<-mean(cur$rFilter[DOY(cur$dateTime)%in%
                              (DOY(eventResponseR$dateTime[i])-window-1):(DOY(eventResponseR$dateTime[i])-1)],na.rm=T)
  curMean<-curMean+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseR$dateTime[i])]
  # docLoad after event - exceeding threshold 
  curR<-0
  curDuration<-0
  curLag<-0
  curBack<-0
  curX<-0
  curMaxDis<-c()
  curThresh<-eventResponseR$threshold[i]+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseR$dateTime[i])] # threshold that is back calculated to real values (Rweighted values )
  nextEvent<-min((length(cur$dateTime)+1),grep(eventResponseR$dateTime[i+1],cur$dateTime),na.rm = T)
  j<-grep(eventResponseR$dateTime[i],cur$dateTime) #variable response start time due to lag in response...
  if(cur$rFilter[j]<eventResponseR$threshold[i]&!is.na(cur$rFilter[j])&j<nextEvent){ # if not immediately exceeding threshold, then find first time when it does exceed 
    while(cur$rFilter[j]<eventResponseR$threshold[i]&!is.na(cur$rFilter[j])&j<nextEvent){
      curLag<-curLag+1
      #       curR<-curR+cur$rFilter[j]
      #       curBack<-curBack+(cur$rQuad[j])
      #       curX<-curX+cur$rFilter[j]
      #       curMaxDis<-rbind(curMaxDis,cur$rFilter[j])
      j<-j+1
    }
  }
  while(cur$rFilter[j]>eventResponseR$threshold[i]&j<nextEvent&!is.na(cur$rFilter[j])){ #using filtered data to find duration of response 
    curR<-curR+cur$rFilter[j]+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[j] #filtered data back-calculated to get real metabolism values 
    curDuration<-curDuration+1
    curBack<-curBack+curMean # background respiration 
    curX<-curX+(cur$rFilter[j]+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[j]-curMean) # excess respiration over background respiration
    curMaxDis<-rbind(curMaxDis,(cur$rFilter[j]+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[j])) # maximum displacemnent 
    j<-j+1
  }
  eventResponseR$Rtotal[i]<-curR
  eventResponseR$duration[i]<-curDuration 
  eventResponseR$lag[i]<-curLag
  eventResponseR$backgroundR[i]<-curBack
  eventResponseR$excessR[i]<-curX
  if(curDuration==0){
    eventResponseR$maxDisplace[i]<-0
  }else{
    eventResponseR$maxDisplace[i]<-max(curMaxDis,na.rm=T)-curMean 
  }
}

eventResponseR


# gpp 
eventResponseGPP<-data.frame(dateTime=eventsNew,threshold=rep(NA,length(eventsNew)),
                             GPPtotal=rep(NA,length(eventsNew)),backgroundGPP=rep(NA,length(eventsNew)),
                             excessGPP=rep(NA,length(eventsNew)),maxDisplace=rep(NA,length(eventsNew)),
                             duration=rep(NA,length(eventsNew)),lag=rep(NA,length(eventsNew)))
window<-5
for(i in 1:length(eventsNew)){
  cur<-temp[strftime(temp$dateTime,'%Y')==strftime(eventsNew[i],'%Y'),] #have to split up data by years 
  eventResponseGPP$threshold[i]<-mean(cur$gppFilter[DOY(cur$dateTime)%in%
                                                      (DOY(eventResponseGPP$dateTime[i])-window-1):(DOY(eventResponseGPP$dateTime[i])-1)],na.rm=T)+
    sd(cur$gppFilter[DOY(cur$dateTime)%in%
                       (DOY(eventResponseGPP$dateTime[i])-window-1):(DOY(eventResponseGPP$dateTime[i])-1)],na.rm=T)*2
  curMean<-mean(cur$gppFilter[DOY(cur$dateTime)%in%
                                (DOY(eventResponseGPP$dateTime[i])-window-1):(DOY(eventResponseGPP$dateTime[i])-1)],na.rm=T)
  curMean<-curMean+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseGPP$dateTime[i])] 
  # docLoad after event - exceeding threshold 
  curGPP<-0
  curDuration<-0
  curLag<-0
  curBack<-0
  curX<-0
  curMaxDis<-c()
  curThresh<-eventResponseGPP$threshold[i]+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseGPP$dateTime[i])] # threshold that is back calculated to real values (Rweighted values )
  nextEvent<-min((length(cur$dateTime)+1),grep(eventResponseGPP$dateTime[i+1],cur$dateTime),na.rm = T)
  j<-grep(eventResponseGPP$dateTime[i],cur$dateTime) #variable response start time due to lag in response...
  if(cur$gppFilter[j]<eventResponseGPP$threshold[i]&!is.na(cur$gppFilter[j])){ # if not immediately exceeding threshold, then find first time when it does exceed 
    while(cur$gppFilter[j]<eventResponseGPP$threshold[i]&j<nextEvent&!is.na(cur$gppFilter[j])){
      curLag<-curLag+1
      #       curGPP<-curGPP+cur$gppFilter[j]
      #       curBack<-curBack+(cur$gppQuad[j]+eventResponseGPP$threshold[i])
      #       curX<-curX+cur$gppFilter[j]
      #       curMaxDis<-rbind(curMaxDis,cur$gppFilter[j])
      j<-j+1
    }
  }
  while(cur$gppFilter[j]>eventResponseGPP$threshold[i]&j<nextEvent&!is.na(cur$gppFilter[j])){ #using filtered data to find duration of response 
    curGPP<-curGPP+cur$gppFilter[j]+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[j] #filtered data back-calculated to get real metabolism values 
    curDuration<-curDuration+1
    curBack<-curBack+curMean
    curX<-curX+(cur$gppFilter[j]+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[j]-curMean) # excess respiration over background respiration
    curMaxDis<-rbind(curMaxDis,(cur$gppFilter[j]+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[j]))
    j<-j+1
  }
  eventResponseGPP$GPPtotal[i]<-curGPP
  eventResponseGPP$duration[i]<-curDuration 
  eventResponseGPP$lag[i]<-curLag
  eventResponseGPP$backgroundGPP[i]<-curBack
  eventResponseGPP$excessGPP[i]<-curX
  if(curDuration==0){
    eventResponseGPP$maxDisplace[i]<-0
  }else{
    eventResponseGPP$maxDisplace[i]<-max(curMaxDis,na.rm=T)-curMean
  }
}

eventResponseGPP


# doing a different NEP response. Using the max duration of a response from GPP and R, calculating excess 
#   R and GPP during that time and using excess GPP - excess R as NEP response 

eventResponseNEP<-data.frame(dateTime=eventsNew,excessNEP=rep(NA,length(eventsNew)),
                             duration=rep(NA,length(eventsNew)),lag=rep(NA,length(eventsNew)))
window<-5
for(i in 1:length(eventsNew)){
  cur<-temp[strftime(temp$dateTime,'%Y')==strftime(eventsNew[i],'%Y'),] #have to split up data by years 
  curMeanGPP<-mean(cur$gppFilter[DOY(cur$dateTime)%in%
                                   (DOY(eventResponseNEP$dateTime[i])-window-1):(DOY(eventResponseNEP$dateTime[i])-1)],na.rm=T)
  curMeanGPP<-curMeanGPP+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseNEP$dateTime[i])] 
  curMeanR<-mean(cur$rFilter[DOY(cur$dateTime)%in%
                               (DOY(eventResponseNEP$dateTime[i])-window-1):(DOY(eventResponseNEP$dateTime[i])-1)],na.rm=T)
  curMeanR<-curMeanR+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[as.POSIXct(cur$dateTime)==as.POSIXct(eventResponseNEP$dateTime[i])] 
  
  if(eventResponseGPP$duration[i]==0&eventResponseR$duration[i]==0){ 
    eventResponseNEP$excessNEP[i]=0
    eventResponseNEP$duration[i]=0
    eventResponseNEP$lag[i]=eventResponseGPP$lag[i]+eventResponseGPP$duration[i]
    next
  }
  responseStart<-min(eventResponseGPP$lag[i],eventResponseR$lag[i]) #start of response is the quickest respnose of either GPP or R (days after event Load)
  if(eventResponseGPP$duration[i]==0&eventResponseR$duration[i]>0){
    responseEnd<-eventResponseR$duration[i]+eventResponseR$lag[i]
  }else if(eventResponseGPP$duration[i]>0&eventResponseR$duration[i]==0){
    responseEnd<-eventResponseGPP$duration[i]+eventResponseGPP$lag[i]
  }else{
    responseEnd<-max((eventResponseGPP$duration[i]+eventResponseGPP$lag[i]),(eventResponseR$lag[i]+eventResponseR$duration[i])) #longest duration of repsonse added on to start = end of response 
  }
  #   if((responseEnd-responseStart)>30){
  #     responseEnd<-min(eventResponseGPP$duration[i],eventResponseR$duration[i])+responseStart # duration is too long, pick smaller of two durations 
  #   }
  eventResponseNEP$lag[i]<-responseStart
  j<-grep(eventResponseNEP$dateTime[i],cur$dateTime) #variable response start time due to lag in response...
  start<-j+responseStart
  end<-j+responseEnd
  eventResponseNEP$duration[i]<-length(na.omit(cur$gppFilter[start:end])) # removing NA metabolism for correct duration 
  excessGPP<-sum(cur$gppFilter[start:end]+min(cur$GPPweight-cur$gppQuad,na.rm=T)+cur$gppQuad[start:end],na.rm=T)-(eventResponseNEP$duration[i]*curMeanGPP)
  excessR<-sum(cur$rFilter[start:end]+min(cur$Rweight-cur$rQuad,na.rm=T)+cur$rQuad[start:end],na.rm=T)-(eventResponseNEP$duration[i]*curMeanR)
  
  eventResponseNEP$excessNEP[i]<-excessGPP-excessR
}

eventResponseNEP


# using window size of 5 days to find a threshold limit on the pulse 
eventLoad<-data.frame(dateTime=eventsNew,threshold=rep(NA,length(eventsNew)),docLoad=rep(NA,length(eventsNew)),
                      tpLoad=rep(NA,length(eventsNew)),
                      maxDisplace=rep(NA,length(eventsNew)),excessLoad=rep(NA,length(eventsNew)),
                      excessLoadTP=rep(NA,length(eventsNew)),
                      duration=rep(NA,length(eventsNew)),lag=rep(NA,length(eventsNew)))
window<-5
for(i in 1:length(eventsNew)){
  cur<-data[strftime(data$dateTime,'%Y')==strftime(eventsNew[i],'%Y'),] #have to split up data by years 
  eventLoad$threshold[i]<-mean(cur$totalAlloC[DOY(cur$dateTime)%in%
                                               (DOY(eventLoad$dateTime[i])-window-1):(DOY(eventLoad$dateTime[i])-1)])+
    sd(cur$totalAlloC[DOY(cur$dateTime)%in%
                       (DOY(eventLoad$dateTime[i])-window-1):(DOY(eventLoad$dateTime[i])-1)])*2
  curMean<-mean(cur$totalAlloC[DOY(cur$dateTime)%in%
                                (DOY(eventLoad$dateTime[i])-window-1):(DOY(eventLoad$dateTime[i])-1)])
  curMeanTP<-mean(cur$totalAlloP[DOY(cur$dateTime)%in%
                                   (DOY(eventLoad$dateTime[i])-window-1):(DOY(eventLoad$dateTime[i])-1)])
  
  # docLoad after event - exceeding threshold 
  curLoad<-0
  curLoadTP<-0
  curDuration<-0
  curLag<-0
  curX<-0
  curMaxDis<-c()
  j<-grep(eventLoad$dateTime[i],cur$dateTime)
  curThresh<-eventLoad$threshold[i] # threshold 
  nextEvent<-min((length(cur$dateTime)+1),grep(eventLoad$dateTime[i+1],cur$dateTime),na.rm = T)
  if(cur$totalAlloC[j]<eventLoad$threshold[i]&j<nextEvent){ # if not immediately exceeding threshold, then find first time when it does exceed 
    while(cur$totalAlloC[j]<eventLoad$threshold[i]&j<nextEvent){
      curLag<-curLag+1
      #       curLoad<-curLoad+cur$totalAlloC[j]
      #       curMaxDis<-rbind(curMaxDis,cur$totalAlloC[j])
      j<-j+1
    }
  }
  while(cur$totalAlloC[j]>eventLoad$threshold[i]&j<(length(cur$dateTime)+1)&j<nextEvent){ 
    curLoad<-curLoad+cur$totalAlloC[j]
    curLoadTP<-curLoadTP+cur$totalAlloP[j]
    curDuration<-curDuration+1
    curMaxDis<-rbind(curMaxDis,cur$totalAlloC[j])
    j<-j+1
  }
  eventLoad$docLoad[i]<-curLoad
  eventLoad$tpLoad[i]<-curLoadTP
  eventLoad$duration[i]<-curDuration 
  eventLoad$lag[i]<-curLag
  eventLoad$excessLoad[i]<-curLoad-(curMean*curDuration)
  eventLoad$excessLoadTP[i]<-curLoadTP-(curMeanTP*curDuration)
  eventLoad$maxDisplace[i]<-max(curMaxDis,na.rm=T)-curMean
}

eventLoad

# mass balance of all events; 2015-11-02 

library(LakeMetabolizer)
zmix<-data.frame()
for(i in 1:length(years)){
  cur<-read.table(file.path('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Data/zMix Data/',years[i],lake,'TEMP_PROF.txt'),
                  stringsAsFactor=F,header=T,sep='\t')
  if(length(cur[1,grep('Outlet',colnames(cur))])>0){
    cur<-cur[,-grep('Outlet',colnames(cur))]
  }
  cols<-colnames(cur)
  cols<-gsub('temp','wtr_',cols)
  colnames(cur)<-cols
  curZ<-ts.thermo.depth(cur,na.rm=T,seasonal=F)
  curZ<-aggregate(curZ$thermo.depth,by = list(as.Date(curZ$datetime)),FUN=max,na.rm=T)
  colnames(curZ)<-c('dateTime','thermo.depth')
  curZ$thermo.depth<-ifelse(curZ$thermo.depth=='-Inf',max(get.offsets(cur)),curZ$thermo.depth) # if INF, set thermo to deepest temp chain depth 
  zmix<-rbind(zmix,curZ)
}

zmix<-zmix[!duplicated(zmix$dateTime),]
zmix$dateTime<-as.character(zmix$dateTime)

data<-merge(data,zmix,by='dateTime',all.x=T)
source('/Users/Jake/Desktop/R functions/layer.volume.R')

# duration based on max NEP or load response across all lakes, except for the event on 2014-06-27 because EL had the 
#  longest response which was 59 days - seems unreasonable 

massBalance<-data.frame(dateTime=c('2013-06-20','2013-07-07','2013-08-26','2014-06-27','2014-08-25','2014-09-10'),
                        duration=c(17,22,20,26,16,6),waterRes=rep(NA,6),waterResOut=rep(NA,6),
                        docRetained=rep(NA,6),docIn=rep(NA,6),docOut=rep(NA,6),nepTotal=rep(NA,6),
                        mixVolume=rep(NA,6))
for(i in 1:length(massBalance$dateTime)){
  # load in data from event 
  start<-as.Date(massBalance$dateTime[i]) #start of event 
  end<-massBalance$duration[i]+start-1 #end event 
  
  cur<-data[as.Date(data$dateTime)>=start&as.Date(data$dateTime)<=end,] #pull out data just for event 
  for(j in 1:length(cur$dateTime)){
    cur$mixVol[j]<-layer.volume(bthA = as.numeric(bathy$area_m2),bthD=as.numeric(bathy$depth_m),top=0,bot=cur$thermo.depth[j]) # volume of water (m3) in mixed layer
  }
  cur$loadTotal<-cur$totalAlloC*lakeVol #converting from mg C m-3 to mg C (using total lake volume to convert back to total load)
  cur$nepTotal<--cur$NEPweight*1000/32*12*cur$mixVol #converting from mg O2 L-1 to mg C and making positive 
  # calculating water residence time 
  # GWin+STin+Precip=Evap+STout+GWout
  cur$gwDisch # m3 
  cur$streamWaterdisch # m3 
  cur$precip #m3 
  cur$evap #m3
  cur$waterHeight_m # m 
  if(lake=='EL'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$overlandDOCdisch
    docOut<-cur$gwDischDOC+(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='MO'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$gwDisch+cur$overland
    waterOut<-cur$evap+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$gwDischDOC+cur$overlandDOCdisch
    docOut<-(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='CR'){
    waterIn<-cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$streamWaterdisch
    docIn<-cur$totalAlloC*lakeVol
    docOut<-cur$gwDischDOC+cur$streamDOCdisch
  }
  waterRes<-lakeVol/mean(waterIn) # m3 / m3 day-1 = days 
  waterResOut<-lakeVol/mean(waterOut)
  docRetained<-(sum(docIn)-sum(docOut))/sum(docIn) # fraction of DOC retained in lake 
  massBalance$waterRes[i]<-waterRes
  massBalance$docRetained[i]<-docRetained
  massBalance$docIn[i]<-sum(docIn)
  massBalance$docOut[i]<-sum(docOut)
  massBalance$nepTotal[i]<-sum(cur$nepTotal)
  massBalance$waterResOut[i]<-waterResOut
  massBalance$mixVolume[i]<-mean(cur$mixVol,na.rm=T)
}

massBalance
massBalance$docInDaily<-massBalance$docIn/massBalance$duration
massBalance$nepDaily<-massBalance$nepTotal/massBalance$duration
plot(massBalance$docRetained~massBalance$waterResOut)
plot(massBalance$docInDaily~massBalance$waterRes)
plot(massBalance$nepDaily~massBalance$docRetained)
plot(massBalance$nepDaily~massBalance$waterRes)

# water residence time using variable duration of events ( not fixed as in mass balance calcs )
eventLoad$waterRes<-rep(NA,length(eventLoad$dateTime))
eventLoad$waterResOut<-rep(NA,length(eventLoad$dateTime))
eventResponseNEP$waterRes<-rep(NA,length(eventResponseNEP$dateTime))
eventResponseNEP$waterResOut<-rep(NA,length(eventResponseNEP$dateTime))

for(i in 1:length(eventLoad$dateTime)){
  # load in data from event 
  start<-as.Date(eventLoad$dateTime[i]) #start of event 
  end<-eventLoad$duration[i]+start-1 #end event 
  
  cur<-data[as.Date(data$dateTime)>=start&as.Date(data$dateTime)<=end,] #pull out data just for event 
  
  if(lake=='EL'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$overlandDOCdisch
    docOut<-cur$gwDischDOC+(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='MO'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$gwDisch+cur$overland
    waterOut<-cur$evap+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$gwDischDOC+cur$overlandDOCdisch
    docOut<-(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='CR'){
    waterIn<-cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$streamWaterdisch
    docIn<-cur$totalAlloC*lakeVol
    docOut<-cur$gwDischDOC+cur$streamDOCdisch
  }
  waterRes<-lakeVol/mean(waterIn) # m3 / m3 day-1 = days 
  waterResOut<-lakeVol/mean(waterOut)
  eventLoad$waterRes[i]<-waterRes
  eventLoad$waterResOut[i]<-waterResOut
  
  ## using duration of NEP response ## 
  start<-as.Date(eventResponseNEP$dateTime[i]) #start of event 
  end<-eventResponseNEP$duration[i]+eventResponseNEP$lag[i]+start-2 #end event 
  
  cur<-data[as.Date(data$dateTime)>=start&as.Date(data$dateTime)<=end,] #pull out data just for event 
  
  if(lake=='EL'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$overlandDOCdisch
    docOut<-cur$gwDischDOC+(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='MO'){
    waterIn<-cur$streamWaterdisch+cur$precip+cur$gwDisch+cur$overland
    waterOut<-cur$evap+cur$outletDischarge
    docIn<-cur$streamDOCdisch+cur$precipDOC+cur$gwDischDOC+cur$overlandDOCdisch
    docOut<-(cur$outletDischarge*cur$doc*1000)
  }else if(lake=='CR'){
    waterIn<-cur$precip+cur$overland
    waterOut<-cur$evap+cur$gwDisch+cur$streamWaterdisch
    docIn<-cur$totalAlloC*lakeVol
    docOut<-cur$gwDischDOC+cur$streamDOCdisch
  }
  waterRes<-lakeVol/mean(waterIn) # m3 / m3 day-1 = days 
  waterResOut<-lakeVol/mean(waterOut)
  eventResponseNEP$waterRes[i]<-waterRes
  eventResponseNEP$waterResOut[i]<-waterResOut
}




# load in all lake data R datafile for figures and analyses 
load('C:/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/R Data/allLakesData_20160818.RData')


# Figure 1 *********************************************


png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_1a.png',  # example of DOC load and pulse detection 
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=6
ylim=c(0.5,3)
par(mar=c(5,6,4,2))
plot(elData$totalAlloC[strftime(elData$dateTime,'%Y')==2013&!is.na(elData$Rweight)]~
       as.POSIXct(elData$dateTime[strftime(elData$dateTime,'%Y')==2013&!is.na(elData$Rweight)]),
     lwd=lwd,xlim=c(as.POSIXct('2013-06-17'),as.POSIXct('2013-09-14')),type='l',xlab='',
     ylab=expression(t-DOC~Load~(mg~C~m^-3~day^-1)),cex.axis=cex,cex.lab=cex,bty='o')
abline(v=as.POSIXct(eventsNew[2]),col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct(eventsNew[2]),y0 = elEventLoad$threshold[2],x1 =as.POSIXct('2013-07-20'),
         y1=elEventLoad$threshold[2],col='gray40',lty=1,lwd=4)
abline(v=as.POSIXct(eventsNew[1]),col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct(eventsNew[1]),y0 = elEventLoad$threshold[1],x1 = as.POSIXct('2013-06-24'),
         y1=elEventLoad$threshold[1],col='gray40',lty=1,lwd=4)
segments(x0=as.POSIXct(eventsNew[3]),y0=-100,x1=as.POSIXct(eventsNew[3]),y1=680,
         col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct(eventsNew[3]),y0 = elEventLoad$threshold[3]-35,x1 = as.POSIXct('2013-08-30'),
         y1=elEventLoad$threshold[3]-35,col='grey40',lty=1,lwd=4)
legend("topright", legend=c("DOC Load","Precipitation Event",'DOC Load Threshold'),
       col=c('black','gray40','grey40'),pt.bg=c('black','gray40','grey40'), 
       ncol=1,lwd=c(lwd,4,4),bty='n',lty=c(1,6,1))
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_1b.png',  # showing raw metab data with sd as filled in error
    res=300, width=7, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=3
ylim=c(0.5,5.5)
par(mar=c(5,6,4,2))
cur<-elData[strftime(elData$dateTime,'%Y')==2013,]
cur<-cur[!is.na(cur$R),]
plot(cur$R~as.POSIXct(cur$dateTime),
     lwd=lwd,bty='o',ylab=expression(R~(mg~O[2]~L^-1~day^-1)),xlab='',ylim=ylim,
     cex.axis=cex,cex.lab=cex,type='l',xlim=c(as.POSIXct('2013-06-17'),as.POSIXct('2013-09-14')))
lower<-elData$R[strftime(elData$dateTime,'%Y')==2013]-elData$sdR[strftime(elData$dateTime,'%Y')==2013]
upper<-elData$R[strftime(elData$dateTime,'%Y')==2013]+elData$sdR[strftime(elData$dateTime,'%Y')==2013]
x<-c(elData$dateTime[strftime(elData$dateTime,'%Y')==2013],rev(elData$dateTime[strftime(elData$dateTime,'%Y')==2013]))
y<-c(upper,rev(lower))
xy<-na.omit(data.frame(x=x,y=y))
polygon(as.POSIXct(xy$x),xy$y,col='gray80',border = 'grey20',lwd=2)
lines(cur$R~as.POSIXct(cur$dateTime),
      lwd=lwd,bty='o',ylab=expression(R~(mg~O[2]~L^-1~day^-1)),xlab='',ylim=ylim,
      cex.axis=cex,cex.lab=cex,type='l',xlim=c(as.POSIXct('2013-06-17'),as.POSIXct('2013-09-14')))
legend("topright", legend=c("R",'Standard Dev'),col=c('black','gray60'),
       pt.bg=c('black','gray60'), ncol=1,lwd=c(lwd,4),bty='n',lty=c(1,1))
dev.off()


png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_1c&d.png', 
    res=300, width=14, height=7, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=6
ylim=c(0.5,3)
par(mar=c(5,6,4,2),mfrow=c(1,2))
plot(elData$Rweight[strftime(elData$dateTime,'%Y')==2013]~as.POSIXct(elData$dateTime[strftime(elData$dateTime,'%Y')==2013]),
     lwd=lwd,bty='o',ylab=expression(Smoothed~R~(mg~O[2]~L^-1~day^-1)),xlab='',
     cex.axis=cex,cex.lab=cex,type='l',xlim=c(as.POSIXct('2013-06-17'),as.POSIXct('2013-09-14')))
lines(elTemp$rQuad[strftime(elTemp$dateTime,'%Y')==2013]~as.POSIXct(elTemp$dateTime[strftime(elTemp$dateTime,'%Y')==2013]),
      lty=6,lwd=4,col='grey40')
legend("topright", legend=c("Smoothed R","Seasonal Trend "),col=c('black','gray40'),
       pt.bg=c('black','gray40'), ncol=1,lwd=c(lwd,4),bty='n',lty=c(1,6))
plot(elTemp$rFilter[strftime(elTemp$dateTime,'%Y')==2013]~as.POSIXct(elTemp$dateTime[strftime(elTemp$dateTime,'%Y')==2013]),
     lwd=lwd,xlim=c(as.POSIXct('2013-06-17'),as.POSIXct('2013-09-13')),type='l',xlab='',
     ylab=expression(Detrended~R~(mg~O[2]~L^-1~day^-1)),cex.axis=cex,cex.lab=cex,bty='o')
abline(v=as.POSIXct(eventsNew[2]),col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct('2013-07-11'),y0 = elEventResponseR$threshold[2],x1 =as.POSIXct('2013-07-29'),
         y1=elEventResponseR$threshold[2],col='grey40',lty=1,lwd=4)
abline(v=as.POSIXct(eventsNew[1]),col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct(eventsNew[1]),y0 = elEventResponseR$threshold[1],x1 = as.POSIXct(eventsNew[2]),
         y1=elEventResponseR$threshold[1],col='grey40',lty=1,lwd=4)
segments(x0=as.POSIXct(eventsNew[3]),y0=-1,x1=as.POSIXct(eventsNew[3]),y1=1.4,
         col='grey40',lwd=4,lty=2)
segments(x0 = as.POSIXct(eventsNew[3]),y0 = elEventResponseR$threshold[3],x1 = as.POSIXct('2013-09-16'),
         y1=elEventResponseR$threshold[3],col='grey40',lty=1,lwd=4)
legend("topright", legend=c("Detrended R","Precipitation Event",'Respiration Threshold'),
       col=c('black','gray40','grey40'),pt.bg=c('black','gray40','grey40'), 
       ncol=1,lwd=c(lwd,4,4),bty='n',lty=c(1,6,1))
dev.off()

# Figure 2 ******************************************************

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_2a_EL.png', 
    res=300, width=15, height=5, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=6
ylim=c(0,max(moData$totalAlloC))
par(mar=c(5,6,4,2))
plot(elData$totalAlloC,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(DOC~Load~(mg~C~m^-3~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
years<-c(2013,2014)
lakes<-c('el')
for(i in 1:length(years)){
  for(j in 1:length(lakes)){
    cur<-eval(parse(text=(paste(lakes[j],'Data',sep=''))))
    cur<-cur[cur$dateTime!='2014-05-10',]
    Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    if(j==1){
      lines(Xindex,cur$totalAlloC[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
            type='l',lwd=lwd,xaxt='n')
    }else if(j==2){
      lines(Xindex,cur$totalAlloC[Yindex],col='grey20',lwd=lwd,lty=6)
    }else if(j==3){
      lines(Xindex,cur$totalAlloC[Yindex],col='gray40',lwd=lwd,lty=1)
    }
  }
}

for(i in 1:length(eventsNew)){
  curDur<-elEventLoad$duration[elEventLoad$dateTime==eventsNew[i]]
  curLag<-elEventLoad$lag[elEventLoad$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur
  x<-c(which(as.Date(elData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  if(as.Date(eventsNew[i])>as.Date('2014-04-01')){
    x<-x-1
  }
  y<-c(max(moData$totalAlloC)+50,max(moData$totalAlloC)+50,max(moData$totalAlloC)+200,max(moData$totalAlloC)+200)
  polygon(x = x,y=y,col = 'gray40',border = 'black')
}
abline(v=which(elData$dateTime%in%eventsNew[1:3]),col='grey40',lwd=2,lty=5)
abline(v=which(elData$dateTime%in%eventsNew[4:6])-1,col='grey40',lwd=2,lty=5)

dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_2b_MO.png', 
    res=300, width=15, height=5, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=6
ylim=c(0,max(moData$totalAlloC))
par(mar=c(5,6,4,2))
plot(moData$totalAlloC,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(DOC~Load~(mg~C~m^-3~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
years<-c(2013,2014)
lakes<-c('mo')
for(i in 1:length(years)){
  for(j in 1:length(lakes)){
    cur<-eval(parse(text=(paste(lakes[j],'Data',sep=''))))
    cur<-cur[cur$dateTime!='2014-05-10',]
    Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    if(j==1){
      lines(Xindex,cur$totalAlloC[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
            type='l',lwd=lwd,xaxt='n')
    }else if(j==2){
      lines(Xindex,cur$totalAlloC[Yindex],col='grey20',lwd=lwd,lty=6)
    }else if(j==3){
      lines(Xindex,cur$totalAlloC[Yindex],col='gray40',lwd=lwd,lty=1)
    }
  }
}

for(i in 1:length(eventsNew)){
  curDur<-moEventLoad$duration[moEventLoad$dateTime==eventsNew[i]]
  curLag<-moEventLoad$lag[moEventLoad$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur
  x<-c(which(as.Date(moData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  if(as.Date(eventsNew[i])>as.Date('2014-04-01')){
    x<-x-1
  }
  y<-c(max(moData$totalAlloC)+50,max(moData$totalAlloC)+50,max(moData$totalAlloC)+200,max(moData$totalAlloC)+200)
  polygon(x = x,y=y,col = 'gray40',border = 'black')
}
abline(v=which(moData$dateTime%in%eventsNew[1:3]),col='grey40',lwd=2,lty=5)
abline(v=which(elData$dateTime%in%eventsNew[4:6])-1,col='grey40',lwd=2,lty=5)

dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_2c_CR.png', 
    res=300, width=15, height=5, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2
lwd=6
ylim=c(0,max(moData$totalAlloC))
par(mar=c(5,6,4,2))
plot(crData$totalAlloC,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(DOC~Load~(mg~C~m^-3~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
years<-c(2013,2014)
lakes<-c('cr')
for(i in 1:length(years)){
  for(j in 1:length(lakes)){
    cur<-eval(parse(text=(paste(lakes[j],'Data',sep=''))))
    cur<-cur[cur$dateTime!='2014-05-10',]
    Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
    if(j==1){
      lines(Xindex,cur$totalAlloC[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
            type='l',lwd=lwd,xaxt='n')
    }else if(j==2){
      lines(Xindex,cur$totalAlloC[Yindex],col='grey20',lwd=lwd,lty=6)
    }else if(j==3){
      lines(Xindex,cur$totalAlloC[Yindex],col='gray40',lwd=lwd,lty=1)
    }
  }
}

for(i in 1:length(eventsNew)){
  curDur<-crEventLoad$duration[crEventLoad$dateTime==eventsNew[i]]
  curLag<-crEventLoad$lag[crEventLoad$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur
  x<-c(which(as.Date(crData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  if(as.Date(eventsNew[i])>as.Date('2014-04-01')){
    x<-x-1
  }
  y<-c(max(moData$totalAlloC)+50,max(moData$totalAlloC)+50,max(moData$totalAlloC)+200,max(moData$totalAlloC)+200)
  polygon(x = x,y=y,col = 'gray40',border = 'black')
}

abline(v=which(crData$dateTime%in%eventsNew[1:3]),col='grey40',lwd=2,lty=5)
abline(v=which(elData$dateTime%in%eventsNew[4:6])-1,col='grey40',lwd=2,lty=5)

dev.off()

# Figure 3 ***********************************************************

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_3a_elMetab.png', # metabolism time series 
    res=300, width=20, height=10, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2.5
lwd=6
ylim=c(-max(moData$Rweight),max(moData$GPPweight))
par(mar=c(5,6,4,2))
plot(elData$Rweight,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(Metabolism~(mg~O[2]~L^-1~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
years<-c(2013,2014)
for(i in 1:length(years)){
  cur<-elData
  cur<-cur[cur$dateTime!='2014-05-10',]
  Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  lines(Xindex,-cur$Rweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n',xlab='',ylab='')
  lines(Xindex,cur$GPPweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
  lines(Xindex,cur$NEPweight[Yindex],ylim=ylim,col='gray60',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
}
for(i in 1:length(eventsNew)){
  curDur<-elEventResponseNEP$duration[elEventResponseNEP$dateTime==eventsNew[i]]
  curLag<-elEventResponseNEP$lag[elEventResponseNEP$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur-2
  x<-c(which(as.Date(elData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  y<-c(ylim[1]+ylim[1]*.01,ylim[1]+ylim[1]*.01,ylim[1]+ylim[1],ylim[1]+ylim[1])
  polygon(x = x,y=y,col = 'gray80',border = 'black')
}
box()
abline(0,0,lty=2,lwd=2,col='black')
abline(v=which(elData$dateTime%in%eventsNew[1:6]),col='grey40',lwd=2,lty=5)

dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_3b_moMetab.png', # metabolism time series 
    res=300, width=20, height=10, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2.5
lwd=6
ylim=c(-max(moData$Rweight),max(moData$GPPweight))
par(mar=c(5,6,4,2))
plot(moData$Rweight,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(Metabolism~(mg~O[2]~L^-1~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
legend("top", legend=c("GPP & R",'NEP'),col=c('black','gray60'),
       pt.bg=c('black','grey60'), ncol=1,cex=cex,bty='n',lwd=lwd,lty=c(1,1))
years<-c(2013,2014)
for(i in 1:length(years)){
  cur<-moData
  cur<-cur[cur$dateTime!='2014-05-10',]
  Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  lines(Xindex,-cur$Rweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n',xlab='',ylab='')
  lines(Xindex,cur$GPPweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
  lines(Xindex,cur$NEPweight[Yindex],ylim=ylim,col='gray60',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
}
for(i in 1:length(eventsNew)){
  curDur<-moEventResponseNEP$duration[moEventResponseNEP$dateTime==eventsNew[i]]
  curLag<-moEventResponseNEP$lag[moEventResponseNEP$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur-2
  x<-c(which(as.Date(moData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  y<-c(ylim[1]+ylim[1]*.01,ylim[1]+ylim[1]*.01,ylim[1]+ylim[1],ylim[1]+ylim[1])
  polygon(x = x,y=y,col = 'gray80',border = 'black')
}
box()
abline(0,0,lty=2,lwd=2,col='black')

abline(v=which(moData$dateTime%in%eventsNew[1:6]),col='grey40',lwd=2,lty=5)

dev.off()


png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_3b_crMetab.png', # metabolism time series 
    res=300, width=20, height=10, units = 'in')
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels

cex=2.5
lwd=6
ylim=c(-max(moData$Rweight),max(moData$GPPweight))
par(mar=c(5,6,4,2))
plot(crData$Rweight,ylim=ylim,xaxt='n',
     lwd=0,xlab='',ylab=expression(Metabolism~(mg~O[2]~L^-1~day^-1)),
     cex.axis=cex,cex.lab=cex,bty='o',cex=0)
axis(1,at=c(50,150),labels=c('2013','2014'),cex.axis=cex)
years<-c(2013,2014)
for(i in 1:length(years)){
  cur<-crData
  cur<-cur[cur$dateTime!='2014-05-10',]
  Yindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  Xindex<-which(strftime(cur$dateTime,'%Y')==years[i])
  lines(Xindex,-cur$Rweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n',xlab='',ylab='')
  lines(Xindex,cur$GPPweight[Yindex],ylim=ylim,col='black',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
  lines(Xindex,cur$NEPweight[Yindex],ylim=ylim,col='gray60',cex.axis=cex,cex.lab=cex,
        type='l',lwd=lwd,xaxt='n')
}
for(i in 1:length(eventsNew)){
  curDur<-crEventResponseNEP$duration[crEventResponseNEP$dateTime==eventsNew[i]]
  curLag<-crEventResponseNEP$lag[crEventResponseNEP$dateTime==eventsNew[i]]
  start<-as.Date(eventsNew[i])+curLag
  end<-as.Date(eventsNew[i])+curLag+curDur-2
  x<-c(which(as.Date(crData$dateTime)%in%c(start,end)))
  x<-c(x,rev(x))
  y<-c(ylim[1]+ylim[1]*.01,ylim[1]+ylim[1]*.01,ylim[1]+ylim[1],ylim[1]+ylim[1])
  polygon(x = x,y=y,col = 'gray80',border = 'black')
}
box()
abline(0,0,lty=2,lwd=2,col='black')

abline(v=which(crData$dateTime%in%eventsNew[1:6]),col='grey40',lwd=2,lty=5)

dev.off()

# Figure 4 *******************************************************

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_4a_load_heterotrophy.png',  
    res=300, width=7, height=7, units = 'in')
cex=2
lwd=6
ylim=c(-1000,18000)
par(mar=c(5,6,4,2))
plot(mean(elEventLoad$excessLoad)~elMeanResOut,ylim=ylim,xlim=c(0,800),pch=16,cex=2,cex.lab=2,cex.axis=2,
     ylab='t-DOC Load or Heterotrophy',xlab='Water Residence Time (days)')
points(mean(moEventLoad$excessLoad)~moMeanResOut,pch=18,cex=2)
points(mean(crEventLoad$excessLoad)~crMeanResOut,pch=17,cex=2)
arrows(elMeanResOut,quantile(elEventLoad$excessLoad,0.05),elMeanResOut,
       quantile(elEventLoad$excessLoad,0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(moMeanResOut,quantile(moEventLoad$excessLoad,0.05),moMeanResOut,
       quantile(moEventLoad$excessLoad,0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(crMeanResOut,quantile(crEventLoad$excessLoad,0.05),crMeanResOut,
       quantile(crEventLoad$excessLoad,0.95),lwd=2,
       length=0.1, angle=90, code=3)
points(-mean(elEventResponseNEP$excessNEPc)~elMeanResOut,pch=16,cex=2,col='grey60')
points(-mean(moEventResponseNEP$excessNEPc)~moMeanResOut,pch=18,cex=2,col='grey60')
points(-mean(crEventResponseNEP$excessNEPc)~crMeanResOut,pch=17,cex=2,col='grey60')
arrows(elMeanResOut,-quantile(elEventResponseNEP$excessNEPc,0.05),elMeanResOut,
       -quantile(elEventResponseNEP$excessNEPc,0.95),lwd=2,
       length=0.1, angle=90, code=3,col='grey60')
arrows(moMeanResOut,-quantile(moEventResponseNEP$excessNEPc,0.05),moMeanResOut,
       -quantile(moEventResponseNEP$excessNEPc,0.95),lwd=2,
       length=0.1, angle=90, code=3,col='grey60')
arrows(crMeanResOut,-quantile(crEventResponseNEP$excessNEPc,0.05),crMeanResOut,
       -quantile(crEventResponseNEP$excessNEPc,0.95),lwd=2,
       length=0.1, angle=90, code=3,col='grey60')
legend("topright", legend=c("t-DOC Load","Heterotrophy"),col=c('black','gray60'),pt.cex = c(2,2),pch=c(16,16),
       pt.bg=c('black','gray60'), ncol=1,bty='n')
dev.off()


png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_4b_fractionMineralized.png',  
    res=300, width=7, height=7, units = 'in')
cex=2
lwd=6
ylim=c(0,1)
par(mar=c(5,6,4,2))
plot(mean(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='EL'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1])~elMeanResOut,
     ylim=ylim,xlim=c(0,800),pch=16,cex=2,cex.lab=2,cex.axis=2,
     ylab='Fraction Mineralized',xlab='Hydrologic Residence Time (days)')
points(mean(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='MO'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1])~moMeanResOut,pch=18,cex=2)
points(mean(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='CR'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1])~crMeanResOut,pch=17,cex=2)
arrows(elMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='EL'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.05),
       elMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='EL'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(moMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='MO'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.05),
       moMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='MO'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(crMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='CR'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.05),
       crMeanResOut,
       quantile(allLakeResponse$fractionProcessed[allLakeResponse$lakeID=='CR'&allLakeResponse$fractionProcessed>0&allLakeResponse$fractionProcessed<1],0.95),lwd=2,
       length=0.1, angle=90, code=3)
dev.off()

png('/Users/Jake/Documents/Jake/MyPapers/Pulses and Lags/Figures/fig_4c_tDOCturnover.png',  
    res=300, width=7, height=7, units = 'in')
cex=2
lwd=6
ylim=c(0.0001,.023)
par(mar=c(5,6,4,2))
plot(mean(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='EL'&allLakeResponse$DOCturnover>0])~elMeanResOut,
     ylim=ylim,xlim=c(0,800),pch=16,cex=2,cex.lab=2,cex.axis=2,
     ylab='t-DOC Turnover',xlab='Hydrologic Residence Time (days)')
points(mean(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='MO'&allLakeResponse$DOCturnover>0])~moMeanResOut,pch=18,cex=2)
points(mean(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='CR'&allLakeResponse$DOCturnover>0])~crMeanResOut,pch=17,cex=2)
arrows(elMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='EL'&allLakeResponse$DOCturnover>0],0.05),
       elMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='EL'&allLakeResponse$DOCturnover>0],0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(moMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='MO'&allLakeResponse$DOCturnover>0],0.05),
       moMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='MO'&allLakeResponse$DOCturnover>0],0.95),lwd=2,
       length=0.1, angle=90, code=3)
arrows(crMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='CR'&allLakeResponse$DOCturnover>0],0.05),
       crMeanResOut,
       quantile(allLakeResponse$DOCturnover[allLakeResponse$lakeID=='CR'&allLakeResponse$DOCturnover>0],0.95),lwd=2,
       length=0.1, angle=90, code=3)
dev.off()


# Table 2 *****************************************************
# percent contribution of inflows to hydrologic and DOC loads 
# hydrologic contribution 
sum(elData$streamWaterdisch)/sum(elData$streamWaterdisch+elData$precip+elData$overland)
sum(elData$precip)/sum(elData$streamWaterdisch+elData$precip+elData$overland)
sum(elData$overland)/sum(elData$streamWaterdisch+elData$precip+elData$overland)

sum(moData$streamWaterdisch)/sum(moData$streamWaterdisch+moData$precip+moData$gwDisch+moData$overland)
sum(moData$precip)/sum(moData$streamWaterdisch+moData$precip+moData$gwDisch+moData$overland)
sum(moData$gwDisch)/sum(moData$streamWaterdisch+moData$precip+moData$gwDisch+moData$overland)
sum(moData$overland)/sum(moData$streamWaterdisch+moData$precip+moData$gwDisch+moData$overland)

sum(crData$overland)/sum(crData$overland+crData$precip)
sum(crData$precip)/sum(crData$precip+crData$overland)

# carbon contribution 
sum(elData$streamDOCdisch)/sum(elData$streamDOCdisch+elData$precipDOC+elData$overlandDOCdisch)
sum(elData$precipDOC)/sum(elData$streamDOCdisch+elData$precipDOC+elData$overlandDOCdisch)
sum(elData$overlandDOCdisch)/sum(elData$streamDOCdisch+elData$precipDOC+elData$overlandDOCdisch)

sum(moData$streamDOCdisch)/sum(moData$streamDOCdisch+moData$precipDOC+moData$gwDischDOC+moData$overlandDOCdisch)
sum(moData$precipDOC)/sum(moData$streamDOCdisch+moData$precipDOC+moData$gwDischDOC+moData$overlandDOCdisch)
sum(moData$gwDischDOC)/sum(moData$streamDOCdisch+moData$precipDOC+moData$gwDischDOC+moData$overlandDOCdisch)
sum(moData$overlandDOCdisch)/sum(moData$streamDOCdisch+moData$precipDOC+moData$gwDischDOC+moData$overlandDOCdisch)

sum(crData$overlandDOCdisch)/sum(crData$overlandDOCdisch+crData$precipDOC)
sum(crData$precipDOC)/sum(crData$precipDOC+crData$overlandDOCdisch)


# Table 3 **********************************************************************
# ccf of t-DOC load 
el_mo_load<-ccf(elData$totalAlloC,moData$totalAlloC,lag.max=20)
cr_mo_load<-ccf(crData$totalAlloC[!is.na(crData$totalAlloC)],moData$totalAlloC[!is.na(crData$totalAlloC)],lag.max=20)
el_cr_load<-ccf(elData$totalAlloC[!is.na(crData$totalAlloC)],crData$totalAlloC[!is.na(crData$totalAlloC)],lag.max=20)
# ccf of ecosystem R 
el_mo_r<-ccf(elData$Rweight,moData$Rweight,lag.max=20)
el_cr_r<-ccf(elData$Rweight,crData$Rweight,lag.max=20)
cr_mo_r<-ccf(crData$Rweight,moData$Rweight,lag.max=20)
# ccf of NEP 
el_mo_nep<-ccf(elData$NEPweight,moData$NEPweight,lag.max=20)
el_cr_nep<-ccf(elData$NEPweight,crData$NEPweight,lag.max=20)
cr_mo_nep<-ccf(crData$NEPweight,moData$NEPweight,lag.max=20)
# ccf of GPP
el_mo_gpp<-ccf(elData$GPPweight,moData$GPPweight,lag.max=20)
el_cr_gpp<-ccf(elData$GPPweight,crData$GPPweight,lag.max=20)
cr_mo_gpp<-ccf(crData$GPPweight,moData$GPPweight,lag.max=20)


# Table 4 *************************************************************
#t-DOC annual average 
mean(moData$totalAlloC)
mean(elData$totalAlloC)
mean(crData$totalAlloC)

# GPP 
mean(moData$GPPweight)*1000*12/32
mean(elData$GPPweight)*1000*12/32
mean(crData$GPPweight)*1000*12/32

# R
mean(moData$Rweight)*1000*12/32
mean(elData$Rweight)*1000*12/32
mean(crData$Rweight)*1000*12/32

# NEP
mean(moData$NEPweight)*1000*12/32
mean(elData$NEPweight)*1000*12/32
mean(crData$NEPweight)*1000*12/32

# HRT 
elMeanResOut<-123734/mean(elData$outletDischarge+elData$evap+elData$gwDisch)
crMeanResOut<-1302506/mean(crData$streamWaterdisch+crData$evap+crData$gwDisch)
moMeanResOut<-142454/mean(moData$streamWaterdisch+moData$precip+moData$gwDisch+moData$overland)

# t-DOC event excess load average and sd 
mean(moEventLoad$excessLoad)
mean(elEventLoad$excessLoad)
mean(crEventLoad$excessLoad)

sd(moEventLoad$excessLoad)
sd(elEventLoad$excessLoad)
sd(crEventLoad$excessLoad)

# excess heterotrophy event-1 average and sd 
mean(moEventResponseNEP$excessNEPc)
mean(elEventResponseNEP$excessNEPc)
mean(crEventResponseNEP$excessNEPc)

sd(moEventResponseNEP$excessNEPc)
sd(elEventResponseNEP$excessNEPc)
sd(crEventResponseNEP$excessNEPc)


#cv of heterotrophy
sd(-elEventResponseNEP$excessNEPc[elEventResponseNEP$flag==0&elEventResponseNEP$excessNEPc<0])/mean(-elEventResponseNEP$excessNEPc[elEventResponseNEP$flag==0&elEventResponseNEP$excessNEPc<0])
sd(-moEventResponseNEP$excessNEPc[moEventResponseNEP$excessNEPc<0])/mean(-moEventResponseNEP$excessNEPc[moEventResponseNEP$excessNEPc<0])
sd(-crEventResponseNEP$excessNEPc[crEventResponseNEP$excessNEPc<0])/mean(-crEventResponseNEP$excessNEPc[crEventResponseNEP$excessNEPc<0])

#cv of t-DOC load 
sd(elEventLoad$excessLoad[elEventResponseNEP$flag==0])/mean(elEventLoad$excessLoad[elEventResponseNEP$flag==0])
sd(moEventLoad$excessLoad)/mean(moEventLoad$excessLoad)
sd(crEventLoad$excessLoad)/mean(crEventLoad$excessLoad)

# extent of time series for t-DOC loads and contribution to t-DOC seasonal flux 
sum(moEventLoad$duration)/length(moData$dateTime)
sum(elEventLoad$duration)/length(elData$dateTime)
sum(crEventLoad$duration)/length(crData$dateTime)

sum(moEventLoad$docLoad)/sum(moData$totalAlloC)
sum(elEventLoad$docLoad)/sum(elData$totalAlloC)
sum(crEventLoad$docLoad)/sum(crData$totalAlloC)


# excess load comparison 
mean(moEventLoad$excessLoad/crEventLoad$excessLoad)
mean(moEventLoad$excessLoad/elEventLoad$excessLoad)
mean(elEventLoad$excessLoad/crEventLoad$excessLoad)

mean(moEventResponseNEP$excessNEP[moEventResponseNEP$excessNEP<0&crEventResponseNEP$excessNEP<0]/
       crEventResponseNEP$excessNEP[moEventResponseNEP$excessNEP<0&crEventResponseNEP$excessNEP<0])
mean(moEventResponseNEP$excessNEP[moEventResponseNEP$excessNEP<0&elEventResponseNEP$excessNEP<0]/
       elEventResponseNEP$excessNEP[moEventResponseNEP$excessNEP<0&elEventResponseNEP$excessNEP<0])
mean(elEventResponseNEP$excessNEP[elEventResponseNEP$excessNEP<0&crEventResponseNEP$excessNEP<0]/
       crEventResponseNEP$excessNEP[elEventResponseNEP$excessNEP<0&crEventResponseNEP$excessNEP<0])


# metabolic response compared to seasonal metabolic responses for each lake 
moSeas2013<-max(moTemp$rQuad[strftime(moTemp$dateTime,'%Y')==2013])-min(moTemp$rQuad[strftime(moTemp$dateTime,'%Y')==2013])
moSeas2014<-max(moTemp$rQuad[strftime(moTemp$dateTime,'%Y')==2014])-min(moTemp$rQuad[strftime(moTemp$dateTime,'%Y')==2014])
mean(c(moEventResponseR$maxDisplace[strftime(moEventResponseR$dateTime,'%Y')==2013]/moSeas2013,moEventResponseR$maxDisplace[strftime(moEventResponseR$dateTime,'%Y')==2014]/moSeas2014))

elSeas2013<-max(elTemp$rQuad[strftime(elTemp$dateTime,'%Y')==2013])-min(elTemp$rQuad[strftime(elTemp$dateTime,'%Y')==2013])
elSeas2014<-max(elTemp$rQuad[strftime(elTemp$dateTime,'%Y')==2014])-min(elTemp$rQuad[strftime(elTemp$dateTime,'%Y')==2014])
mean(c(elEventResponseR$maxDisplace[strftime(elEventResponseR$dateTime,'%Y')==2013]/elSeas2013,elEventResponseR$maxDisplace[strftime(elEventResponseR$dateTime,'%Y')==2014]/elSeas2014))

crSeas2013<-max(crTemp$rQuad[strftime(crTemp$dateTime,'%Y')==2013])-min(crTemp$rQuad[strftime(crTemp$dateTime,'%Y')==2013])
crSeas2014<-max(crTemp$rQuad[strftime(crTemp$dateTime,'%Y')==2014])-min(crTemp$rQuad[strftime(crTemp$dateTime,'%Y')==2014])
mean(c(crEventResponseR$maxDisplace[strftime(crEventResponseR$dateTime,'%Y')==2013]/crSeas2013,crEventResponseR$maxDisplace[strftime(crEventResponseR$dateTime,'%Y')==2014]/crSeas2014))

moSeas2013<-max(moTemp$gppQuad[strftime(moTemp$dateTime,'%Y')==2013])-min(moTemp$gppQuad[strftime(moTemp$dateTime,'%Y')==2013])
moSeas2014<-max(moTemp$gppQuad[strftime(moTemp$dateTime,'%Y')==2014])-min(moTemp$gppQuad[strftime(moTemp$dateTime,'%Y')==2014])
mean(c(moEventResponseGPP$maxDisplace[strftime(moEventResponseGPP$dateTime,'%Y')==2013]/moSeas2013,moEventResponseGPP$maxDisplace[strftime(moEventResponseGPP$dateTime,'%Y')==2014]/moSeas2014))

elSeas2013<-max(elTemp$gppQuad[strftime(elTemp$dateTime,'%Y')==2013])-min(elTemp$gppQuad[strftime(elTemp$dateTime,'%Y')==2013])
elSeas2014<-max(elTemp$gppQuad[strftime(elTemp$dateTime,'%Y')==2014])-min(elTemp$gppQuad[strftime(elTemp$dateTime,'%Y')==2014])
mean(c(elEventResponseGPP$maxDisplace[strftime(elEventResponseGPP$dateTime,'%Y')==2013]/elSeas2013,elEventResponseGPP$maxDisplace[strftime(elEventResponseGPP$dateTime,'%Y')==2014]/elSeas2014))

crSeas2013<-max(crTemp$gppQuad[strftime(crTemp$dateTime,'%Y')==2013])-min(crTemp$gppQuad[strftime(crTemp$dateTime,'%Y')==2013])
crSeas2014<-max(crTemp$gppQuad[strftime(crTemp$dateTime,'%Y')==2014])-min(crTemp$gppQuad[strftime(crTemp$dateTime,'%Y')==2014])
mean(c(crEventResponseGPP$maxDisplace[strftime(crEventResponseGPP$dateTime,'%Y')==2013]/crSeas2013,crEventResponseGPP$maxDisplace[strftime(crEventResponseGPP$dateTime,'%Y')==2014]/crSeas2014))



