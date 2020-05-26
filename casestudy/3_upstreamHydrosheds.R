

# Part 3: Calculate upstream measurements

# import packages

require(sf)
require(data.table)

# read the data

setwd("../data/casestudy")
connections <- readRDS("stationConnectionsHydrosheds.rds")
stations <- readRDS("hydroshedsMeasurementStations.rds")

# set up an data frame with all required columns

measurementFrame <- data.frame(SYCode = NA_real_,
                               IBI = NA_real_,
                               ICI = NA_real_,
                               IDraAre = NA_real_,
                               ILatitu = NA_real_,
                               ILongit = NA_real_,
                               IElev = NA_real_,
                               ISlope = NA_real_,
                               IChanne = NA_real_,
                               ICover = NA_real_,
                               IGraMet = NA_real_,
                               IPool = NA_real_,
                               IRiffle = NA_real_,
                               IRipari = NA_real_,
                               ISubstr = NA_real_,
                               IBOD5 = NA_real_,
                               ICOD = NA_real_,
                               ISpeCon = NA_real_,
                               ICaCO3 = NA_real_,
                               IP = NA_real_,
                               ITKN = NA_real_,
                               IpH = NA_real_,
                               ITDS = NA_real_,
                               ITSS = NA_real_,
                               IPAFInd = NA_real_,
                               UDraAre = NA_real_,
                               ULatitu = NA_real_,
                               ULongit = NA_real_,
                               UElev = NA_real_,
                               USlope = NA_real_,
                               UChanne = NA_real_,
                               UCover = NA_real_,
                               UGraMet = NA_real_,
                               UPool = NA_real_,
                               URiffle = NA_real_,
                               URipari = NA_real_,
                               USubstr = NA_real_,
                               UBOD5 = NA_real_,
                               UCOD = NA_real_,
                               USpeCon = NA_real_,
                               UCaCO3 = NA_real_,
                               UP = NA_real_,
                               UTKN = NA_real_,
                               UpH = NA_real_,
                               UTDS = NA_real_,
                               UTSS = NA_real_,
                               UPAFInd = NA_real_)

# progress bar

pb <- txtProgressBar(min = 1, max = length(stations$geometry))

# loop over all stations

for (i in 1:length(stations$geometry)) {
  
  # update progress bar and get stations
  
  setTxtProgressBar(pb, i)
  station <- stations[i,]
  
  # get all upstream stations and loop over them
  upstreamStations <- connections$upstreamStations[connections$stationID == station$SYCode]
  if(length(upstreamStations) > 0) {
      for (j in 1:length(upstreamStations)) {
        
      # for each upstream station, one row in the data frame, later we average the rows
        
      upstreamStation <- stations[stations$SYCode == upstreamStations[j],]
      newRow = data.frame(SYCode = station$SYCode,
                          IBI = station$IBI,
                          ICI = station$ICI,
                          IDraAre = station$DraAre,
                          ILatitu = station$Latitu,
                          ILongit = station$Longit,
                          IElev = station$Elev,
                          ISlope = station$Slope,
                          IChanne = station$Channe,
                          ICover = station$Cover,
                          IGraMet = station$GraMet,
                          IPool = station$Pool,
                          IRiffle = station$Riffle,
                          IRipari = station$Ripari,
                          ISubstr = station$Substr,
                          IBOD5 = station$BOD5,
                          ICOD = station$COD,
                          ISpeCon = station$SpeCon,
                          ICaCO3 = station$CaCO3,
                          IP = station$P,
                          ITKN = station$TKN,
                          IpH = station$pH,
                          ITDS = station$TDS,
                          ITSS = station$TSS,
                          IPAFInd = station$PAFInd,
                          UDraAre = upstreamStation$DraAre,
                          ULatitu = upstreamStation$Latitu,
                          ULongit = upstreamStation$Longit,
                          UElev = upstreamStation$Elev,
                          USlope = upstreamStation$Slope,
                          UChanne = upstreamStation$Channe,
                          UCover = upstreamStation$Cover,
                          UGraMet = upstreamStation$GraMet,
                          UPool = upstreamStation$Pool,
                          URiffle = upstreamStation$Riffle,
                          URipari = upstreamStation$Ripari,
                          USubstr = upstreamStation$Substr,
                          UBOD5 = upstreamStation$BOD5,
                          UCOD = upstreamStation$COD,
                          USpeCon = upstreamStation$SpeCon,
                          UCaCO3 = upstreamStation$CaCO3,
                          UP = upstreamStation$P,
                          UTKN = upstreamStation$TKN,
                          UpH = upstreamStation$pH,
                          UTDS = upstreamStation$TDS,
                          UTSS = upstreamStation$TSS,
                          UPAFInd = upstreamStation$PAFInd)
      measurementFrame <- rbind(measurementFrame, newRow)
    }
  }
}

# scrap the first row

measurementFrame <- measurementFrame[2:nrow(measurementFrame),]

# average for multiple upstream stations

measurementFrame <- aggregate(measurementFrame,by=list(ID=measurementFrame$SYCode),data=measurementFrame,FUN=mean, na.rm = TRUE)

# remove empty variable

measurementFrame <- within(measurementFrame, rm(SYCode))

# save the dataset

saveRDS(measurementFrame, "datasetHydrosheds.rds")