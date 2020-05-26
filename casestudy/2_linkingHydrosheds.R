
# Part 2: Link the measurement stations to the upstream ones


# this function finds for a specific station all closest stations upstream

recursiveLinking <- function(stationID, segments, stations, segmentID) {
  
  # first, if this is the segment of the station, we have to apply special rules
  if (segmentID == stations$segment[stations$SYCode == stationID]) {
    # if there is only one station (the station we want to link something to) in the segment, then we do not have to do anything
    # thus, we only look at the segment if there is more than one station in it
    if (segments$Associated[segmentID] > 1) {
      
      # get all stations and order by elevation (as the water flows down, this is equivalent to the order the stations appear in the segment)
      inSegment <- stations[which(stations$segment == segmentID),]
      inSegment <- inSegment[order(inSegment$Elev, decreasing = TRUE),]
      
      # find the starting segment
      stationPos <- which(inSegment$SYCode == stationID)
      
      # if the starting segment is the first, we have to look further upstream. Otherwise, take it and return it
      if (stationPos > 1) {
        
        upstreamStations <- inSegment$SYCode[stationPos - 1]
        upstreamFrame <- data.frame(stationID, upstreamStations)
        return(upstreamFrame)
      }
    }
  } else {
    
    # if this is any other segment, we have to check whether there are any stations in it
    if (segments$Associated[segmentID] > 0) {
      
      # if so, order them, and take the last station in the segment
      inSegment <- stations[which(stations$segment == segmentID),]
      inSegment <- inSegment[order(inSegment$Elev, decreasing = FALSE),]
      
      # get the station and return it
      upstreamStations <- inSegment$SYCode[1]
      upstreamFrame <- data.frame(stationID, upstreamStations)
      return(upstreamFrame)
    }
  }
  # if neither of the first two conditions triggered, we have to look further upstream 
  # there is a maximum of three connections to any segment, so we go over them and call the recursive functions
  if (is.na(segments$connection1[segmentID]) == FALSE) {
    upstreamStation1 <- recursiveLinking(stationID = stationID, segments = segments, stations = stations, segments$connection1[segmentID])
    if (is.na(segments$connection2[segmentID]) == FALSE) {
      upstreamStation2 <- recursiveLinking(stationID = stationID, segments = segments, stations = stations, segments$connection2[segmentID])
      if (is.na(segments$connection3[segmentID]) == FALSE) {
        upstreamStation3 <- recursiveLinking(stationID = stationID, segments = segments, stations = stations, segments$connection3[segmentID])
        # depending on the number of connections, we then combine all the stations and return them
        upstreamStations <- rbind(upstreamStation1, upstreamStation2, upstreamStation3)
        return(upstreamStations)
      } else {
        upstreamStations <- rbind(upstreamStation1, upstreamStation2)
        return(upstreamStations)
      }
    } else {
      upstreamStations <- upstreamStation1
      return(upstreamStations)
    }
  } else {
    upstreamStations <- NA
    upstreamFrame <- data.frame(stationID, upstreamStations)
    return(upstreamFrame)
  }
}


# import packages

require(sf)
require(data.table)

# read the data
setwd("../data/casestudy") #data directory
segments <- readRDS("hydroshedsRiverSegments.rds") #read data
stations <- readRDS("hydroshedsMeasurementStations.rds")

# generate an empty data frame to append to

upstreamStations <- NA
stationID <- NA
segmentFrame <- data.frame(stationID, upstreamStations)

# progress bar

pb <- txtProgressBar(min = 1, max = length(stations$geometry))

# loop over all stations

for (i in 1:length(stations$geometry)) {
  
  # update progress bar
  
  setTxtProgressBar(pb, i)
  
  # get the ID and the segment for the station
  
  stationID <- stations$SYCode[i]
  segmentID <- stations$segment[i]
  
  # call the function
  
  upstreamFrame <- recursiveLinking(stationID = stationID, stations = stations, segments = segments, segmentID = segmentID)
  
  # store the output
  
  segmentFrame <- rbind(segmentFrame, upstreamFrame)
}

# scrap the first row of the data frame

connections <- segmentFrame[complete.cases(segmentFrame),]

# save the data

saveRDS(connections, "stationConnectionsHydrosheds.rds")