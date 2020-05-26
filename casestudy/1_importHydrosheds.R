
# Part 1: Open flow1k and find the connection between the segments #####


##### Possible improvements:
##### include a buffer for the part that is cut out from river network
##### convert to Cartesian coordinates before cropping


require(sf)
library(rgdal)
require(readxl)
require(data.table)
library(elevatr)

setwd("../data/casestudy/na_riv") #data directory
hydrosheds <- read_sf(dsn = ".", layer = "na_riv_15s") #read hydrosheds
setwd("..") #data directory

# read the data, name the columns consistently, only keep those that are interesting and convert to numeric

datapoints_Fish <- as.data.frame(read_excel("Ohio_data_correction_April2017_stand.xlsx", sheet = "Fish")) # read measurement stations for fish
datapoints_Fish <- datapoints_Fish[3:nrow(datapoints_Fish),] # drop first two rows (mean and sd)
datapoints_Fish <- datapoints_Fish[,c("SYCode", "DraAre", "Latitu", "Longit", "Slope_C", "Channe_C",
                                      "Cover_C", "GraMet_C", "Pool_C", "Riffle_C", "Ripari_C", "Substr_C",
                                      "BOD5", "COD", "SpeCon_median", "CaCO3_median", "P_median",
                                      "TKN_median", "pH_median", "TDS", "TSS", "PAFInd_avg", "IBI")]
names(datapoints_Fish) <- c("SYCode", "DraAre", "Latitu", "Longit", "Slope", "Channe", "Cover", "GraMet", "Pool", "Riffle", "Ripari", "Substr", "BOD5", "COD", "SpeCon", "CaCO3", "P", "TKN", "pH", "TDS", "TSS", "PAFInd", "IBI")
datapoints_Fish <- transform(datapoints_Fish,
                             DraAre = as.numeric(DraAre),
                             LatituTemp = as.numeric(Latitu), 
                             LongitTemp = as.numeric(Longit),
                             Latitu = as.numeric(Latitu), 
                             Longit = as.numeric(Longit),
                             Slope = as.numeric(Slope),
                             Channe = as.numeric(Channe),
                             Cover = as.numeric(Cover),
                             GraMet = as.numeric(GraMet),
                             Pool = as.numeric(Pool),
                             Riffle = as.numeric(Riffle),
                             Ripari = as.numeric(Ripari),
                             Substr = as.numeric(Substr),
                             BOD5 = as.numeric(BOD5),
                             COD = as.numeric(COD),
                             SpeCon = as.numeric(SpeCon),
                             CaCO3 = as.numeric(CaCO3),
                             P = as.numeric(P),
                             TKN = as.numeric(TKN),
                             pH = as.numeric(pH),
                             TDS = as.numeric(TDS),
                             TSS = as.numeric(TSS),
                             PAFInd = as.numeric(PAFInd),
                             IBI = as.numeric(IBI))

datapoints_Inv <- as.data.frame(read_excel("Ohio_data_correction_April2017_stand.xlsx", sheet = "Invertebrates 1")) # read measurement stations for invertebrates
names(datapoints_Inv) <- as.character(datapoints_Inv[3,]) # rename the columns
datapoints_Inv <- datapoints_Inv[6:nrow(datapoints_Inv),3:ncol(datapoints_Inv)] #drop the uninformative rows and columns
datapoints_Inv <- datapoints_Inv[,c("Site Code", "ICI")]
names(datapoints_Inv) <- c("SYCode", "ICI")
datapoints_Inv <- transform(datapoints_Inv, ICI = as.numeric(ICI))

# combine fish and invertebrate tables

datapoints <- merge(datapoints_Fish, datapoints_Inv, by = "SYCode", all = TRUE)
data <- as.data.table(datapoints)

# remove incomplete cases

data <- na.omit(data, cols = c("SYCode", "DraAre", "LatituTemp", "LongitTemp", "Latitu", "Longit", "Slope", "Channe", "Cover", "GraMet", "Pool", "Riffle", "Ripari", "Substr", "BOD5", "COD", "SpeCon", "CaCO3", "P", "TKN", "pH", "TDS", "TSS", "PAFInd"))

datapoints <- st_as_sf(data, coords = c("LongitTemp", "LatituTemp")) # convert to simple features
st_crs(datapoints) <- 4326 # set crs to WGS84

# set coordinates

coordinates <- data.frame(datapoints$Longit, datapoints$Latitu)
crs <- st_crs(datapoints)

# get the elevation at the measurement locations

elevationPoints <- get_elev_point(coordinates, prj = crs$proj4string, src = "aws")$elevation
datapoints$Elev <- elevationPoints

# transform to cartesian coordinates

hydrosheds <- st_transform(hydrosheds, crs = 32717) # convert to cartesian coordinates (UTM 17S)
datapoints <- st_transform(datapoints, crs = 32717) # convert measurement stations to cartesian coordinates, too

# from hydrosheds, extract the part that belongs to the measurement region

box <- st_bbox(datapoints) + c(-10000,-10000,10000,10000)
hydrosheds_cropped <- st_crop(x = hydrosheds, box) # cut the part out where there are actually measurement stations

# generate a plot of locations and hydrosheds

pdf('../../plots/casestudy/locations.pdf')
plot(hydrosheds_cropped$geometry, col = 'blue', xlim = c(box[[1]], box[[3]]), ylim = c(box[[2]], box[[4]]))#, axes = TRUE)
par(new = TRUE)
plot(datapoints$geometry, pch = 4, cex = 0.5, xlim = c(box[[1]], box[[3]]), ylim = c(box[[2]], box[[4]]))#, axes = TRUE)
dev.off()
# empty data table for river segments

frame <- data.table(matrix(NA_real_, nrow = 100000, ncol = 5)) # 100000 is just a large number so that we don't run out of space; data.table for runtime
names(frame) <- c("StartX", "StartY", "EndX", "EndY", "Segment") # rename columns, segments tells us which the original segment was in river data to reassemble smaller segments later

counter <- 1 # counter to store all line segments in data.table

pb <- txtProgressBar(min = 1, max = length(hydrosheds_cropped$geometry)) # progress bar

# loop over all simple feature structures (lines and multilines) and convert them to single lines; store in dataframe

for (i in 1:length(hydrosheds_cropped$geometry)) {
  setTxtProgressBar(pb, i)
  coordinates <- st_coordinates(hydrosheds_cropped$geometry[i])
  for (j in 1:(nrow(coordinates)-1)) {
    frame[counter, StartX := coordinates[j,1]]
    frame[counter, StartY := coordinates[j,2]]
    frame[counter, EndX := coordinates[j+1,1]]
    frame[counter, EndY := coordinates[j+1,2]]
    frame[counter, Segment := i]
    counter <- counter + 1
  }
}

frame <- frame[complete.cases(frame),] # remove the rows that were not needed

pb <- txtProgressBar(min = 1, max = nrow(frame)) # new progress bar

# loop over all single line segments and find their connections (which segment follows which segment); store these connections in data.table (downstream -> upstream)

for (i in 1:nrow(frame)) {
  setTxtProgressBar(pb, i)
  x <- frame$StartX[i]
  y <- frame$StartY[i]
  connection <- which(frame$EndX == x & frame$EndY == y)
  if(length(connection) > 0) {
    for (j in 1:length(connection)) {
      frame[i, paste0('connection', toString(j)) := connection[j]]
    }
  }
}

pb <- txtProgressBar(min = 1, max = length(hydrosheds_cropped$geometry)) # new progress bar

# loop over all line segments, link them to "old" (larger) segments again and find the connection to other "old" segments

for (i in 1:length(hydrosheds_cropped$geometry)) {
  
  setTxtProgressBar(pb, i)
  
  segment_frame <- frame[frame$Segment == i]
  connected_segments <- c(segment_frame$connection1, segment_frame$connection2, segment_frame$connection3)
  connected_segments <- connected_segments[!is.na(connected_segments)]
  connected_segments <- frame$Segment[connected_segments]
  connected_segments <- unique(connected_segments)
  connected_segments <- connected_segments[which(connected_segments != i)]
  
  if(length(connected_segments) > 0) {
    for (j in 1:length(connected_segments)) {
      hydrosheds_cropped[i, paste0('connection', toString(j))] <- connected_segments[j]
    }
  }
}

hydrosheds_cropped["Associated"] <- 0 # variable for the number of measurements stations per segment (1 if there is one station in the segment, 2 for two, ...)

pb <- txtProgressBar(min = 1, max = nrow(datapoints)) # new progress bar

# loop over all measurement stations and connect them to a single segment (the one with the closest distance)

for (i in 1:nrow(datapoints)) {
  
  setTxtProgressBar(pb, i)

  segment <- datapoints$geometry[i]
  nearest_points <- st_nearest_points(hydrosheds_cropped, segment) # find nearest points in all segments
  distances <- st_distance(hydrosheds_cropped, segment) # find distances to all nearest points
  projection <- st_cast(nearest_points[which.min(distances)], "POINT")[1] # projection on line segment
  distance_to_start <- st_distance(projection, st_cast(hydrosheds_cropped$geometry[which.min(distances)], "POINT")[1]) # distance to the start of the segment (used for the order of measurement station)

  datapoints$segment[i] <- which.min(distances) # save all information
  datapoints$distance[i] <- min(distances)
  datapoints$to_start[i] <- distance_to_start
  
  hydrosheds_cropped$Associated[which.min(distances)] <- hydrosheds_cropped$Associated[which.min(distances)] + 1
  
}

saveRDS(hydrosheds_cropped, file = "hydroshedsRiverSegments.rds")
saveRDS(datapoints, file = "hydroshedsMeasurementStations.rds")
