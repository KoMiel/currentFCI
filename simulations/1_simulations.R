
# import libraries

library(rjson)
library(foreach)
library(doParallel)

registerDoParallel()

# function that generates a random tree (used for river network)

randomNetwork <- function(totalNodes, splitProb) {
  
  # generate a data frame for the tree structure
  tree <- data.frame(matrix(nrow = totalNodes, ncol = 3))
  
  # set the names of the columns
  names(tree) <- c("Node", "Parent", "Level")
  
  # first measurement station
  tree[1,] <- c(1, NA, 1)
  
  # set up measurement station index and level
  nodeInd <- 2
  level <- 1
  
  # repeat until enough stations are found
  while(nodeInd <= totalNodes) {
    
    # get all stations on the current level
    current <- tree[which(tree$Level == level),]
    
    # loop over stations
    for (node in 1:nrow(current)){
      
      # randomly decide whether there is a split in the branch
      binom <- rbinom(n = 1, size = 1, prob = splitProb)

      # add one or two childs in the tree
      for (childs in 1:(1 + binom)) {
        tree[nodeInd,] <- c(nodeInd, current[node,1], level + 1)
        nodeInd <- nodeInd + 1
      }
    }
    
    # increase the level
    level <- level + 1
  }
  
  # return the tree
  tree <- tree[1:totalNodes,]
  return(tree)
}

# function that generates a random causal structure

randomCausalStructure <- function(denseness, nIn, nOut, nMisc, pLat) {
  
  # generate variable names
  if (nIn > 0) {
    namesUp <- paste0("U", seq(1:nIn))
    namesIn <- paste0("I", seq(1:nIn))
  } else {
    namesUp <- vector()
    namesIn <- vector()
  }
  if (nOut > 0) {
    namesOut <- paste0("O", seq(1:nOut))
  } else {
    namesOut <- vector()
  }
  if (nMisc > 0) {
    namesMisc <- paste0("R", seq(1:nMisc))
  } else {
    namesMisc <- vector()
  }
  
  namesLat <- list()

  # order the names
  names <- sample(c(namesUp, namesOut, namesIn, namesMisc))
  
  # generate an empty matrix for the structure
  mat <- matrix(ncol = length(names), nrow = length(names), 0)
  
  # set row and column names
  rownames(mat) <- names
  colnames(mat) <- names
  
  latCount <- 1
  
  nVar <- nrow(mat)
  
  # loop over the matrix
  for (i in 1:nVar) {
    for (j in 1:nVar) {
      
      # lower triangle matrix is set to 0 (so we do not simulate any confounders here)
      if (j <= i) {
        next
      
      # randomly sample whether a connection in the lower triangle exists or not
      } else {
        con <- rbinom(prob = denseness, n = 1, size = 1)
        if (con == 1) {
          mat[j, i] <- 1
        } else if (!(rownames(mat)[i] %in% namesUp) | !(rownames(mat)[j] %in% namesUp)) {
          lat <- rbinom(prob = pLat, n = 1, size = 1)
          if (lat == 1) {
            mat <- rbind(mat, 0)
            rownames(mat)[nrow(mat)] <- paste0("L", latCount)
            mat <- cbind(mat, 0)
            colnames(mat)[ncol(mat)] <- paste0("L", latCount)
            namesLat <- c(unlist(namesLat), paste0("L", latCount))
            names <- c(paste0("L", latCount), names)
            latCount <- latCount + 1  
            mat[i, nrow(mat)] <- 1
            mat[j, ncol(mat)] <- 1
          }
        }
      }
    }
  }
  
  # set connections up-up and out-up to 0, as well as connections in -> up to 0
  mat[namesUp, namesUp] <- 0
  mat[namesUp, namesOut] <- 0
  mat[namesOut, namesUp] <- 0
  mat[namesUp, namesIn] <- 0
  mat[namesOut, namesIn] <- 0
  mat[namesUp, namesMisc] <- 0
  mat[namesOut, namesMisc] <- 0

  # put latent variables on top
  mat <- mat[names, names]
  
  # loop over entries of matrix again
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      
      # if connection exists
      if (mat[i, j] == 1) {
        
        # sample a random value from (-1,-0.1) and (0.1,1)
        sign <- c(-1, 1)
        sign <- sample(x = sign, size = 1, replace = TRUE)
        link <- runif(n = 1, min = 0.1, max = 1)
        mat[i, j] <- sign * link
      }
    }
  }

  # set the connections for instream to upstream partner to the maximum (1)
  for (name in namesIn) {
    partner <- paste0("U", substr(name, start = 2, stop = nchar(name)))
    mat[name, partner] <- 1
  }
  
  # divide values by colsums to ensure that values do not explode in the simulation
  mat <- mat/max(colSums(abs(mat)))
  
  # return causal structure
  return(list(mat = mat, names= names, namesIn = namesIn, namesOut = namesOut, namesMisc = namesMisc, namesUp = namesUp, namesLat = namesLat))
}

# function that generates measurement data from station network and causal structure

simulateData <- function(network, structure, namesIn) {
  
  # generate an empty data frame for data
  data <- data.frame(matrix(nrow = nrow(network), ncol = ncol(structure) + 1, 0))
  
  # get variable names
  names <- colnames(structure)
  
  # set names of data frame
  names(data) <- c("Station", names)
  
  # generate upstream names from instream ones
  namesUp <- paste0("U", substr(namesIn, start = 2, stop = nchar(namesIn)))
  
  # order the network by level (highest level up)
  network <- network[order(network$Level, decreasing = TRUE),]
  
  # get all names besides upstream ones
  namesRest <- setdiff(names, namesUp)

  # loop over all measurement stations
  for (i in 1:nrow(network)) {
    
    # set the station
    data[i, "Station"] <- network$Node[i]
    
    # if there are no upstream stations
    if (is.na(network$Node[network$Parent == network$Node[i]])){
      
      # sample the upstream variables from normal distribution
      for (name in namesUp) {
        data[i, name] <- rnorm(1)
      }
      
      # if not, then select the upstream stations and get the data
    } else {
      upstreamStations <- network$Node[network$Parent == network$Node[i]]
      upstreamStations <- upstreamStations[!is.na(upstreamStations)]
      dataUpstream <- data[data$Station %in% upstreamStations,]

      # loop over all upstream variables
      for (name in namesUp) {

        # take the partner values, averaged over all stations
        partner <- paste0("I", substr(name, start = 2, stop = nchar(name)))
        data[i, name] <- mean(dataUpstream[,partner])
      }
    }
    
    # loop over all other variables
    for (name in namesRest) {
      
      # sample from normal distribution
      temp <- rnorm(1)

      # get the parents
      parents <- names(structure[name,][structure[name,] != 0])

      # loop over parents and use their causal effect
      for (j in parents) {
        temp <- temp + structure[name,j] * data[i,j]
      }
      
      # write values into data frame
      data[i, name] <- temp
    }
  }
  
  # eliminate the station ID, as it is not useful
  data <- data[,2:ncol(data)]
  
  # return data frame
  return(data)
}

# directory for data

dir.create(path = "../data/", showWarnings = FALSE)
dir.create(path = "../data/simulations/", showWarnings = FALSE)

# set working directory

setwd("../data/simulations/")

# read settings file

f <- "../../settings.json"
settings <- fromJSON(file = f)

# generate a random seed and use it

seed <- runif(1, min = 0, max = 1000000)
set.seed(seed)

# save the random seed 

randomSeed <- c(seed, '\n')
filename <- "seed.txt"
sink(file = filename, append = FALSE)
cat(randomSeed)
sink()

# read parameters
nRuns <- settings$nRuns
nDatas <- settings$nDatas
densenesses <- settings$densenesses
splitProb <- settings$splitProb
nIn <- settings$nIn
nOut <- settings$nOut
nMisc <- settings$nMisc
pLat <- settings$pLat

# loop over conditions and generate data (in parallel)

foreach (denseness = densenesses) %dopar% {
  dir.create(path = paste0(denseness, "/"), showWarnings = FALSE)
  
  foreach (nData = nDatas) %dopar% {
    dir.create(path = paste0(denseness, "/", nData, "/"), showWarnings = FALSE)
    
    foreach (i = 1:nRuns) %dopar% {
      cat(paste(denseness, nData, i, "\n"))
      
      # generate a random network
      network <- randomNetwork(totalNodes = nData, splitProb = splitProb)
      
      # generate a random causal structure
      structure <- randomCausalStructure(denseness = denseness, nIn = nIn, nOut = nOut, nMisc = nMisc, pLat = 0.01)

      # generate observations
      observations <- simulateData(network = network, structure = structure$mat, namesIn = structure$namesIn)
    
      # combine to data object
      data <- list(network, structure, observations)
      
      # save the data
      filename <- paste0(denseness, "/", nData, "/data_", i, ".rdata")
      save(data, file = filename)
    }
  }
}
