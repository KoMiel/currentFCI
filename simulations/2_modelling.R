
# import libraries

library(rjson)
library(pcalg)
library(foreach)
library(doParallel)

# parallel computing

registerDoParallel()


# function that converts causal structure to a dag representation

convertStructure <- function(structure) {
  
  # every value that is not 0 is a cause (= 2)
  structure[structure != 0] <- 2
  
  # for everything else, we start with no edge
  structure[structure == 0] <- 0
  
  # loop over the entire matrix
  for (i in 1:(ncol(structure))) {
    for (j in 1:(ncol(structure))) {
      
      # for each cause, we have to add a tail mark (= 3) to the other side of the edge
      if (structure[i, j] == 2) {
        structure[j, i] <- 3
      }
    }
  }
  
  # transpose the structure
  structure <- t(structure)
  
  # return the structure
  return(structure)
}

# function that eliminates latent variables from the dataset

latentVariables <- function(data, namesLat) {
  
  # get all variable names
  names <- names(data)
  
  # get all other variable names and the corresponding observations
  names <- setdiff(names, namesLat)
  data <- data[, names]
 
  # return the data 
  return(data = data)
}

# function that marginalizes latent nodes from a causal graph

marginalize <- function(structure, nodes) {
  
  for (node in nodes) {
    
    # get all parents and childs of the node
    childs <- names(which(structure[node, ] == 2))
    
    # expand the childs, eliminate duplicates
    confounder <- expand.grid(childs, childs, stringsAsFactors = FALSE)
    confounder <- confounder[confounder[,1] != confounder[,2],]
    confounder <- confounder[confounder[,1] < confounder[,2],]
    
    # set the childs as confounded
    for (row in 1:nrow(confounder)) {
      structure[confounder[row, 1], confounder[row, 2]] <- 2
      structure[confounder[row, 2], confounder[row, 1]] <- 2
    }
  }
  
  # eliminate latent variables
  structure <- structure[(length(nodes)+1):nrow(structure), (length(nodes)+1):ncol(structure)]
  
  # return the causal graph
  return(structure)
}


# set working directory

setwd("..")

# import source file

source("currentFCI.R")

# read settings file

f <- "settings.json"
settings <- fromJSON(file = f)

# generate a random seed and use it

seed <- runif(1, min = 0, max = 1000000)
set.seed(seed)

# save the random seed 

randomSeed <- c(seed, '\n')
filename <- "models/simulations/seed.txt"
sink(file = filename, append = FALSE)
cat(randomSeed)
sink()

# read parameters
nRuns <- settings$nRuns
nDatas <- settings$nDatas
densenesses <- settings$densenesses
p <- settings$p
alphas <- settings$alphas

# loop over conditions

foreach (denseness = densenesses) %dopar% {
  dir.create(path = paste0('models/simulations/', denseness, "/"), showWarnings = FALSE)
  
  foreach (nData = nDatas) %dopar% {
    dir.create(path = paste0('models/simulations/', denseness, "/", nData, "/"), showWarnings = FALSE)
    
    foreach (alpha = alphas) %dopar% {
      dir.create(path = paste0('models/simulations/', denseness, "/", nData, "/", alpha, "/"), showWarnings = FALSE)
      
      foreach (i = 1:nRuns) %dopar% {
        cat(paste(denseness, nData, alpha, i, "\n"))
        # load the data
        filename <- paste0("data/simulations/", denseness, "/", nData, "/data_", i, ".rdata")
        load(file = filename)
        network <- data[[1]]
        structure <- data[[2]]
        observations <- data[[3]]
        
        # convert adjacency matrix to directed acyclic graph representation
        dag <- convertStructure(structure$mat)
        
        # eliminate latent variables and marginalize them out in graph
        latObservations <- latentVariables(data = observations, namesLat = structure$namesLat)
        latDag <- marginalize(structure = dag, nodes = structure$namesLat)
        
        # generate sufficient statistics
        mat <- as.matrix(latObservations)
        suffStat <- list(C = cor(mat), n = nrow(mat))
        
        # fit a causal model with standard fci
        fciModelFull <- fci(m.max = 3, suffStat, indepTest = gaussCItest, alpha = alpha, labels = rownames(suffStat$C), doPdsep = TRUE, conservative = TRUE, verbose = FALSE)        

        # find the positions in the data matrix belonging to out of stream and misc variables
        posOut <- which(rownames(suffStat$C) %in% structure$namesOut)
        posMisc <- which(rownames(suffStat$C) %in% structure$namesMisc)
        
        # generate a 2 column matrix with partner variabls positions (upstream-instream pairs)
        posPartners <- matrix(nrow = length(structure$namesUp), ncol = 2)
        j <- 1
        for (name in structure$namesIn) {
          partner <- paste0("U", substr(name, start = 2, stop = nchar(name)))
          posPartners[j,1] <- which(rownames(suffStat$C) %in% partner)
          posPartners[j,2] <- which(rownames(suffStat$C) %in% name)
          j <- j + 1
        }      
  
        # fit a causal model with current FCI
        currentFciModel <- currentFci(m.max = 3, suffStat, indepTest = gaussCItest, alpha = alpha, labels = rownames(suffStat$C), doPdsep = TRUE, conservative = TRUE, verbose = FALSE, partners = posPartners, posOut = posOut, posRem = posMisc)

        # eliminate out of stream and upstream variables from data
        inObservations <- latObservations[,c(structure$namesIn, structure$namesMisc)]
        inMat <- as.matrix(inObservations)
        suffStat <- list(C = cor(inMat), n = nrow(inMat))
        
        # fit a causal model with standard FCI on limited data set
        fciModelPart <- fci(m.max = 3, suffStat, indepTest = gaussCItest, alpha = alpha, labels = rownames(suffStat$C), doPdsep = TRUE, conservative = TRUE, verbose = FALSE)  
        
        results <- list(latDag, currentFciModel, fciModelFull, fciModelPart)
        filename <- paste0('models/simulations/', denseness, "/", nData, "/", alpha, "/model_", i, ".rdata")
        save(results, file = filename)
      }
    }
  }
}
