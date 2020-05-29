
# import libraries

library(rjson)

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
filename <- "models/casestudy/seed.txt"
sink(file = filename, append = FALSE)
cat(randomSeed)
sink()

# read variable names and combine them

namesIn <- settings$namesIn
namesUp <- settings$namesUp
namesOut <- settings$namesOut
namesMisc <- settings$namesMisc
names <- c(namesIn, namesUp, namesOut, namesMisc)
names <- names[order(names)]

# read data

data <- readRDS("data/casestudy/datasetHydrosheds.rds")

# add IBI to variable names

namesFish <- c(names, "IBI")

# remove everything but the requested variables

dataFish <- data[,namesFish]
dataFish <- dataFish[complete.cases(dataFish),]

# generate sufficient statistics

mat <- as.matrix(dataFish)
suffStat <- list(C = cor(mat), n = nrow(mat))

# get positions of out of stream, miscellaneous and biotic variables

posOut <- which(rownames(suffStat$C) %in% namesOut)
posMisc <- which(rownames(suffStat$C) %in% c(namesMisc, "IBI"))
posBio <- which(rownames(suffStat$C) %in% "IBI")

# get counterpart pairs of upstream and instream variables

partners <- matrix(nrow = length(namesIn), ncol = 2)
j <- 1

for (name in namesIn) {
  partner <- paste0("U", substr(name, start = 2, stop = nchar(name)))
  partners[j,1] <- which(rownames(suffStat$C) %in% partner)
  partners[j,2] <- which(rownames(suffStat$C) %in% name)
  j <- j + 1
}

# calculate ranks for all variables
for(col in 1:ncol(dataFish)) {
    dataFish[, col] <- rank(dataFish[, col])
}

# perform bootstrapping
fishModel <- bootstrap(data = dataFish, nBoot = 1000, pBoot = 1, indepTest = gaussCItest, alpha = 0.01,
                 partners = partners, posOut = posOut, posMisc = posMisc,
                 doPdsep = TRUE, conservative = TRUE, m.max = 3, posBio = posBio, correlation = "spearman")

# save the model
filename <- paste0('models/casestudy/fishModel.rdata')
save(fishModel, file = filename)

# add ICI to variable names

namesInv <- c(names, "ICI")

# remove everything but the requested variables

dataInv <- dataInv[,namesInv]
dataInv <- dataInv[complete.cases(dataInv),]

# generate sufficient statistics

mat <- as.matrix(dataInv)
suffStat <- list(C = cor(mat), n = nrow(mat))

# get positions of out of stream and miscellaneous variables

posOut <- which(rownames(suffStat$C) %in% namesOut)
posMisc <- which(rownames(suffStat$C) %in% c(namesMisc, "ICI"))
posBio <- which(rownames(suffStat$C) %in% "ICI")

# get counterpart pairs of upstream and instream variables

partners <- matrix(nrow = length(namesIn), ncol = 2)
j <- 1

for (name in namesIn) {
  partner <- paste0("U", substr(name, start = 2, stop = nchar(name)))
  partners[j,1] <- which(rownames(suffStat$C) %in% partner)
  partners[j,2] <- which(rownames(suffStat$C) %in% name)
  j <- j + 1
}

# calculate ranks for all variables
for(col in 1:ncol(dataInv)) {
    dataInv[, col] <- rank(dataInv[, col])
}

# perform bootstrapping
invModel <- bootstrap(data = dataInv, nBoot = 1000, pBoot = 1, indepTest = gaussCItest, alpha = 0.01,
                 partners = partners, posOut = posOut, posMisc = posMisc,
                 doPdsep = TRUE, conservative = TRUE, m.max = 3, posBio = posBio, correlation  = "spearman")

# save the model
filename <- paste0('models/casestudy/invModel.rdata')
save(invModel, file = filename)
