
# import libraries

library(rjson)
library(pcalg)


# function that calculates the confusion matrix for two causal adjacency matrices
# structure 1 is the ESTIMATED matrix
# structure 2 is the TRUE matrix

confusionMatrix <- function(structure1, structure2, partners, posMisc) {
  
  # set up counts for edges and edge marks (rows: estimated matrix, cols: true matrix)
  edges <- matrix(nrow = 2, ncol = 2, 0)
  rownames(edges) <- c("present", "absent")
  colnames(edges) <- c("present", "absent")
  marks <- matrix(nrow = 3, ncol = 2, 0)
  rownames(marks) <- c("arrow", "circle", "tail")
  colnames(marks) <- c("arrow", "tail")
  
  # loop over rows and columns of the matrices, counting edges
  for (rowname in rownames(structure1)) {
    for (colname in colnames(structure1)) {
      
      # check cases: no edge in either matrix
      if (structure1[rowname, colname] == 0 & structure2[rowname, colname] == 0) {
        edges["absent", "absent"] <- edges["absent", "absent"] + 1  
        
      # only an edge in matrix 1
      } else if (structure1[rowname, colname] != 0 & structure2[rowname, colname] == 0) {
        edges["present", "absent"] <- edges["present", "absent"] + 1
        
      # only an edge in matrix 2
      } else if (structure1[rowname, colname] == 0 & structure2[rowname, colname] != 0) {
        edges["absent", "present"] <- edges["absent", "present"] + 1
        
      # edges in both matrices
      } else {
        edges["present", "present"] <- edges["present", "present"] + 1
      }
    }
  }
  
  # correction for counting edges on diagonal
  edges["absent", "absent"] <- edges["absent", "absent"] - nrow(structure1)

  # correction for counting edges double
  edges <- edges/2
  
  # loop again over rows and columns of the matrices, this time counting edge marks
  for (rowname in rownames(structure1)) {
    for (colname in colnames(structure1)) {
      
      # ensure that there is an edge in both matrices
      if (structure1[rowname, colname] != 0 & structure2[rowname, colname] != 0) {
        
        # both arrowheads
        if (structure1[rowname, colname] == "2" & structure2[rowname, colname] == "2") {
          marks["arrow", "arrow"] <- marks["arrow", "arrow"] + 1
          
        # arrowhead in 1, tail in 2
        } else if (structure1[rowname, colname] == "2" & structure2[rowname, colname] == "3") {
          marks["arrow", "tail"] <- marks["arrow", "tail"] + 1
          
        # tail in 1, arrowhead in 2
        } else if (structure1[rowname, colname] == "3" & structure2[rowname, colname] == "2") {
          marks["tail", "arrow"] <- marks["tail", "arrow"] + 1
          
        # tail in 1, tail in 2
        } else if (structure1[rowname, colname] == "3" & structure2[rowname, colname] == "3") {
          marks["tail", "tail"] <- marks["tail", "tail"] + 1
          
        # circle in 1, arrowhead in 2
        } else if (structure1[rowname, colname] == "1" & structure2[rowname, colname] == "2") {
          marks["circle", "arrow"] <- marks["circle", "arrow"] + 1
          
        # circle in 1, tail in 2
        } else if (structure1[rowname, colname] == "1" & structure2[rowname, colname] == "3") {
          marks["circle", "tail"] <- marks["circle", "tail"] + 1
        }
      }
    }
  }
  
  # return edges and edge marks
  return(list(edges = edges, marks = marks))
}



# set working directory

setwd("..")

# import source file

source("currentFCI.R")

# read settings file

f <- "settings.json"
settings <- fromJSON(file = f)

# read parameters
nRuns <- settings$nRuns
nDatas <- settings$nDatas
densenesses <- settings$densenesses
alphas <- settings$alphas

# generate directories
dir.create(path = 'results/', showWarnings = FALSE)

for (noise in 1:6) {

dir.create(path = paste0('results/simulations', noise), showWarnings = FALSE)

# generate an empty data frame and set names
df <- data.frame(matrix(nrow = length(densenesses) * length(nDatas) * length(alphas), ncol = 3+30))
names(df) <- c("denseness", "nData", "alpha", rep(c("E_pp", "E_ap", "E_pa", "E_aa", "M_aa", "M_ca", "M_ta", "M_at", "M_ct", "M_tt"), 3))

# generate a second empty data frame and set names
dfMetrics <- data.frame(matrix(nrow = length(densenesses) * length(nDatas) * length(alphas), ncol = 3+18))
names(dfMetrics) <- c("denseness", "nData", "alpha", rep(c("E_prec", "E_reca", "M_a_prec", "M_a_reca", "M_t_prec", "M_t_reca"), 3))

# initialize a counter
counter <- 1

# loop over conditions
for (denseness in densenesses) {

  for (nData in nDatas) {

    for (alpha in alphas) {

      # vector full of 0s
      conf <- integer(30)
      
      for (i in 1:nRuns) {
        cat(paste(denseness, nData, alpha, i, "\n"))
        
        # load the data
        filename <- paste0('models/simulations', noise, '/', denseness, "/", nData, "/", alpha, "/model_", i, ".rdata")
        load(file = filename)
        
        # extract all dags/pags
        trueDag <- results[[1]]
        currentPag <- attributes(results[[2]])$amat
        fullFciPag <- attributes(results[[3]])$amat
        partFciPag <- attributes(results[[4]])$amat
        
        # only take the part of the larger dags/pags that is to be evaluated
        trueDag <- trueDag[rownames(partFciPag), colnames(partFciPag)]
        currentPag <- currentPag[rownames(partFciPag), colnames(partFciPag)]
        fullFciPag <- fullFciPag[rownames(partFciPag), colnames(partFciPag)]
        
        # calculate the confusion matrices
        currentConf <- confusionMatrix(currentPag, trueDag)
        fullFciConf <- confusionMatrix(fullFciPag, trueDag)
        partFciConf <- confusionMatrix(partFciPag, trueDag)
        
        # combine everything into one vector and sum up
        conf <- conf + c(as.vector(currentConf$edges), as.vector(currentConf$marks), as.vector(fullFciConf$edges), as.vector(fullFciConf$marks), as.vector(partFciConf$edges), as.vector(partFciConf$marks))
      }
      
      df[counter, ] <- c(denseness, nData, alpha, conf)

      # calculate metrics
      dfMetrics[counter, 1:3] <- c(denseness, nData, alpha)
      dfMetrics[counter, 4] <- conf[1]/(conf[1] + conf[3])
      dfMetrics[counter, 5] <- conf[1]/(conf[1] + conf[2])
      dfMetrics[counter, 6] <- conf[5]/(conf[5] + conf[8])
      dfMetrics[counter, 7] <- conf[5]/(conf[5] + conf[6] + conf[7])
      dfMetrics[counter, 8] <- conf[10]/(conf[7] + conf[10])
      dfMetrics[counter, 9] <- conf[10]/(conf[8] + conf[9] + conf[10])
      dfMetrics[counter, 4 + 6*1] <- conf[1 + 10*1]/(conf[1 + 10*1] + conf[3 + 10*1])
      dfMetrics[counter, 5 + 6*1] <- conf[1 + 10*1]/(conf[1 + 10*1] + conf[2 + 10*1])
      dfMetrics[counter, 6 + 6*1] <- conf[5 + 10*1]/(conf[5 + 10*1] + conf[8 + 10*1])
      dfMetrics[counter, 7 + 6*1] <- conf[5 + 10*1]/(conf[5 + 10*1] + conf[6 + 10*1] + conf[7 + 10*1])
      dfMetrics[counter, 8 + 6*1] <- conf[10 + 10*1]/(conf[7 + 10*1] + conf[10 + 10*1])
      dfMetrics[counter, 9 + 6*1] <- conf[10 + 10*1]/(conf[8 + 10*1] + conf[9 + 10*1] + conf[10 + 10*1])
      dfMetrics[counter, 4 + 6*2] <- conf[1 + 10*2]/(conf[1 + 10*2] + conf[3 + 10*2])
      dfMetrics[counter, 5 + 6*2] <- conf[1 + 10*2]/(conf[1 + 10*2] + conf[2 + 10*2])
      dfMetrics[counter, 6 + 6*2] <- conf[5 + 10*2]/(conf[5 + 10*2] + conf[8 + 10*2])
      dfMetrics[counter, 7 + 6*2] <- conf[5 + 10*2]/(conf[5 + 10*2] + conf[6 + 10*2] + conf[7 + 10*2])
      dfMetrics[counter, 8 + 6*2] <- conf[10 + 10*2]/(conf[7 + 10*2] + conf[10 + 10*2])
      dfMetrics[counter, 9 + 6*2] <- conf[10 + 10*2]/(conf[8 + 10*2] + conf[9 + 10*2] + conf[10 + 10*2])
      counter <- counter + 1
    }
  }
}

# save the results
save(df, file = paste0('results/simulations', noise, '/confusion.rdata'))
save(dfMetrics, file = paste0('results/simulations', noise, '/metrics.rdata'))

}
