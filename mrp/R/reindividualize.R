# One might need to take the poll data structure that is constructed by the mrp package, modify it, and then want to run
# again.  The problem is that mrp cannot take cell-structured poll data.  This function, convertCellToIndividual, 
# will take such a data.frame and create a new individual-level data.frame that can be used for a subsequent mrp. 
# An example of usage is given.






createObservationCopies <- function (aggregate.data, index.var = "Freq", retain.freq = FALSE) 
{
  output <- NULL
  for (i in 1:nrow(aggregate.data)) {
    if (retain.freq) {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[, 
                                                                   which(names(aggregate.data) == index.var)][i]), 
                                             ])
    }
    else {
      output <- rbind(output, aggregate.data[rep(i, aggregate.data[, 
                                                                   which(names(aggregate.data) == index.var)][i]), 
                                             ][, -which(names(aggregate.data) == index.var)])
    }
  }
  data.frame(output, row.names = 1:nrow(output))
}


convertCellToIndividual <- function( cell.version){
  cell.version$response.yes <- round(cell.version$response.yes)
  cell.version$response.no <- round(cell.version$response.no)
  
  temp.yes <- createObservationCopies(cell.version, index.var="response.yes", retain.freq = TRUE)
  temp.no <- createObservationCopies(cell.version, index.var="response.no", retain.freq = TRUE)
  
  temp.yes$response <- 1
  #temp.yes <- subset(temp.yes, select = -c(response.no)  )
  
  
  temp.no$response <- 0
  #  temp.no <- subset(temp.no, select = -c(response.yes)  )
  
  
  temp <- rbind(temp.yes, temp.no)
  temp <- subset(temp, select = -c(N, design.effect.cell,ybar.w, response.yes, response.no, finalrow))
  temp$response <- factor(ifelse(temp$response ==1, "yes", "no"))
  temp
}


## usage
## let mrp.testing be an mrp object

## this flattens the population slot
#flatpop <- as.data.frame.table(mrp.testing@population)
## this joins the cell level data with the population information including frequencies. 
#jrl <- join(mrp.testing@data,flatpop)
## this could now be modified as needed
## the below cbinds the fitted values which could be used for modification
#jrl <- cbind(jrl, fitted(mrp.testing@multilevelModel))




## Below is how to convert cell level data to individual
# new.individualized.polldata <- convertCellToIndividual(jrl)
