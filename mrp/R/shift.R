###  goal:  take  a vector of opinion levels -- say, state opinion levels --- and figure out what they would look like if projected to a different national mean
##### function that takes an opinion vector across states, a weighting vector of state population percentage of national total, and a projected national mean and gives you a new opinion vector across states with a shift on the logit scale that leads to that national mean
###   It could be across age levels or another poststatified breakdown.  Or even the full set of cells.  So long as th weight.vec sms to one.

shifted.op.vec <- function(old.op.vec, weight.vec, nat.target){
  new_shifted_national <- function(x){   ((invlogit(logit(old.op.vec)  + x )  %*%  weight.vec   ) - nat.target)^2
  }
  needed.shift <-  optimize(f = new_shifted_national, c(-100,100))$minimum
  new.op.vec  <- invlogit(logit(old.op.vec)  + needed.shift )
  new.natop <- new.op.vec %*% weight.vec
  new.op.vec
}

#example of usage
## options(digits =2)
## old.op.vec <- .07* c(1:10)
## weight.vec <-  runif(10)
## weight.vec <-  weight.vec/ sum(weight.vec)
## old.natop <- old.op.vec %*%  weight.vec

## old.op.vec
## old.op.vec %*% weight.vec

## shifted.op.vec(old.op.vec,weight.vec, .5)
## shifted.op.vec(old.op.vec,weight.vec, .5) %*% weight.vec

## shifted.op.vec(old.op.vec,weight.vec, .75)
## shifted.op.vec(old.op.vec,weight.vec, .75) %*% weight.vec



##### NEW SHIFT FUNCTION for state vote total.
##### only one margin now.

## shift(mrp,
## turnoutData (vector or char name of data col),
## ~shiftvar)

## p1 = collapse across population by shiftvar
## delta = apply (p1, shiftfun, turnoutdata)
##  do this optimization and get a bunch of
##  shiftresults
## p2 = sweep delta+population array
## return a full-dimension shifted array
