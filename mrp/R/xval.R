.xval <- function(object, formula, folds, loss.type, ...){
    ## create a list of length folds that holds different partitions
    require("parallel")
    K <- folds
    M <- object
    if(missing(formula)) {
        fm <- M@formula
    }
    else {
        fm <- update.formula(M@formula, formula)
    }

    fm.terms <- terms(fm)
    cl.labels <- attr(fm.terms, "term.labels")
    cl.labels <- gsub("1 \\| ", "", cl.labels)
    cl.fm <- paste("response~", cl.labels[1], sep="")
    for(i in 2:length(cl.labels))
        cl.fm <- paste(cl.fm, "+", cl.labels[i], sep="")
    cl.fm <- as.formula(cl.fm)

    cl.terms <- attr(terms(cl.fm), "term.labels")
    cl.terms <- grep(":", cl.terms, invert=T, value=T)
    index <- NULL
    for(i in cl.terms)
        index <- paste(index, M@data[, i], sep=":")

    if(missing(loss.type)) loss.type <- "log"
    ##            require(doMC, quietly=T)
    response <- M@data[,c("response.yes", "response.no")];
    response <- ceiling(response) # annoying floating point rounding errors
    partition <- array(0, dim=c(2, 2, nrow(response), K))
    for(i in 1:nrow(response)) {
        part <- sample(1:K, sum(response[i, ]), replace=T);
        full <- rep(c(1,0), response[i,])
        for(j in 1:K){
            partition[,1,i,j] <- c(sum(full[part!=j]==1), sum(full[part!=j]==0));
            partition[,2,i,j] <- c(sum(full[part==j]==1), sum(full[part==j]==0));
        }
    }

    listofpartition <- lapply(1:K, function(k){
        newdata <- M@data
        newdata$response.yes <- partition[1,1,,k];
        newdata$response.no <- partition[2,1,,k];
        testdata <- data.frame(t(partition[,2,,k]));
        colnames(testdata) <- c("response.yes", "response.no")
        list(training=newdata, testing=testdata)
    })

    if(loss.type=="log") {
        k <- 1 ## bound below but this silences 'note' from CHECK
        loss <- parallel::mclapply(1:K, function(k){
            response <- as.matrix((listofpartition[[k]]$training)[, c("response.yes", "response.no")])
            attr(fm, ".Environment") <- environment() ## crucial!! Environment of formula!!
            foo <- blmer(fm, data=listofpartition[[k]]$training, family=quasibinomial, ...)
            yhat <- try(subset(fitted(foo), rowSums(listofpartition[[k]]$testing)>0),
                        silent=TRUE)
            S <- listofpartition[[k]]$testing
            S <- subset(S, rowSums(S)>0)
                                        #    loss <- na.omit(loss)
            pred <- yhat
            logloss <- -(S$response.yes*log(yhat) + S$response.no*log(1-yhat))
            n <- rowSums(S)
            loss <- data.frame(pred, logloss, n)
            return(loss)
        })
        cl.loss <- parallel::mclapply(1:K, function(k) {
            response <- as.matrix((listofpartition[[k]]$training)[, c("response.yes", "response.no")])
            attr(cl.fm, ".Environment") <- environment()
            foo <- glm(cl.fm, data=listofpartition[[k]]$training, family=quasibinomial)
            yhat <- try(subset(fitted(foo), rowSums(listofpartition[[k]]$testing)>0),
                        silent=TRUE)
            S <- listofpartition[[k]]$testing
            S <- subset(S, rowSums(S)>0)
                                        #    loss <- na.omit(loss)
            pred <- yhat
            logloss <- -(S$response.yes*log(yhat) + S$response.no*log(1-yhat))
            n <- rowSums(S)
            loss <- data.frame(pred, logloss, n)
            return(loss)
        })
        ## lb.loss <-  foreach(k = 1:K, .verbose=FALSE) %dopar% {

        ##   S <- listofpartition[[k]]$testing
        ##   S2 <- tapply(S$response.no, index, sum)
        ##   S1 <-  tapply(S$response.yes, index, sum)
        ##   S <- data.frame(response.yes=S1, response.no=S2)
        ##   S <- subset(S, rowSums(S)>0)
        ##                           #    loss <- na.omit(loss)
        ##   yhat <- S$response.yes/(S$response.no + S$response.yes)
        ##   yhat <- pmin(0.999, yhat)
        ##   yhat <- pmax(0.001, yhat)
        ##   pred <- yhat
        ##   logloss <- -(S$response.yes*log(yhat) + S$response.no*log(1-yhat))
        ##   n <- rowSums(S)
        ##   loss <- data.frame(pred, logloss, n)
        ##   return(loss)
        ## }
    }

    mat <-
        sapply(loss, function(l){ # fold
            sum(l$logloss)
        })
    aa <- sum(mat)

    cl.mat <-
        sapply(cl.loss, function(l){ # fold
            sum(l$logloss)
        })
    cl.aa <- sum(cl.mat)

    ## lb.mat <-
    ##   sapply(lb.loss, function(l){ # fold
    ##     sum(l$logloss)
    ##   })
    ## lb.aa <- sum(lb.mat)

    response <- as.matrix(M@data[, c("response.yes", "response.no")])
    attr(cl.fm, ".Environment") <- environment()
    insampleM <- glm(cl.fm, data=M@data,  family=binomial)
    yhat <- fitted(insampleM)
    LB <- -sum(M@data$response.yes*log(yhat) + M@data$response.no*log(1-yhat))
    return(list(MulRes=aa, MulFormula=fm, ClaRes=cl.aa,  ClaFormula=cl.fm, LB=LB))
}
##' Run k-fold cross validations on mrp model
##'
##' Run a (binomial) multilevel regression in survey data for later
##' poststratification.
##'
##' @rdname xval-methods
##' @export
##' @name xval
##' @param object A \code{mrp} object.
##' @param formula A formula specification for the multilevel model on which
##' the k-folds cross validation is run. The default is the formula of the
##' \code{mrp} object. The left-hand side should always be
##' \sQuote{\code{response}}. For convenience, the formula is handled by
##' \code{update.formula} so that \code{.} indicates the current formula
##' contents on either side of the \code{~}, e.g., \code{.~.+newVar}. The
##' initial default formula is constructed as just an intercept term for each
##' of the variables in the main formula specification
##' (\code{(1|way1)+(1|way2)} etc.)
##' @param folds Number of folds for cross valiation. Default is 4.
##' @param loss.type Type of loss measure used in cross validation. Currently
##' \dQuote{log} is supported.
##' @param \dots Additional arguments to be passed to the multilevel regression
##' step, which uses \code{\link[blme]{bglmer}} by default.
##' @return a list with named elements:
##' MulRes, MulFormula, ClaRes, ClaFormula, and LB
##' WEI-- DOC THIS RETURN PLEASE
##' @docType methods
setGeneric("xval", function(object, formula, folds, loss.type, ...) {standardGeneric("xval")})
##' @rdname xval-methods
##' @aliases xval,mrp-method
setMethod(f="xval",
          signature=signature(object="mrp"),
          definition=.xval
          )

########## CAN WE KILL THIS unused fn ?
##########
## lb <- function(object){
##   M <- object@data
##   resp <- ddply(M, .(state, income), function(x) {
##     if(is.null(dim(x)))
##       return(c(0,0))
##     else
##       apply(cbind(x$response.yes, x$response.no), 2, sum)
##   }
##                 )
##   resp1 <- ddply(resp, .(state, income), function(x) c(x$V1/(x$V1+x$V2+1), (x$V1+x$V2)))
##   with(resp1, {-sum(((V1*log(V1+.0001)+(1-V1)*log(1-V1+0.0001)))*V2)})
## }
