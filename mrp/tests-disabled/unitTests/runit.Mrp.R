## Use case: happy path

## 1) create object with data.
##    - implicit is the validity check
## 2) run multilevel regression (mr)
##    - check that multilevel regression is done correctly
##    - this means running the "standard" model
## 3) poststratify
##    - check that average values match up
##    - verify that we can graph every thing properly.
##
createFakeData <- function (n=500) {
    response <- rep(c(1,0), n)
    var1 <- factor (rep (state.abb, length.out=n))
    var2 <- factor (rep (1:3, length.out=n))
    var3 <- factor (rep (1:5, length.out=n))
    weight <- rep (1, length.out=n)
    
    return (data.frame (response=response, var1=var1, var2=var2, var3=var3, weight=weight))
}

createFakePopulation <- function () {
    population <- expand.grid (state.abb, factor (1:3), factor(1:5))
    population <- cbind (population, rep(1, nrow(population))) 
    names (population) <- c("var1", "var2", "var3", "population")
    
    return (population)
}


test.creation <- function () {
    fakeData <- createFakeData()
    fakePopulation <- createFakePopulation()
    obj <- mrp (formula = response ~ var1 + var2 + var3,
                poll=fakeData,
                poll.weights="weight",
                population=fakePopulation,
                use="population")
    
    checkEqualsNumeric (3, getNumberWays(obj)["poll"])
    checkTrue (all (getPopulation (obj)@.Data == 1))
}

# Current version of the code does not allow for non 0/1 responses.
test.creation.failure <- function () {
    fakeData <- createFakeData()
    # This is what we are testing: non 0/1 response
    fakeData$response <- factor (fakeData$response, labels=c("no", "yes"))
    
    checkException (mrp (formula = response ~ var1 + var2 + var3,
            poll=fakeData))
}


test.multilevelRegression <- function () {
    fakeData <- createFakeData()
    fakePopulation <- createFakePopulation()
    
    obj <- mrp (formula = response ~ var1 + var2 + var3,
            poll=fakeData,
            poll.weights="weight",
            population=fakePopulation,
            use="population")
    
    checkEqualsNumeric (0.5, poststratify (obj), tolerance=1e-5)
    checkEqualsNumeric (rep (0.5, nlevels (fakeData$var2)*nlevels(fakeData$var3)), poststratify (obj, ~var2+var3), tolerance=1e-5)
    
    checkEqualsNumeric (fitted (obj@multilevelModel), poststratify (obj, ~var1+var2+var3), tolerance=1e-5)
}
