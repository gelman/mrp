\name{mrp.census}
\alias{mrp.census}
\docType{data}
\title{
US Census Public Use Microdata Sample (PUMS) for Poststratification
}
\description{
  Census data from 2000 with updated weights for estimates in 2004 and
  2008. \bold{Names and levels of \sQuote{poll} data must match these for poststratification.}
}
\usage{data(mrp.census)}
\format{
  A data frame with 1263600 observations on the following 14 variables.
  \describe{
    \item{\code{state}}{a factor with levels \code{AL} \code{AK} \code{AZ} \code{AR} \code{CA} \code{CO} \code{CT} \code{DE} \code{DC} \code{FL} \code{GA} \code{HI} \code{ID} \code{IL} \code{IN} \code{IA} \code{KS} \code{KY} \code{LA} \code{ME} \code{MD} \code{MA} \code{MI} \code{MN} \code{MS} \code{MO} \code{MT} \code{NE} \code{NV} \code{NH} \code{NJ} \code{NM} \code{NY} \code{NC} \code{ND} \code{OH} \code{OK} \code{OR} \code{PA} \code{RI} \code{SC} \code{SD} \code{TN} \code{TX} \code{UT} \code{VT} \code{VA} \code{WA} \code{WV} \code{WI} \code{WY}}
    \item{\code{race}}{a factor with levels \code{White} \code{Black} \code{Hispanic} \code{Other}}
    \item{\code{income}}{a factor with levels \code{less than $20,000} \code{$20,000-$40,000} \code{$40,000-$75,000} \code{$75,000-$150,000} \code{$150,000+}}
    \item{\code{age}}{a factor with levels \code{18-29} \code{30-44} \code{45-64} \code{65+}}
    \item{\code{sex}}{a factor with levels \code{Male} \code{Female}}
    \item{\code{education}}{a factor with levels \code{less than a high school education} \code{high school graduate} \code{some college} \code{college graduate} \code{postgraduate}}
    \item{\code{married}}{a factor with levels \code{Married} \code{Single}}
    \item{\code{kid}}{a factor with levels \code{has kids} \code{does not have kids}}
    \item{\code{raw2000}}{a numeric vector}
    \item{\code{weighted2000}}{a numeric vector}
    \item{\code{raw2004}}{a numeric vector}
    \item{\code{weighted2004}}{a numeric vector}
    \item{\code{raw2008}}{a numeric vector}
    \item{\code{weighted2008}}{a numeric vector}
  }
}
% \details{
  % ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  \url{http://www.census.gov/main/www/pums.html}, or see \url{http://dataferrett.census.gov/}
}
%\references{
%}
\examples{
data(mrp.census)
}
\keyword{datasets}
