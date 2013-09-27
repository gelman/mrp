\name{dpsp}
\alias{dpsp}
\docType{data}
\title{
Demographically Purged State Predictor (DPSP)
}
\description{
The average desired state level intercept shift across a wide set of policies.
}
\usage{data(dpsp)}
\format{
  A data frame with 51 observations on the following 2 variables.
  \describe{
    \item{\code{state}}{a factor with levels \code{AK} \code{AL} \code{AR} \code{AZ} \code{CA} \code{CO} \code{CT} \code{DC} \code{DE} \code{FL} \code{GA} \code{HI} \code{IA} \code{ID} \code{IL} \code{IN} \code{KS} \code{KY} \code{LA} \code{MA} \code{MD} \code{ME} \code{MI} \code{MN} \code{MO} \code{MS} \code{MT} \code{NC} \code{ND} \code{NE} \code{NH} \code{NJ} \code{NM} \code{NV} \code{NY} \code{OH} \code{OK} \code{OR} \code{PA} \code{RI} \code{SC} \code{SD} \code{TN} \code{TX} \code{UT} \code{VA} \code{VT} \code{WA} \code{WI} \code{WV} \code{WY}}
    \item{\code{DPSP.20130104}}{a numeric vector}
  }
}
\details{
DPSP stands for demographically purged state predictor. The state-level predictors we use, such as
presidential vote or ideology, are usually things we think are correlated with the actual state-level true
values -- which are, after all, connected strongly to demographics -- rather than being directly correlated to
the state level random effects, which are the corrections to a purely demographic model. These intercept
shifts are to be the corrections to whatever the demographic and other variables would produce. Therefore, it
might be odd to use a model for them that assumes that the linear relationship between ultimate state opinion
and presidential vote is the same as the linear relationship between state-level corrections and presidential
vote. We constructed the DPSP we use herein from the full set of 39 survey sets in Lax and Phillips (2009a).
Very simply, this DPSP is the state random effects vector found by excluding all state level predictors and
running a somewhat standard model otherwise. This state random effects vector (from a model with 200K
observations across many survey questions) is the average desired state level intercept shift across a wide
set of policies. We find that DPSP does at least weakly better than other state predictors most of the time
and, indeed, when used, reduces the variation in state random effects in MRP applications, showing that
more state variation (at the level of corrections to demographic effects) is explained by DPSP than other
state-level predictors.
}
\source{
  \href{http://www.columbia.edu/~jrl2124/mrp2.pdf}{How Should We Estimate Sub-National Opinion
Using MRP?: Preliminary Findings and Recommendations}
}
\keyword{datasets}
