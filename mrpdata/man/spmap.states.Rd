\name{spmap.states}
\alias{spmap.states}
\docType{data}
\title{U.S. States and Territories map data}
\description{
A US map with state abbreviations and FIPS codes suitable for plotting
poststratified survey results.}
\usage{data(spmap.states)}
\format{
  A data frame with 56 observations on the following 5 variables.
  \describe{
    \item{\code{STATE}}{a factor with levels \code{AK} \code{AL} \code{AR} \code{AS} \code{AZ} \code{CA} \code{CO} \code{CT} \code{DC} \code{DE} \code{FL} \code{GA} \code{GU} \code{HI} \code{IA} \code{ID} \code{IL} \code{IN} \code{KS} \code{KY} \code{LA} \code{MA} \code{MD} \code{ME} \code{MI} \code{MN} \code{MO} \code{MS} \code{MT} \code{NC} \code{ND} \code{NE} \code{NH} \code{NJ} \code{NM} \code{NV} \code{NY} \code{OH} \code{OK} \code{OR} \code{PA} \code{PR} \code{RI} \code{SC} \code{SD} \code{TN} \code{TX} \code{UM} \code{UT} \code{VA} \code{VI} \code{VT} \code{WA} \code{WI} \code{WV} \code{WY}}
    \item{\code{NAME}}{a factor with levels \code{Alabama} \code{Alaska} \code{American Samoa} \code{Arizona} \code{Arkansas} \code{California} \code{Colorado} \code{Connecticut} \code{Delaware} \code{District of Columbia} \code{Florida} \code{Georgia} \code{Hawaii} \code{Idaho} \code{Illinois} \code{Indiana} \code{Iowa} \code{Kansas} \code{Kentucky} \code{Louisiana} \code{Maine} \code{Maryland} \code{Massachusetts} \code{Michigan} \code{Minnesota} \code{Mississippi} \code{Missouri} \code{Montana} \code{Nebraska} \code{Nevada} \code{New Hampshire} \code{New Jersey} \code{New Mexico} \code{New York} \code{North Carolina} \code{North Dakota} \code{Ohio} \code{Oklahoma} \code{Oregon} \code{Pennsylvania} \code{Puerto Rico} \code{Rhode Island} \code{South Carolina} \code{South Dakota} \code{Tennessee} \code{Texas} \code{Utah} \code{Vermont} \code{Virgin Islands} \code{Virginia} \code{Washington} \code{West Virginia} \code{Wisconsin} \code{Wyoming}}
    \item{\code{FIPS}}{a factor with levels \code{00} \code{01} \code{02} \code{04} \code{05} \code{06} \code{08} \code{09} \code{10} \code{11} \code{12} \code{13} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22} \code{23} \code{24} \code{25} \code{26} \code{27} \code{28} \code{29} \code{30} \code{31} \code{32} \code{33} \code{34} \code{35} \code{36} \code{37} \code{38} \code{39} \code{40} \code{41} \code{42} \code{44} \code{45} \code{46} \code{47} \code{48} \code{49} \code{50} \code{51} \code{53} \code{54} \code{55} \code{56} \code{60} \code{72} \code{78}}
    \item{\code{LON}}{a numeric vector}
    \item{\code{LAT}}{a numeric vector}
  }
}
\details{
  A \link[sp]{SpatialPolygonsDataFrame}; the \sQuote{data} slot is
  described above. The polygons have been projected using an Albers
  Equal-Area Conic projection with parameters equal to those suggested
  by USGS for the contiguous US: (A PROJ4 specification of \code{+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +ellps=WGS84}).
}
\source{
\url{http://www.weather.gov/geodata/catalog/national/html/us_state.htm}
}
\references{
  \url{http://nationalatlas.gov/atlasftp.html}; for projection see
  \url{http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html}
  and to apply other projections, see \link[rgdal]{spTransform}.
  and the Spatial CRAN task view (\url{http://cran.r-project.org/web/views/Spatial.html}) }
\examples{
data(spmap.states)
}
\keyword{datasets}
