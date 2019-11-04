\name{quake}

\docType{data}

\alias{quake}
\alias{plate}

\title{Geographical locations of earthquakes}

\description{
The \code{quake} data set contains the geographical locations 
of severe earthquakes in the years 100 and 2016 inclusive. The
\code{plate} data set contains the geographical locations of the tectonic
plate boundaries.    
}

\usage{data(quake)
data(plate)}

\format{--For \code{quake}, a matrix with 5871 rows and  5 columns.
  Each row corresponds to an earthquake.
  The first column is the year (negative years indicate B.C.E.),
  the second is the longitude (decimal degrees),
  the third is the latitude (decimal degrees),
  the fourth is the depth beneath the Earth's surface (km),
  the fifth is a flag for the location inside the circum-Pacific belt
  (aka Pacific Ring of Fire). 

  --For \code{plate}, a matrix with 3 columns and 6276 rows.
  Each row corresponds to an location of the tectonic plate boundaries.
  The first is the longitude,
  the second is the latitude,
  the third is the label of the tectonic plate.} 

\source{
Bird, P. (2003) An updated digital model of plate boundaries,
\emph{Geochemistry, Geophysics, Geosystems} \bold{4(3)}, 1-52. 1027.
%Data set accessed 2016-03-24 from
%\url{http://peterbird.name/publications/2003_PB2002/2003_PB2002.htm}.
  
NGDC/WDS (2017) Global significant earthquake database, National
  Geophysical Data Center, NOAA, doi:10.7289/V5TD9V7K.
National Geophysical Data Center/World Data Service. Accessed 2017-03-30.
}


\keyword{datasets}

