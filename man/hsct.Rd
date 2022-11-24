\name{hsct}
\docType{data}
\alias{hsct}
\title{Haematopoietic stem cell transplant}

\description{
  This data set contains the haematopoietic stem cell transplant (HSCT)
  measurements obtained by a flow cytometer from mouse subjects. A flow
  cytometer measures the spectra of fluorescent signals from biological
  cell samples to study their properties.    
}

\usage{data(hsct)}

\format{A matrix with 39128 rows and 6 columns.
  The first column is the FITC-CD45.1 fluorescence (0-1023), 
  the second is the PE-Ly65/Mac1 fluorescence (0-1023),
  the third is the PI-LiveDead fluorescence (0-1023),
  the fourth is the APC-CD45.2 fluorescence (0-1023),
  the fifth is the class label of the cell type (1, 2, 3, 4, 5),
  the sixth the mouse subject number (5, 6, 9, 12).
}

\source{
Aghaeepour, N., Finak, G., The FlowCAP Consortium, The DREAM Consortium,
Hoos, H., Mosmann, T. R., Brinkman, R., Gottardo, R. & Scheuermann, R. H. (2013) Critical assessment of automated flow cytometry data analysis techniques, \emph{Nature Methods} \bold{10}, 228-238.
}

\keyword{datasets}

