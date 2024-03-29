% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individualshrubfuel.R
\name{individualshrubfuel}
\alias{individualshrubfuel}
\alias{individualshrubarea}
\alias{individualshrubvolume}
\title{Shrub fuel biomass}
\usage{
individualshrubfuel(x, type = "total", allometric = TRUE,
  excludeSSP = TRUE, equationNames = FALSE, var = FALSE,
  agg = "none", customParams = NULL, outside = "warning",
  na.rm = TRUE)

individualshrubarea(x, excludeSSP = TRUE, var = FALSE,
  customParams = NULL, outside = "warning", na.rm = TRUE)

individualshrubvolume(x, agg = "none", na.rm = TRUE)
}
\arguments{
\item{x}{data frame with columns 'plot', 'species', 'H' (height in cm), 'D1' and 'D2' (in cm). If 'D2' is ommitted then
shrub crowns are assumed to be circular (i.e. D2 = D1). Column 'group' may be used to specify a functional group to be used when
the species has not a specific allometry (see details).}

\item{type}{'total'  (total fuel) or 'fine' (fine fuels)}

\item{allometric}{wether to use allometric equations or bulk density estimates}

\item{excludeSSP}{excludes subspecies information for species matching}

\item{equationNames}{a flag to indicate that the name of species or group allometry or bulk density used is desired}

\item{var}{a flag to indicate that variance of estimates is desired}

\item{agg}{aggregation of results. Either 'none', 'species', 'speciesplot', 'plotspecies' or 'plot'}

\item{customParams}{custom allometry parameter table (for species not in default params)}

\item{outside}{Treament of values outside calibration range: either 'warning' (to prompt a warning) or 'missing' (prompt a warning and return NA)}

\item{na.rm}{whether to exclude missing values when aggregating biomass}
}
\value{
a vector of dry weight (biomass) in kg
}
\description{
Calculates dry weight (biomass, in kg) of total or fine fuels corresponding to individual shrub data
}
\details{
The function determines the allometry to be applied using the following rules, sequentially:
\enumerate{
  \item{If the species name is included in the list of species with parameter values, it takes the parameters from that species.}
  \item{If the species name is a synonym for a species included the list of species with parameter values, it takes the parameters from that species (and gives a warning).}
  \item{If the user has specified a column 'group' and the value is included in the list of valid groups, it takes the parameters from that group.}
  \item{If the species is listed within the checklist in 'species_groups', it takes the parameters from the group specified in that file.}
  \item{Otherwise, it gives a warning and no allometry is applied.}
}
}
\examples{
plot = c(1,1,2,2,2)
species = c("Erica arborea","Genista cinerea", "Erica arborea", "Chamaerops humilis", "Genista balansae")
H = c(60,70,200,100,10)
D1 = c(60,40,100, 100,25)
D2 = D1
x = data.frame(plot, species, H, D1, D2, stringsAsFactors=FALSE)

individualshrubfuel(x)
}
