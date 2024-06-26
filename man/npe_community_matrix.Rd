% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{npe_community_matrix}
\alias{npe_community_matrix}
\title{Create a species abundance or occurrence matrix}
\usage{
npe_community_matrix(
  x,
  scale = "plot",
  trace_cover = 0.5,
  timescale = "annual",
  input = "neon_div_object",
  binary = FALSE
)
}
\arguments{
\item{x}{Input object. See input argument help for more details.}

\item{scale}{what level of aggregation? This can be "1m", "10m", "100m",
"plot", which is the default, or "site".}

\item{trace_cover}{cover value for subplots where only occupancy was recorded}

\item{timescale}{what temporal resolution? can be "subannual", which is really
only applicable at sites where there are multiple bouts per year, "annual" or
"all", which dissolves together the entire time series.}

\item{input}{by default, longform dataframe is calculated from the diversity
object and then converted to a community matrix, set this option to "lf"
to use a longform data frame that was created separately (and perhaps modified).
Another option is input = "divStack", which is using the output from the
divStack function in the neonPlants package. Using a premade longform data
frame or a divStack output will use the spatial and temporal scale of that
input data separately}

\item{binary}{should the matrix be converted from percent cover to binary?}

\item{neon_div_object}{the raw diversity data downloaded using
neonPlantEcology::download_plant_div() or the function
neonUtilities::loadByProduct() with the dpID arguement set to "DP1.10058.001".}
}
\value{
a data frame with each row a site aggregated at the spatial and
temporal scales chosen by the user. Each column is a single species, and cell
values can be either cover (a value between 0 and 100) or occurrence (1 or 0)
}
\description{
npe_community_matrix creates a wide matrix of species cover or binary (presence/absence)
values with the plot/subplot/year as rownames. This is useful for the vegan
package, hence the name.
}
\examples{
data("D14")
comm <- npe_community_matrix(D14)

}
