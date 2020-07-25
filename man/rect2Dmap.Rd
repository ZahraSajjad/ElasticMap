\name{rect2Dmap}
\alias{rect2Dmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Map Geometry
%%  ~~function to do ... ~~
}
\description{The function that creates two dimensional maps geometry with rectangular grid. contains two argument which are number of rows and number of columns.
}
\usage{
rect2Dmap(rows, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rows}{number of rows of Map Geometry}
  \item{cols}{number of columns of Map Geometry}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list that describe the feature of map geometry such as:
 \describe{
   \item{Internal coordinates:}{is the set of coordinates for each node in the map defined coordinates}
   \item{Map Links}{is the set of map edges. Each edge is the fragment of straight line which connect the nearest nodes}
   \item{Map Ribs}{is set of three adjacent nodes which belonged to one straight line in the internal coordinates. For triangular two dimension map there are three directions of ribs}
   \item{Map Faces}{is set of three nodes. We use this variable for 2d maps to project data points to the faces of map.}
   \item{Map Pcs}{is principal components of map that is null. calculated at initialization step with function (\code{\link{preprocessDataInit}}) }
   \item{Mapped coordinates}{is the set of coordinates of nodes in the data space that is null. These coordinates are initially defined by initializing procedure (\code{\link{init}}) and then adjusted by the map fitting. Procedure of map fitting is external with respect to the map and can be provided by SOM or EM fitting process.}
   \item{Map dispersion}{is dispersion measure for PQSQ approach. To calculate disp we have to calculate distance from each data point to the nearest initial node and then take maximum of these distances. Firstly it is null and then it is calculated at}
 }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zahra Rezaei
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Define dimension of map
rows <- 20
cols <- 20
map <- rect2Dmap(rows,cols)


map.internal
map$map.ribs
map$map.faces

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


