\name{initialize_map}
\alias{initialize_map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Map Initializing
%%  ~~function to do ... ~~
}

\description{
Function of map initialization. This function perform required data preprocessing (see functions \code{\link{preprocessDataInit}} and \code{\link{preprocessData}}). This function defines an initial mapped coordinates.
}

\usage{
init(map, data, type, reduce)
}

\arguments{
  \item{map}{The map is created by Map Geometry Functions (\code{\link{rect2Dmap}})}
  \item{data}{The data is n-by-m matrix with m data points and n coordinates for each point (each column is one data point)}
  \item{type}{The type of initialization: principal component initialization (Pci) or Random.}
  \item{reduce}{is integer. If reduce is positive and is less than n then specified number of the first principal components are used. If 'reduce' is zero and m>n then the first n-1 principal components is used. If 'reduce' is positive and is greater or equal to n or 'reduce' is zero and n>m then dimensionality reduction is not performed. If reduce is negative then reduce PCs are calculated but dimensionality reduction is not performed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ zahra Rezaei
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rect2Dmap}}, \code{\link{preprocessDataInit}} , \code{\link{preprocessData}}
}
\examples{
data <- BreastCancer_Data
map <- rect2Dmap(10,10)
Data_Map <- init(map, data,"Pci")

# return principal component of map and map coordinate in the data space

Data_Map$map$map.Pcs
Data_map$map$map.mapped


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
