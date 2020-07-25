\name{EM}
\alias{EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Elastic Map
%%  ~~function to do ... ~~
}
\description{EM is function which fit parameters of elastic map. The main idea of this approach is to search map as solution of optimization problem.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EM(Data_Map, data, varargin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data_Map}{is the list of preprocessed data and initialized map variable that are obtaind with \code{\link{init}} function}
  \item{data}{The data is n-by-m matrix with m data points and n coordinates for each point (each column is one data point)}
  \item{varargin}{is a list of fitting options that can be one of following strings:

  \describe{

\item{type:}{ can be hard, medium or soft. "hard" is hard map with bend = 1 and stretch = 1. "medium" is flexible map with bend = 0.7 and stretch = 0.7 and "soft" is soft map with bend = 0.5 and stretch = 0.5}

\item{stretch:}{ is a positive numeric value which is the value of stretching modulo, "bend" is a positive numeric value which is the value of the bending modulo.}

\item{weights:}{ is n-by-1 vector of weights for data points. Weights must be nonnegative.}

\item{Number_of_intervals:}{ specifies the number of intervals to automatic interval calculation. Default value is 5. Ignored if "potential" is not specified.}

\item{intshrinkage:}{ is fraction of maximal distance from data points to original map which is used for intervals shrinkage. Default value is 1 (no shrinkage). Ignored if "potential" is not specified.}

\item{potential:}{ is majorant function for PQSQ. L2 distance without shrinkage is used if 'potential' is not specified.}

}
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{it is a List same as Data_Map variable that map.mapped is changed because of Elastic Map training.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zahra Reazaei
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Data Points
data <- BreastCancer_Data

# Create Map Geometry
map <- rect2Dmap(10,10)

# initialize Map
Data_Map <- init(map, data,"Pci")

# Elastic Map Training
varargin <- list("stretch", 0.01,"bend", 0.1)
Data_Map <- EM(Data_Map, data, varargin)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
