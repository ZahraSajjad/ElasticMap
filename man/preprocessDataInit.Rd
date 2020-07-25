\name{preprocessDataInit}
\alias{preprocessDataInit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data Preprocessing
%%  ~~function to do ... ~~
}
\description{the function to perform data preprocessing . to calculate map principal components (Orthonormal Matrix)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
preprocessDataInit(map,data,reduce)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{a list that was created by Map Geometry Functions \code{\link{rect2Dmap}}}
  \item{data}{The data is n-by-m matrix with m data points and n coordinates for each point (each column is one data point)}
  \item{reduce}{is integer. If reduce is positive and is less than n then specified number of the first principal components are used. If 'reduce' is zero and m>n then the first n-1 principal components is used. If 'reduce' is positive and is greater or equal to n or 'reduce' is zero and n>m then dimensionality reduction is not performed. If reduce is negative then reduce PCs are calculated but dimensionality reduction is not performed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with two variable:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{map}{update map that was created before with map geometry function, add value to map Pcs }
  \item{data}{is the result of data preprocessing after we calculate principal components of map. we call \code{\link{preprocessData}}} function to calculation.
%% ...
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
\code{\link{init}},\code{\link{preprocessDataInit}},\code{\link{preprocessData}},


}
\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
