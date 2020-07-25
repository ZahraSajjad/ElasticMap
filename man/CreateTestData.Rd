\name{CreateTestData}
\alias{CreateTestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Test Data Points
}
\description{This function creates test data point by use of uniform distribution and dispersion, and then transform by the rotation matrix along coordinate axes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CreateTestData(dims, N, disper, angles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dims}{ dimension of data points
%%     ~~Describe \code{dims} here~~
}
  \item{N}{ Number of data points
%%     ~~Describe \code{N} here~~
}
  \item{disper}{ dispersion vector for each dimension
%%     ~~Describe \code{disper} here~~
}
  \item{angles}{ angles vector for each dimension
%%     ~~Describe \code{angles} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ dims-by-N matrix with N data points and dims coordinates for each point (each column is one data point)
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
Zahra Rezaei}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dims <- 100  # dimension of data points
N <- 5000     # Number of data points

# Define dispersion of data points
# disper is the vector with size dims.
# for each dimension we define different dispersion.
# we can either same amount of dispersion for each dimension.

disper1 <- rep(0.2,dims)
disper2 <- c(rep(0.2,30),rep(-0.2,40),rep(0.3,30))

# Define angles of data points
# angles is the vector of degrees with size dims
# angles vector creates transformation matrix
# definition of angles has same condition as dispersion.

angles1 <- rep(30,dims)
angles2 <- c(rep(30,30),rep(45,40),rep(-15,30))


data1 <- CreateTestData(dims, N, disper1, angles1)
data2 <- CreateTestData(dims, N, disper2, angles2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
