fit_map = function(Data_Map, data, varargin){

  map <- Data_Map$map

  if (map$map.preproc){

  data <- Data_Map$data

  }
  n <- nrow(data)

  m <- ncol(data)

  # Default values of customisable variables

  constStretching <- 0.7

  constBending <- 0.7

  weightss <- NULL

  func <- NULL

  intervals <- NULL

  nInt <- 5

  delta <- 1

  for(i in seq(1,length(varargin),2))  {

    if (varargin[i] == 'type') {

      if (varargin[i+1]=='hard')

      {constStretching <- 1;constBending <- 1}

      else if (varargin[i+1]=='medium')

      {constStretching <- 0.7;constBending <- 0.7}

      else if (varargin[i+1]=='soft')

      {constStretching <- 0.5;constBending <- 0.5}

      else   {stop('Incorrect value for type argument')}

    }

    if (varargin[i] == 'stretch') { tmp <- as.numeric(varargin[i+1]) ;constStretching <- tmp }

    if (varargin[i] == 'bend') { tmp <- as.numeric(varargin[i+1]);constBending <- tmp }

    if (varargin[i] == 'weights') { weights <- varargin[i + 1] }

    if (varargin[i] == 'intervals') { intervals <- varargin[i + 1] }

    if (varargin[i] == 'Number_of_intervals') { nInt <- unlist(varargin[i + 1]) }

    else if (varargin[i] == 'intshrinkage') { delta <- unlist(varargin[i + 1]) }

    else if (varargin[i] == 'potential') { func <- unlist(varargin[i + 1]) }

  }

  if (is.null(weightss))  {

    weightss <- rep(1, n)

  } else { weightss <- c(weightss) }

  # Define total weights

  TotalWeight <- sum(weightss)

  weigh <- weightss

  pFunc <- NULL

  if (!is.null(func)) {

    if (is.null(intervals)) {

      # Function has to create intervals by automatic way
      # nInt must be positive integer scalar

      if (is.complex(nInt) || !is.finite(nInt) || nInt < 1) {

        stop ('Incorrect value of "number_of_intervals"
             argumentIt must be positive integer scalar')
      } else {nInt = floor(nInt) }

      # delta has to be positive real scalar

      if (is.complex(delta) || !is.finite(delta) || delta < 0) {

        stop('Incorrect value of "intshrinkage" argument
             It must be positive real scalar') }

      pr <- nInt - 1

      # Intervals is the product of row and
      # maximal coefficient multiplied by delta

      intervals <- (map$map.disp * delta) * ((0:pr) / pr) ^ 2

      potentialFunction.intervals <- c(intervals, Inf)

      potentialFunction.sqint <- potentialFunction.intervals ^ 2

      # Get dimensions of intervals

      p <- length(intervals)

      # Preallocate memory

      A <- rep(0,p)

      B <- rep(0,p)

      if (func == 'L2' ){ pxk <- intervals^2 }

      if (func =='L1_5') {pxk <- abs(intervals)^1.5 }

      if (func == 'L1') {pxk <- abs(intervals)}

      if (func == 'LLog') {pxk <- log1p(intervals)}

      if (func == 'LSqrt') {pxk <- sqrt(intervals)}

      sxk <- intervals^2

      A [1:p-1] <- (pxk[1:p-1] - pxk[2:p])/(sxk[1:p-1] - sxk[2:p])

      B [1:p-1] <- (pxk[2:p]*sxk[1:p-1] - pxk[1:p-1]*sxk[2:p])/(sxk[1:p-1] - sxk[2:p])

      B[p] <- pxk[p]

      potentialFunction.B <- B

      potentialFunction.A <- A

      pFunc <- list(potentialFunction.sqint = potentialFunction.sqint,
                    potentialFunction.intervals = potentialFunction.intervals,
                    potentialFunction.A = potentialFunction.A,
                    potentialFunction.B = potentialFunction.B )


    } else {

      # intervals must contains non nerative values in ascending order.
      # The first value must be zero.

      if (intervals[1]!=0||!all(is.finite(intervals))||
          any((intervals[2:length(intervals)]-intervals[1:length(intervals)-1])<=0))

      { stop('Incorrect values in argument intervals: intervals must
             contains finite non negative values in ascending order.
             The first value must be zero.') }

      pFunc.intervals <- c(t(c(intervals)), Inf)

      # Get dimensions of intervals

      p <- length(intervals)

      # Preallocate memory

      A <- rep(0,p)

      B <- rep(0,p)

      if (func == 'L2' ){ pxk <- intervals^2 }

      if (func =='L1_5') {pxk <- abs(intervals)^1.5 }

      if (func == 'L1') {pxk <- abs(intervals)}

      if (func == 'LLog') {pxk <- log1p(intervals)}

      if (func == 'LSqrt') {pxk <- sqrt(intervals)}

      sxk <- intervals^2

      A [1:p-1] <- (pxk[1:p-1] - pxk[2:p])/(sxk[1:p-1] - sxk[2:p])

      B [1:p-1] <- (pxk[2:p]*sxk[1:p-1] - pxk[1:p-1]*sxk[2:p])/(sxk[1:p-1] - sxk[2:p])

      B[p] <- pxk[p]

      pFunc.A <- A
      pFunc.B <- B

    }
  }

  # Get initial state of nodes

  nodes <- map$map.mapped

  N <- nrow(nodes)

  # Form matrices B and C

  tmp <- map$map.links

  B <- diag(histc(c(tmp), 1:N)$cnt)

  tmp <- accumarray(tmp, rep(1, nrow(tmp)), c(N, N))

  B = B - tmp - t(tmp)

  tmp <- map$map.ribs

  C1 <- histc(c(tmp[, 1], tmp[, 3]), 1:N)$cnt

  C2 <- rep(0, N)

  C2[1:N-1] <- accumarray(tmp[, 2], rep(4, length(tmp[, 2])))

  C <- diag(C1 + C2)

  w <- accumarray(tmp[, c(1,3)], rep(1, length(tmp[, 1])), c(N, N))

  tmp <- accumarray(rbind(tmp[, 1:2], tmp[, 2:3]), rep(2, 2*length(tmp[, 1])), c(N, N))

  C = C + w + t(w) - tmp - t(tmp)

  # Start iterative process

  epoch <- 1   # Number of iteration

  ass <- rep(0,n)   # Initial associations. It is impossible combination

  qInd <- rep(0,n)

  # Get initial modulo

  stretch <- constStretching
  bend <- constBending

  while (TRUE) {

    # Save old associations and q indices.
    oldAss <- ass;
    oldQInd <- qInd;

    # Find new associations
    x1 <- rowSums(data^2)
    x2 <- t(rowSums(nodes^2))

    x1 <- matrix(x1)
    x1 <- kronecker(matrix(1, 1, dim(x2)[2]), x1)

    x2 <- kronecker(matrix(1, dim(x1)[1],1), x2)

    dist1 <- bsxfun('+', x1, x2) - 2*(data%*%t(nodes))

    dist <- apply(dist1, 1, min)

    k <- nrow(dist1)

#    klas <- rep(0, k)

    klas <- apply(dist1, 1,which.min)

  #  for (i in 1:k) {klas[i] = which.min(dist1[i, ])}

    ass <- klas

    # Find indeces for PQSQ if required

    if (!is.null(pFunc)){

      qInd <- histc(dist, pFunc$potentialFunction.sqint)$bin

      weigh <- weightss * t(pFunc$potentialFunction.A[qInd])

    }

    # If nothing is changed then we have end of epoch
    if (all(oldAss == ass) && all(oldQInd == qInd)){

      epoch <- epoch + 1

      tmp <- constStretching

      tmp1 <- constBending

      if (tmp == 0 && tmp1 == 0) {break}

      if (abs(tmp - stretch) + abs(tmp1 - bend) == 0){break}

      stretch <- tmp;
      bend <- tmp1;
    }

    # Form matrix A
    # For further robust and so on we consider possibility of zeros in
    # ass and create dummy element

    ass <- ass + 1;

    # Calculate number of points for each node

    tmp <- accumarray(ass, weigh, N + 1)

    # Normalise and remove dummy element

    NodeClusterRelativeSize <- tmp[2:length(tmp)] / TotalWeight

    # Create centroids

    NodeClusterCenters <- matrix(0, N + 1, m)

    for (k in 1:m) {

      NodeClusterCenters[, k] <- accumarray(ass, data[, k] * weigh, N + 1) / TotalWeight

    }

    # Remove dummy element

    NodeClusterCenters <- NodeClusterCenters[2:nrow(NodeClusterCenters), ]

    # Form SLAE

    SLAUMatrix <- diag(NodeClusterRelativeSize) + stretch * B + bend * C

    nodes <- solve(SLAUMatrix, NodeClusterCenters,tol=.Machine$double.eps)

    # Restore ass

    ass <- ass - 1

  }

  Data_Map$map$map.mapped <- nodes

  return(Data_Map)


}
