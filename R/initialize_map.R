initialize_map = function(map, data, type, reduce){

  n <- nrow(data)

  m <- ncol(data)

  reduce <- floor(reduce)

  if (reduce >= m || (reduce == 0 && n > m)){

    # Calculate 3 PCs
    k <- 3

    if (k > m){ k = m }

    map$map.preproc <- FALSE

  }else if (reduce < 0){

    k <- -reduce

    if (k > m){ k = m }

    map$map.preproc <- FALSE

  }else {

    # Define required number of PCs
    k <- n - 1

    if (reduce > 0 && reduce > k) {k = reduce}

    map$map.preproc <- TRUE

  }

  # Search required number of PCs

  map.mean <- colMeans(data)

  MeanRep <- kronecker(matrix(1, n, 1), t(map.mean))

  Element <- bsxfun('-', data, MeanRep)

  tmp <- svdr(Element, k)

  D <- tmp$d

  V <- tmp$v

  Dsort <- sort(D, index.return = TRUE, decreasing = TRUE)

  ind <- Dsort$ix

  V <- V[, ind]

  # Standardise Direction of PCs

  ind <- (diag(V)<0)

  V[,ind] <- -V[, ind]

  # Preprocess Data

  if (map$map.preproc == TRUE){

    data <- bsxfun('-', data, MeanRep) %*% V

  }

  # Store Result

  map$map.mean <- map.mean

  map$map.Pcs <- V

  if (type == "random")  {

    # Random generation
    # Calculate intervals for each coordinates

    mini <- apply(data, 2, min)

    maxi <- apply(data, 2, max) - mini

    # Generate random coordinates
    data <- rand(nrow(map$map.internal), ncol(data))

    # Scale coordinates
    data <- bsxfun('*', data, kronecker(matrix(1, nrow(data), 1), t(maxi)));

    # shift coordinates
    map$map.mapped <- bsxfun('+', data, kronecker(matrix(1, nrow(data), 1), t(mini)))
  }

  if (type == "Pci")  {

    # Principal component initialization

    # Calculate mean and PCs

    if (map$map.preproc) {

      # Data were preprocessed

      V <- diag(1, ncol(data), map$map.dimension)

      tmp <- data[, 1:map$map.dimension]

      MeanDat <- zeros(1, size(data, 2))

    } else {

      # Get required number of Pcs

      MeanDat <- map$map.mean

      V <- map$map.Pcs[, 1:map$map.dimension]

      tmp <- data %*% V

    }

    # Calculate mean and dispersion along each PCs

    mini <- apply(tmp, 2, min)

    maxi <- apply(tmp, 2, max)

    meant <- colMeans(tmp)

    M1 <- matrix(c(meant-mini, maxi-meant), nrow = 2, byrow = TRUE)

    disper <- apply(M1, 2, min)

    # Calculate mean and dispersion along internal coordinates

    minI <- apply(map$map.internal, 2, min)

    maxI <- apply(map$map.internal, 2, max)

    meanI <- colSums(map$map.internal) / dim(map$map.internal)[1]

    M2 <- matrix(c(meanI-minI, maxI-meanI), nrow = 2, byrow = TRUE)

    disP <- apply(M2, 2, min)

    # Auxiliary calculations

    V <- bsxfun('*', V, kronecker(matrix(1, nrow(V), 1), t(disper / disP)))

    # final values

    meanI <- kronecker(matrix(1, dim(map$map.internal)[1], 1), t(as.matrix(meanI)))

    MeanDat <- kronecker(matrix(1, dim(map$map.internal)[1],1), (as.matrix(MeanDat)))

    map$map.mapped <- bsxfun('+', bsxfun('-', map$map.internal, meanI) %*% t(V), MeanDat)

  }

  # calculate map disp

  x1 = rowSums(data^2)

  x2 = t(rowSums(map$map.mapped^2))

  x1 = matrix(x1)

  x1 = kronecker(matrix(1, 1, dim(x2)[2]), x1)

  x2 = kronecker(matrix(1, dim(x1)[1], 1), x2)

  dist1 = bsxfun('+', x1, x2) - 2*(data%*%t(map$map.mapped))

  dist = apply(dist1, 1, min)

  map$map.disp <- sqrt(max(dist))

  Data_Map <- list(data = data, map = map)

  return(Data_Map)

}
