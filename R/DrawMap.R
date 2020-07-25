DrawMap = function(Data_Map, data, projectType, varargin){

  map <- Data_Map$map

  if (map$map.preproc){

    data <- Data_Map$data

  }

  N <- nrow(data)

  dims <- ncol(data)

  # Parse varargin
  Source <- NULL

  smooth <- 0.15


  for (i in 1:length(varargin)){

    if (names(varargin)[i] == 'classes'){ classes <- varargin$classes }
    if (names(varargin)[i] == 'markColour'){ markColour <- varargin$markColour }
    if (names(varargin)[i] == 'coloring'){ Source <- varargin$coloring }
    if (names(varargin)[i] == 'ColorMap'){ colMap <- varargin$ColorMap }
    if (names(varargin)[i] == 'flatColoring'){ flatColoring <- varargin$flatColoring }
    if (names(varargin)[i] == 'Surf'){ Surf <- varargin$Surf }

  }
  if (isempty(classes)) {
    classes <- ones(size(data, 1), 1)
    cls <- 1
  } else { cls <- unique(classes) }
  nCls <- length(cls)

  # Get map internal
  intern <- map$map.internal

  # Get map links
  links <- map$map.links

  # Get map faces
  faces <- map$map.faces

  # Get map mapped
  mapped <- map$map.mapped

  kind <-'internal'

  # Project data onto map
  dataP <- NULL

  if (!isempty(data)) {

    cType = strcmpi('mapped', kind)

    points <- data

    if (projectType == 0) {

      # Projection to the nearest node
      # Search the nearest node

      x1 <- rowSums(points^2)
      x2 <- t(rowSums(mapped^2))

      x1 <- matrix(x1)
      x1 = kronecker(matrix(1, 1, dim(x2)[2]), x1)

      x2 = kronecker(matrix(1, dim(x1)[1], 1), x2)

      dist1 <- bsxfun('+', x1, x2)-2*(points%*%t(mapped))

      dist <- apply(dist1,1,min)

      k <- nrow(dist1)

      num <- rep(0,k)

      for (i in 1:k) {num[i] <- which.min(dist1[i, ])}

      if (cType) { dataP <- mapped[num, ]

      } else { dataP <- intern[num, ] }

    }

    if (projectType == 1) {

      # projection onto nearest edge

      # Get array of edges end

      V2 <- t(mapped[links[, 2], ])

      # Form matrix of edges directions

      dir <- t(mapped[links[, 1], ]) - V2

      # Calculate squared length of edge directions

      len <- (colSums(dir^2))


      # Calculate projections length (l in documentation, matrix analogue of (2))

      pr <- bsxfun('-', points%*%dir, kronecker(matrix(1, nrow(points), 1), t(colSums(V2*dir))))

      # Copy projections to normalize (l* in documentation)

      prn <- bsxfun('/', pr, kronecker(matrix(1, nrow(points), 1), t(len)))

      # Non negativity

      prn[prn<0] <- 0

      #Cut too long projections (it is the same as step 3 of algorithm in documentation)

      prn[prn>1] <- 1

      # Calculate distances:

      x <- kronecker(matrix(1, 1, dim(links)[1]), rowSums(points^2))

      y <- kronecker(matrix(1, 1, N), as.matrix(colSums(V2^2)))

      dist <- bsxfun('+', x, t(y)) - 2*points%*%V2 +
        prn*( bsxfun('*', prn, kronecker(matrix(1, N, 1), t(len)))-2*pr)

      # Select the nearest edge

      edge <- apply(dist, 1, which.min)

      dist <- apply(dist, 1, min)

      # form index to find length of projections

      ind <- t(1:N) + (edge - 1)*N

      if (cType) {

        dataP <- bsxfun('*', kronecker(matrix(1, 1, ncol(mapped)), (1-prn[ind])), mapped[links[edge, 2], ])
        + bsxfun('*', kronecker(matrix(1, 1, ncol(mapped)), prn[ind]), mapped[links[edge, 1], ])

      } else {

        c11 <- kronecker(matrix(1,1,2), (1-prn[ind]))
        c12 <- kronecker(matrix(1,1,2), prn[ind])

        c21 <- intern[links[edge, 2], ]
        c22 <- intern[links[edge, 1], ]

        dataP <- bsxfun('*', c11, c21) + bsxfun('*', c12, c22)}

    }

    if (projectType == 2) {


      #form auxiliary vectors

      Y2 <- t(mapped[faces[, 3], ])
      Y20 <- Y2-t(mapped[faces[, 1], ])
      Y21 <- Y2-t(mapped[faces[, 2], ])
      Y10 <- t(mapped[faces[, 2], ] - mapped[faces[, 1], ])
      Y20Y20 <- colSums(Y20 ^ 2)
      Y21Y20 <- colSums(Y20 * Y21)
      Y21Y21 <- colSums(Y21 ^ 2)

      # Calculate projections

      A20 <- bsxfun('-', kronecker(matrix(1, nrow(points), 1), t(as.matrix(colSums(Y2*Y20)))), points%*%Y20)

      A21 <- bsxfun('-', kronecker(matrix(1, nrow(points), 1), t(as.matrix(colSums(Y2*Y21)))), points%*%Y21)

      A10 <- bsxfun('/', bsxfun('-', (A20-A21), kronecker(matrix(1, nrow(points), 1), t(as.matrix(colSums(Y21*Y10))))),
                   kronecker(matrix(1, nrow(points), 1), t(as.matrix(colSums(Y10^2)))))

      tmp <- Y20Y20*Y21Y21 - Y21Y20^2

      A0 <- bsxfun('/', bsxfun('*', A20,kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y21Y21))))
                  -bsxfun('*', A21, kronecker(matrix(1, nrow(points),1),t(as.matrix(Y21Y20)))),
                  kronecker(matrix(1, nrow(points),1), t(as.matrix(tmp))))
      A1 <- bsxfun('/',bsxfun('*', A21,kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y20Y20))))
                  -bsxfun('*',A20, kronecker(matrix(1, nrow(points),1), t(as.matrix(Y21Y20)))),
                  kronecker(matrix(1, nrow(points), 1), t(as.matrix(tmp))))

      A20 <- bsxfun('/', A20, kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y20Y20))))

      A21 <- bsxfun('/', A21, kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y21Y21))))

      # Normalize projections
      A20N <- A20
      A20N[A20N<0] <- 0
      A20N[A20N>1] <- 1
      A21N <- A21
      A21N[A21N<0] <- 0
      A21N[A21N>1] <- 1
      A10[A10<0] <- 0
      A10[A10>1] <- 1
      tmp <- A0<0
      A0[tmp] <- 0
      A1[tmp] <- A21N[tmp]
      tmp <- A1<0
      A0[tmp] <- A20N[tmp]
      A1[tmp] <- 0
      tmp <- (1-(A0+A1))<0
      A0[tmp] <- A10[tmp]
      A1[tmp] <- 1-A10[tmp]

      # Calculate distances

      dist1 <- bsxfun('+', kronecker(matrix(1, 1, nrow(faces)), rowSums(points^2)),
                     kronecker(matrix(1, nrow(points), 1), t(colSums(Y2^2)))) - 2*points%*%Y2

      dist2 <- bsxfun('*', A0, kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y20Y20))))*(A0-2*A20)

      dist3 <- bsxfun('*', A1, kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y21Y21))))*(A1-2*A21)

      dist4 <- 2 * bsxfun('*', A0*A1, kronecker(matrix(1, nrow(points), 1), t(as.matrix(Y21Y20))))

      dist <- dist1 + dist2 + dist3 + dist4

      tmp <- apply(dist, 1, which.min)

      dist <- apply(dist, 1, min)

      ind <- t(1:N) + (tmp - 1)*N

      if (cType){

        dataP <- bsxfun('*', A0[ind], mapped[faces[tmp, 1], ])
        + bsxfun('*', A1[ind], mapped[faces[tmp, 2], ])
        + bsxfun('*', 1-A0[ind]-A1[ind], mapped[faces[tmp, 3], ])

      } else {

        c1 <- bsxfun('*', kronecker(matrix(1,1,2), A0[ind]), intern[faces[tmp, 1], ])
        c2 <- bsxfun('*', kronecker(matrix(1,1,2), A1[ind]), intern[faces[tmp, 2], ])
        c3 <- bsxfun('*', kronecker(matrix(1,1,2), (1-A0[ind]-A1[ind])), intern[faces[tmp, 3], ])

        dataP <- c1 + c2 + c3

      }


    }
  }

  if (!isempty(Source)) {

    d <- ncol(mapped)

    # Extract grid

    for(i in 1:nrow(faces)) { faces[i, ] <- sort(faces[i, ]) }

    # Form list of all edges in grid
    edges <- rbind(faces[, 1:2], faces[, 2:3], faces[, c(1, 3)])

    # Search unique values
    edges <- as.matrix(distinct(tibble(edges)))

    # Step 2. Form list of nodes in new node list
    nN <- nrow(mapped)
    nE <- nrow(edges)
    maps <- rbind(mapped, (mapped[edges[, 1], ] + mapped[edges[, 2], ]) / 2)
    inter <- rbind(intern, (intern[edges[, 1], ] + intern[edges[, 2], ]) / 2)

    # Form list of indexes for nodes
    ind <- matrix(0, nE + nN,nE + nN)
    siz <- dim(ind)

    indL <- (edges[,1]) + ((edges[,2]) - 1)*siz[1]
    ind[indL] <- (nN + 1):(nN + nE)

    # Step 3. form list of six nodes for each face
    face <- cbind(faces, ind[(faces[, 1]) + ((faces[, 2]) - 1)*siz[1]],
                 ind[(faces[, 2]) + ((faces[, 3]) - 1)*siz[1]],
                 ind[(faces[, 1]) + ((faces[, 3]) - 1)*siz[1]])

    # Step 4. Form final list of triangles
    grid <- rbind(face[, c(2, 4, 5)], face[, c(1, 4, 6)],
                 face[, c(3, 5, 6)], face[, c(4, 5, 6)])

    gridReady <- grid
    nodeMap <- maps
    nodeInt <- inter

    if (is.numeric(Source)) {

      if (length(Source) == 1) {

        # It is number of coordinate or PCs
        Source <- round(Source)

        if (Source > 0 ) {

          # It is coordinate
          # Get data in original space

         temp <- bsxfun( '+', nodeMap %*% t(map$map.PCs), kronecker(matrix(1, nrow(nodeMap), 1),
                              matrix(map$map.means, nrow = 1, byrow = TRUE)))

          t <- mapped[, Source]

           if (Source > ncol(temp)) {f <- temp[, Source]}
        }

        if (Source > 0 ) {

          if (map$map.preproc) {

            # Get required coordinate

            f <- nodeMap[, -source]

          }else {

            # Calculate projection on required PC

            f <- bsxfun('-', nodeMap, map.means)* map$map.PCs[, -source]
          }
         }

      }

      if (is.array(Source))    {

        Source <- c(Source)
        if (length(Source) != nrow(data)){
          stop('Number of elements in vector "coloring" ', nrow(Source),
               ' must coincides with number of data points ', nrow(data))
        }

        smooth <- -1 / (smooth * mean(diag(var(dataP))));

        # Calculate distances from each node to each datapoint
        d1 <- kronecker(matrix(1,1,nrow(nodeInt)),as.matrix(rowSums(dataP^2)))
        d2 <- kronecker(matrix(1,nrow(dataP),1),as.matrix(t(rowSums(nodeInt^2))))

        dist <- bsxfun('+',d1,d2) - 2 * (dataP %*% t(nodeInt))

        # Calclulate RBF in each point
        tmp <- bsxfun('*', exp(dist * smooth), kronecker(matrix(1,1,ncol(dist)),Source))

        # Calculate result
        res <- t(colSums(tmp))
        mins <- min(res)
        res <- (res - mins) / (max(res) - mins)

        f <- res
      }

      }

    if (is.character(Source) && strcmp(Source,'density')){
    # Calculate variances for all attributes
    smooth <- -1 / (smooth * mean(diag(var(dataP))));

    # Calculate distances from each node to each datapoint
    d1 <- kronecker(matrix(1, 1, nrow(nodeInt)), as.matrix(rowSums(dataP^2)))
    d2 <- kronecker(matrix(1, nrow(dataP), 1), as.matrix(t(rowSums(nodeInt^2))))

    dist <- bsxfun('+',d1,d2) - 2 * (dataP %*% t(nodeInt))

    # Calclulate RBF in each point
    tmp <- bsxfun('*', exp(dist * smooth), kronecker(matrix(1, 1, ncol(dist)), matrix(1, nrow(data), 1)))

    # Calculate result
    res <- t(colSums(tmp))
    mins <- min(res)
    res <- (res - mins) / (max(res) - mins)

    # Draw colouring

    f <- res

    }

    m <- cbind(nodeInt,t(f))

    zmean <- apply((gridReady), MARGIN = 1, function(row){mean(m[row,3])})

    facecolor <- colour_ramp(viridis_pal()(8))( rescale( x = zmean ) )

    if (flatColoring) {

      # draw 2d surface

      if (colMap == 'continuous'){

        p <- plot_ly(x = intern[,1], y = intern[,2],type = 'scatter',

                    mode = 'markers', marker = list(size = 3, color = 'red'),

                    name = 'Map nodes')%>%

          add_trace(x = nodeInt[,1], y = nodeInt[,2], z = array(f/mean(f)),
                    type = 'heatmap',colors = facecolor, name = 'Heat Map Plot')%>%
          layout(title = ' Map Colouring ',xaxis = list(title = 'x'),
                 yaxis = list(title = 'y'))
      }

      if (colMap == 'discrete'){

        p <- plot_ly(x = intern[,1], y = intern[,2],type = 'scatter',

                    mode = 'markers',marker = list(size = 3,color = 'red'),

                    name = 'Map nodes')%>%
          add_trace(x = nodeInt[,1], y = nodeInt[,2], z = array(f/mean(f)),
                    type = 'contour', name = 'Contour Plot')%>%
          layout(title = ' Map Colouring ',xaxis = list(title = 'x'),
                 yaxis = list(title = 'y'))

      }

      for (i in 1:nCls){

        p <- add_trace(p,
                      x = dataP[which(classes == i),1],
                      y = dataP[which(classes == i),2],
                      type = 'scatter', mode = 'markers',
                      marker = list(color = markColour[i], size = 7, symbol = 201),
                      name = paste('Class',toString(i)))

      }

      return(p)

    }

    if (!flatColoring) {

      if (Surf == "Internal"){

      p <- plot_ly(x = nodeInt[,1], y = nodeInt[,2], z = array(f),type = 'mesh3d',
                  i = gridReady[,1]-1,
                  j = gridReady[,2]-1,
                  k = gridReady[,3]-1,
                  facecolor = facecolor)
      }

      if (Surf == "External"){

        p <- plot_ly(x = nodeMap[,1], y = nodeMap[,2], z = nodeMap[,3],
                    type = 'mesh3d',
                    i = gridReady[,1]-1,
                    j = gridReady[,2]-1,
                    k = gridReady[,3]-1,
                    facecolor = facecolor,
                    name = '3D surface')

        for (i in 1:nCls){

        p <- add_trace(p,
                         x = data[which(classes == i),1],
                         y = data[which(classes == i),2],
                         z = data[which(classes == i),3],
                         type = 'scatter3d',mode = 'markers',
                         marker = list(color = markColour[i], size = 6, symbol = 201),
                         name = paste('Class',toString(i)))
        }
      }

      return(p)

    }
  }

  if (isempty(Source)) {

  data <- dataP

  # Draw Map

  nrows <- map$map.sizes[1]
  ncols <- map$map.sizes[2]

  Vmap_xy <- matrix(NaN,(nrows*ncols+(ncols-1)),2)
  Hmap_xy <- matrix(NaN,1,2)


  for (i in 1:ncols){

    k1 <- 1 + (i-1)*(nrows+1)

    k2 <- 1 + (i-1)*(nrows)

    Vmap_xy[k1:(k1 + (nrows - 1)),] <- intern[k2:(k2 + (nrows - 1)), ]

  }

  for (i in 1:nrows){
    Hmap_xy <- rbind(Hmap_xy,
                    matrix(c(1:ncols, rep(i,ncols)),
                           ncol = 2, nrow = ncols), matrix(NaN,1,2))
  }

  map_xy <- rbind(Vmap_xy,Hmap_xy)

  # Project data to map

  if (projectType == 0) {

    lims <- matrix(0, 1, 4)

    lims[c(1, 3)] <- min(intern) - 0.5

    lims[c(2, 4)] <- max(intern) + 0.5

    plot(NULL, xlim = lims[1:2], ylim = lims[3:4], ylab="y", xlab="x")

    L <- map$map.sizes[1]

    for(i in 1:L) {

      points2D(intern[(L*i-L+1):(L*i), 1],
               intern[(L*i-L+1):(L*i), 2],
               xlim = lims[1:2], ylim = lims[3:4],
               type ='l', col = "red", add = TRUE)

    }

    for(i in 1:L) {

      points2D(intern[seq(i, (i + L*(L-1)),L), 1],
               intern[seq(i, (i + L*(L-1)),L), 2],
               type ='l',col = "black",add = TRUE)

    }

    if (!is.null(data)) {

      d <- 10*data[, 1] + data[, 2]

      # Search unique points

      d1 <- uniq(d)

      dat <- data[d1$m, ]

      ic <- d1$n

      count <- histc(ic, 1:length(d1$b))$cnt

      ma <- max(count)

      if (nCls==1) {   # No classes

        plot(dat[,1], dat[,2],
             pch = 19, col = 'green',cex = 1.35*count,
             xlim = c(1,10), ylim = c(1,10), xlab = 'X', ylab = 'Y')

      }

      nDat <- nrow(dat)

      Proportions <- matrix(0, nDat, nCls)

      for (k in 1:nCls){

        ind <- classes == cls[k]

        Proportions[,k] <- histc(ic[ind], 1:nDat)$cnt

      }

      for (k in 1:nDat){

          # Transform counts to proportions

          prop <- Proportions[k,] / sum(Proportions[k,])

          # Starts from zero angle

          angle <- 0

          # Draw chart

          if (length(prop) == 1){

            # One class chart

            for (i in 1:length(Proportions[k,])){

              color <- markColour[i,]

              t <- seq(0,2*pi,0.05)

              xp <- count[k]/(ma*2)*cos(t) + dat[k, 1]

              yp <- count[k]/(ma*2)*sin(t) + dat[k,2]

              polygon2D(xp, yp, col = color, xlim = lims[1:2], ylim = lims[3:4])

            }

          } else {

            # Several classes case

            for (i in 1:length(Proportions[k,])){

              th1 <- angle

              th2 <- angle + prop[i]*2*pi

              color <- markColour[i]

              t <- linspace(th1,th2)

              xp <- count[k]/(ma*2)*cos(t) + dat[k, 1]

              yp <- count[k]/(ma*2)*sin(t) + dat[k,2]

              xp <- c(xp,dat[k, 1],xp[1])

              yp <- c(yp,dat[k,2],yp[1])

              polygon2D(xp, yp, col = color, add = TRUE, xlim = lims[1:2], ylim = lims[3:4])

              angle <- th2

            }

          }
      }

    }

  }

  if (projectType == 1) {

    p <- plot_ly(x = map_xy[, 1],
                y = map_xy[, 2],
                type = 'scatter', mode ='lines',
                marker = list(color = 'red'),
                name = 'Map frame')

    for (i in 1:nCls){

    p <-  add_trace(p,
                     x = data[which(classes == i),1],
                     y = data[which(classes == i),2],
                     type = 'scatter',mode = 'markers',
                     marker = list(color = markColour[i], size = 8, symbol = 201),
                     name = paste('Class',toString(i)))%>%
        layout(title = "Projection of data points to the nearest edge",
               xaxis = list (title = "x"),
               yaxis = list (title = "y"))
    }

     return(p)

  }

  if (projectType == 2) {

     p <- plot_ly (x = map_xy[, 1],
                  y = map_xy[, 2],
                  type = 'scatter', mode ='lines',
                  marker = list(color = 'red'),
                  name = 'Map frame')

    for (i in 1:nCls){

     p <- add_trace(p,
                     x = data[which(classes == i),1],
                     y = data[which(classes == i),2],
                     type = 'scatter',mode = 'markers',
                     marker = list(color = markColour[i], size = 7, symbol = 201),
                     name = paste('Class',toString(i)))%>%
        layout(title = "Projection of data points to the face",
               xaxis = list (title = "x"),
               yaxis = list (title = "y"))

    }

    return(p)

  }

  }

  }







