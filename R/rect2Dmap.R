rect2Dmap = function(rows,cols){

   map.sizes = c(rows,cols)  #Store size

   map.dimension = length(map.sizes)

   N = rows*cols

   a1 = kronecker(matrix(1,rows,1),t(1:cols))

   a2 = kronecker(matrix(1,1,cols),1:rows)

   map.internal = cbind(c(a1),c(a2))

   A = seq(1,N)

   dim(A) = c(rows,cols)

   B = A [1:(dim(A)[1]-1),]

   C = A [2:dim(A)[1],]

   D = A[,1:(dim(A)[2]-1)]

   E = A[,2:(dim(A)[2])]

   link1 = matrix(c(c(B),c(D)))

   link2 = matrix(c(c(C),c(E)))

   map.links = cbind(link1,link2)


   B1 = A [1:(dim(A)[1]-2),]

   B2 = A [2:(dim(A)[1]-1),]

   B3 = A [3:(dim(A)[1]),]

   C1 = A [,1:(dim(A)[2]-2)]

   C2 = A [,2:(dim(A)[2]-1)]

   C3 = A [,3:(dim(A)[2])]

   ribs1 = matrix(c(c(B1),c(C1)))

   ribs2 = matrix(c(c(B2),c(C2)))

   ribs3 = matrix(c(c(B3),c(C3)))

   map.ribs = cbind(ribs1,ribs2,ribs3)


   B1 = A [1:(dim(A)[1]-1),1:(dim(A)[2]-1)]

   B2 = A [2:(dim(A)[1]),1:(dim(A)[2]-1)]

   B3 = A [1:(dim(A)[1]-1),2:(dim(A)[2])]

   C1 = A [1:(dim(A)[1]-1),2:(dim(A)[2])]

   C2 = A [2:(dim(A)[1]),1:(dim(A)[2]-1)]

   C3 = A [2:(dim(A)[1]),2:(dim(A)[2])]

   face1 = matrix(c(c(B1),c(C1)))

   face2 = matrix(c(c(B2),c(C2)))

   face3 = matrix(c(c(B3),c(C3)))

   map.faces = cbind(face1,face2,face3)

   map = list(map.sizes = map.sizes,
              map.dimension = map.dimension,
              map.internal = map.internal,
              map.links = map.links,
             map.ribs = map.ribs,
             map.faces = map.faces,
             map.Pcs = NULL,
             map.mean = NULL,
             map.preproc = NULL,
             map.mapped = NULL,
             map.disp = NULL)

   return(map)

}
