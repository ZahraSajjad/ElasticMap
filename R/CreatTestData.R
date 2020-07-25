CreateTestData = function(dims, N, disper, angles){

  library(pracma)

  data = matrix( runif( N*dims ), N)

  disper = kronecker(matrix(1,nrow(data),1),t(disper))

  data = bsxfun("*",data,disper)

  K = length(angles)

  if (K > dims-1 ){ K = dims-1}

  for (i in 1:K){
  c = cosd(angles[i])
  s = sind(angles[i])
  R = matrix(c(c,s,-s,c),nrow = 2, byrow = TRUE)
  data[,i:(i+1)] = data[,i:(i+1)]%*%R
  }

  return(t(data))
}
