#' Riemman Sums of Location Vectors
#'
#' @param x \eqn{\tau \times 1} local vector of coefficient estimates from FDR in stage 1
#'
#' @return Riemman sum of the given location
#' @export

riemSum = function(x) {
  # Uses the Simpson's Rule (numerical method).
  # https://link.springer.com/chapter/10.1007/978-1-4612-4974-0_26
  nn = length(x)
  if(length(x) %% 2 == 0){
    xlast = x[length(x)]
    x = x[-length(x)]
  }

  n = length(x)-1

  xValues = seq(from=1, to=n, by=2)

  sum = list()
  for(i in 1:length(xValues)){
    n_sub = xValues[[i]]-1
    n = xValues[[i]]
    n_add = xValues[[i]]+1

    v1 = x[[n_sub+1]]
    v2 = x[[n+1]]
    v3 = x[[n_add+1]]

    s = (1/3)*(v1+4*v2+v3)
    sum = append(sum, s)
  }
  sum = unlist(sum)

  auc = sum(sum)
  ##jh edit: add trapezoid for last two data points to result
  if(nn %% 2 == 0){
    auc = auc + (x[length(x)] + xlast)/2
  }


  return(auc)

}
