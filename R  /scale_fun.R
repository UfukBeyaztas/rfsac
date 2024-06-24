# scale function
scale_fun <- function(data, method.type){
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  if(method.type == "robust"){
    if(p == 1){
      cen.data <- apply(data, 2, median)
    }else{
      cen.data <- l1median(data)
    }
  }else if(method.type == "classical"){
    cen.data <- apply(data, 2, mean)
  }
  
  scl.data <- (data - matrix(cen.data, nrow = n, ncol = p, byrow = TRUE))
  
  return(scl.data)
}
