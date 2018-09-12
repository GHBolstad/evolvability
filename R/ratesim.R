#' Simulating evolutionary rate model
#' 
#' \code{ratesim} Simulating evolutionary rate model
#' 
#' \code{ratesim} Simulates different evolutionary rates models where a trait x evolves according to a Brownian motion or geometric Bronian motion 
#' process (constant evolutionary rate) while the rate of a second trait, y, is a linear function of x.  
#' 
#' @param tree a \code{\link{phylo}} object. Must be ultrametric and scaled to unit length.
#' @param starv_x star value for x
#' @param sigma_x evolutionary rate of x  
#' @param a parameter of the evolutionary rate of y (see details)
#' @param b parameter of the evolutionary rate of y (see details)
#' @param model either Brownaian motion model of x "predictor_BM", or geometric Brownian motion model of x "predictor_geometricBM"
#' 
#' 
#' @return \code{ratesim} returns a dataset of x and y values at the tips.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export


ratesim <- function(tree, startv_x, sigma_x, a, b, model = "predictor_BM"){
  
  EDGE <- cbind(tree$edge, round(tree$edge.length, 3))
  EDGE[,3][EDGE[,3]==0] <- 0.001 #adding one timestep to the edges that was rounded down to zero
  EDGE <- EDGE[order(tree$edge[,1]),]
  x_evo <- list()
  y_evo <- list()
  
  # Simulating x-values
  for(i in 1:nrow(EDGE)){
    # If the ancestor is among the descendants, the startingvalue is updated
    if(EDGE[i,1] %in% EDGE[1:i,2]) x0 <- rev(x_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1]
    else x0 <- startv_x
    x_evo[[i]] <- x0 + cumsum(rnorm(n = EDGE[i,3]*1000, mean = 0, sd = sqrt(1/1000)*sigma_x))
  }

  if(model == "predictor_BM"){ # finding a k value to avoid negative roots
    k <- abs(b*min(sapply(x_evo, min)))+1
  } 
  
  # Simulating y-values
  for(i in 1:nrow(EDGE)){
    # If the ancestor is among the descendants, the startingvalue is updated
    if(EDGE[i,1] %in% EDGE[1:i,2]) y0 <- rev(y_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1]
    else y0 <- 0
    #if(model == "predictor_BM")   y_evo[[i]] <- y0 + cumsum(sqrt(a + b*x_evo[[i]]) * rnorm(n = EDGE[i,3]*1000, mean = 0, sd = 1/sqrt(1000))) 
    if(model == "predictor_BM")  y_evo[[i]] <- y0 + cumsum(sqrt(k*a + b*x_evo[[i]]) * rnorm(n = EDGE[i,3]*1000, mean = 0, sd = sqrt(1/1000*abs(a + b*x_evo[[i]])/(k*a + b*x_evo[[i]]))))
    # !!!!! must double check that the algorithm for this trick actually works (that the correct x value is used) !!!!!
    if(model == "predictor_geometricBM") y_evo[[i]] <- y0 + cumsum(sqrt(a + b*exp(x_evo[[i]])) * rnorm(n = EDGE[i,3]*1000, mean = 0, sd = 1/sqrt(1000))) 
  }
  
  EDGE <- cbind(EDGE, sapply(x_evo, function(x) rev(x)[1]), sapply(y_evo, function(x) rev(x)[1]))
  DATA <- EDGE[EDGE[,2]<EDGE[1,1],c(2,4,5)]
  DATA <- DATA[order(DATA[,1]),]
  DATA <- as.data.frame(DATA)
  colnames(DATA) <- c("species", "x", "y")
  if(model == "predictor_geometricBM") DATA$x <- exp(DATA$x)
  return(DATA)
}


