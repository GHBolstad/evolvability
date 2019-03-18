#' Simulating evolutionary rate model
#' 
#' \code{simulate_rate} Simulating evolutionary rate model
#' 
#' \code{simulate_rate} Simulates different evolutionary rates models where a trait x evolves according to a Brownian motion or geometric Bronian motion 
#' process (constant evolutionary rate) while the rate of a second trait, y, is a linear function of x.  
#' 
#' @param tree a \code{\link{phylo}} object. Must be ultrametric and scaled to unit length.
#' @param starv_x starting value for x
#' @param sigma_x evolutionary rate of x  
#' @param a parameter of the evolutionary rate of y (see details)
#' @param b parameter of the evolutionary rate of y (see details)
#' @param x optional fixed values of x (only for the "residual_rate" model), 
#' must equal number of tips in the pylogeny, must be the same order as the tip labels
#' @param model either Brownaian motion model of x "predictor_BM", 
#' geometric Brownian motion model of x "predictor_geometricBM", or "residual_rate".
#' 
#' 
#' @return \code{simulate_rate} returns a dataset of x and y values at the tips.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export

simulate_rate <- function(tree, startv_x = NULL, sigma_x = NULL, a, b, x = NULL, model = "predictor_BM"){
  
  # include error messages so that ERROR if
  # Trees with negative edge lengths
  # start_x is negative for model == "predictor_geomtricBM"
  
  if(model == "predictor_BM" | model == "predictor_geometricBM"){
    EDGE <- cbind(tree$edge, round(tree$edge.length, 3))
    EDGE[,3][EDGE[,3]==0] <- 0.001 #adding one timestep to the edges that was rounded down to zero
    EDGE <- EDGE[order(tree$edge[,1]),]
    x_evo <- list()
    y_evo <- list()
    
    # Simulating x-values
    for(i in 1:nrow(EDGE)){
      if(EDGE[i,1] %in% EDGE[1:i,2]) x0 <- rev(x_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1]  # If the ancestor already has a starting value
      else{
        if(model == "predictor_BM") x0 <- startv_x
        if(model == "predictor_geometricBM") x0 <- log(startv_x)
      }
      x_evo[[i]] <- x0 + cumsum(rnorm(n = EDGE[i,3]*1000, mean = 0, sd = sqrt(1/1000)*sigma_x))
    }
  
    # percentage negative square roots
    if(model == "predictor_BM"){
      percent_negative <- length(which((a + b*unlist(x_evo))<0))/length(unlist(x_evo))*100
      if(percent_negative>0) warning(paste("Number of negative a + bx is ", round(percent_negative, 0), "%. The term a + bx is set to zero for these iterations", sep = ""))
    }
      
    # Simulating y-values
    for(i in 1:nrow(EDGE)){
      if(EDGE[i,1] %in% EDGE[1:i,2]) y0 <- rev(y_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1] # If the ancestor already has a starting value
      else y0 <- 0
      if(model == "predictor_BM"){
        r <- a + b*x_evo[[i]]
        r[r<0] <- 0
        y_evo[[i]] <- y0 + cumsum(sqrt(r) * rnorm(n = EDGE[i,3]*1000, mean = 0, sd = 1/sqrt(1000))) 
      }   
      if(model == "predictor_geometricBM")  y_evo[[i]] <- y0 + cumsum(sqrt(a + b*exp(x_evo[[i]])) * rnorm(n = EDGE[i,3]*1000, mean = 0, sd = 1/sqrt(1000))) 
    }
    
    EDGE <- cbind(EDGE, sapply(x_evo, function(x) rev(x)[1]), sapply(y_evo, function(x) rev(x)[1]))
    DATA <- EDGE[EDGE[,2]<EDGE[1,1],c(2,4,5)]
    DATA <- DATA[order(DATA[,1]),]
    DATA <- as.data.frame(DATA)
    colnames(DATA) <- c("species", "x", "y")
    DATA$species <- tree$tip.label[DATA$species]
    if(model == "predictor_geometricBM") DATA$x <- exp(DATA$x)
  }
  
  if(model == "residual_rate"){
    n <- length(tree$tip.label)
    A <-ape::vcv(tree)
    if(is.null(x)) x <- c(startv_x + t(chol(A))%*%rnorm(n, 0, sigma_x))
    V <- A + a*diag(n) + diag(a + b*x)
    if(length(diag(V)[diag(V)<0])>0) stop("Negative values in the variance matrix of y. Change a, b or x.")
    y <- t(chol(V))%*%rnorm(n)
    DATA <- data.frame(species = tree$tip.label, x = x, y = y)
  }
  
  return(DATA)
}

