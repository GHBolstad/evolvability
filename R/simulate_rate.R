#' Simulating evolutionary rate model
#' 
#' \code{simulate_rate} Simulating evolutionary rate model
#' 
#' \code{simulate_rate} Simulates different evolutionary rates models where a trait x evolves according to a Brownian motion or geometric Bronian motion 
#' process (constant evolutionary rate) while the rate of a second trait, y, is a linear function of x.  
#' 
#' @param tree a \code{\link{phylo}} object. Must be ultrametric and scaled to unit length.
#' @param starv_x starting value for x.
#' @param sigma_x evolutionary rate of x.
#' @param a parameter of the evolutionary rate of y (see details).
#' @param b parameter of the evolutionary rate of y (see details).
#' @param sigma_y evolutionary rate for the macroevolution of y (Brownian motion process) used in the "recent_evolution" model.
#' @param x optional fixed values of x (only for the "recent_evol" model), 
#' must equal number of tips in the pylogeny, must correspond to the order of the tip labels
#' @param model either Brownaian motion model of x "predictor_BM", 
#' geometric Brownian motion model of x "predictor_gBM", or "recent_evol".
#' 
#' @return \code{simulate_rate} returns a dataset of x and y values for the tips.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export

simulate_rate <- function(tree, startv_x = NULL, sigma_x = NULL, a, b, sigma_y = NULL, x = NULL, model = "predictor_BM"){
  
  # include error messages so that ERROR if
  # Trees with negative edge lengths
  # start_x is negative for model == "predictor_gBM"
  
  #### predictor_BM or predictor_gBM ####
  if(model == "predictor_BM" | model == "predictor_gBM"){
    EDGE <- cbind(tree$edge, round(tree$edge.length, 3))
    EDGE[,3][EDGE[,3]==0] <- 0.001 #adding one timestep to the edges that was rounded down to zero
    EDGE <- EDGE[order(tree$edge[,1]),]
    x_evo <- list()
    y_evo <- list()
    
    # Simulating x-values
    for(i in 1:nrow(EDGE)){
      if(EDGE[i,1] %in% EDGE[1:i,2]){
        x0 <- rev(x_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1]  # If the ancestor already has a starting value
        }else{
          if(model == "predictor_BM")  x0 <- startv_x
          if(model == "predictor_gBM") x0 <- log(startv_x)
        }
      x_evo[[i]] <- x0 + cumsum(rnorm(n = EDGE[i,3]*1000, mean = 0, sd = sqrt(1/1000)*sigma_x))
    }

    if(a == "scaleme"){
      if(model == "predictor_BM")  a <- -min(b*unlist(x_evo))
      if(model == "predictor_gBM") a <- -min(b*exp(unlist(x_evo)))
    } 
    
    # percentage of negative square roots
    if(model == "predictor_BM")  percent_negative <- length(which((a + b*unlist(x_evo))<0))/length(unlist(x_evo))*100
    if(model == "predictor_gBM") percent_negative <- length(which((a + b*exp(unlist(x_evo)))<0))/length(exp(unlist(x_evo)))*100
    if(percent_negative>0.01) warning(paste("Number of negative a + bx is ", round(percent_negative, 1), "%. The term a + bx is set to zero for these values of x", sep = ""))
    
    # Simulating y-values
    for(i in 1:nrow(EDGE)){
      if(EDGE[i,1] %in% EDGE[1:i,2]){
        y0 <- rev(y_evo[[which(EDGE[1:i,2] == EDGE[i,1])]])[1] # If the ancestor already has a starting value
        }else{
          y0 <- 0
          }
      if(model == "predictor_BM")  r <- a + b*x_evo[[i]]
      if(model == "predictor_gBM") r <- a + b*exp(x_evo[[i]])
      r[r<0] <- 0
      y_evo[[i]] <- y0 + cumsum(sqrt(r) * rnorm(n = length(r), mean = 0, sd = 1/sqrt(1000))) 
    }
    

    EDGE <- cbind(EDGE, sapply(x_evo, function(x) rev(x)[1]), sapply(y_evo, function(x) rev(x)[1]))
    DATA <- EDGE[EDGE[,2]<EDGE[1,1],c(2,4,5)]
    DATA <- DATA[order(DATA[,1]),]
    DATA <- as.data.frame(DATA)
    colnames(DATA) <- c("species", "x", "y")
    DATA$species <- tree$tip.label[DATA$species]
    if(model == "predictor_gBM") DATA$x <- exp(DATA$x)
  }
  
  #### recent_evol ####
  if(model == "recent_evol"){
    n <- length(tree$tip.label)
    A <- ape::vcv(tree)
    sigma_y <- sigma_y[1]
    if(is.null(x)) x <- rnorm(n, startv_x, sigma_x)
    if(is.null(sigma_y)) stop("Specify a positive sigma_y") 
    if(sigma_y <= 0)    stop("Specify a positive sigma_y") 
    Vmacro <- c(sigma_y^2)*A
    diag_Vmicro <- c(a + b*x)
    diag_Vmicro[diag_Vmicro<0] <- 0  #Ensures that there are no negative variances
    Vmicro <- diag(x = diag_Vmicro)
    V <- Vmacro + Vmicro
    #if(length(diag(V)[diag(V)<0])>0) stop("Negative values along the diagonal of the variance matrix of y. Change a, b or x.")
    y <- t(chol(V))%*%rnorm(n)
    DATA <- data.frame(species = tree$tip.label, x = x, y = y)
  }
  return(DATA)
}



