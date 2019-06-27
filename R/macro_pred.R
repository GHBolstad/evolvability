#' Macroevolutionary predictions
#' 
#' \code{macro_pred} Macroevolutionary predictions
#' 
#' \code{macro_pred} Gives a vector of macroevolutionary predictions for each 
#' species based on the other species given the phylogeny and a phylogenetic
#' variance matrix.
#' 
#' @param y vector of species means
#' @param V phylogenetic variance matrix, must have same order as \code{y}
#' @param useLFO excludes the focal species when calculating the corresponding species' 
#' mean. The correct way is to use TRUE, but in practice it has little effect and FALSE 
#' will speed up the model fit.
#' 
#' @return \code{macro_pred} returns a of macroevolutionary preditions at the 
#' tips.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export

macro_pred <- function(y, V, useLFO = TRUE){
  Vinv <- solve(V)
  if(useLFO){# here the mean is calculated leaving out the focal species
    X <- matrix(rep(1, length(y)-1), ncol = 1) # design matrix
    y_mean <- as.matrix(
      apply(cbind(y), 1, 
            function(x){
              j <- which(c(y)==x)
              new_y <- y[-j]
              new_Vinv <- Vinv[-j,-j]
              solve(t(X)%*%new_Vinv%*%X)%*%t(X)%*%new_Vinv%*%new_y
            }
            )
      )
  } else {# not leaving out the focal species
    X <- matrix(rep(1, length(y)), ncol = 1)
    y_mean <- c(solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv%*%y)
  }
  y_predicted <- y_mean + (V - solve(diag(x = diag(Vinv))))%*%Vinv%*%(y-y_mean)
}

