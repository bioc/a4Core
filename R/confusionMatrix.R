#' Generic function to produce a confusion matrix (related to a classification problem)
#' @param x object (usually a model fit object) that contains all information 
#' needed to produce the confusion matrix. 
#' @param ... further arguments for a specific method
#' @return A confusion matrix
#' @author Tobias Verbeke
#' @keywords models
#' @export
confusionMatrix <- function(x, ...){ # common to a4Classif (pamClass) and nlcv
  UseMethod("confusionMatrix")
}



