#' Methods for topTable
#' 
#' Methods for topTable. topTable extracts the top n most important features
#' for a given classification or regression procedure 
#' @section Methods:
#' \describe{glmnet and lognet
#' \itemize{
#' \item{fit = "glmnet", n = "numeric"}{glmnet objects are produced by \code{lassoClass} (a4Classif) or 
#' \code{lassoReg} (a4Base)}
#' \item{fit = "lognet", n = "numeric"}{lognet objects are produced by \code{lassoClass} (a4Classif) or 
#'   \code{lassoReg} (a4Base)}
#' \item{fit = "elnet", n = "numeric"}{elnet objects are produced by \code{lassoClass} (a4Classif) or 
#'   \code{lassoReg} (a4Base)}
#' }
#' }
#' @name topTable-methods
#' @docType methods
#' @aliases topTable-methods
#' @aliases topTable,glmnet-method
#' @aliases topTable,lognet-method
#' @aliases topTable,elnet-method
#' @param fit object resulting from a classification or regression procedure
#' @param n number of features that one wants to extract from a table that
#' ranks all features according to their importance in the classification
#' or regression model; defaults to 10 for limma objects
#' @keywords methods manip
NULL

#' S4 Generic for obtaining a top table 
#' 
#' a top table is a rectangular object (e.g. data frame) which
#' lists the top n most relevant variables
#' @param fit object for which to obtain a top table, generally a fit object for a given model class
#' @param n  number of features (variables) to list in the top table, ranked by importance
#' @param ... further arguments for specific methods
#' @return Top table with top n relevant variable.
#' @usage topTable(fit, n, ...)
#' @author Tobias Verbeke
#' @exportMethod topTable
#' @importFrom methods setGeneric
setGeneric("topTable", function(fit, n, ...){ # common to nlcv and (at least) a4Classif
      standardGeneric("topTable")
    })

#' @import glmnet
#' @importFrom methods setOldClass
setOldClass("glmnet")

#' @import glmnet
#' @importFrom methods setOldClass
setOldClass("elnet")

#' @import glmnet
#' @importFrom methods setOldClass
setOldClass("lognet") # from glmnet 1.5.3 onwards

#' @importFrom stats coef
glmnetUtil <- function(fit, n){
  summary.output <- summary(fit)
  coef.output <- coef(fit) # extract coefficients at a single value of lambda
  last.coef.output <- coef.output[, ncol(coef.output)]
  selProbeSets <- last.coef.output[which(last.coef.output != 0)]
  selProbeSetsGeneSymbol <- fit$featureData[names(selProbeSets), "SYMBOL"]
  
  selGenesOutput <- cbind.data.frame(selProbeSetsGeneSymbol, selProbeSets)
  rownames(selGenesOutput) <- names(selProbeSets)
  colnames(selGenesOutput) <- c('Gene','Coefficient')
  # remove the estimate of the intercept (typically, but not always, the first row)
  exclIntercept <- which(rownames(selGenesOutput)%in%c('','(Intercept)'))
  selGenesOutput <- selGenesOutput[-exclIntercept,]
  
  numberSelGenes <- nrow(selGenesOutput)
  topList <- selGenesOutput[order(abs(selGenesOutput[,2]),decreasing=TRUE),][seq_len(min(n, numberSelGenes)),] # first row is the estimate of the intercept.
  retval <- list(topList = topList, numberSelGenes = numberSelGenes, n = n)
  return(retval)
}

#' @export
#' @importFrom methods setMethod
setMethod("topTable",
    "glmnet",
    function(fit, n){
      res <- glmnetUtil(fit = fit, n = n)
      class(res) <- "topTableGlmnet"
      return(res)
    }
)

#' @export
#' @importFrom methods setMethod
setMethod("topTable",
    "lognet",
    function(fit, n){
      res <- glmnetUtil(fit = fit, n = n)
      class(res) <- "topTableLognet"
      return(res)
    }
)

#' @export
#' @importFrom methods setMethod
setMethod("topTable",
    "elnet",
    function(fit, n){
      res <- glmnetUtil(fit = fit, n = n)
      class(res) <- "topTableElnet"
      return(res)
    }
)


#' @export
print.topTableGlmnet <- function(x,  ...){
  cat("The lasso selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}

#' @export
print.topTableLognet <- function(x,  ...){
  cat("The lasso selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}

#' @export
print.topTableElnet <- function(x,  ...){
  cat("The lasso selected ", x$numberSelGenes, " genes. The top ", x$n, " genes are:\n\n", sep = "")
  print(x$topList, ...)
}
