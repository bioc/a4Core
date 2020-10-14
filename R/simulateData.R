# simulate data for package testing and demonstration purposes
# nCols: number of samples
# nRows: number of features (genes)
# nEffectRows: number of differentially expressed features
# nNoEffectCols: number of samples for which the profile 
# of a differentially expressed feature will be
# set similar to the other class 


#' Simulate Data for Package Testing and Demonstration Purposes
#' 
#' Simulate Data for Package Testing and Demonstration Purposes
#' @param nCols number of samples; currently this should be an even number
#' @param nRows number of features (genes)
#' @param nEffectRows number of differentially expressed features
#' @param nNoEffectCols number of samples for which the profile
#'  of a differentially expressed feature will be set similar 
#'  to the other class
#' @param betweenClassDifference Average mean difference between the two classes 
#' to simulate a certain signal in the features for which an effect was introduced;
#' the default is set to 1
#' @param withinClassSd Within class standard deviation used to add a certain noise 
#' level to the features for which an effect was introduced; the default standard
#' deviation is set to 0.5
#' @return object of class ExpressionSet with the characteristics specified 
#' @note The simulation assumes the variances are equal between the two classes.
#' Heterogeneity could easily be introduced in the simulation if this would
#' be requested by the users.
#' @examples 
#' someEset <- simulateData(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5)
#' someEset
#' @usage simulateData(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5,
#'  betweenClassDifference = 1, withinClassSd = 0.5)
#' @author W. Talloen and T. Verbeke
#' @keywords manip
#' @importFrom stats rnorm
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @export
simulateData <- function(nCols = 40, nRows = 1000, nEffectRows = 5, nNoEffectCols = 5,
    betweenClassDifference = 1, withinClassSd = 0.5){
  ## response
  if (nCols %% 2 == 1) stop("'nCols' should be even")
  if (nNoEffectCols < 0) stop("'nNoEffectCols' should be positive (or zero)")
  # TODO: only works when even number of nCols !! TV (cf. infra as well)
  
  yData <- as.factor(rep(c("A", "B"), each = nCols / 2))
  
  ## no effect
  xData <- matrix(rnorm(nRows * nCols, mean = 8.73), ncol = nCols)
  colnames(xData) <- paste("Sample", seq_len(nCols), sep = "")
  rownames(xData) <- paste("Gene", seq_len(nRows), sep = ".") # change TV to prevent 
  # possible bug in new MLInterfaces version
  
  ## add effect to a certain number of rows
  xDataEffect <- matrix(0, nrow = nEffectRows, ncol = nCols)
  
  if (nNoEffectCols == 0){
    
    if (nEffectRows > 0) {
      
      ## create signal  
      for (i in seq(nEffectRows)){
        
        xDataEffect[i, ] <- rnorm(nCols, mean = 0, sd = withinClassSd) + 
            betweenClassDifference * as.numeric(yData)
      }
      xData[seq(nEffectRows), ] <- xDataEffect
      
    }  
  } else {
    ## create misbehaving samples
    
    yDataAnti <- factor(levels(yData)[(as.numeric(yData) %% 2) + 1], levels = levels(yData))
    yDataNoEffectCols <- yData
    yDataNoEffectCols[seq_len(nNoEffectCols)] <- yDataAnti[seq_len(nNoEffectCols)]
    
    if (nEffectRows > 0) {
      ## create signal  
      for (i in seq(nEffectRows)){
        
        xDataEffect[i, ] <- rnorm(nCols, mean = 0, sd = withinClassSd) + 
            betweenClassDifference * as.numeric(yDataNoEffectCols)
      }
      xData[seq(nEffectRows), ] <- xDataEffect
    }
  }  
  
  ### turn data into an ExpressionSet object
  info <- data.frame(yData, row.names = colnames(xData))
  colnames(info) <- "type"
  labels <- data.frame(type = "type")
  pd <- AnnotatedDataFrame(data = info, varMetadata = labels)
  # new("ExpressionSet", exprs = as.matrix(xData), phenoData = pd)
  eset.all <- ExpressionSet(assayData = as.matrix(xData), phenoData = pd) 
  # "ExpressionSet")  # will generate a lot of warnings
  return(eset.all)
}




