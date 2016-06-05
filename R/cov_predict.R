single_predict <- function(rit,instance){
  class_names <- names(rit[["Class_priors"]])
  nb_class <- length(class_names)
  
  probs <- vector(mode="numeric",length=nb_class) + log10(rit[["Class_priors"]])
  names(probs) <- class_names
  
  model <- rit[["Model"]]
  isCat <- rit[["cat_attr"]]
  
  for(elem in model){
    if(has_inter(elem[["interaction"]],instance,isCat)){
      probs <- probs + log10(elem[["sm"]])
    }
  }
  
  # Compute argmax
  response <- names(probs)[which.max(probs)]
  response
}

#' @title Classification Rule for Coverage-based Random Intersection Trees.
#' @description Applies a basic \code{argmax} rule in order to classify new instances.
#'
#' @return A response vector for the \code{testset} instances
#'
#' @param rit A model produced by \code{cov_RIT}
#' @param testset A dataframe containing the instances to classify
#' 
#' @references Ballarini Robin. Random intersection trees for genomic data analysis. Ecole polytechnique de Louvain, Université catholique de Louvain, 2016. Prom. : Dupont, Pierre.
#' @export
#'
cov_predict <- function(rit,testset){    
  # Predict
  response <- vector(mode="character",length=nrow(testset))
  
  for(i in 1:nrow(testset)){
    response[[i]] <- single_predict(rit,testset[i,])
  }
  response
}