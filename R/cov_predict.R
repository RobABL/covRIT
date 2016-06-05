single_predict <- function(rit,instance,count){
  class_names <- names(rit[["Class_priors"]])
  nb_class <- length(class_names)
  
  probs <- vector(mode="numeric",length=nb_class) + log10(rit[["Class_priors"]])
  names(probs) <- class_names
  
  model <- rit[["Model"]]
  isCat <- rit[["cat_attr"]]
  
  count_has <- 0
  for(elem in model){
    if(has_inter(elem[["interaction"]],instance,isCat)){
      probs <- probs + log10(elem[["sm"]])
      count_has <- count_has + 1 
    }
  }
  
  # Compute argmax
  response <- names(probs)[which.max(probs)]
  if(count){
    return(list(response=response,count_has=count_has/length(model)))
  }
  else{
    return(response)
  }
}

cov_predict <- function(rit,testset,count=FALSE){    
  # Predict
  if(count)
    response <- vector(mode="list",length=nrow(testset))
  else
    response <- vector(mode="character",length=nrow(testset))
  
  for(i in 1:nrow(testset)){
    response[[i]] <- single_predict(rit,testset[i,],count)
  }
  
  if(count){
    r <- sapply(response,function(i){
      i[["response"]]
    })
    c <- sapply(response,function(i){
      i[["count_has"]]
    })
    return(list(response=r,count_has=sum(c)))
  }
  else{
    return(response)
  }
}