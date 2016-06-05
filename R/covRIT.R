stopParameters <- function(data,classes,epsilon,epsilon_cat,n_trees,depth,min_inter_sz,branch,split_nb){
  # check parameters
  if(!is.data.frame(data))
    stop("Data parameter must be a dataframe.")
  if(length(classes) != nrow(data))
    stop("Data and response vector must have the same number of rows.")
  if(nrow(data) == 0)
    stop("Data cannot be empty.")
  if(any(is.na(classes)))
    stop("Response vector cannot contain missing values.")
  if(epsilon < 0 || epsilon > 1)
    stop("epsilon hp should be between 0 and 1.")
  if(epsilon_cat < 0 || epsilon_cat > 1)
    stop("epsilon_cat hp should be between 0 and 1.")
  if(depth < 0)
    stop("depth should be positive.")
  if(n_trees < 1)
    stop("n_trees should be at least 1.")
  if(min_inter_sz < 2)
    stop("min_inter_sz should be at least 2.")
  if(branch < 1)
    stop("branch should be at least 1.")
  if(split_nb < 1)
    stop("split_nb should be at least 1.")
}

cov_RIT <- function(data,classes,theta,epsilon,epsilon_cat,n_trees=100L,depth=10L,min_inter_sz=2L,branch=5,split_nb=1,es=TRUE){
  
  stopParameters(data,classes,epsilon,epsilon_cat,n_trees,depth,min_inter_sz,branch,split_nb)
  
  # Class information
  classes <- factor(classes)
  class_names <- levels(classes)
  nb_class <- nlevels(classes)
  
  # Check parameter theta
  if(missing(theta)){
    theta <- vector(mode="numeric",length=nb_class) - 1
    names(theta) <- class_names
  }
  else{
    if(length(theta) == 1){ # theta is scalar
      theta <- vector(mode="numeric",length=nb_class) + theta
      names(theta) <- class_names
    }
    else if(!(all(class_names %in% names(theta)))){
      stop("Named vector theta doesn't contain prevalence thresholds for some classes.")
    }
  }
  

  d <- sapply(data,function(attr){
    if(is.ordered(attr))
      as.numeric(attr)
    else
      attr
  })
  data <- as.data.frame(d)
  isCat <- sapply(data,is.factor) # Logical vector: if each feature/attr is categorical or not
  
  # Split dataset per class and order theta vector
  datas <- vector("list",length=nb_class)
  names(datas) <- class_names
  o_theta <- vector(mode="numeric",length=nb_class)
  names(o_theta) <- class_names
  for(cls in class_names){
    datas[[cls]] <- data[which(classes==cls),]
    o_theta[[cls]] <- theta[[cls]]
  }
  
  # Compute spans
  spans <- sapply(data,function(attr){
    if(is.factor(attr))
      nlevels(attr)
    else
      max(attr) - min(attr)
  })
  
  # Call main algo
  all_leaves <- cpp_cov_RIT(datas, o_theta, isCat, epsilon, epsilon_cat, spans, n_trees, depth, branch, split_nb,min_inter_sz,es)
  
  # Compute class prior probabilities
  priors <- data.table(classes)[,.N,keyby=classes]
  priors[,"N"] <- priors[,N] / length(classes)
  priors_vec <- vector(mode="numeric",length=nb_class)
  names(priors_vec) <- class_names
  for(i in 1:nrow(priors)){
    priors_vec[[priors[i,classes]]] <- priors[i,N]
  }
  
  # Smoothing
  cls_sums <- sapply(datas,nrow)
  sm_leaves <- lapply(all_leaves,function(inter){
    nb <- inter[["prevalence"]]*cls_sums
    nb_neg <- {1-inter[["prevalence"]]}*cls_sums
    list(interaction=inter[["interaction"]],prevalence=inter[["prevalence"]],sm={nb+1}/{nb+nb_neg+2},sm_neg={nb_neg+1}/{nb+nb_neg+2})
  })
  
  list(Model=sm_leaves,Class_priors=priors_vec,cat_attr=isCat)
}
