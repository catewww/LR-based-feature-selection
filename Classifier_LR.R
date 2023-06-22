#Function that performs model parameter values estimation
#estimates also prediction for observations classification purposes for binary classification problems
#estimates accuracy, log-likelihood and BIC values for the individual model

ClassifierLR <- function(LR, anno, epsilon = 0.1, alpha = 0.001, max_iter = 50)

{
  #define sigmoid function which will return 0-1 probabilities for each cell of being class 1
  sigmoid <- function(z)
  {
    return(1/(1+exp(-(z))))
  }
  
  #epochs -> number of repeated runs of the algorithm
  #alpha -> learning rate
  #main logistic regression algorithm to find coefficients basing on training dataset
  theta <- matrix(ncol = ncol(LR)+1, nrow = 1) #matrix where all parameters will be placed
  theta[1:ncol(theta)] = 0 #set initial theta values to 0
  prediction <- matrix(ncol = 1, nrow = nrow(LR)) #initialize matrix with all predictions for each cell within each epoch
  accuracy <- matrix() #matrix with accuracy estimates for each epoch
  annotation <- prediction #matrix with cells annotations basing on probabilities in prediction matrix
  L <- matrix()
  stop = 0 #stop iterations if difference between likelihoods is small enough
  i = 1
  #repeat all algorithm till difference between likelihoods will be less than epsilon or till number of iterations reach 10
  while (stop != 1)
  {
    #repeat for each cell
    for (j in 1 : nrow(prediction))
    {
      #calculate probability
      z = theta[1] 
      for (p in 1 : (ncol(theta)-1))
        z = z + (theta[p+1] * LR[j, p])
      
      #estimate prediction
      prediction[j, 1] = sigmoid(z)
      
      #calculate new theta parameters estimates for each feature - cross entropy as cost function
      for (t in 1 : ncol(theta))
      {
        if (t == 1)
          theta[t] = theta[t] - alpha * sum((prediction[j, 1] - anno[j]) * 1)
        if (t != 1)
          theta[t] = theta[t] - alpha * sum((prediction[j, 1] - anno[j]) * LR[j, t-1])
      }
    }
    
    #accuracy
    for (a in 1 : nrow(annotation))
    {
      #change probability to annotation 0 or 1
      if (prediction[a, 1] < 0.5)
        annotation[a, 1] = 0
      if (prediction[a, 1] >= 0.5)
        annotation[a, 1] = 1
    }

    #number of correct predictions
    N1_correct = as.numeric(length(which(annotation[which(annotation[,1] == '0'), 1] == anno[which(annotation[,1] == '0')]) == T))
    N2_correct = as.numeric(length(which(annotation[which(annotation[,1] == '1'), 1] == anno[which(annotation[,1] == '1')]) == T))
    N1 = length(which(annotation[,1] == '0'))
    N2 = length(which(annotation[,1] == '1'))
    acc1 = N1_correct/N1
    acc2 = N2_correct/N2
    if (is.na(acc1))
      acc1 = 0
    if (is.na(acc2))
      acc2 = 0
    
    accuracy[i] = ((acc1 + acc2) / 2) * 100
    
  
    #calculate likelihood difference
    L[i] = sum(prediction[which(anno == '1'), 1]) + sum(1 - prediction[which(anno == '0'), 1])
    if (i >= 2)
    {
      dif = abs(L[i] - L[i-1])
      if (dif < epsilon)#stop if difference is small enough
        stop = 1
    }
    i = i + 1
    if (i > max_iter)
    {
      stop = 1
      #cat(sprintf("%s did not converge\n", colnames(LR)[ncol(LR)]))
    }
    
  }
  
  #plot accuracy per epoch for training set
  #tiff(paste(colnames(LR)[ncol(LR)], "accuracy.tiff", sep = "_"), units="in", width=10, height=10, res=300)
  #plot(1:(i-1), accuracy, 'l', ylim = c(50, 100), main = paste('Accuracy per epoch', colnames(LR)), xlab = 'Epoch', ylab = 'Accuracy [%]', 
  #     col = 'red')
  #dev.off()
  
  #save accuracy values
  #write.table(accuracy, paste(colnames(LR)[ncol(LR)], "accuracyPERepoch.txt"), col.names = F, row.names = F)
  #save parameters
  #write.table(theta, paste(colnames(LR)[ncol(LR)], "parameters.txt"), col.names = F, row.names = F)
  #save likelihood values for each learning epoch (till dif < epsilon reached)
  #write.table(L, paste(colnames(LR)[ncol(LR)], "LikelihoodPERepoch.txt"), col.names = F, row.names = F)
  
  #Log-Likelihood for each model
  ind = which(prediction[,1] == 1)
  if (is.null(ind) == F)
    prediction[ind] = 0.99999
  LL = sum(log(prediction[which(anno == '1'), 1])) + sum(log(1 - prediction[which(anno == '0'), 1]))
  #calculate BIC for each model
  BIC = (ncol(LR)+1)*log(nrow(prediction)) - 2*LL
  
  result <- list("BIC" = BIC, "Likelihood" = LL, "Accuracy" = accuracy, "Parameters" = theta)
  return(result)


  
}
