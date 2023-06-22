#Function that performs feature selection utilizing wrapper method with build in logistic regression with log-likelihood estimation procedures
#BIC - choose better model form N-feature models
#BF - choose better model comparing N-1 with N-feature models
#stop if max number of features reached or if N-1-feature model is better
#create a loop over k-algorithm-runs and obtain k individual models created with 
#accuracy, BIC, BF, feature name and model parameter values saved as .txt files
#works with ClassifierLR function utilized inside


#Full Logistic regression classifier
FullLRclassifier <- function(LR, anno, epsilon1 = 0.1, max_iter1 = 100, alpha1 = 0.001, draw, ngenes=ncol(LR))
{

  Mk = ncol(LR) #number of models that have to be made in each step of choosing best model with attachment method
  BIC_compare <- matrix() #keep all BIC values for each step (BIC values are means from oll cells at particular step)
  BIC_min <- matrix() #save BIC values for best models on each step of models building
  step = 1 #number of genes in particular model
  stop = 0 #criterion to stop adding more features to model
  LL <- matrix() #matrix to keep Likelihood values
  LL_best <- matrix() #matrix with likelihoods of the best models after each step
  BF <- matrix()
  
  while (stop != 1)
  {
    for (M in 1 : (Mk-step+1)) #at each level one feature less Mk-step+1
    {
      #first step - single models
      if (step == 1)
      {
        LR_M = as.data.frame(LR[, M]) #pass only one gene expression to parameters estimate
        colnames(LR_M) = colnames(LR)[M]
        #estimate parameters and calculate BIC for each model (each cell has its own model with only one gene at the very beginning)
        BIC = ClassifierLR(LR = LR_M, anno, epsilon = epsilon1, max_iter = max_iter1, alpha = alpha1)
        
        BIC_compare[M] = BIC$BIC
        #choose the smallest BIC from first step
        BIC_compare <- as.matrix(BIC_compare)
        BIC_min = as.numeric(which(BIC_compare == min(BIC_compare)))
        #find the name of gene with the best BIC
        gene_best = colnames(LR)[BIC_min]
        #keep log-likelihood
        LL[M] = BIC$Likelihood
      }
      
      #next steps - more features in the model
      else
      {
        #next step - multiple features models - on each next step add one gene more to the best model till BIC is decreasing
        if (M == 1)
        {
          BIC_compare <- matrix()
          #save best gene in LR_M matrix to add more genes in further steps to this best from previous step
          LR = LR[, -which(colnames(LR) == gene_best)]
          LR_M = cbind(LR_M, LR[, M])
          colnames(LR_M)[ncol(LR_M)] = colnames(LR)[M]
        }
        else
        {
          LR_M[,step] = LR[, M]
          colnames(LR_M)[step] = colnames(LR)[M]
        }
        
        #estimate parameters and calculate BIC for each model (each cell has its own model with only one gene at the very beginning)
        BIC = ClassifierLR(LR = LR_M, anno, epsilon = epsilon1, max_iter = max_iter1, alpha = alpha1)
        BIC_compare[M] = BIC$BIC
        #choose the smallest BIC
        BIC_min = as.numeric(which(BIC_compare == min(BIC_compare)))
        #if there are two models with 
        #find the name of gene with the best BIC
        gene_best = colnames(LR)[BIC_min]
        #keep likelihood
        #calculate mean likelihood among all cells
        LL[M] = BIC$Likelihood
      }
      cat(sprintf("| "))
      
    }
    if (step == 1)
    {
      LR_M = as.data.frame(LR[, BIC_min])
      colnames(LR_M) = colnames(LR)[BIC_min]
    }
    if (step >= 2)
    {
      LR_M = as.data.frame(cbind(as.data.frame(LR_M[, 1:(step-1)]), LR[, BIC_min]))
      colnames(LR_M)[step] = colnames(LR)[BIC_min]
    }
    
    cat(sprintf("step %s\n", step))
    
    #keep best model likelihood
    LL_best[step] = LL[BIC_min]
    cat(sprintf("Log-Likelihood %s\n", LL_best[step]))
    cat(sprintf("BIC %s\n", BIC_compare[BIC_min]))
    cat(sprintf("Gene %s\n", gene_best))
    
    
    #save gene name that is included in the final model
    write.table(gene_best, paste(draw, 'MODEL_gene_step', step, '.txt', sep = ""), quote = F)
    #write.table(BIC_compare[BIC_min], paste('MODEL_BIC_step', step, '.txt', sep = ""), quote = F)
    #write.table(LL_best[step], paste('MODEL_LL_step', step, '.txt', sep = ""), quote = F)
    
    #save Likelihood of the best model on particular step in LL matrix returned by ClassifierLR function
    #compare Bayes factors on different steps of the model (level of model's complication - number of parameters/features)
    if (step >= 2)
    {
      #calculate Bayes Factor (index of Bayes Factor is one less than step because BF is calculated for step 2 and more)
      BF[step-1] = exp(LL_best[step] - LL_best[step-1])
      cat(sprintf("Bayes Factor %s\n", BF[step-1]))
      #write.table(BF[step-1], paste('MODEL_BF_step', step, '.txt', sep = ""), quote = F)
      
      #if next level model is worse basing on BF than previous one stop building more complicated models
      if (BF[step-1] < (10^1.5))
      {
        stop = 1
        cat(sprintf("STOPPED at step %s because of Bayes Factor=%s choosing OLD model at this level\n", step, BF[step-1]))
      }
    }
    
    step = step + 1#next gene to be added to model building
    
    #stop if no more genes are to be added in further steps
    if (step == length(ngenes))
    {
      stop = 1
      cat(sprintf("STOPPED at step %s because of no more genes are to be added to the model\nFinal model build of all provided features", step))
    }
    
    #save results from the last best model
    if (stop != 1)
    {
      write.table(BIC_compare[BIC_min], paste(draw, '_MODEL_BIC.txt', sep = ""), quote = F)
      write.table(LL_best[step-1], paste(draw, '_MODEL_LL.txt', sep = ""), quote = F)
      write.table(BIC$Accuracy[length(BIC$Accuracy)], paste(draw, '_MODEL_Acc.txt', sep = ""), quote = F)
      write.table(BIC$Parameters, paste(draw, '_MODEL_Parameters.txt', sep = ""), quote = F)
    }
  }
  
}