calLK <- function(fam.cancer.data, penetrance.all){
  invasive.cut <- 100
  # Aim: calculate the likelihood, Pr(D|G) for each individual
  # Author: Gang Peng
  #
  # input: 
  # fam.cancer.data: family data
  # penetrance.all: penetrance in the population for male and female
  #
  # output:
  # likelihood
  
  num.individual <- length(fam.cancer.data$id)
  lik <- matrix(NA, nrow=num.individual, ncol=3)
  
  # penetrance length, oldest age in the penetrance table
  length.pene <- nrow(penetrance.all$fMX)
  
  for(i in 1:num.individual){
    if(fam.cancer.data$num.cancer[i]==0){
      if(is.na(fam.cancer.data$gender[i])){
        pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
        lik[i,] <- lkNoneAffect(pene.mix, fam.cancer.data$age[i])
      }
      else{        
        if(fam.cancer.data$gender[i]==1){
          lik[i,] <- lkNoneAffect(penetrance.all$fMX, fam.cancer.data$age[i])
        }
        else{
          lik[i,] <- lkNoneAffect(penetrance.all$fFX, fam.cancer.data$age[i])
        }
      }
    }
    else{
      flag.benign <- TRUE
      for(j in 1:fam.cancer.data$num.cancer[i]){
        if(fam.cancer.data$cancer.info[[i]]$cancer.type[j] < invasive.cut){
          flag.benign <- FALSE
          break
        }
      }
      
      if(flag.benign){
        if(is.na(fam.cancer.data$gender[i])){
          pene.mix <- (penetrance.all$fMX + penetrance.all$fFX)/2
          lik[i,] <- lkNoneAffect(pene.mix, fam.cancer.data$age[i])
        }
        else{        
          if(fam.cancer.data$gender[i]==1){
            lik[i,] <- lkNoneAffect(penetrance.all$fMX, fam.cancer.data$age[i])
          }
          else{
            lik[i,] <- lkNoneAffect(penetrance.all$fFX, fam.cancer.data$age[i])
          }
        }
      }
      else
      {
        age.tmp <- min(fam.cancer.data$cancer.info[[i]]$diag.age)
        if(age.tmp < 1){
          age.tmp = 1
        }
        if(age.tmp > length.pene){
          age.tmp = length.pene
        }
        if(is.na(fam.cancer.data$gender[i])){
          lik[i,] <- (penetrance.all$fMX[age.tmp,]+penetrance.all$fFX[age.tmp,])/2
        }
        else{
          if(fam.cancer.data$gender[i]==1){
            lik[i,] <- penetrance.all$fMX[age.tmp,]
          }
          else{
            lik[i,] <- penetrance.all$fFX[age.tmp,]
          }
        }
      }
      
      
    }
  }
  
  return(lik)
}


lkNoneAffect <- function(penetrance, age){
  # calculate the likelihood for unaffected sample
  # Args:
  #   age: sample's age when he/she still doesn't have disease
  #
  # Returns:
  #   likelihood
  
  length.pene <- nrow(penetrance)
  if(age==0)
  {
    age <- 1
  }
  if(age > length.pene)
  {
    age <- length.pene
  }
  rlt <- rep(1,3)
  rlt[1] <- 1-sum(penetrance[1:age,1])
  rlt[2] <- 1-sum(penetrance[1:age,2])
  rlt[3] <- 1-sum(penetrance[1:age,3])
  #   for(i in 1:age){
  #     rlt[1] <- rlt[1]*(1-penetrance[i,1])
  #     rlt[2] <- rlt[2]*(1-penetrance[i,2])
  #     rlt[3] <- rlt[3]*(1-penetrance[i,3])
  #   }
  rlt
}