#Exploring SNPs score simulation
#After http://www.r-bloggers.com/instrumental-variables-simulation/

#Set working directory
setwd("Z:/Data/Fotios/Simulation")

#Load used libraries
library(MASS)
library(ggplot2)

mr_corr_sim<-function(n_sample,score_corr,met_corr){
  
  #Generate x1 (metabolite 1) and x2 (Metabolite 2)
  x1StarAndx2Star <- mvrnorm(n_sample, c(20, 15), matrix(c(1, met_corr, met_corr, 1), 2, 2))
  
  #Summary of generated metabolites
  summary(x1StarAndx2Star)
  
  #Define the two metabolites
  x1Star <- x1StarAndx2Star[, 1]
  x2Star <- x1StarAndx2Star[, 2]
  
  #Generate scores
  z1Andz2 <- mvrnorm(n_sample, c(1, 1), matrix(c(1, score_corr, score_corr, 1), 2, 2))
  
  #Generate genetic scores z1 and z2 
  z1<-z1Andz2[,1]
  z2<-z1Andz2[,2]
  
  #Generate x1 and x2
  x1 <- x1Star+z1
  x2 <- x2Star+z2
  
  #Tesing the correlations
  #cor(x1Star,x2Star)
  #cor(x1,x2)
  #cor(x1,z2)
  #cor(x2,z1)
  #cor(z1,z2)
  
  #Generate the risk factor
  y<-1+x1+x2+rnorm(n_sample, 0, 0.5)
  
  #Estimate full model
  full<-lm(y ~ x1 + x2)
  full_x1<-full$coefficients[2]
  
  #Estimate model without x2
  single<-lm(y ~ x1)
  single_x1<-single$coefficients[2]
  
  #Run IV two-stage IV
  x1Hat <- lm(x1 ~ z1)$fitted.values
  mrx1<-lm(y ~ x1Hat)
  mrx1_x1<-mrx1$coefficients[2]
  
  #Build result list
  out_list<-list(n_sample,score_corr,met_corr,full_x1,single_x1,mrx1_x1)
  
  #Retun list
  return(out_list)
}

#Test run
#mr_corr_sim(1000,0,0.5)

#Values to be considered
values_cons<-seq(0,1,0.1)

#Initialise matrix
mr_score_corr<-matrix(NA, nrow = length(values_cons), ncol = 6)

#Loop though values
for(i in 1:length(values_cons)){
  #Run mr_corr_sim function
  mr_score_corr[i,]<-unlist(mr_corr_sim(10000,values_cons[i],0.5))
  
}

#Results to dataframe
mr_score_dat<-as.data.frame(mr_score_corr)

#Name columns
names(mr_score_dat)<-c("N_of_sample","Gen_score_corr","Met_corr","Full_model","One_var_model","MR_model")

#Generate difference bwetween true and MR estimates
mr_score_dat$True_MR_per<-100*(-mr_score_dat$Full_model+mr_score_dat$MR_model)/mr_score_dat$Full_model

#Model between score correlation and bias oercentage
bias_mode<-lm(mr_score_dat$True_MR_per ~ mr_score_dat$Gen_score_corr)

# Add points
ggplot(data=mr_score_dat, aes(x=Gen_score_corr, y=True_MR_per, group=1)) +
  geom_abline(intercept = coef(bias_mode)[1], slope = coef(bias_mode)[2],
              colour="red", linetype="solid", size=1.5) +
  geom_point(colour="red", size=4, shape=21, fill="white") +
  labs(x = "Correlation between gene scores",
       y = "Percentage bias in MR estimate %") +
  theme_bw() +
  theme(axis.text = element_text(size=12))






