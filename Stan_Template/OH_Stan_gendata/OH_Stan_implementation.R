library(rethinking)
library(rstan)
library(tidybayes)
library(tidyverse)
source("gen_OH_data.R")
set.seed("23627632")

#Generate random data
N_id = 2
dataList <- gen_OH_data(N_id = N_id, group = "both", N_trials_test_sess = 20) #2 mixed and 2 uniform group IDs, no test for now
str(dataList)

dataList$log_dur

init_fun <- function() {
  list(
    logit_phi = matrix(0.0, nrow = dataList$N_id, ncol = dataList$N_sess),
    log_lambda = matrix(0.0, dataList$N_id, dataList$N_sess),
    V_p = matrix(0.0, dataList$N_id, dataList$N_food),
    z_ID = matrix(0.0, 3, dataList$N_id),
    sigma_ID = rep(1.0, 3),
    Rho_ID = diag(1.0, 3)
  )
}

####Run Model####
options(mc.cores = parallel::detectCores()-2)

nIter     = 2000
nChains   = 1
nWarmup   = floor(nIter/2)
nThin     = 1


modelFile = 'OH_stan.stan'

fit = stan(modelFile, 
           data    = dataList, 
           chains  = nChains,
           iter    = nIter,
           warmup  = nWarmup,
           thin    = nThin,
           #init = init_fun
           ) 

fit

#Extract posterior samples
post <- rstan::extract(fit, pars = c("phi","lambda","alpha","beta","alpha_T","Prob_F"))

#Print mean alpha_T for ID1 for each container and both food types
for (i in 1:11){
  a = mean(post$alpha[,1,i,1]) #[all samples, ID 1, session i, food type 1]
  b = mean(post$alpha[,1,i,2]) #[all samples, ID 1, session i, food type 2]
  print(paste(round(a, digits = 2),round(b, digits = 2)))
}

#Create alpha_T data frame for plotting
plot <- data.frame(alpha_T = post$alpha_T[,,4][,224])

#Plot alpha T 
ggplot(plot, aes(x = alpha_T))+
  geom_density()+
  scale_x_continuous(limits = c(0,1))+
  theme_minimal()

#Animate evolution of alpha_T for ID1 across all trials
for (i in 1:120){#Trials for ID 1
  plot <- data.frame(alpha = post$alpha_T[,,4][,i])
  
  p <- ggplot(plot, aes(x = alpha))+
    geom_density()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,100))+
    theme_minimal()
  
  plot(p)
  
  Sys.sleep(1)
}

#Animate evolution of Prob_F for ID1 in test
for (i in 101:120){#ID1 test trials 
  plot <- data.frame(pf = post$Prob_F[,,4][,i])
  
  p <- ggplot(plot, aes(x = pf))+
    geom_density()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,100))+
    theme_minimal()
  
  plot(p)
  
  Sys.sleep(1)
}

#Animate evolution of Prob_F for ID2 in test
for (i in 221:240){#ID2 test trials
  plot <- data.frame(pf = post$Prob_F[,,4][,i])
  
  p <- ggplot(plot, aes(x = pf))+
    geom_density()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,100))+
    theme_minimal()
  
  plot(p)
  
  Sys.sleep(1)
}


dim(post$Prob_F) # 1000 iterations, 480 trials, 5 columns

iterations = 1000
N = dataList$N
ids = 1:dataList$N_id

#Flatten the array and create a data frame
df_pF <- data.frame(
  iteration = rep(1:iterations, times = N),
  trial = rep(1:N, each = iterations),
  id = rep(ids, each = iterations),
  pF1 = as.vector(post$Prob_F[,,4]),
  pF2 = as.vector(post$Prob_F[,,5]) 
)

pF_ID1 <- df_pF[df_pF$id == 1,]
plot <- pF_ID1[pF_ID1$trial == 1,]

ggplot(plot, aes(x = pF1)) +
  geom_density(alpha = 0.5)+
  scale_x_continuous(limits = c(0,1))+
  theme_minimal()

