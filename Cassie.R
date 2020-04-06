#simulation function  
popsim <- function(N.init,lamb,se.lamb,sd.lamb,se.sd.lamb,sims,yrs){
  N <- matrix(NA,nrow <- yrs,ncol <- sims) #set up results storage
  #pull from sampling distributions for lambda and stochasticity in lambda (sd.lamb)
  r.lamb <- r.sd.lamb <- rep(NA,sims)
  r.lamb <- rnorm(sims,lamb,se.lamb) 
  r.sd.lamb <- rnorm(sims,sd.lamb,se.sd.lamb)
  N[1,] <- N.init #initialize for all sims
  for(k in 1:sims){
    #based on values from sampling distributions, simulate over years with stochasticity
    for(y in 2:yrs){
      N[y,k] <- rpois(1,N[(y-1),k]*rnorm(1,r.lamb[k],r.sd.lamb[k]))
    } #yrs 
  } #sims
return(N)
} #function

#function for plotting and reporting key results 
plot.fun <- function(N){
  yrs <- nrow(N)
  sims <- ncol(N)
  N.mean <- round(apply(N,1,mean), dig = 0)
  #plot each simulation 
  plot(x = 1:yrs,y = N[,1], type = 'l', col = "black",ylim = c(0,1000))
  for(k in 2:sims){
    lines(x = 1:yrs,y = N[,k], col = "black")
  }
  lines(x = 1:yrs,y = N.mean, col = "red") #plot mean across simulations 
  #calculate geometric mean lambda
  geom.lamb <- rep(NA,sims)
  for(k in 1:sims){
    geom.lamb[k] <- (N[yrs,k]/N[1,k])^(1/(yrs-1))
  }
  lamb.mean <- mean(geom.lamb)
  #calculate extinction probability 
  extinct.prob <- ifelse(length(which(N[yrs,]==0))==0,0,length(which(N[yrs,]==0))/sims)
  #bundle and name output 
  output <- c(round(lamb.mean,dig=5),round(extinct.prob,dig=5))
  names(output) <- c("Mean Lambda", "Extinction Probability")
  return(output)
}

#simulate 
sims <- 10000 #simulations to run - do enough to stabilize output from repeated runs of the function 
yrs <- 20 #years to simulate forward 

#    popsim(N.init, lamb, se.lamb, sd.lamb, se.sd.lamb, sims, yrs)

#here we start where sd.lamb and se.lamb are equal
N <- popsim(100,    1.01, 0.01,    0.01,    0.001,       sims, yrs)
plot.fun(N)

#notice what happens when we increase the sd.lamb - the variability over years
#simply by increasing the variability in the system, we increase extinction probability
#but it doesn't increase it that much, because for any simulation, r.lambda is very likely to be >1
N <- popsim(100,    1.01, 0.01,    0.1,    0.001,       sims, yrs)
plot.fun(N)

#now notice what happens when we increase the se.lamb - our uncertainty in mean lambda
#in this case it has a much bigger effect *with the same increase*, because now, for any simulation, r.lambda has a meaningful chance of being <1  
N <- popsim(100,    1.01, 0.1,    0.01,    0.001,       sims, yrs)
plot.fun(N)

#lets now see how we would get data like these to estimate these parameters
#lets say that we have d.yrs of data from each of d.reps different sites
#where we think that each of the different sites are subject to the same stochasticity
#we can think of site as a replicate of the same process, but with sampling uncertainty
#i.e., a good year in one site is a good year in another
#lets generate such data and then extract the observed lambdas 
#them model them in a random effects model 

#note that I would prefer to do all this in JAGS or similar because it seems more transparent to me! 

set.seed(1) #set the random seed so we can get the same output each time
d.yrs <- 20 #number of years
d.reps <- 100 #number of sites 
lamb <- 1.01 #the mean lambda across all years (and sites)
sd.lamb <- 0.05 #the stochasticity in lambda across years
se.lamb <- 0.01 #sampling uncertainty from rep to rep
r.lamb <- rnorm(d.yrs,lamb,sd.lamb) #draw true lambda for each year 
N.data <- matrix(NA,nrow=d.yrs,ncol=d.reps) #hold data 
N.data[1,] <- 100 #initialize
for(k in 1:d.reps){
  for(y in 2:d.yrs){ #calculate N 
   N.data[y,k] <- N.data[y-1,k]*rnorm(1,r.lamb[y],se.lamb)
  }#sims
}#yrs

#now extract the realized lambdas 
est.lamb <- matrix(NA,nrow=(d.yrs-1),ncol=d.reps)
for(y in 1:(d.yrs-1)){
  est.lamb[y,] <- N.data[y+1,]/N.data[y,]
}
#set data up to go into the random effects model 
vec.lamb <- as.vector(est.lamb)
year <- rep(1:(d.yrs-1),d.reps)
#run the random effects model 
model <- lmer( vec.lamb ~ 1 |year)

#look at an overall summary
summary(model)
#estimate of lamb (mean lambda)
summary(model)$coefficients[1]
#estimate of se lambda (sampling uncertainty in lambda over replicates)
summary(model)$coefficients[2]
#estimate of sd lambda (variation in lambda over time)
sqrt(as.data.frame(VarCorr(model))$vcov[1])

