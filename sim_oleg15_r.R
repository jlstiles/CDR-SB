library(simcausal)
library(mvtnorm)
library(lme4)
t.end <- 6
options(simcausal.verbose=FALSE)
B=5000

rmnom = function(n, p1,p2,p3) {
  mm = rmultinom(n, prob = c(p1,p2,p3), size = 1)
  as.factor(unlist(lapply(1:n, FUN = function(x) which(mm[,x] == 1))))
}

rnormT = function(n, mean, sd) {
  pmax(rnorm(n, mean, sd), 0)
}

# first initialize
D <- DAG.empty()
D <- D +
  node("male", t = 0, distr = "rbern",
       prob = 0.5) +
  node("APOE", t = 0, distr = "rcat.factor",
       probs = c(.4,.4,.2)) +
  node("age", t = 0, distr = "rnorm",
       mean = 70, sd = 7.5) + 
  node("ed", t = 0, distr = "rnorm", mean = 15, sd = 4) +
  node("vol", t = 0, distr = "rnorm", mean = 3000, sd = 370) +
  node("comed", t = 0, distr = "rbern", prob = 0.5) + 
  node("Y", t = 0, distr = "runif", params = list(min = 1, max = 4))+
  node("A", t = 0, distr = "rbern", prob = 0.5) + 
  node("C", t = 0, distr = "rbern", prob = 0) 



A_rate = (.334 + -.06*.5 + .03*.4 + .07*.2)/4
A_rate
# define variables based on past values in time
D <- D + 
  node("Y", t = 1:t.end, distr = "rnormT",
       mean = Y[0] + (.334 - .06*male[0] + 0.03*(APOE[0] == 2) + 
                        0.07*(APOE[0] == 3) + (age[0] - 70)/7.5*(.01) +
                        .02*(ed[0] - 15)/4 + (vol[0] - 3000)/370*(-.04) - (3/5)*.0825*A[0])*t, sd = 2, EFU = FALSE) +
  node("A", t = 1:t.end, distr = "rbern",
       prob = {A[0]}) + 
  node("C", t = 1:t.end, distr = "rbern",
       prob = {exp(-4 + .189*Y[t-1])},
       EFU = TRUE)


lDAG <- set.DAG(D)

act_A <-c(node("A", t = 0:6, distr = "rbern", prob = setA),
          node("C", t = 0:6, distr = "rbern", prob = 0))

lDAG = lDAG + action("A1", nodes = act_A, setA = 1) + action("A0", nodes = act_A, setA = 0) 

time = 0:6

beta_s = function(n, shape1, shape2, fac, addon) {
  (rbeta(n, shape1, shape2)*fac+addon)
}

tester = beta_s(1e6, 2,5,4,addon=0)
hist(tester)
mean(tester)
addon = - mean(tester) 
rate_test = beta_s(1e5, 2,5,4,addon=addon)
hist(rate_test)
mean(rate_test)
sd(rate_test)

rateA_test = beta_s(1e6, 2,2,1/3,0)
addon1 = -mean(rateA_test)
rateA_test = beta_s(1e5, 2,2,1/3,addon=addon1)
hist(rateA_test)
mean(rateA_test)
sd(rateA_test)


rho = .5
sigma = .8
rho = .5
sigma = 1
mm = rho^time
mm = c(mm[order(mm)], mm[2:length(mm)])
G = vapply(time, FUN = function(x) mm[(length(time)-x):(2*length(time)-x-1)], 
           FUN.VALUE = rep(1, length(time)))
G = G*sigma^2
G[,1] = G[1,]=0

library(parallel)

test = mclapply(1:B, FUN = function(x) {
  n=1600
  df_alz <- sim(DAG = lDAG, n=n, wide = FALSE)
  
  rateA = beta_s(n, 2,2,1/3,addon=addon1)

  df_alz$rateA = unlist(lapply(1:n, FUN = function(x) {
    df_alz[df_alz$ID==x,]$t*rateA[x]
  }))
  
  rate = beta_s(n, 2,5,4,addon=addon)
  # rate = rnorm(n, 0,5)
  df_alz$rate = unlist(lapply(1:n, FUN = function(x) {
    df_alz[df_alz$ID==x,]$t*rate[x]
  }))
  
  noise = rmvnorm(n,  rep(0,length(time)), G)
  noise = as.vector(t(noise))
  
  df_alz$noise = unlist(lapply(1:n, FUN = function(x) {
    nn = sum(df_alz$ID==x)
    noise[1:nn +(length(time))*x]
  }))
  
  df_alz = df_alz[!is.na(df_alz$Y),]
  df_alz$Y = pmax(df_alz$Y + df_alz$noise + df_alz$rate + df_alz$rateA, 0)

  test_alz <- lmer(Y ~ A:t +  t + (t | ID),data= df_alz)
  ss = summary(test_alz)
  return(list(ss = ss$coefficients, cover = abs(ss$coefficients[3,3])>=1.96))
  }, mc.cores = getOption("mc.cores", 24L))


save(test, file = "test15_r.RData")

