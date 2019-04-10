library(simcausal)
t.end <- 6
options(simcausal.verbose=FALSE)
B=5
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
  node("N", t = 0, distr = "rbern",prob = 0.5) +
  node("A", t = 0, distr = "rbern", prob = 0.5) + 
  node("C", t = 0, distr = "rbern", prob = 0) + 
  node("Y", t = 0, distr = "runif", params = list(min = 1, max = 4))


A_rate = (.334 + -.06*.5 + .03*.4 + .07*.2)/4
A_rate
# define variables based on past values in time
D <- D + 
  node("A", t = 1:t.end, distr = "rbern",
       prob = {A[0]}) + 
  node("N", t = 1:t.end, distr = "rbern",
       prob = 0) +
  node("C", t = 1:t.end, distr = "rbern",
       prob = {exp(-4 + .189*Y[t-1])},
       EFU = TRUE) +
  node("Y", t = 1:t.end, distr = "rnormT",
       mean = Y[0] + (.334 + N[0]*0 - .06*male[0] + 0.03*(APOE[0] == 2) + 
                        0.07*(APOE[0] == 3) + (age[0] - 70)/7.5*(.01) +
                        .02*(ed[0] - 15)/4 + (vol[0] - 3000)/370*(-.04) -.0825*A[0])*t, sd = 2, EFU = FALSE)

lDAG <- set.DAG(D)

act_A <-c(node("A", t = 0:6, distr = "rbern", prob = setA),
          node("C", t = 0:6, distr = "rbern", prob = 0))

lDAG = lDAG + action("A1", nodes = act_A, setA = 1) + action("A0", nodes = act_A, setA = 0) 

time = 0:6
# Y1_t = lapply(time, FUN = function(tt) {
#     Dt <- set.targetE(lDAG, outcome = "Y", t = tt, param = "A1")
#     return(eval.target(Dt, n = 100000)$res)
#   })
# 
# Y0_t = lapply(time, FUN = function(tt) {
#   Dt <- set.targetE(lDAG, outcome = "Y", t = tt, param = "A0")
#   return(eval.target(Dt, n = 100000)$res)
# })
# 
# 
# MSM_df = data.frame(Y = c(unlist(Y1_t), unlist(Y0_t)), t = rep(time, 2), A = c(rep(1,7), rep(0,7)))
# MSM_df
# 
# MSM_fit = glm(formula = formula("Y~ t + A:t"), family = gaussian, data = MSM_df)
# coef(MSM_fit)[3]*4 + coef(MSM_fit)[2]
#  
# correlation terms in correlation matrix, auto-regressive

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

  cutout = sample(1:1600, 400)
  temp = lapply(1:1600, FUN = function(x) ifelse(x %in% cutout, 1, 0))
  temp = unlist(lapply(temp, FUN = function(x) rep(x, length(time))))
  temp[1:50]
  
  df_alz = lapply(1:n,FUN = function(x) {
      df = subset(df_alz, ID==x)
      if (nrow(df) <= 5) return(df) else {
      if (temp[df_alz$ID==x][1] == 1) return(df[1:5, ]) else return(df)
    }
  })
  
  df_alz = do.call(rbind, df_alz)
  test_alz <- lmer(Y ~ A:t +  t + (t | ID),data= df_alz)
  ss = summary(test_alz)
  return(list(ss = ss$coefficients, cover = abs(ss$coefficients[3,3])>=1.96))
  }, mc.cores = getOption("mc.cores", 24L))

save(test, file = "test_i.RData")

