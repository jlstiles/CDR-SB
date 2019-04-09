library(CDRSB)

library(CDRSB)

# These variables are predefined so as to be nicely explained
# time points, these could be 6 month intervals
time = 0:6

# correlation terms in correlation matrix, auto-regressive
rho = .5
sigma = .8
rho = .5
mm = rho^time
mm = c(mm[order(mm)], mm[2:length(mm)])
G = vapply(time, FUN = function(x) mm[(length(time)-x):(2*length(time)-x-1)], 
           FUN.VALUE = rep(1, length(time)))
G = G*sigma^2

vars = list(W1 = list(dist = rbinom, params = list(size = 1, prob = .5)),
            W2 = list(dist = rbinom, params = list(size = 1, prob = .5)),
            W3 = list(dist = rmultinom, params = list(size = 1, prob = c(.4, .4, 0.2))),
            W4 = list(dist = rnorm, params = list(mean = 15, sd = 4)),
            W5 = list(dist = rnorm, params = list(mean = 70, sd = 7.5)),
            W6 = list(dist = rnorm, params = list(mean = 3000, sd = 370))
)


beta_s = function(n, shape1, shape2, fac, addon) {
  (rbeta(n, shape1, shape2)*fac+addon)
}

random_effects = list(Y_0 = list(dist = runif, params = list(min = 1, max = 4)),
                      rate = list(dist = beta_s, params = list(shape1 = 2, 
                                                               shape2 = 3, fac = 1/1.2,addon = 0)),
                      # rate_A = list(dist = rnorm, params = list(mean = .05,
                      #                                            sd = .161))
                      rate_A = list(dist = beta_s, params = list(shape1 = 3,
                                                                 shape2 = 5, fac = -.5, addon = .122 ))
)


effects = list(f1 = function(W) {
  L = length(W[[1]])
  re = rnorm(mean = .06, sd = .2, n = L)
  e = unlist(lapply(1:L, FUN = function(x) W$W1[x]*time*re[x]))
  return(e)
},
f2 = function(W) {
  L = length(W[[1]])
  re = rnorm(mean = -.06, sd = .2, n = L)
  e = unlist(lapply(1:L, FUN = function(x) W$W2[x]*time*re[x]))
  return(e)
},
f3 = function(W) {
  L = length(W[[1]])
  re = rmvnorm(mean = c(-.01,.03,.07), 
               sigma = matrix(c(.05^2,0,0 ,0,.04^2,0,0,0, .1^2),nrow = 3), n = L)
  e = unlist(lapply(1:L, FUN = function(x) sum(W$W3[,x]*re[x,])*time))
  return(e)
},
f4 = function(W) {
  L = length(W[[1]])
  re = rbeta(shape1 = 2, shape2 = 1, n = L) - 2/3 + .03
  e = unlist(lapply(1:L, FUN = function(x) (W$W4[x]-15)/4*time*re[x]))
  return(e)
},
f5 = function(W) {
  L = length(W[[1]])
  re = rbeta(shape1 = 2, shape2 = 1, n = L) - 2/3 + .01
  e = unlist(lapply(1:L, FUN = function(x) (W$W5[x]-70)/7.5*time*re[x]*.1))
  return(e)
},
f6 = function(W) {
  L = length(W[[1]])
  re = rbeta(shape1 = 2, shape2 = 1, n = L) - 2/3 - .04
  e = unlist(lapply(1:L, FUN = function(x) (W$W6[x]-3000)/370*time*re[x]))
}
)

library(parallel)
test = lapply(1:5000, FUN = function(x) {
  df_alz = sim_CDR(n, time, G, W_ind = c(1,2,3,4,5,6), effect_ind = c(1,2,3,4,5,6), 
                   vars, effects, main_effects)
  df_alz1 = ddply(df_alz, "ID", .fun = function(x) {
    probs = pmin(1,exp(-3 + .189*x$Y_t)[1:(length(time)-1)])
    cens = rbinom((length(time)-1), 1, probs)
    probs
    if (any(cens==1)) {
      C = min(which(cens==1))+1
      return(x[1:C,])
    } else return(x)
  })
  test_alz <- lmer(Y_t ~ A:t +  t + (t | ID), df_alz1)
  ss = summary(test_alz)
  return(list(ss = ss, cover = abs(ss$coefficients[3,3])>=1.96))
},  mc.cores = getOption("mc.cores", 24L))

save(test, file = "test1.RData")
