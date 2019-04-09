library(boot)
############## Simulating Alzheimer's data

#' @title sim_CDR
#' @description Simulating Alzheimer's data
#' @param n sample size
#' @param time time intervals, as in, 0:5
#' @param G correlation matrix
#' @param vars patient characteristics
#' @param effects effects placed on the patient characteristics 
#' @param main_effects specifying for Y_0, rate and rate_A (treatment reduction of rate)
#' @export
sim_CDR = function(n, time, G, W_ind, effect_ind, vars, effects, main_effects) {
  
  # baseline outcome
  # undebug(effect)
  Y_0 = do.call(effect, append(main_effects$Y_0, list(n = n)))
  
  # baseline rate (control arm)
  rate = do.call(effect, append(main_effects$rate, list(n = n)))

  Y_t = lapply(1:n, FUN = function(x) Y_0[x]+time*rate[x])
  Y_t = unlist(Y_t)
  # rate for the treated
  rate_A = do.call(effect, append(main_effects$rate_A, list(n = n)))
  
  # randomly assign treatment
  A = lapply(1:n, FUN = function(x) rep(rbinom(1, 1, 0.5), length(time)))
  effect_A = lapply(1:n, FUN = function(x) A[[x]]*(time*rate_A[x]))
  A = unlist(A)
  
  noise = rmvnorm(n,  rep(0,length(time)), G)
  noise = as.vector(t(noise))
  # noise = t(vapply(1:n, FUN = function(x) runif(length(time)), FUN.VALUE = rep(1,length(time))))
  
  if (!is.null(W_ind)) {
  W = lapply(vars[W_ind], FUN = function(x) {
      temp = do.call(effect, append(list(dist=x$dist,params = x$params), list(n=n)))
    })
  if (length(W)==1) {
    othereffects = lapply(effects[effect_ind], FUN = function(eff) eff(W))
    } else {
    othereffects = as.data.frame(lapply(effects[effect_ind], FUN = function(eff) eff(W)))
    othereffects = rowSums(othereffects)
  }


  W = lapply(W, FUN = function(x) {
    if (is.matrix(x)) {
      temp = t(x) 
      temp = vapply(1:ncol(temp), FUN = function(c) {
        tt= unlist(lapply(temp[,c], FUN = function(a) rep(a, length(time))))
        return(tt)
        },FUN.VALUE = rep(1:n, length(time)))
      } else {
      temp = unlist(lapply(x, FUN = function(a) rep(a, length(time))))
    }
    return(temp)
  }
  )
  if (length(W)==1) W = as.data.frame(W) else W = as.data.frame(do.call(cbind, W))
  Y_t = Y_t + noise + othereffects + unlist(effect_A)
  Y_t = pmax(0,Y_t)
  } else {
    Y_t = Y_t + noise + unlist(effect_A)
    # Y_t = pmax(0,Y_t)
  }
  ID = unlist(lapply(1:n, FUN = function(x) rep(x,length(time))))
  t = rep(time, n)

  
  # MAR = lapply(1:n, FUN = function(x) {
  #   cens = rbinom(length(time)-1, 1, c(.97,.97^2, .97^3, .97^4))
  #   if (any(cens==0)) {
  #     tt = min(which(cens==0))+1
  #     aa = rep(1,length(time))
  #     aa[tt:length(time)] = 0
  #     } else {
  #       aa = rep(1,length(time))
  #     }
  #   return(aa)
  #   })
  
  # MAR = lapply(1:n, FUN = function(x) {
  #   cens = rbinom(length(time)-1, 1, c(.97,.97^2, .97^3, .97^4))
  #   return(c(1,cens))
  # })
  # 
  # MAR = unlist(MAR)
  
  if (!is.null(W_ind)) df = cbind(ID, W, A, t, Y_t) else {
    df = cbind(ID, A, t, Y_t)
  }
  return(as.data.frame(df))
}

#' @export
effect = function(n, dist, params) {
  params = append(params, list(n = n))
  do.call(dist, params)
}

