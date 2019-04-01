
sim_CDR = function(n, d, G, CDR_0, effects, time) {
  
  
  # n=1000
  CDR_0 = effect(append(CDR_0, n = n))

  noise = rmvnorm(n,  rep(0,length(time)), do.call(cbind,Mcor))
  # head(noise)
  noise = cbind(rep(0,n), noise)
  # as.vector(t(noise))[1:10]
  
  rate = rnorm(n, randomeff$rate, randomeff$sig)
  # rate[rate < 0] = 0
  CDR_t = lapply(1:n, FUN = function(x) c(CDR_0[x], CDR_0[x]+1:5*rate[x]))
  
  rate_CDR = rnorm(n, randomeff$rate_CDR, randomeff$sig_CDR)
  rate_CDR = pmin(0, CDR_0-4)*rate_CDR 
  c_CDR = unlist(lapply(1:n, FUN = function(x) c(0, 1:5*rate_CDR[x])))
  
  rate_ed = rnorm(n, randomeff$rate_ed, randomeff$sig_ed)
  ed = lapply(1:n, FUN = function(x) rep(rbinom(1, 1, 0.3),6))
  c_ed = unlist(lapply(1:n, FUN = function(x) ed[[x]]*c(0,1:5*rate_ed[x])))
  
  rate_gender = rnorm(n, randomeff$rate_gender, randomeff$sig_gender)
  gender = lapply(1:n, FUN = function(x) rep(rbinom(1, 1, 0.5),6))
  c_gender = unlist(lapply(1:n, FUN = function(x) gender[[x]]*c(0, 1:5*rate_gender[x])))
  
  
  A = lapply(1:n, FUN = function(x) rep(rbinom(1, 1, 0.5),6))
  rate_A = rnorm(n, randomeff$rate_A, randomeff$sig_A)
  effect_A = lapply(1:n, FUN = function(x) A[[x]]*c(0, 1:5*rate_A[x]))
  A = unlist(A)
  CDR_i = unlist(CDR_t) + 
    unlist(c_CDR) +
    unlist(c_gender) +
    unlist(effect_A) +
    c_ed +
    as.vector(t(noise))
  
  
  ID = unlist(lapply(1:n, FUN = function(x) rep(x,6)))
  CDR_0 = unlist(lapply(CDR_0, FUN = function(x) rep(x,6)))
  df_alz = data.frame(ID=ID, A = A, gender = unlist(gender), CDR_i = CDR_i, CDR_0 = CDR_0, ed = unlist(ed), time = time)
  df_alz
}

#' @export
effect = function(dist, n, params) {
  params = append(params, list(n = n))
  do.call(dist, params)
}

