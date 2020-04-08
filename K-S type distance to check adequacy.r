library(copula)

rm(list = ls())

TIMES = 173 #simulate how many times?
time = 0 #count we already had simulated how many times

clay_s = array(0, dim = TIMES)
gumbel_s = array(0, dim = TIMES)
frank_s = array(0, dim = TIMES)
indep_s = array(0, dim = TIMES)

NUM = 100

while(time < TIMES){
  w1 = rweibull(NUM,1.79,220)
  w2 = rweibull(NUM,1.53,180)
  w = cbind(w1,w2)

  m = pobs(as.matrix(w))

  clay_model = claytonCopula(dim = 2)
  clay_fit = fitCopula(clay_model, m, method = 'itau')
  clay_alpha = coef(clay_fit)

  gumbel_model = gumbelCopula(dim = 2)
  gumbel_fit = fitCopula(gumbel_model, m, method = 'itau')
  gumbel_alpha = coef(gumbel_fit)

  frank_model = frankCopula(dim = 2)
  frank_fit = fitCopula(frank_model, m, method = 'itau')
  frank_alpha = coef(frank_fit)

  q = array(0, NUM)

  for(i in 1:NUM){
    temp = 0
    for(j in 1:NUM){
      if((w[j,1] < w[i,1]) && (w[j,2] < w[i,2])){
        temp = temp + 1
      }
    }
    q[temp+1] = q[temp+1] + 1
  }

  K_hat = cumsum(q) / NUM

  K <- function(t){
    t - t * log(t) #independence
  }

  clay_K <- function(t){
    t - (t^(1+clay_alpha)-t) / clay_alpha #clayton
  }

  gumbel_K <- function(t){
    t - (t * log(t)) / gumbel_alpha #gumbel
  }

  frank_K <- function(t){
    t - (exp(t*frank_alpha)-1) / frank_alpha * log((exp(-t*frank_alpha)-1)/(exp(-frank_alpha)-1)) #frank
  }

  indep_ans = 0
  clay_ans = 0
  gumbel_ans = 0
  frank_ans = 0

  i = 1
  while(i < NUM){
    indep_ans = indep_ans + (K_hat[i+1] - K(i/NUM))^2 * (K_hat[i+1] - K_hat[i])
    clay_ans = clay_ans + (K_hat[i+1] - clay_K(i/NUM))^2 * (K_hat[i+1] - K_hat[i])
    gumbel_ans = gumbel_ans + (K_hat[i+1] - gumbel_K(i/NUM))^2 * (K_hat[i+1] - K_hat[i])
    frank_ans = frank_ans + (K_hat[i+1] - frank_K(i/NUM))^2 * (K_hat[i+1] - K_hat[i])
    
    i = i + 1
  }

  clay_s[time+1] = clay_ans
  gumbel_s[time+1] = gumbel_ans
  frank_s[time+1] = frank_ans
  indep_s[time+1] = indep_ans

  time = time + 1

}

clay_s_cdf = ecdf(clay_s)
gumbel_s_cdf = ecdf(gumbel_s)
frank_s_cdf = ecdf(frank_s)
indep_s_cdf = ecdf(indep_s)

clay_s_cdf(0.001)
gumbel_s_cdf(0.001)
frank_s_cdf(0.001)
indep_s_cdf(0.001)