weibull_F1 <- function(x){
  1-exp(-(x/b_1)^a_1)
}
weibull_F2 <- function(y){
  1-exp(-(y/b_2)^a_2)
}

Gumbel_Cop <- function(x,y){
  F_1 <- weibull_F1(x)
  F_2 <- weibull_F2(y)
  a <- (-log(F_1))^gumbel_alpha
  b <- (-log(F_2))^gumbel_alpha
  return(exp(-(a+b)^(1/gumbel_alpha)))
}

clay_Cop <- function(x,y){
  F_1 <- weibull_F1(x)
  F_2 <- weibull_F2(y)
  return((F_1^(-clay_alpha)+F_2^(-clay_alpha)-1)^(-1/clay_alpha))
}

frank_Cop <- function(x,y){
  F_1 <- weibull_F1(x)
  F_2 <- weibull_F2(y)
  a <- (exp(-frank_alpha*F_1)-1)
  b <- (exp(-frank_alpha*F_2)-1)
  c <- (exp(-frank_alpha)-1)
  return((-1/frank_alpha)*log(1+(a*b)/c))
}

indep_Cop <-function(x,y){
  F_1 <- weibull_F1(x)
  F_2 <- weibull_F2(y)
  return(F_1*F_2)
}