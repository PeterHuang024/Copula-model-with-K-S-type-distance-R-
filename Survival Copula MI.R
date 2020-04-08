library(readr)
rm(list = ls())
Sdata<-read_csv('DM_MI_X_time.csv')
MI_MI <- cbind(Sdata[,1],Sdata[,3])

#MI_Stroke <- cbind(Sdata[,1],Sdata[5])
time <- 0
TIME <- 100
k <- 0
clay_s = array(0, dim = TIME)
gumbel_s = array(0, dim = TIME)
frank_s = array(0, dim = TIME)
indep_s = array(0, dim = TIME)
NUM = 100

while (time < TIME){
  while (k < TIME) {
    a <- 100*k+1
    b <- 100*k+100
MI_MIx <- MI_MI[a:b,]

m_1 = pobs(as.matrix(MI_MIx))

NUM = 100

clay_model = claytonCopula(dim = 2)
clay_fit = fitCopula(clay_model, m_1, method = 'itau')
clay_alpha = coef(clay_fit)

gumbel_model = gumbelCopula(dim = 2)
gumbel_fit = fitCopula(gumbel_model, m_1, method = 'itau')
gumbel_alpha = coef(gumbel_fit)

frank_model = frankCopula(dim = 2)
frank_fit = fitCopula(frank_model, m_1, method = 'itau')
frank_alpha = coef(frank_fit)

q = array(0, NUM)

for(i in 1:NUM){
  temp = 0
  for(j in 1:NUM){
    if((MI_MIx[j,1] < MI_MIx[i,1]) && (MI_MIx[j,2] < MI_MIx[i,2])){
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
k = k+1
  }
  
}

clay_s_cdf = ecdf(clay_s)
gumbel_s_cdf = ecdf(gumbel_s)
frank_s_cdf = ecdf(frank_s)
indep_s_cdf = ecdf(indep_s)

clay_s_cdf(0.001)
gumbel_s_cdf(0.001)
frank_s_cdf(0.001)
indep_s_cdf(0.001)
# choose the best copula
MI_choose <- cbind(clay_s, gumbel_s, frank_s, indep_s)
Family = c()
for (i in 1:100) {
  a = MI_choose[i,]
  min = min(a)
  if(min == a[1]){
    copula = 'clay'
  }
  else if(min == a[2]){
    copula = 'gumbel'
  }
  else if(min == a[3]){
    copula = 'frank'
  }
  else if(min == a[4]){
    copula = 'indep'
  }
  Family <- rbind(Family, copula)
}
table(Family)
#CLAY 8, FRANK 23, GUMBLE 55, INDEP 14
#p-value: 0.21, 0.67, 0.53, 0.16
##relation
#theat=1.158776
#tau = 0.1370204
######################
#survival prob
library(EnvStats)
eweibull(MI_MI[,1])#SHAPE =0.9440789, SCALE = 135.8916674
eweibull(MI_MI[,2])#SHAPE = 1.248839, SCALE = 707.640737
m_MI = pobs(as.matrix(MI_MI))

gumbel_model = gumbelCopula(dim = 2)
gumbel_fitMI = fitCopula(gumbel_model, m_MI, method = 'itau')
gumbel_alpha = coef(gumbel_fitMI)


a_1 = 0.9441
b_1 = 135.89
a_2 = 1.2488
b_2 = 707.64
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
MIMI_prob = Gumbel_Cop(MI_MI[,1],MI_MI[,2])

MIMI_cdf_x = weibull_F1(MI_MI[,1])
MIMI_cdf_y = weibull_F2(MI_MI[,2])


one = rep(1,31412)
MIMI_Sprob = one - MIMI_cdf_x - MIMI_cdf_y + MIMI_prob

#plot surface
MIMI_ALL = cbind(MI_MI, MIMI_Sprob)
x <- c(MIMI_ALL[1:1000,1])
y <- c(MIMI_ALL[1:1000,2])
z <- c(MIMI_ALL[1:1000,3])
df <- data.frame(x=x,y=y,z=z)
library(ggplot2)
library(plotly)  
library(magrittr)
p <- plot_ly(data = df, x=~x,y=~y, z=~z, type = "contour", colorscale='Jet') %>%
  layout(title = 'DM -> MI -> MI',
         xaxis = list(title = 'DM -> MI (days)',zeroline = TRUE,range=c(0,600)),
         yaxis = list(title = 'MI -> MI (days)',range=c(0,2000)))
p

MIMI_ALL2 = cbind(MI_MI, MIMI_prob)
x <- c(MIMI_ALL2[1:1000,1])
y <- c(MIMI_ALL2[1:1000,2])
z <- c(MIMI_ALL2[1:1000,3])
df <- data.frame(x=x,y=y,z=z)
library(ggplot2)
library(plotly)  
library(magrittr)
p6 <- plot_ly(data = df, x=~x,y=~y, z=~z, type = "contour", colorscale='Jet') %>%
  layout(title = 'DM -> MI -> MI',
         xaxis = list(title = 'DM -> MI (days)',zeroline = TRUE,range=c(0,365)),
         yaxis = list(title = 'MI -> MI (days)',range=c(0,1460)))
plot(MI_MI[1:1000,1],MI_MI[1:1000,2],xlim = c(0,365), ylim = c(0,1460),
     main="DM -> MI -> MI",xlab="DM -> MI (days)", ylab="MI -> MI (days)")
clip(0,100,0,250)
abline(v=100,col="purple",lwd = 5)
clip(0,100,0,250)
abline(h=250,col="purple",lwd = 5)#18.4%
clip(0,200,0,500)
abline(v=200,col="yellow",lwd = 5)
clip(0,200,0,500)
abline(h=500,col="yellow",lwd = 5)#45.3%
clip(0,300,0,1000)
abline(v=300,col="green",lwd = 5)
clip(0,300,0,1000)
abline(h=1000,col="green",lwd = 5)#70.4%