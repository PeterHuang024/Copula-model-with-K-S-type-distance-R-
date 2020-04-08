MI_Stroke <- cbind(Sdata[,1],Sdata[,5])
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
    MI_Strokex <- MI_Stroke[a:b,]
    
    m_1 = pobs(as.matrix(MI_Strokex))
    
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
        if((MI_Strokex[j,1] < MI_Strokex[i,1]) && (MI_Strokex[j,2] < MI_Strokex[i,2])){
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
Stroke_choose <- cbind(clay_s, gumbel_s, frank_s, indep_s)
Family = c()
for (i in 1:100) {
  a = Stroke_choose[i,]
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
# CLAY 9 , FRANK 25, GUMBEL 48, INDEP,18
#p-value = 0.19, 0.63, 0.57, 0.22
#####relation
#theat=1.180821
#tau = 0.1531316
######################
library(EnvStats)
eweibull(MI_Stroke[,1])#SHAPE =0.9440789, SCALE = 135.8916674
eweibull(MI_Stroke[,2])#SHAPE = 1.157, SCALE = 2179.975374
m_Stroke = pobs(as.matrix(MI_Stroke))

gumbel_model = gumbelCopula(dim = 2)
gumbel_fitStroke = fitCopula(gumbel_model, m_Stroke, method = 'itau')
gumbel_alpha = coef(gumbel_fitStroke)

a_1 = 0.9441
b_1 = 135.89
a_2 = 1.157
b_2 = 2179.975374
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
MIStroke_prob = Gumbel_Cop(MI_Stroke[,1],MI_Stroke[,2])

MIStroke_cdf_x = weibull_F1(MI_Stroke[,1])
MIStroke_cdf_y = weibull_F2(MI_Stroke[,2])

one = rep(1,31412)
MIStroke_Sprob = one - MIStroke_cdf_x - MIStroke_cdf_y + MIStroke_prob
#plot surface
MIStroke_ALL = cbind(MI_Stroke, MIStroke_Sprob)
x <- c(MIStroke_ALL[1:1000,1])
y <- c(MIStroke_ALL[1:1000,2])
z <- c(MIStroke_ALL[1:1000,3])
df <- data.frame(x=x,y=y,z=z)
library(ggplot2)
library(plotly)  
library(magrittr)
p3 <- plot_ly(data=df, x=~x,y=~y, z=~z,  type = "contour", colorscale='Jet') %>%
  layout(title = 'DM -> MI -> Stroke',
         xaxis = list(title = 'DM -> MI (days)',zeroline = TRUE,range=c(0,600)),
         yaxis = list(title = 'MI -> Stroke (days)',range=c(0,6000)))
p3

MIStroke_ALL2 = cbind(MI_Stroke, MIStroke_prob)
x <- c(MIStroke_ALL2[1:1000,1])
y <- c(MIStroke_ALL2[1:1000,2])
z <- c(MIStroke_ALL2[1:1000,3])
df <- data.frame(x=x,y=y,z=z)
library(ggplot2)
library(plotly)  
library(magrittr)
p4 <- plot_ly(data=df, x=~x,y=~y, z=~z,  type = "contour", colorscale='Jet') %>%
  layout(title = 'DM -> MI -> Stroke',
         xaxis = list(title = 'DM -> MI (days)',zeroline = TRUE,range=c(0,365)),
         yaxis = list(title = 'MI -> Stroke (days)',range=c(0,1460)))
plot(MI_Stroke[1:1000,1],MI_Stroke[1:1000,2],xlim = c(0,365), ylim = c(0,1460),
     main="DM -> MI -> Stroke",xlab="DM -> MI (days)", ylab="MI -> Stroke (days)")
clip(0,100,0,250)
abline(v=100,col="purple",lwd = 5)
clip(0,100,0,250)
abline(h=250,col="purple",lwd = 5)#3.9%
clip(0,200,0,500)
abline(v=200,col="yellow",lwd = 5)
clip(0,200,0,500)
abline(h=500,col="yellow",lwd = 5)#12.9%
clip(0,300,0,1000)
abline(v=300,col="green",lwd = 5)
clip(0,300,0,1000)
abline(h=1000,col="green",lwd = 5)#31%