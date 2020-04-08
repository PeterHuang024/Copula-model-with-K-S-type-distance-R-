Sur_risk = cbind(MIMI_Sprob,MICHF_Sprob,MIStroke_Sprob)
Route = c()
for (i in 1:31412) {
  a = Sur_risk[i,]
  max = max(a)
  if(max == a[1]){
    disease = 'MI_MI'
  }
  else if(max == a[2]){
    disease = 'MI_CHF'
  }
  else if(max == a[3]){
    disease = 'MI_Stroke'
  }
 
  
  Route <- rbind(Route, disease)
}
table(Route[,1])
#MI_CHF 10453, MI_MI 10611, MI_Stroke 10348