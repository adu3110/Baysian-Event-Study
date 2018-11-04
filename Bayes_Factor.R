
 #Bayes Factor====================================================================
  Mu_i=Eta=Sig2_i=matrix(0,n,1)
  Psi=matrix(0,1,1)
  
  for(i in 1:n)
  {
    
    Uoi=Alpha[i]+ Beta[i]*Ret_m[i,180+t]
    Eta[i]=(SSEi+ (Koi/(1+Koi))*(Ret_i[i,180+t]- (Alpha[i]+Beta[i]*Ret_m[i,180+t]))^2)/2
    Sig2_i[i]=Eta[i]/((Ei-1)/2)
    Mu_i[i]=(Koi*Uoi + Ret_i[i,221])/(Koi+1)
  }
  Psi=sum((Ei-1)/(2*Eta))
  
  BFstar1[t]=exp(-1/2 * (Zsig^2/Psi))
  
  sum=0
  for(i in 1:n)
  {
    x= (Ret_i[i,180+t]-Mu_i[i]) * Eta[i]^(-0.5) *
      (gamma(Ei/2)/gamma((Ei-1)/2))
    sum=sum+x
  }
  
  ES= (1/n)* sum
  BFstar2[t]= exp(-n/2 * ES^2)
}


days=0:30
df=data.frame(days,PvalLRrho,PvalMLRrho,PvalSLRrho,PvalMSLRrho)
colnames(df)=c("days","LR","Modified LR","SLR","Modified SLR")
df=round(df,5)
write.csv(df,file="4.csv")
