
n=112 #multi securities
Ei=175
#########
e=matrix(0,n,175)
r=res=matrix(0,n,175+31)

tP=tCross=tBW=tBMP=
  PvalP=PvalCross=PvalBW=PvalBMP=
  tMM=tLR=tMLR=tSLR=tMSLR=
  PvalMM=PvalLR=PvalMLR=PvalSLR=PvalMSLR=matrix(0,31,1)

tLRrho=tMLRrho=tSLRrho=tMSLRrho=
  PvalLRrho=PvalMLRrho=PvalSLRrho=PvalMSLRrho=matrix(0,31,1)

BFstar1=BFstar2=matrix(0,31,1)

AR=SAR=C=CCum=matrix(0,n,31)
Alpha=Beta=matrix(0,n,1)

for(i in 1:112)
{
  m=lm(R[i,1:175]~Rm[i,1:175])
  Alpha[i]=m$coef[1]
  Beta[i]=m$coef[2]
  e[i,]=m$residuals
  
  
  for(t  in 1:31)
  {
    
    AR[i,t]= R[i,180+t] - (Alpha[i] + Beta[i]*Rm[i,180+t])
    
    #estimation period
    Rmmean=mean(Rm[i,1:175])
    
    C[i,t] = (1+ (1/Ei) + ((Rm[i,180+t]-Rmmean)^2/(sum((Rm[i,1:175]-Rmmean)^2))))
    CCum[i,t] =   (1+ (t/Ei) +((sum(Rm[i,(180+1):(180+t)])- (t*Rmmean))^2)/(t*sum((Rm[i,1:175]-Rmmean)^2)))
    SAR[i,t]= AR[i,t]/(sd(e[i,])*sqrt(C[i,t]))
  }# for tau
}# for i in 1:n stocks

for(t in 1:31)
{
  #BW=================================================
  sum=0
  for(days in 1:175)
  {
    x=(mean(e[,days])-mean(e))^2
    sum= sum + x
  }
  denomBW=sqrt(sum/174)
  
  tBW[t]=mean(AR[,t])/denomBW
  PvalBW[t]=  2*(1-pt(abs(tBW[t]),174))
  
  tBWCum[t]= sum(tBW[1:t])/sqrt(t)  #coz denom is const
  PvalBWCum[t]=  2*(1-pt(abs(tBWCum[t]),174))
  
  
  #Patell====================================================================================
  
  tP[t]= sum(SAR[,t])/sqrt(n*(173/171))
  PvalP[t]=2*(1-pnorm(abs(tP[t])))
  
  tPCum[t]= sum(tP[1:t])/sqrt(t)
  PvalPCum[t]=2*(1-pnorm(abs(tP[t])))
  
  
  #BMP=====================================================================================
  tBMP[t]=sum(SAR[,t])/sqrt(var(SAR[,t])*n)
  PvalBMP[t]= 2*(1-pt(abs(tBMP[t]),n-1))
  
  
  tBMPCum[t]=sum(tBMP[1:t])/sqrt(t)
  PvalBMPCum[t]=2*(1-pt(abs(tBMPCum[t]),n-1))
  #=================================================================
  #MM=================
  #MM denom
  sum=0
  
  for(i in 1:n)
  {
    x=var(e[i,])*C[i,t]
    sum=sum+x
  }
  
  tMM[t]= n*mean(AR[,t])/sqrt(sum)
  PvalMM[t]=2*(1-pnorm(abs(tMM[t])))
  
  
  #=====================================================================================
  
  #Corrado===================
  
  #==================Do not touch coz coded from original paper
  
  for(i in 1:n)
  {
    res[i,1:175]=e[i,]
    res[i,176:(175+t)]=AR[i,t]
    r[i,]=rank(res[i,])
  }
  
  #Corr denom
  x=matrix(0,175+t,1)
  for(days in 1:(175+t))
  {
    x[days]= (1/n)*sum(r[,days]-mean(1:(Ei+1)))
  }
  
  denom=sqrt((1/(175+t))* sum(x^2))
  
  Corr[t]= (1/n)* sum(r[,175+t]- mean(1:(Ei+1)))/ denom
  CorrPval[t]=2*(1-pnorm(abs(Corr[t])))
  
  #===================Do not touch
  
  rCum = resCum = matrix(0,n,175+t)
  for(i in 1:n)
  {
    resCum[i,1:175]=e[i,]
    resCum[i,176:(175+t)]=AR[i,1:t]
    rCum[i,]=rank(resCum[i,])
  }
  sumrank=matrix(0,n,1)
  for(l in 1:t)
  {
    x = rCum[,175+l]
    sumrank = sumrank + x
  }
  
  CorrCum[t]= sqrt(1/n)* sum((sumrank - t * mean(1:(Ei+t+1)) ) / (
    sqrt(t*Ei*(t+Ei+1)/12)))
  CorrPvalCum[t]=2*(1-pnorm(abs(CorrCum[t])))
  
  
  #============================================================================================
  
  #=============== do not touch
  #Corrado & Zivney=================
  for(i in 1:n)
  {
    res[i,1:175]=e[i,]/sd(e[i,])
    res[i,175+t]=(AR[i,t]/(sd(e[i,])*C[i,t]))/(sd(AR[,t]/sd(AR[,t])))
    r[i,]=rank(res[i,])
  }
  
  CorrZ[t]=sqrt(1/n)* sum( (r[,175+t] -mean(1:(Ei+1))) /sd(1:(Ei+1)))
  CorrZPval[t]=2*(1-pnorm(abs(CorrZ[t])))
  
  #================ do not touch
  
  rCum = resCum = matrix(0,n,175+t)
  
  for(i in 1:n)
  {
    resCum[i,1:175]=e[i,]/sd(e[i,])
    #resCum[i,1:200]=e[i,]
    for(l in 1:t)
    {
      resCum[i,175+l]=AR[i,l]/(sd(e[i,])*C[i,l])/(sd(AR[,l]/sd(AR[,l])))
      #resCum[i,200+l]=SAR[i,l]
    }
    
    rCum[i,]=rank(resCum[i,])
  }
  
  sumrank=matrix(0,n,1)
  for(l in 1:t)
  {
    x = rCum[,175+l]
    sumrank = sumrank + x
  }
  
  CorrZCum[t]= sqrt(1/n)* sum((sumrank - t *
                                 mean(1:(Ei+t+1)))/(sqrt(t*Ei*(t+Ei+1)/12)))
  CorrZPvalCum[t]=2*(1-pnorm(abs(CorrZCum[t])))
  
  
  
  #=======================================================================
  
  #sign===========
  x=matrix(0,n,1)
  sign=matrix(0,n,31)
  for(l in 1:31)
  {
    for(i in 1:n)
    {
      x[i]=median(e[i,])
      if(AR[i,l]>x[i])
      {
        sign[i,l]=1
      }
      
      if(AR[i,l]<=x[i])
      {
        sign[i,l]=0
      }
    }
  }
  
  #Sign[t]= sqrt(n)^-1 * sum((sign[,t]- 0.5)/sqrt(0.5*0.5))
  Sign[t]= sqrt(n) *  (sum(sign[,t])/n - 0.5)/(0.5)
  SignPval[t]=2*(1-pnorm(abs(Sign[t])))
  
  
  
  x=sumsign=matrix(0,n,1)
  for(l in 1:t)
  {
    x=sign[,l]
    sumsign=sumsign+x
  }
  SignCum[t]=sqrt(n)^-1 * sum((sumsign- t*0.5)/sqrt(t*0.5*0.5))
  SignPvalCum[t]=2*(1-pnorm(abs(SignCum[t])))
  
  
  #==================================================================================
  #cross==========================
  tCross[t]=t.test(AR[,t])$statistic
  PvalCross[t]=t.test(AR[,t])$p.value
  
  sum=x=matrix(0,n,1)
  for(l in 1:t)
  {
    x=AR[,l]
    sum=sum+x
  }
  tCrossCum[t]=t.test(sum)$statistic
  PvalCrossCum[t]=t.test(sum)$p.value
  
  
  
  
}

#df=data.frame(PvalBW,PvalP,PvalBMP,PvalCross,PvalLR,PvalSLR,PvalMLR,PvalMSLR,PvalMM,
#              CorrPval,CorrZPval,SignPval)
#days=0:30
#df=data.frame(days,PvalBW,PvalP,PvalBMP,PvalCross,PvalMM,PvalLR,PvalMLR,PvalSLR,PvalMSLR)
#colnames(df)=c("days","Brown and Warner","Patell","BMP", "Cross sectional","Method of Moments","LR","Modified LR","SLR","Modified SLR")
#df=round(df,5)
#write.csv(df,file="9.csv")

#days=0:30
#df=data.frame(days,BFstar1,BFstar2)
#colnames(df)=c("days","BF*1","BF*2")
#df=round(df,5)
#write.csv(df,file="2.csv")

days=0:30
df=data.frame(days,PvalLRrho,PvalMLRrho,PvalSLRrho,PvalMSLRrho)
colnames(df)=c("days","LR","Modified LR","SLR","Modified SLR")
df=round(df,5)
write.csv(df,file="4.csv")
