#Exercise 11
#Hannah Gordon
#Due: November 13, 2020

#Imagine a cancer cell in a tumor that spontaneously exhibited a mutation that confers drug resistance. The
#mutation does not have any positive or negative effects on growth rate of that sub-population when the
#cancer drug is absent. However, when the cancer drug is present the mutant sub-population grows at 50% of
#its growth rate in the absence of the drug and the non-mutant sub-population declines rapidly. The model
#we will use to represent the growth of the two sub-populations is this:
#  Nt+1 = Nt + rN Nt(1 ???(Nt+Mt)/K)
#  Mt+1 = Mt + rMMt(1 ???(Nt+Mt)/K)
#Assume in the absence of the cancer drug that the cells grow at a rate of 0.1 per day (rN = rM = 0.1) and the
#carrying capapcity (K) of the tumor is one million cells. The mutation of a single cell occurred early in the
#tumor growth and when it occurred there were 100 total cells in the tumor. Drug treatment of non-mutant
#cells results in a negative growth rate of -0.1.

#Generate a script that simulates growth of the two sub-populations in the tumor to equilibrium
#followed by drug treatment. Plot your results using a line graph.

#set working directory
setwd("/Users/Public/Documents/Biocomputing2020_Tutorial13")

#define parameters
#r is growth rate, same for both normal and mutants, and k is carrying capacity
r=0.1
K=1000000
#rDN is growth rate of normal/wild-type cells with drug, rDM is growth rate of mutant cells with drug
rDN=-0.1
rDM=0.05

#set timesteps to simulate
times<-1:600

#place to store output
output<-matrix(data=NA,nrow=length(times),ncol=3)
output[,1]=times
#because mutation occured at 100 cells, set one column to start at 1 (mutant) and the next to 99(normal)
output[1,2]=1
output[1,3]=99

#create a for loop to simulate the populations to equilibrium
for(i in times[-1]){ 
  if(times[i]<=175){ 
    output[i,2]=output[(i-1),2]+r*output[(i-1),2]*(1-(output[(i-1),2]+output[(i-1),3])/K)
    output[i,3]=output[(i-1),3]+r*output[(i-1),3]*(1-(output[(i-1),2]+output[(i-1),3])/K)
    
  }else{ 
    #drug introduced after 175 timesteps
    output[i,2]=output[(i-1),2]+rDM*output[(i-1),2]*(1-(output[(i-1),2]+output[(i-1),3])/K)
    output[i,3]=output[(i-1),3]+rDN*output[(i-1),3]*(1-(output[(i-1),2]+output[(i-1),3])/K)
    
  }
}
#Output the data
outputData<-data.frame(time=output[,1], Mutants=output[,2], Wildtypes=output[,3])

#Make sure the for loop worked
tail(outputData)

#Plot the data with a line graph
library(ggplot2)
ggplot()+
  geom_line(data = outputData, aes(x=time,y=Mutants),color='steelblue')+ #mutants as steelblue
  geom_line(data = outputData, aes(x=time,y=Wildtypes),color='maroon')+ #wildtype cells as maroon
  xlab('Time')+
  ylab('Number of cells')

#done :)
