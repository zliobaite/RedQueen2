# 2019 11 04 I.Zliobaite
# logistic: proportional prize competition

pop1 <- 0.01
pop2 <- 1 - pop1
pop11 <- pop1

alfa <- 1.01
cik <- 1000

pop <- c(pop1,pop2,pop11)
  
for (sk in 1:cik){
  
  beta <- 1/ (alfa*pop1 + 1 - pop1)
  
  pop1 <- pop1*alfa*beta
  pop2 <- pop2*beta
  pop11 <- pop11*alfa
  
  pop <-rbind(pop,c(pop1,pop2,pop11))
}

pdf('fig_Fig1A.pdf',width = 4.5, height = 4)
plot(c(0:cik),pop[,1],xlim = c(1,cik),ylim = c(0,1), xlab = 'Time step', ylab = 'Relative population size',type = "l", lwd = 3, main = 'Expanding species (A)')
lines(c(0:cik),pop[,3],lty = 2 )
legend('bottomright',lty= c(2,1), lwd = c(1,3), legend = c('Exponentail growth','Proportional prize competition'),cex = 0.8, bty = 'n')
dev.off()

pdf('fig_Fig1B.pdf',width = 4.5, height = 4)
plot(c(0:cik),pop[,2],xlim = c(1,cik),ylim = c(0,1), xlab = 'Time step', ylab = 'Relative population size',type = "l", lwd = 3, main = 'Declining species (B)')
dev.off()
