# 2019 10 30 I.Zliobaite

# experimental plot -- exponential growth, logistic growth, Lotka growth, Red Queen

W1 <- 1
S1 <- 1.3
S2 <- 1
E1 <- 1
E2 <- 999
EE <- E1 + E2

a_Lotka <- S1 - 1 - 0.25

cik <- 70

E_exp <- E1
E_log <- E1
E_Lotka <- E1
E_RQ <- E1



results_all <- c()

for (sk in 1:cik){
  E_exp <- E_exp*S1  
  E_log <- E_log + (S1-1)*E_log*(1 - E_log/EE)
  E_Lotka <- E_Lotka + (S1-1)*E_Lotka*(1 - E_Lotka/EE)*(1 - a_Lotka)
  E_RQ <- E_RQ*S1/(S1*E_RQ/EE + 1 - E_RQ/EE)
  results_all <- rbind(results_all,cbind(sk,E_exp/EE,E_log/EE,E_Lotka/EE,E_RQ/EE))
}
colnames(results_all) <- c('sk','E_exp','E_log','E_Lotka','E_RQ')

pdf('fig_FigA1.pdf',height = 4, width = 4.5)
plot(NA,NA,xlim = c(1,cik),ylim = c(0,1), xlab = 'Time step', ylab = 'Relative population size')
lines(results_all[,'E_exp'],col='grey', lwd = 2)
lines(results_all[,'E_log'],lty=5, lwd = 2)
lines(results_all[,'E_Lotka'],lty=3,col='black', lwd = 2)
lines(results_all[,'E_RQ'],lty=1, lwd = 2)
legend('bottomright',lwd = 2, col = c('gray','black','black','black'),lty = c(1,5,1,3),legend = c('Exponential growth (rmax = 0.3)','Logistic growth','Proportional prize','Lotka-Voltera competitive (a = 0.05)'),cex = 0.7, bty = 'n')
dev.off()
