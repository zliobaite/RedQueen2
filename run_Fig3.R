# 2019 10 19 I.Zliobaite
# hat patterns?
# small and large

resources <- 10000
N <- 10000
prob <- 0.75

#body mass range
bmr <- c(1,10,100,1000)
#bmr <- c(1,1,1,1)
all_all <- c()
for (sk in 1:length(bmr)){
  count_now <- round(bmr[sk]^(-0.75)*100000)
  all_all <- c(all_all,rep(bmr[sk],count_now))
}

#body_masses <- sample(all_all,N,replace = TRUE)
body_masses <- sample(bmr,N,replace = TRUE)

metabolic_rate <- body_masses^0.75
metabolic_rate <- N*metabolic_rate/sum(metabolic_rate)

abundances <- rep(0, N)
max_abundances <- rep(resources,N)/metabolic_rate
speed <- rep(0, N)
time_eating <- 1

set.seed(1985)

rec_abundances <- c()
rec_biomass <- c()
rec_alfa <- c()
rec_ene <- c()
rec_bio <- c()
for (sk in 1:N){
  print(sk)
  if (sk==1){
    new_spec <- 1  
  }else{
    #if (runif(1, 0, 1)>prob){
      new_spec <- rnorm(1,1,0.1)
      #new_spec <- 1.1  
    #}else{
      #new_spec <- -1
    #}
  }
  
  
  abundances[sk] <- max(metabolic_rate)/metabolic_rate[sk]
  #abundances[sk] <- metabolic_rate[sk]/min(metabolic_rate)
  #abundances[sk] <- 10
  #if (sk==1)
  #  speed[sk] <- 1
  #}else{
    speed[sk] <- max(new_spec,0)
  #}
    rec_alfa <- rbind(rec_alfa,new_spec)
  
  
  time_eating <- resources/sum(abundances*speed*metabolic_rate)
  #print(time_eating)
  energy_acquired_individual <- time_eating*speed*metabolic_rate
  #print(energy_acquired_individual)
  energy_acquired_population <- energy_acquired_individual*abundances
  biomass_population <- body_masses*abundances
  #print(energy_acquired_population)
  abundances <- energy_acquired_population/metabolic_rate  
  #print(abundances)

  ind1 <- which(body_masses==1)
  ind2 <- which(body_masses==10)
  ind3 <- which(body_masses==100)
  ind4 <- which(body_masses==1000)
  
  rec_ene <- rbind(rec_ene,cbind(sum(energy_acquired_population[ind1]),sum(energy_acquired_population[ind2]),sum(energy_acquired_population[ind3]),sum(energy_acquired_population[ind4])))
  #rec_ene <- log(rec_ene + 0.1)
  rec_bio <- rbind(rec_bio,cbind(sum(biomass_population[ind1]),sum(biomass_population[ind2]),sum(biomass_population[ind3]),sum(biomass_population[ind4])))
  ind <- which(energy_acquired_population<max(metabolic_rate))
  #ind <- which(energy_acquired_population<1)
  abundances[ind] <- 0
  speed[ind] <- 0
  rec_abundances <- rbind(rec_abundances,abundances/max_abundances)
  #rec_abundances <- rbind(rec_abundances,energy_acquired_population/resources)
  rec_biomass <- rbind(rec_biomass,abundances*body_masses)
  
  speed <- speed*time_eating
}

rec_abundances <- as.matrix(rec_abundances)

write.table(rec_abundances,file = 'ab_all.csv',row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE )

cols <- c('#e66101','#fdb863','#b2abd2','#5e3c99')
#cols_grey <- c('black','grey30','grey60','grey90')
cols_grey <- c('black','black','black','black')
lty_grey <- c(1,3,5,6)
lwd_grey <- c(2,2,2,2)

pdf('fig_Fig2_col.pdf',width = 15, height = 4)
plot(NA,NA,xlim = c(1,N),ylim = c(0,1), xlab = 'Time step', ylab = 'Relative abundance')
for (sk in 1:N ){
  if (body_masses[sk] == 1){
    cc <- 1
  }else{
    if (body_masses[sk] == 10){
      cc <- 2
    }else{
      if (body_masses[sk] == 100){
      cc <- 3
        }else{ cc <- 4}
    }
  }
  #print(cc)
  #print(max(rec_abundances[1:N,sk]))
  lines(rec_abundances[1:N,sk],col = cols[cc],lwd = 3)
}
dev.off()

pdf('fig_Fig2_black.pdf',width = 12, height = 3.5)
plot(NA,NA,xlim = c(1,N),ylim = c(0,1), xlab = 'Time step', ylab = 'Relative abundance')
for (sk in 1:N ){
  if (body_masses[sk] == 1){
    cc <- 1
  }else{
    if (body_masses[sk] == 10){
      cc <- 2
    }else{
      if (body_masses[sk] == 100){
        cc <- 3
      }else{ cc <- 4}
    }
  }
  #print(cc)
  #print(max(rec_abundances[1:N,sk]))
  lines(rec_abundances[1:N,sk],lwd = lwd_grey[cc], col = cols_grey[cc], lty = lty_grey[cc])
}
dev.off()

data_hats <- c()
for (sk in 1:N){
  ind <- which(rec_abundances[,sk]>0)
  if (length(ind)>0){
    #data_hats <- rbind(data_hats,c(sk,body_masses[sk],rec_alfa[sk],length(ind),max(rec_abundances[,sk])))
    data_hats <- rbind(data_hats,c(sk,body_masses[sk],rec_alfa[sk],length(ind),log(max(rec_abundances[,sk]))))  
  }
}
colnames(data_hats) <- c('sk','M','a','lng','mx')
write.table(data_hats,file = 'data_hats.csv',row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE )

library(ggplot2)

plot_wd <- 4.5
plot_ht <- 3.5

data_hats <- as.data.frame(data_hats)
pdf('fig_Fig3A_box.pdf',height = plot_ht,width = plot_wd)
boxplot(lng~M,data=data_hats, xlab="Body mass class", ylab="Species duration")
dev.off()

data_hats$M <- as.factor(data_hats$M)
pdf('fig_Fig3A.pdf',height = plot_ht,width = plot_wd)
p <- ggplot(data_hats, aes(x=M, y=lng)) + geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05 ) + labs(x="Body mass class", y = "Species duration", title = '(A)') + theme_minimal()
print(p)
dev.off()

pdf('fig_Fig3B_box.pdf',height = plot_ht,width = plot_wd)
boxplot(mx~M,data=data_hats, xlab="Body mass class", ylab="Log. peak relative abundance")
dev.off()

pdf('fig_sFig3B.pdf',height = plot_ht,width = plot_wd)
p <- ggplot(data_hats, aes(x=M, y=mx)) + geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05 ) + labs(x="Body mass class", y = "Log. peak relative abundance", title = '(B)') + theme_minimal()
print(p)
dev.off()

rec_eve_flat <- cbind(rep(1, N),rec_ene[,1])
rec_eve_flat <- rbind(rec_eve_flat,cbind(rep(10, N),rec_ene[,2]))
rec_eve_flat <- rbind(rec_eve_flat,cbind(rep(100, N),rec_ene[,3]))
rec_eve_flat <- rbind(rec_eve_flat,cbind(rep(1000, N),rec_ene[,4]))
colnames(rec_eve_flat) <- c('M','E')
rec_eve_flat <- as.data.frame(rec_eve_flat)
rec_eve_flat$M <- as.factor(rec_eve_flat$M)

pdf('fig_Fig3C_box.pdf',height = plot_ht,width = plot_wd)
boxplot(E~M,data=rec_eve_flat, xlab="Body mass class", ylab="Energy controlled")
dev.off()

pdf('fig_Fig3C.pdf',height = plot_ht,width = plot_wd)
p <- ggplot(rec_eve_flat, aes(x=M, y=E)) + geom_violin(trim=FALSE)  + labs(x="Body mass class", y = "Energy controlled", title = '(C)') + theme_minimal() + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1 )
print(p)
dev.off()

rec_bio_flat <- cbind(rep(1, N),rec_bio[,1])
rec_bio_flat <- rbind(rec_bio_flat,cbind(rep(10, N),rec_bio[,2]))
rec_bio_flat <- rbind(rec_bio_flat,cbind(rep(100, N),rec_bio[,3]))
rec_bio_flat <- rbind(rec_bio_flat,cbind(rep(1000, N),rec_bio[,4]))
colnames(rec_bio_flat) <- c('M','Mtot')
rec_bio_flat <- as.data.frame(rec_bio_flat)
rec_bio_flat$M <- as.factor(rec_bio_flat$M)

pdf('fig_Fig3D_box.pdf',height = plot_ht,width = plot_wd)
boxplot(Mtot~M,data=rec_bio_flat, xlab="Body mass class", ylab="Biomass")
dev.off()

pdf('fig_Fig3D.pdf',height = plot_ht,width = plot_wd)
p <- ggplot(rec_bio_flat, aes(x=M, y=Mtot)) + geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1 ) + labs(x="Body mass class", y = "Biomass", title = '(D)') + theme_minimal()
print(p)
dev.off()


