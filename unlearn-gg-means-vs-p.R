library(ggplot2)
library(tikzDevice)

delta <- c("1_0","0_5")[1]

for(iter1 in 1:3){
  for(iter2 in 1:2){

mval <- c("01","05","10")[iter1]
errmet <- c("ged","frgt")[iter2]

astr <- "C:/Users/arnab/Desktop/R-codes/unlearning-results/Output_del"

readfile <- gsub(" ", "", paste(astr, delta,"_","m",mval,".csv"))
writefile <- gsub(" ", "", paste(errmet,"_p_del", delta,"_","m",mval,".png"))
  

df_full <- read.csv(readfile, header = TRUE)

df_full <- as.matrix(df_full)
dimnames(df_full) <- NULL
df_full <- df_full[,-1]

dfmat <- colMeans(df_full)

eps_i = 3  #chooses eps = 1

eps_val = seq(0.25,1.5,0.25)[eps_i]

pvec <- rep(NA, 6)
l2errvec <- rep(NA, 6)
inerrmat <- matrix(NA, 6, 2)
gedmat <- matrix(NA, 6, 2)

for(k in 1:6){
  
  datavec = dfmat[(32*(k-1)+1):(32*k)]
  
  pvec[k] <- (datavec[1])
  l2errvec[k] <- (datavec[2])
  
  gedmat[k,] <- (matrix(datavec[-c(1,2)],ncol=6)[2:3,eps_i])
  
  inerrmat[k,] <- (matrix(datavec[-c(1,2)],ncol=6)[4:5,eps_i])
}

lm(log(l2errvec)~log(pvec))
lm(log(gedmat[,1])~log(pvec))
lm(log(gedmat[,2])~log(pvec))
lm(log(inerrmat[,1])~log(pvec))
lm(log(inerrmat[,2])~log(pvec))


ged_df <- data.frame(rbind(cbind(gedmat[,1],pvec),cbind(gedmat[,2],pvec)))

ged_df$noise <- c(rep("Laplace",length(pvec)),rep("Normal",length(pvec)))

colnames(ged_df)[1:2] <- c("errs","p")

inerr_df <- data.frame(rbind(cbind(inerrmat[,1],pvec),cbind(inerrmat[,2],pvec)))

inerr_df$noise <- c(rep("Laplace",length(pvec)),rep("Normal",length(pvec)))

colnames(inerr_df)[1:2] <- c("errs","p")

if(errmet == "ged"){
  ggsave(
    writefile,
  ggplot(ged_df, aes(x = p, y = errs, color = noise)) +
    geom_point() +
    geom_line() +
    labs(x = expression("Dimension"), y = "Generalization error divergence", fill = "Perturbation")+
    ggtitle(bquote("GED vs. p for n = p, m ="~ .(as.numeric(mval)) * " and " ~ epsilon * " = "~ .(eps_val)))+
         #title=paste("GED vs. p for m=", as.numeric(mval),"and",expression(~epsilon),"=",as.numeric(eps_val)))+
    scale_y_continuous(trans='log10',labels = label_number(accuracy = 0.001)) +
    scale_x_continuous(trans='log10') +
    theme_minimal() +
    #scale_x_discrete(labels = pvec)+#as.factor(seq(0.25,1.5,0.25))) +
    scale_fill_discrete(name= "Perturbation", labels= c("Laplace","Gaussian"))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    width = 5,  height = 3,
    dpi = 400
  )
}

if(errmet == "frgt"){
  ggsave(
    writefile,
    ggplot(inerr_df, aes(x = p, y = errs, color = noise)) +
      geom_point() +
      geom_line() +
      labs(x = expression("Dimension"), y = "Error on unlearned set", fill = "Perturbation")+
      ggtitle(bquote("UED vs. p for n = p, m ="~ .(as.numeric(mval)) * " and " ~ epsilon * " = "~ .(eps_val)))+
      scale_y_continuous(trans='log10',labels = label_number(accuracy = 0.001)) +
      scale_x_continuous(trans='log10') +
      theme_minimal() +
      #scale_x_discrete(labels = pvec)+#as.factor(seq(0.25,1.5,0.25))) +
      scale_fill_discrete(name= "Perturbation", labels= c("Laplace","Gaussian"))+
      theme(legend.position = "none")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    width = 5,  height = 3,
    dpi = 400
  )
}
}}
