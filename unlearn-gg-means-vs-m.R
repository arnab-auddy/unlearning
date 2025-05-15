library(ggplot2)
library(tikzDevice)

delta <- c("1_0","0_5")[1]

for(iter1 in 1:3){
for(iter2 in 1:2){
  
pval <- c("0500","1500","2500")[iter1]
errmet <- c("ged","frgt")[iter2]

astr <- "C:/Users/arnab/Desktop/R-codes/unlearning-results/Output_mvar_del"

readfile <- gsub(" ", "", paste(astr, delta,"_","p",pval,".csv"))
writefile <- gsub(" ", "", paste(errmet,"_m_del", delta,"_","p",pval,".png"))


df_full <- read.csv(readfile, header = TRUE)

df_full <- as.matrix(df_full)
dimnames(df_full) <- NULL
df_full <- df_full[,-1]

dfmat <- colMeans(df_full)

eps_i = 3  #chooses eps = 1

mvec <- rep(NA, 6)
l2errvec <- rep(NA, 6)
inerrmat <- matrix(NA, 6, 2)
gedmat <- matrix(NA, 6, 2)

for(k in 1:6){
  
  datavec = dfmat[(32*(k-1)+1):(32*k)]
  
  mvec[k] <- (datavec[1])
  l2errvec[k] <- (datavec[2])
  
  gedmat[k,] <- (matrix(datavec[-c(1,2)],ncol=6)[2:3,eps_i])
  
  inerrmat[k,] <- (matrix(datavec[-c(1,2)],ncol=6)[4:5,eps_i])
}

lm(log(l2errvec)~log(mvec))
lm(log(gedmat[,1])~log(mvec))
lm(log(gedmat[,2])~log(mvec))
lm(log(inerrmat[,1])~log(mvec))
lm(inerrmat[,2]~(mvec))


ged_df <- data.frame(rbind(cbind(gedmat[,1],mvec),cbind(gedmat[,2],mvec)))

ged_df$noise <- c(rep("Laplace",length(mvec)),rep("Gaussian",length(mvec)))

colnames(ged_df)[1:2] <- c("errs","m")

inerr_df <- data.frame(rbind(cbind(inerrmat[,1],mvec),cbind(inerrmat[,2],mvec)))

inerr_df$noise <- c(rep("Laplace",length(mvec)),rep("Gaussian",length(mvec)))

colnames(inerr_df)[1:2] <- c("errs","m")

if(errmet == "ged"){
  ggsave(
    writefile,
    ggplot(ged_df, aes(x = m, y = errs, color = noise)) +
      geom_point() +
      geom_line() +
      labs(x = expression("Unlearning size (m)"), y = "Generalization error divergence", fill = "Perturbation") +
      ggtitle(bquote("GED vs. m for n = p ="~ .(as.numeric(pval)) * " and " ~ epsilon * " = "~ .(eps_val)))+
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
    ggplot(inerr_df, aes(x = m, y = errs, color = noise)) +
      geom_point() +
      geom_line() +
      labs(x = expression("Unlearning size (m)"), y = "Error on unlearned set", fill = "Perturbation") +
      ggtitle(bquote("UED vs. m for n = p ="~ .(as.numeric(pval)) * " and " ~ epsilon * " = "~ .(eps_val)))+
      scale_y_continuous(trans='log10',labels = label_number(accuracy = 0.001)) +
      scale_x_continuous(trans='log10') +
      theme_minimal() +
      #scale_x_discrete(labels = mvec)+#as.factor(seq(0.25,1.5,0.25))) +
      scale_fill_discrete(name= "Perturbation", labels= c("Laplace","Gaussian"))+
      theme(legend.position = "none")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    width = 5,  height = 3,
    dpi = 400
  )
}
}}

