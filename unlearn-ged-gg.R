library(ggplot2)
library(tikzDevice)
library(scales)

delta <- c("1_0","0_5")[1]

errmet <- "ged"

for(iter1 in 1:3){
  
  mval <- c("01","05","10")[iter1]

astr <- "C:/Users/arnab/Desktop/R-codes/unlearning-results/Output_del"

readfile <- gsub(" ", "", paste(astr, delta,"_","m",mval,".csv"))
writefile <- gsub(" ", "", paste(errmet,"_eps_del", delta,"_","m",mval,".png"))

dfmat <- read.csv(readfile, header = TRUE)

dfmat <- as.matrix(dfmat)
dimnames(dfmat) <- NULL
dfmat <- dfmat[,-1]

df <- matrix(NA,nrow(dfmat),12)

for(i in 1:nrow(dfmat)){
  
   k = 3
   datavec = (dfmat[i,])[(32*(k-1)+1):(32*k)]
  
  #datavec = sampler(300)
  df[i,] = c(matrix(datavec[-c(1,2)],ncol=6)[c(2,3),])
}

epsmat <- c(rbind(seq(0.25,1.5,0.25),seq(0.25,1.5,0.25)))

longdf <- data.frame(c(t(df)))
longdf$eps <- rep(epsmat, nrow(df))
longdf$noise <- rep(c("Laplace","Normal"),nrow(df)*6)

colnames(longdf)[1] <- "errs"

ggsave(
  writefile,
ggplot(longdf, aes(x = factor(eps), y = errs, fill = factor(noise))) +
  geom_boxplot(outlier.shape=NA, fatten = NULL) +
  labs(x = expression("Certifiability Parameter"~epsilon), y = "Generalization error divergence", fill = "noise") +
  ggtitle(bquote("GED vs."~epsilon *" for n = p = 1255 and m ="~ .(as.numeric(mval))))+
  scale_y_continuous(trans='log2',labels = label_number(accuracy = 0.001)) +
  theme_minimal() +
  scale_x_discrete(labels = as.factor(seq(0.25,1.5,0.25))) +
  scale_fill_discrete(name= "noise", labels= c("Laplace","Gaussian"))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)), # Tilt X-axis labels for readability,
width = 5,
height = 3,
dpi = 400
)
}