### library used ###

library(ggplot2)
library("lmerTest", lib.loc="~/R/win-library/3.5")
library(cowplot)

############# study of global polymorphism ##########

# evolution of pi, Ho, MAF in region of 25kb around S locus and in control region
 ###### test unilateral ######
 ####### use csv file "VCF_analysis_mean" of each pop generate by python pipeline ###########
 ####### one analysis by population, one synthetic figure by analysis*specie ###########


Summary=NULL # synthesis table

 #halleri

pop="Japan"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x3=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  ylab("pi of all sites")+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF>=max3)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


pop="Nivelle"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x7=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of all sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


pop="Mortagne"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=20) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x11=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of all sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF>=max3)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("*", "**", "***", "**","**", "**"), c(0.35,0.52, 0.33,0.52,0.34,0.52), c(0.92,0.92, 0.59,0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.03, 0.03, 0), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.2, 0.25)+
  draw_plot(x11,0, 0, 0.8, 0.25) +
  draw_plot(x12, .8, 0,0.2, 0.25)+
  draw_plot_label(c("**","***", "**","**", "***","**"), c(0.34,0.51, 0.34,0.52,0.33,0.52),c(0.92,0.92, 0.59,0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.03, 0.03, 0), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

 #lyrara

pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x3=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  ylab("pi of all sites")+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF>=max3)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

pop="Spiterstullen"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x7=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of all sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF>=max3)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

pop="N.America"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$Ho)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$Ho>=max)
unilateral_test_0.975=subset(mlm3,mlm3$Ho>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$Ho>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$Ho/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("Ho",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=Ho)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of all sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=12) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=Ho)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))

x=sort(mlm2$pi)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi>=max3)

effet=mlm3$pi/med
index=rep("pi",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x11=ggplot(data=mlm3, aes(x=as.factor(region), y=pi)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of all sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


x=sort(mlm2$MAF)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF>=max3)

effet=mlm3$MAF/med
index=rep("MAF",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("**", "*", "***","***"), c(0.34,0.54,0.33,0.51), c(0.92, 0.59, 0.26,0.26), size = c(20,20,20,20),fontface = c("italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.02), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.2, 0.25)+
  draw_plot(x11,0, 0, 0.8, 0.25) +
  draw_plot(x12, .8, 0,0.2, 0.25)+
  draw_plot_label(c("**","*","*", "*", "***","**"), c(0.34,0.54,0.35,0.54,0.33,0.52), c(0.92,0.92, 0.59, 0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.02), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

colnames(Summary)=c("Pop","index","control_max_(med)","S_flanking_regions","effect")
write.table(Summary, "global_polymorphism_comparison.csv", sep=";", quote= FALSE)

# evolution of proportion, pi, Ho, MAF of polymorphic sites in region of 25kb around S locus and in control region
###### test unilateral ######
####### use csv file "VCF_analysis_distri" of each pop generate by python pipeline ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########

Summary=NULL

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)

Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")
pop="Plech"

#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x2=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#Ho polymorphic site
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x3=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x5=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x6=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x19=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x20=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")
pop="Spiterstulen"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x7=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  ylab("prop. polymorphic")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x9=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  ylab("Ho mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x10=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x11=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("pi mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x21=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("MAF mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x22=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm
pop="N.America"
mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x13=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x14=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x15=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x16=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x17=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x18=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x23=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x24=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))



ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.3, 0.25)+
  draw_plot(x13,0, 0, 0.8, 0.25) +
  draw_plot(x14, .8, 0,0.3, 0.25)+
  draw_plot_label(c("***","**","**", "***", "***","**"), c(0.33,0.525,0.345,0.51,0.33,0.525), c(0.92,0.92, 0.59, 0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x9,0, 0.33, 0.8, 0.25) +
  draw_plot(x10, .8, 0.33,0.3, 0.25)+
  draw_plot(x15,0, 0, 0.8, 0.25) +
  draw_plot(x16, .8, 0,0.3, 0.25)+
  draw_plot_label(c("*"), c(0.535), c(0.26), size = c(20),fontface = c("italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x5,0, 0.66, 0.8, 0.25) +
  draw_plot(x6, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x11,0, 0.33, 0.8, 0.25) +
  draw_plot(x12, .8, 0.33,0.3, 0.25)+
  draw_plot(x17,0, 0, 0.8, 0.25) +
  draw_plot(x18, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.06, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x19,0, 0.66, 0.8, 0.25) +
  draw_plot(x20, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x21,0, 0.33, 0.8, 0.25) +
  draw_plot(x22, .8, 0.33,0.3, 0.25)+
  draw_plot(x23,0, 0, 0.8, 0.25) +
  draw_plot(x24, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


#halleri

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")
pop="Japan"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x19=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x20=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x21=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x22=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x23=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x24=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x1=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x2=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")
pop="Nivelle"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x25=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  ylab("prop. polymorphic")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x26=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x27=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  ylab("Ho mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x28=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x29=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("pi mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x30=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x3=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("MAF mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm
pop="Mortagne"
mlm2=subset(mlm,mlm$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x31=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x32=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x33=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x34=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x35=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x36=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x5=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x6=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))




ggdraw() +
  draw_plot(x19,0, 0.66, 0.8, 0.25) +
  draw_plot(x20, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x25,0, 0.33, 0.8, 0.25) +
  draw_plot(x26, .8, 0.33,0.3, 0.25)+
  draw_plot(x31,0, 0, 0.8, 0.25) +
  draw_plot(x32, .8, 0,0.3, 0.25)+
  draw_plot_label(c("***","**","***", "**", "***","**"), c(0.32,0.52,0.32,0.52,0.32,0.52), c(0.92,0.92, 0.59, 0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x21,0, 0.66, 0.8, 0.25) +
  draw_plot(x22, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x27,0, 0.33, 0.8, 0.25) +
  draw_plot(x28, .8, 0.33,0.3, 0.25)+
  draw_plot(x33,0, 0, 0.8, 0.25) +
  draw_plot(x34, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x23,0, 0.66, 0.8, 0.25) +
  draw_plot(x24, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x29,0, 0.33, 0.8, 0.25) +
  draw_plot(x30, .8, 0.33,0.3, 0.25)+
  draw_plot(x35,0, 0, 0.8, 0.25) +
  draw_plot(x36, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x3,0, 0.33, 0.8, 0.25) +
  draw_plot(x4, .8, 0.33,0.3, 0.25)+
  draw_plot(x5,0, 0, 0.8, 0.25) +
  draw_plot(x6, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


colnames(Summary)=c("Pop","index","control_max_(med)","S_flanking_regions","effect")
write.table(Summary, "global_polymorphism_comparison_polymorphic_site.csv", sep=";", quote= FALSE)


# evolution of pi, Ho, MAF in S flanking region with distance at the S locus
###### test GLM ######
####### use csv file "VCF_analysis_distri" of each pop generate by python pipeline ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########



pop="Japan"
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a


pop="Nivelle"
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a



pop="Mortagne"
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a



pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a



pop="Spiterstulen"
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a


pop="N.America"
Slocus=NULL
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
#region S locus
mlm3=subset(mlm,mlm$region=="25kb")
DIST=mlm3$dist+9376729
DIST=DIST-9376728
Slocus=cbind(mlm3,DIST)
mlm3=subset(mlm,mlm$region=="50kb")
DIST=mlm3$dist+9401730
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="75kb")
DIST=mlm3$dist+9426731
DIST=DIST-9376728
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-75kb")
DIST=mlm3$dist+9264458
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-50kb")
DIST=mlm3$dist+9289459
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)
mlm3=subset(mlm,mlm$region=="-25kb")
DIST=mlm3$dist+9314460
DIST=9339461-DIST
mlm3=cbind(mlm3,DIST)
Slocus=rbind(mlm3,Slocus)



a=summary(glm(Slocus$Ho~Slocus$DIST))
a
a=summary(glm(Slocus$pi~Slocus$DIST))
a
a=summary(glm(Slocus$MAF~Slocus$DIST))
a



############# study of divergence ##########

###### test unilateral ######
####### use csv file "comparaison_alignement_lyrata_thaliana"  ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########


mlm4=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm4,mlm4$divergence_pourcentage<=0.3)
mlm=subset(mlm,mlm$cover_pourcentage>=80)

mlm2=subset(mlm,mlm$type=="control")
mlm3=subset(mlm4,mlm4$type!="control")

divergence_mean_control=tapply(mlm2$divergence_pourcentage, mlm2$region, mean)
divergence_mean_S=tapply(mlm3$divergence_pourcentage, mlm3$region, mean)

region=c("-75kb","-50kb","-25kb","25kb","50kb","75kb")
mean_divergence=c(divergence_mean_S[3],divergence_mean_S[2],divergence_mean_S[1],divergence_mean_S[4],divergence_mean_S[5],divergence_mean_S[6])
mlm3=cbind(region,mean_divergence)
mlm3=as.data.frame(mlm3)

x=sort(divergence_mean_control)
max=x[0.975*length(x)]
min=x[0.025*length(x)]
med=median(x)
divergence_mean_control=as.data.frame(divergence_mean_control)

x1=ggplot(data=mlm3, aes(x=region, y=as.numeric(mean_divergence))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.25)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("mean divergence by gene")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "black"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))

x2=ggplot(divergence_mean_control, aes(x=divergence_mean_control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.25)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))

ggdraw() +
  draw_plot(x1,0, 0, 0.8, 1) +
  draw_plot(x2, .8, 0,0.2, 1)

############# study of ratio pi0F/pi4F ##########

###### test unilateral ######
####### use csv file "VCF_analysis_mean"  ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########


Summary=NULL
#halleri

pop="Japan"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(plot.title = element_text( size=25, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))



pop="Nivelle"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="Mortagne"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20,color="white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.02, 0.02, 0), c(1,0.99, 0.66, 0.34), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


#lyrara

pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(plot.title = element_text( size=25, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="Spiterstullen"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="N.America"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$X0f.4f)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X0f.4f>=max
mlm3$X0f.4f>=max2
mlm3$X0f.4f>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
effet=mlm3$X0f.4f/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("pi0F/pi4F",length(mlm3[,1]))
val=format(mlm3$X0f.4f,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=X0f.4f)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Pi0f/Pi4f")+
  theme(axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,color="white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=X0f.4f)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("*"), c(0.22), c(0.26), size = c(20),fontface = c("italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.05, 0, 0.02), c(1,0.99, 0.66, 0.34), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

write.table(Summary, "pi_ratio.csv", sep=";", quote= FALSE)


############# study of tajima D ##########

###### test bilateral ######
####### use csv file "gene_control_tajima_"  ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########


Summary=NULL
#halleri

pop="Japan"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)


x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)



x1=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x2=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="Nivelle"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)

x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)



x5=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab("TajimaD")+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "black"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x6=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))




pop="Mortagne"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)


x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)



x9=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "black"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x10=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20, colour = "black"))



ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.04, 0.02, 0), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


#lyrara

pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)

x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)



x1=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x2=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="Spiterstullen"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)


x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)



x5=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab("TajimaD")+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "black"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x6=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


pop="N.America"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


control=tapply(mlm2$Dtajima, mlm2$region, mean)
dtajima=tapply(mlm3$Dtajima, mlm3$region, mean)

Region=c("-25kb","-50kb","-75kb","25kb","50kb","75kb")
Slocus=cbind(Region,dtajima)
Slocus=as.data.frame(Slocus)
rm(Region)
rm(dtajima)

x=sort(as.numeric(control))
max=as.numeric(x[0.975*length(x)])
min=as.numeric(x[0.025*length(x)])
med=median(as.numeric(x))

Pop=rep(pop,length(Slocus[,1]))
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(Slocus[,1]))
effet=as.numeric(Slocus$dtajima)/med
effet=format(effet, scientific=TRUE, digits=3)
index=rep("tajima",length(Slocus[,1]))
val=format(Slocus$dtajima,scientific=TRUE, digits=3)
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x9=ggplot(data=Slocus, aes(x=Region, y=as.numeric(dtajima))) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=0,col="grey")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  geom_hline(yintercept=min,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(-3,4)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Dtajima")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "black"),axis.text.y = element_text(size=20,colour="black"))

control=as.data.frame(control)

x10=ggplot(control, aes(x=control)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(-3,4)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=0,col="grey")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  geom_vline(xintercept=min,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20, colour = "black"))




ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.05, 0, 0.02), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


colnames(Summary)=c("pop","indice","control","region","effet")
write.table(Summary, "Dtajima.csv", sep=";", quote= FALSE)


############# study of 0-fold polymorphism ##########

# evolution of pi, Ho, MAF in region of 25kb around S locus and in control region
###### test unilateral ######
####### use csv file "VCF_analysis_mean" of each pop generate by python pipeline ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########


Summary=NULL # synthesis table

#halleri

pop="Japan"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")



x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Ho")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi0fold>=max3)

effet=mlm3$pi0fold/med
index=rep("pi0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x3=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  ylab("pi")+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


pop="Nivelle"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of 0-fold degenerate sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi0fold>=max3)

effet=mlm3$pi0fold/med
index=rep("pi0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x7=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of 0-fold degenerate sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


pop="Mortagne"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=20) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))

x=sort(mlm2$pi0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi0fold>=max3)

effet=mlm3$pi0fold/med
index=rep("pi0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x11=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c( "*", "*","**", "*"), c( 0.36,0.54,0.345,0.54), c( 0.59,0.59, 0.26,0.26), size = c(20,20,20,20),fontface = c("italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.03, 0.03, 0), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.2, 0.25)+
  draw_plot(x11,0, 0, 0.8, 0.25) +
  draw_plot(x12, .8, 0,0.2, 0.25)+
  draw_plot_label(c("*", "*", "*","**", "*","*"), c(0.54,0.36,0.54,0.345,0.54,0.27), c(0.92, 0.59,0.59, 0.26,0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.03, 0.03, 0), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

#lyrara

pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylab("Ho")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x2=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x3=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  ylab("pi")+
  #ggtitle(paste("A.halleri \n",pop)) +
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))


x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

pop="Spiterstullen"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x5=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho of 0-fold degenerate sites")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=0) ,axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))


x6=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$pi0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi0fold>=max3)

effet=mlm3$pi0fold/med
index=rep("pi0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x7=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi of 0-fold degenerate sites")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=0, colour = "white"))

x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

pop="N.America"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


x=sort(mlm2$ho_of)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$ho_of>=max)
unilateral_test_0.975=subset(mlm3,mlm3$ho_of>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$ho_of>=max3)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$ho_of/med
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$ho_of,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
index=rep("ho_of",length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x9=ggplot(data=mlm3, aes(x=as.factor(region), y=ho_of)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0,0.0128)+
  #ggtitle(pop) +
  ylab("Ho")+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x =element_text(size=12) ,axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))


x10=ggplot(mlm2, aes(x=ho_of)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0128)+
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))

x=sort(mlm2$pi0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$pi0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$pi0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$pi0fold>=max3)

effet=mlm3$pi0fold/med
index=rep("pi0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x11=ggplot(data=mlm3, aes(x=as.factor(region), y=pi0fold)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  xlab(NULL)+
  theme_classic()+
  #ggtitle(pop) +
  ylab("pi ")+
  ylim(0,0.0144)+
  theme(plot.title = element_text( size=20, hjust=0.5),axis.title.x = element_text(size=0),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi0fold)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.0144)+
  ylim(0,15)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.y = element_text(size=0, colour = "white"),axis.text.x = element_text(size=20,colour="black"))


x=sort(mlm2$MAF_0fold)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
med=median(x)

unilateral_test_0.95=subset(mlm3,mlm3$MAF_0fold>=max)
unilateral_test_0.975=subset(mlm3,mlm3$MAF_0fold>=max2)
unilateral_test_0.99=subset(mlm3,mlm3$MAF_0fold>=max3)

effet=mlm3$MAF_0fold/med
index=rep("MAF_0fold",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_0fold,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x5,0, 0.33, 0.8, 0.25) +
  draw_plot(x6, .8, 0.33,0.2, 0.25)+
  draw_plot(x9,0, 0, 0.8, 0.25) +
  draw_plot(x10, .8, 0,0.2, 0.25)+
  draw_plot_label(c( "*", "**","**", "**"), c( 0.36,0.53,0.62,0.53), c(0.92, 0.59,0.59, 0.26), size = c(20, 20,20,20),fontface = c( "italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.02), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.2, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.2, 0.25)+
  draw_plot(x11,0, 0, 0.8, 0.25) +
  draw_plot(x12, .8, 0,0.2, 0.25)+
  draw_plot_label(c( "**", "*","*","*", "*", "*"), c( 0.34,0.27,0.54,0.63,0.36,0.54), c(0.92,0.92, 0.59,0.59, 0.26, 0.26), size = c(20, 20,20,20,20,20),fontface = c( "italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.02), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

colnames(Summary)=c("Pop","index","control_max_(med)","S_flanking_regions","effect")
write.table(Summary, "0fold_polymorphism_comparison.csv", sep=";", quote= FALSE)

# evolution of proportion, pi, Ho, MAF of polymorphic sites in region of 25kb around S locus and in control region
###### test unilateral ######
####### use csv file "VCF_analysis_distri" of each pop generate by python pipeline ###########
####### one analysis by population, one synthetic figure by analysis*specie ###########

Summary=NULL

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")


Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

#region control
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

pop="Plech"

#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x1=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x2=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#Ho polymorphic site
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x3=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x5=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x6=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x19=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x20=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm

#region control
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

pop="Spiterstulen"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x7=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  ylab("prop. polymorphic")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x8=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x9=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  ylab("Ho mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x10=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x11=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("pi mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x12=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x21=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("MAF mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x22=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=mlm
pop="N.America"
#region control
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x13=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x14=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x15=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x16=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x17=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x18=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x23=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x24=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))



ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x7,0, 0.33, 0.8, 0.25) +
  draw_plot(x8, .8, 0.33,0.3, 0.25)+
  draw_plot(x13,0, 0, 0.8, 0.25) +
  draw_plot(x14, .8, 0,0.3, 0.25)+
  draw_plot_label(c( "**", "*","***","***", "***", "**", "*"), c( 0.33,0.26,0.23,0.5,0.6, 0.34,0.53), c(0.92,0.92, 0.59,0.59,0.59, 0.26, 0.26), size = c(20, 20,20,20,20,20,20),fontface = c( "italic","italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x3,0, 0.66, 0.8, 0.25) +
  draw_plot(x4, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x9,0, 0.33, 0.8, 0.25) +
  draw_plot(x10, .8, 0.33,0.3, 0.25)+
  draw_plot(x15,0, 0, 0.8, 0.25) +
  draw_plot(x16, .8, 0,0.3, 0.25)+
  draw_plot_label(c("*"), c(0.535), c(0.26), size = c(20),fontface = c("italic")) +
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x5,0, 0.66, 0.8, 0.25) +
  draw_plot(x6, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x11,0, 0.33, 0.8, 0.25) +
  draw_plot(x12, .8, 0.33,0.3, 0.25)+
  draw_plot(x17,0, 0, 0.8, 0.25) +
  draw_plot(x18, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.06, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x19,0, 0.66, 0.8, 0.25) +
  draw_plot(x20, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x21,0, 0.33, 0.8, 0.25) +
  draw_plot(x22, .8, 0.33,0.3, 0.25)+
  draw_plot(x23,0, 0, 0.8, 0.25) +
  draw_plot(x24, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.lyrata", "Plech", "Spiterstulen","N.America"), c(0.4, 0.07, 0, 0.03), c(1,0.97, 0.64, 0.31),  size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


#halleri

mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

pop="Japan"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x19=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x20=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x21=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x22=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x23=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x24=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x1=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x2=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(plot.title = element_text( size=12),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")

pop="Nivelle"
#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x25=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  ylab("prop. polymorphic")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x26=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0,colour = "white"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x27=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  ylab("Ho mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x28=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x29=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("pi mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x30=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x3=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  ylab("MAF mean pol")+
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=20),axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x4=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=0, colour = "white"),axis.text.y = element_text(size=0, colour = "white"))


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=subset(mlm,mlm$degenerate=="0fold")
Region=as.factor(mlm$region)
Region=levels(Region)
Table=NULL

for(i in 1:length(Region))
{
  mlm2=subset(mlm,mlm$region==Region[i])
  nb_considered=length(mlm2[,1])
  mlm2=subset(mlm2,mlm2$pi!=0)
  nb_pol=length(mlm2[,1])/nb_considered
  Ho_pol=mean(mlm2$Ho)
  MAF_pol=mean(mlm2$MAF)
  pi_pol=mean(mlm2$pi)
  x=cbind(Region[i],nb_pol,Ho_pol,MAF_pol,pi_pol)
  Table=rbind(Table,x)
}

colnames(Table)=c("region","%_pol","Ho_pol","MAF_pol","pi_pol")
write.table(Table, "evolution_pol_pos.csv", sep=";", quote= FALSE)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
mlm2=subset(mlm,mlm$region!="-25kb")
mlm2=subset(mlm2,mlm2$region!="-50kb")
mlm2=subset(mlm2,mlm2$region!="-75kb")
mlm2=subset(mlm2,mlm2$region!="25kb")
mlm2=subset(mlm2,mlm2$region!="50kb")
mlm2=subset(mlm2,mlm2$region!="75kb")
#region S locus
mlm3=subset(mlm,mlm$region=="-75kb"|mlm$region=="-50kb"|mlm$region=="-25kb"|mlm$region=="25kb"|mlm$region=="50kb"|mlm$region=="75kb")


#proportion of site
x=sort(mlm2$X._pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$X._pol>=max
mlm3$X._pol>=max2
mlm3$X._pol>=max3
med=median(x)

Pop=rep(pop,length(mlm3[,1]))
effet=mlm3$X._pol/med
index=rep("prop_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$X._pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)

x31=ggplot(data=mlm3, aes(x=region, y=X._pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.06)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25,colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x32=ggplot(mlm2, aes(x=X._pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.06)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#Ho polymorphic
x=sort(mlm2$Ho_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$Ho_pol>=max
mlm3$Ho_pol>=max2
mlm3$Ho_pol>=max3
med=median(x)

effet=mlm3$Ho_pol/med
index=rep("Ho_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$Ho_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x33=ggplot(data=mlm3, aes(x=region, y=Ho_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.5)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x34=ggplot(mlm2, aes(x=Ho_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.5)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0, colour = "white"))

#pi polymorphic
x=sort(mlm2$pi_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$pi_pol>=max
mlm3$pi_pol>=max2
mlm3$pi_pol>=max3
med=median(x)

effet=mlm3$pi_pol/med
index=rep("pi_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$pi_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)




x35=ggplot(data=mlm3, aes(x=region, y=pi_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x36=ggplot(mlm2, aes(x=pi_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+ 
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))

#MAF polymorphic
x=sort(mlm2$MAF_pol)
max=x[0.95*length(x)]
max2=x[0.975*length(x)]
max3=x[0.99*length(x)]
mlm3$region
mlm3$MAF_pol>=max
mlm3$MAF_pol>=max2
mlm3$MAF_pol>=max3
med=median(x)

effet=mlm3$MAF_pol/med
index=rep("MAF_site_pol",length(mlm3[,1]))
effet=format(effet, scientific=TRUE, digits=3)
val=format(mlm3$MAF_pol,scientific=TRUE, digits=3)
MED=format(med, scientific=TRUE, digits=3)
MAX=format(max, scientific=TRUE, digits=3)
med2=rep(paste(MAX,"(",MED,")"),length(mlm3[,1]))
summ=cbind(Pop,index,med2,val,effet)
Summary=rbind(Summary,summ)


x5=ggplot(data=mlm3, aes(x=region, y=MAF_pol)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limits=c("-75kb","-50kb","-25kb","Slocus","25kb","50kb","75kb"))+
  geom_hline(yintercept=med,col="black")+
  geom_hline(yintercept=max,col="black",linetype="dashed")+
  theme_classic()+
  xlab(NULL)+
  ylim(0,0.45)+ 
  theme(axis.title.x = element_text(size=0, colour = "white"),axis.title.y = element_text(size=25, colour = "white"),axis.text.x = element_text(size=20, colour = "white"),axis.text.y = element_text(size=20,colour="black"))

x6=ggplot(mlm2, aes(x=MAF_pol)) + 
  geom_histogram(fill="white", col="black")+
  coord_flip()+
  xlim(0,0.45)+ 
  ylim(0,20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  geom_vline(xintercept=med,col="black")+
  geom_vline(xintercept=max,col="black",linetype="dashed")+
  theme(axis.text.x = element_text(size=20,colour="black"),axis.text.y = element_text(size=0,colour = "white"))




ggdraw() +
  draw_plot(x19,0, 0.66, 0.8, 0.25) +
  draw_plot(x20, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x25,0, 0.33, 0.8, 0.25) +
  draw_plot(x26, .8, 0.33,0.3, 0.25)+
  draw_plot(x31,0, 0, 0.8, 0.25) +
  draw_plot(x32, .8, 0,0.3, 0.25)+
  draw_plot_label(c("**", "**", "**","*","***", "*"), c(0.33,0.52,0.33,0.53,0.32,0.53), c(0.92,0.92, 0.59,0.59, 0.26,0.26), size = c(20,20,20,20,20,20),fontface = c("italic","italic","italic","italic","italic","italic")) +
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x21,0, 0.66, 0.8, 0.25) +
  draw_plot(x22, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x27,0, 0.33, 0.8, 0.25) +
  draw_plot(x28, .8, 0.33,0.3, 0.25)+
  draw_plot(x33,0, 0, 0.8, 0.25) +
  draw_plot(x34, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x23,0, 0.66, 0.8, 0.25) +
  draw_plot(x24, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x29,0, 0.33, 0.8, 0.25) +
  draw_plot(x30, .8, 0.33,0.3, 0.25)+
  draw_plot(x35,0, 0, 0.8, 0.25) +
  draw_plot(x36, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 

ggdraw() +
  draw_plot(x1,0, 0.66, 0.8, 0.25) +
  draw_plot(x2, .8, 0.66,0.3, 0.25)+ 
  draw_plot(x3,0, 0.33, 0.8, 0.25) +
  draw_plot(x4, .8, 0.33,0.3, 0.25)+
  draw_plot(x5,0, 0, 0.8, 0.25) +
  draw_plot(x6, .8, 0,0.3, 0.25)+
  draw_plot_label(c("A.halleri", "Japan", "Nivelle","Mortagne"), c(0.4, 0.07, 0.05, 0.03), c(1,0.97, 0.64, 0.31), size = c(30,25,25,25),fontface = c("italic","bold","bold","bold")) 


colnames(Summary)=c("Pop","index","control_max_(med)","S_flanking_regions","effect")
write.table(Summary, "0fold_polymorphism_comparison_polymorphic_site.csv", sep=";", quote= FALSE)


############# study of 0-fold/del mutations in heterozygous state in individual scale ##########


###### GLM ######
####### use csv file "vcf_analysis_ind_" of each pop generate by python pipeline ###########
####### one analysis by population, one synthetic figure by analysis ###########

Table=NULL
Summary=NULL

pop="Japan"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)

a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
pop="Nivelle"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)
a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
pop="Mortagne"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)
a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
pop="Plech"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)
a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
pop="Spiterstulen"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)
a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
pop="N.America"
mlm=read.csv(file.choose(), sep=";",dec=".",h=T)
Pop=rep(pop,length(mlm[,1]))
mlm=cbind(mlm,pop)
fold0=mlm$nb_mut_0fold_het/((mlm$nb_0fold/2))
del=mlm$nb_del_SIFT_het/((mlm$nb_cds))
mlm=cbind(mlm,fold0)
mlm=cbind(mlm,del)
a=summary(glm(mlm$fold0~mlm$dist))
a=a$coefficients
a
x=cbind("0fold_het",a[2],a[8],pop)
Summary=rbind(Summary,x)
a=summary(glm(mlm$del~mlm$dist))
a=a$coefficients
a
x=cbind("del_het_prop_mut",a[2],a[8],pop)
Summary=rbind(Summary,x)
Table=rbind(Table,mlm)
write.table(Table, "analyse_ind.csv", sep=";", quote= FALSE)
write.table(Summary, "Summary_analyse_ind.csv", sep=";", quote= FALSE)


mlm=read.csv(file.choose(), sep=";",dec=".",h=T)



ggplot(Table, aes (dist,fold0, color=pop))+
  geom_point()+
  geom_smooth(method=glm) +
  theme_classic()+xlab("Distance to S locus (b)")+ylab("Proportion of heterozygous 0-fold mutations")

ggplot(Table, aes (dist,del, color=pop))+
  geom_point()+
  geom_smooth(method=glm) +
  theme_classic()+xlab("Distance to S locus (b)")+ylab("Proportion of heterozygous del mutations")
