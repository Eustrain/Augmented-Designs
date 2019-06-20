
##***********Augmented Desing************###

library(plyr)


dat <- read.csv("data.csv",head=T)
head(dat)
str(dat)

dat$block <- as.factor(dat$block)
dat$Trat <- as.factor(dat$Trat)

Aug_dis(dat$block,dat$Trat,dat$y,dat)

Aug_dis<- function(block,Trat,y,dat) {
  

####################### analysis of variance 

fm1 <-lm(y~Trat+block,data = dat)
av <- anova(fm1)
 
#####################################################
###Calculate the  frequency for the varieties and checks
trats<- data.frame(table(dat$Trat))
vari = as.vector (trats$Var1[trats$Freq == length(3)])
chec = as.vector (trats$Var1[trats$Freq ==3])
######################################################
variedades <- filter(dat, Trat %in% vari)
checks <- filter(dat, Trat %in% chec)
##################################################

cdata <- ddply(checks, c("Trat"),summarise,
su= sum(y))
######################################################
nblock <- length(levels(dat$block))

######################################################
#####  Sums sq and Mean sq for checks
c <- length(cdata$Trat)
cFC <- sum(checks$y)^2/(nblock*c)         
css <- sum(cdata$su^2)/nblock-cFC       
cms <- css/(c-1) ####### 
#######################
######Sums sq and Mean sq for  varieties
vcf <- sum(variedades$y)^2/ length(variedades$y)
vss <- sum(variedades$y^2)-vcf
vms <- vss/(length(variedades$y)-1)

#############################################
#############################################
######Sums sq and Mean sq for  varieties  vs checks
cvss= av$`Sum Sq`[1]-css-vss
cvms <- cvss/1
############################################
############################################
###### F calculated for  varieties, checks and   varieties  vs checks
fcc <-cms/av$`Mean Sq`[3]
fcv <- vms/av$`Mean Sq`[3]
fccv <- cvms/av$`Mean Sq`[3]
############################################
############################################
fv <- c("Test","Checks","TratvsChec")
ss <- c(vss, css,cvss)
ms <- c(vms,cms,cvms)
f <- c(fcv,fcc,fccv)
df<- c((length(variedades$y)-1),(c-1),1)
#####################################################
fval <-NULL

for (i in 1:3) {
  
  
  pr <- pf(f[i],df[i],av$Df[3], lower.tail=F)#####
  
  fval[[i]] <-pr 
  
}
########################################
##########ANOVA
anv <- data.frame("."=fv,"Df"=df,"Sum sq"=ss,"Mean sq"=ms
                  ,"F value"=f,"Pr(>F)"=fval)


########################################################
####Standard Error of mean difference#######

#################between 2 checks
SEC <- sqrt(2 * av$`Mean Sq`[3]/nblock)    
#####SEC <-( sqrt(2* av$`Mean Sq`[3]))/sqrt((nblock)) 
########################between 2 vars
SEV <-sqrt( (2* av$`Mean Sq`[3]))
########SEV <- sqrt( (2* av$`Mean Sq`[3]))

################between between any tw~ entries of the same block
SE2 <- sqrt(av$`Mean Sq`[3]*1+(1/c))

#############between means of a chec+ and a test variety)
SE3 <- sqrt(av$`Mean Sq`[3]*(1+1/nblock+1/c+1/(nblock*c)))

##################between a vars
Err_st <-c(SEC,SEV,SE2,SE3 )
er <- c("between 2 checks means","between 2 Test"," two entries of the same block","Test Treatment and a Control Treatment")
DF <- data.frame("."=er,"Std. Error of Diff"=Err_st)
###########################
result_anv <- list(anv,av)
result_est <- list(DF)
re<- list(result_anv,result_est)
print(re)


}

#blk <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3)
#y1 <- c(92, 79, 87, 81, 96, 89, 82, 79, 81, 81, 91, 79, 78, 83, 77, 78, 78,
#trt <- c(1, 2, 3, 4, 7, 11, 12, 1, 2, 3, 4, 5, 9, 1, 2, 3, 4, 8, 6, 10)
#70, 75, 74)
#y2 <- c(258, 224, 238, 278, 347, 300, 289, 260, 220, 237, 227, 281, 311, 250,
#240, 268, 287, 226, 395, 450)

#####Example 2
dat<- data.frame(blk,trt,y1)
str(dat)
names(dat) <- c("block","Trat","y")
dat$block <- as.factor(dat$block)
dat$Trat <- as.factor(dat$Trat)
Aug_dis(dat$block,dat$Trat,dat$y,dat)
######Sharma, J.R. (1998) Statistical and Biometrical Techniques in Plant Breeding
###### Aravind, J., Mukesh Sankar, S., Wankhede, D. P., and Kaur, V.
#(2019).  augmentedRCBD: Analysis of Augmented Randomised
#Complete Block Designs. R package version 0.1.0.9000,
#https://aravind-j.github.io/augmentedRCBD/https://cran.r-project.org/package=augmentedRCBD.






