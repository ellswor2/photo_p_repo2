# Clean up workspace
rm(list = ls())
# load library: need to be done each time you open R Studio
library (devtools)
library(plantecophys)
library(lubridate)
library(doBy)
library(plyr)
library(scales)
library(lme4)

# read in the data
mydata <- read.csv("C:/Manuscripts21-25/Ells_P-in-photosynthesis2021_REV/Ellsworth_Photosynthesis_N_P/00_raw_data/Ellsworth_NCOMMS_Figure3_fulldata.csv")
View(mydata)
#set working directory so all files go in same place
setwd("C:/Manuscripts21-25/Ells_P-in-photosynthesis2021_REV/Ellsworth_Photosynthesis_N_P/Graphs_n_Figures")

ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from,to,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}

Acibothfits<-mydata
P_status_set <- c("plum", "purple4")   # subst "orchid" for "skyblue2"
Acibothfits$Plim_status<-as.factor(Acibothfits$Plim_status)
Acibothfits2<-Acibothfits
# READY
Acibothfits2$Continent <-factor(Acibothfits2$Continent)

# We also need to concatenate Continent with P_status in order to create the correct legend
Acibothfits2$Pstat.Cont <- as.factor(paste(Acibothfits2$Plim_status, Acibothfits2$Continent, sep = "_"))

# set up colours
orchblkfun <- colorRampPalette(c("purple4","plum"))
col2<-palette(orchblkfun(4))

#FINAL FIGURE 3
# Cut the dataset into 4 equal N/P classes # log makes it a little more even
Acibothfits2$PPbin<-cut(log(Acibothfits2$Pmass), breaks=4, include.lowest=TRUE, labels=c("Low","Med1", "Med2", "High"), na.rm=TRUE)
# want really to cut the dataset into 3 equal-population factors so do this
Acibothfits2$PPbin<-cut(Acibothfits2$Pmass,quantile(Acibothfits2$Pmass,(0:4)/4))
summary(Acibothfits2$PPbin)   # here the bins are (27.4,105]  (105,141]  (141,459]
# Rename levels of a  
levels(Acibothfits2$PPbin) <- c("low","med1","med2","high")

Acibothfits3 <- subset(Acibothfits2, PPbin == "low" | PPbin == "high")
Acibothfits3.low <- subset(Acibothfits2, PPbin == "low")
Acibothfits3.hi <- subset(Acibothfits2, PPbin == "high")

mean(Acibothfits3.low$Pmass)
mean(Acibothfits3.hi$Pmass)

# Ok so let's test if there is a slope difference in Jmax-Vcmax with NPratio
# is there an interaction between NPratio class and the slopes?
Regr.Vc_Jm.Pmass <- lm(Jmax~Vcmax + Pmass, data=Acibothfits3)
summary(Regr.Vc_Jm.Pmass)  # so Pmass is significant in this model, P = 0.00136

Regr.Vc_Jm.twoP <- lm(Jmax~Vcmax * PPbin, data=Acibothfits3)
summary(Regr.Vc_Jm.twoP)
Regr.Vc_Jm.twoP2 <- lm(Jmax~Vcmax + PPbin, data=Acibothfits3)
summary(Regr.Vc_Jm.twoP2)
Regr.Vc_Jm.2ponly <- lm(Jmax~Vcmax, data=Acibothfits3)
summary(Regr.Vc_Jm.2ponly)

anova(Regr.Vc_Jm.2ponly, Regr.Vc_Jm.twoP) #P = 0.0007231, compares slopes and intercepts vs no slope or intercept
anova(Regr.Vc_Jm.twoP2, Regr.Vc_Jm.twoP) # P = 0.035, compares slopes & same intercept

Regr.Vc_Jm.lowP<- lm(Jmax~Vcmax, data=Acibothfits3.low)
Regr.Vc_Jm.hiP<- lm(Jmax~Vcmax, data=Acibothfits3.hi)
confint(Regr.Vc_Jm.lowP)
confint(Regr.Vc_Jm.hiP)
#############
# Make Figure3
#############
col3<-palette(orchblkfun(4))

#par(pty="s")
#png(filename="FINAL_Fig3_Jmax-Vcmax_regression_by_Pconct.png", width=600, height=600, bg = "white")
par(mar=c(5,6,5,1)+.1)
plot(Jmax ~ Vcmax, data=Acibothfits2, xaxt="none", ylim=c(0, 250), pch=19, cex=1.7, cex.lab=1.8, 
     cex.axis=1.6, xlab=expression(italic(V)[cmax]~(mu~mol~ m^{-2}~s^{-1})), ylab=expression(italic(J)[max]~(mu~mol~ m^{-2}~s^{-1})), 
     col=alpha(col3, 0.7))
axis(1, seq(0,125,25), cex.axis=1.6)
# use the ablinepiece function
ablinepiece(Regr.Vc_Jm.lowP, from=min(Acibothfits3$Vcmax), to =max(Acibothfits3$Vcmax), lwd = 3, col = "purple4")
ablinepiece(Regr.Vc_Jm.hiP, from=min(Acibothfits3$Vcmax), to=max(Acibothfits3$Vcmax), lwd = 3, col = "plum")

my.expressions<-expression(paste(P[mass]," of ","1.78 mg ",g^-1), paste(P[mass]," of ","0.44 mg ",g^-1, sep=""))
#my.expressions<-c(expression(paste(P[mass], "1.8 mg/g"), paste(P[mass], "0.45 mg/g")))

legend("topleft", bty='n',legend=my.expressions,
       col = c("plum", "purple4"),
       lwd = 3, cex = 1.7)  # increase legend size

dev.off()

# END

