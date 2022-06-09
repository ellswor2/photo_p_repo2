# This makes Figure 1 Ellsworth et al. NCOMMS manuscript
# and associated analyses (Tables 1 and 2)
# v. 07.06.2022

# Clean up workspace
rm(list = ls())
#.libPaths("C:/R-cran/R-4.1.1/library")

# load library
library (devtools)
library(plantecophys)
library(lubridate)
library(doBy)
library(plyr)
library(scales)
# read in data and view to check
mydata <- read.csv("C:/Manuscripts21-25/Ells_P-in-photosynthesis2021_REV/Ellsworth_Photosynthesis_N_P/00_raw_data/Ellsworth_NCOMMS_Figure1and2_fulldata.csv")
View(mydata)
#set working directory so all files go in same place
setwd("C:/Manuscripts21-25/Ells_P-in-photosynthesis2021_REV/Ellsworth_Photosynthesis_N_P/Graphs_n_Figures")
# Check that massbased gas exchange is correct by recomputing these
mydata$Photomass<-mydata$Photo/mydata$LMA*1000
mydata$Jmaxmass<-mydata$Jmax/mydata$LMA*1000
mydata$Vcmaxmass<-mydata$Vcmax/mydata$LMA*1000

mydata$Nmass<-mydata$Perc_N*10
mydata$Pmass<-mydata$Perc_P*10 #this converts %P to mg/g P

Acibothfits<-mydata
Acibothfits$NPratio<-Acibothfits$Perc_N/Acibothfits$Perc_P
Acibothfits<-Acibothfits[complete.cases(Acibothfits[ , c(21,22,26)]),]

# Set Plim_status and then subset
# I previously used 0.105 as an arbitrary cutoff but 0.0915 is the median so gives equal sample sizes
x <- c("Low P", "Mod. P")
Acibothfits$Plim_status <- ifelse(Acibothfits$Perc_P< 0.0915, x[1], ifelse(Acibothfits$Perc_P >= 0.0915, x[2], NA))
Acinew.Pnonlim<-subset(Acibothfits, Acibothfits$Plim_status == "Mod. P")
Acinew.Plim<-subset(Acibothfits, Acibothfits$Plim_status == "Low P")

# Different colours for our P_status classses
P_status_set <- c("orchid", "black")   # subst "orchid" for "skyblue2"
Acibothfits$Plim_status<-as.factor(Acibothfits$Plim_status)
Acibothfits2<-Acibothfits
# READY to set up Figure 1
Acibothfits2$Continent <-factor(Acibothfits2$Continent)
# We also need to concatenate Continent with P_status in order to create the correct legend
Acibothfits2$Pstat.Cont <- as.factor(paste(Acibothfits2$Plim_status, Acibothfits2$Continent, sep = "_"))
# so if I'm correct, then we would expect pch15 to correspond with Africa (filled square), pch=25 (filled downward triangle) with Asia, pch = 17 (Filled uptriangle) with Australia and pch=19 (filled circle)with S_Amer

### Main relationship area-based
Regr.Vcmaxmass_NxP.overall <- lm(log(Vcmaxmass)~log(Nmass)+log(Pmass)+log(Nmass)*log(Pmass), data=Acibothfits2, na.action=na.exclude)
Regr.Vcmax_NaxPa.overall <- lm(log(Vcmax)~log(Narea)+log(Parea)+log(Narea)*log(Parea), data=Acibothfits2, na.action=na.exclude)
Regr.Jmaxmass_NxP.overall <- lm(log(Jmaxmass)~log(Nmass)+log(Pmass)+log(Nmass)*log(Pmass), data=Acibothfits2, na.action=na.exclude)
Regr.Jmax_NaxPa.overall <- lm(log(Jmax)~log(Narea)+log(Parea)+log(Narea)*log(Parea), data=Acibothfits2, na.action=na.exclude)

### FIG_1 GRAPH mass-based Anet
Amass_range <-c(5, 500)
Nma_range<-c(2, 55)   # c(3e+06, 6e+08)
Vcmass_range <-c(50, 2000)
Jmass_range <-c(100, 3500)

png(file="Ellsworth_Figure1_Vcmass-logN_Jmass-N_Pcontrast_4panel.png",  width=600, height=550)

op1 <-par(mfrow=c(2,2), mar=c(5,5,4,1), cex.lab=1.5)

# here want to define symbols 
pchs <- c(15,25,17,19,15,25,17,19)  #pch=6 becomes pch=25 (closed)
allcols <- c("orchid","orchid","orchid","orchid","black","black","black","black")

# Create 'Acinew.Pnonlim2' that goes together with Acibothfits2 
Acinew.Pnonlim2 <-subset(Acibothfits2, Acibothfits2$Plim_status == "Mod. P")
Acinew.Plim2 <-subset(Acibothfits2, Acibothfits2$Plim_status == "Low P") #check if this is correct!!

# Then you just need to add a 2 to the data= below
logmod.Amass_N <- lm(log(Photomass)~log(Nmass), data=Acinew.Pnonlim2)
logypredAmass <- predict(logmod.Amass_N)
logmod.Amass_Npl <- lm(log(Photomass)~log(Nmass), data=Acinew.Plim2)
logypredAmasspl <- predict(logmod.Amass_Npl)  

#Plotting PANEL A of Figure 1 using Vcmax
logmod.Vcmaxmass_N2 <- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Pnonlim2) #removed , na.action=na.exclude
logypredVcmaxmass_N2 <- predict(logmod.Vcmaxmass_N2)
regr.Vcmaxmass_N2pnonlim <- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Pnonlim2, na.action=na.exclude)

newx <- seq(min(Acinew.Plim2$Nmass), max(Acinew.Plim2$Nmass), length.out=100)
newx2 <- seq(min(Acinew.Pnonlim2$Nmass), max(Acinew.Pnonlim2$Nmass), length.out=100)

plot(Vcmaxmass ~ Nmass, data=Acibothfits2, log = "xy", pch=c(15,25,17,19)[Continent], cex=1.1,
     col=alpha(P_status_set[Acibothfits$Plim_status], 0.4), cex.lab=1.6, cex.axis=1.4, bty="l", 
     xaxs="i",  yaxs="i",# = list(type = "log"), yaxis = list(type = "log"),
     ylab=expression(italic(V)[cmax_mass]~(nmol~ g^{-1}~s^{-1})),
     xlab=expression("Leaf"~N[mass]~(mg~ g^{-1})),
     xlim=Nma_range, ylim=Vcmass_range)  

logmod.Vcmaxmass_N2pl <- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Plim2) # , na.action=na.exclude
regr.Vcmaxmass_N2pl <- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Plim2, na.action=na.exclude) 

logypredVcmaxmassN2pl <- predict(logmod.Vcmaxmass_N2pl)
# full model to test for differences
regr.Vcmass.Nm.all <- lm(log(Vcmaxmass) ~ log(Nmass)*Plim_status, data=Acibothfits2, na.action=na.exclude) 

# diagnostics plots
#plot(regr.Jm.Pa.nonlim) #this produces 4 plots; can do par(mfrow = c(2,2)) prior to the plot cmd
summary(logmod.Vcmaxmass_N2pl)
#done this part, now this one non-Plim #FIX FIX FIX
regr.Vcmaxmass_Nnlim <- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Pnonlim2, na.action=na.exclude) 

lines(exp(logypredVcmaxmass_N2)~exp(logmod.Vcmaxmass_N2$model$`log(Nmass)`), col="black", lty=5, lwd=3)

# now the P-nonlim confidence interval
preds4 <- predict.lm(regr.Vcmaxmass_Nnlim, newdata = data.frame(Nmass=newx2), interval = "confidence")
predna4 <- predict.lm(regr.Vcmaxmass_Nnlim, newdata = data.frame(Nmass=newx2))

# add transparent fill for the 95% CI
greytrans <- rgb(190, 190, 190, 187, maxColorValue=255)
plumtrans1 <- rgb(221, 160, 221, 187, maxColorValue=255)  # this is plum colour, 127 for half-transparent
violtrans <- rgb(104, 34, 139, 107, maxColorValue=255)

polygon(c(rev(newx2), newx2), c(rev(exp(preds4[ ,3])), exp(preds4[ ,2])), col = greytrans, border = NA)

# this one Plim
lines(exp(logypredVcmaxmassN2pl)~exp(logmod.Vcmaxmass_N2pl$model$`log(Nmass)`), col="orchid", lty=5, lwd=3)
# now predicts + conf interval
preds3 <- predict.lm(regr.Vcmaxmass_N2pl, newdata = data.frame(Nmass=newx), interval = "confidence")
predna3 <- predict.lm(regr.Vcmaxmass_N2pl, newdata = data.frame(Nmass=newx))
# add transparent fill for the 95% CI
polygon(c(rev(newx), newx), c(rev(exp(preds3[ ,3])), exp(preds3[ ,2])), col = plumtrans1, border = NA)

#overlay lines again
lines(exp(logypredVcmaxmassN2pl)~exp(logmod.Vcmaxmass_N2pl$model$`log(Nmass)`), col="orchid", lty=5, lwd=3)
lines(exp(logypredVcmaxmass_N2)~exp(logmod.Vcmaxmass_N2$model$`log(Nmass)`), col="black", lty=5, lwd=3)

legend(x="topleft", pch=c(15,17), legend=c("Low P", "Moderate P"), col=c("orchid","black"),
       cex=1.0, bty="n", lwd=1.5)

mtext("(a)", 3, line=0.5, adj=1, cex=1.4)

# third panel
logmod.Jmaxmass_N2 <- lm(log(Jmaxmass)~log(Nmass), data=Acinew.Pnonlim2, na.action=na.exclude)
logypredJmaxmassN2 <- predict(logmod.Jmaxmass_N2)
# So both of the regressions below produce NaNs which are a problem
# Adding a minimum quantity to the argument of the logarithm function resolved the issue: log(response_variable + 0.0001234)
regr.Jmaxmass_Nnonlim <- lm(log(Jmaxmass)~log(Nmass), data=Acinew.Pnonlim2, na.action=na.exclude) 
regr.Jmaxmass.Nm.all <- lm(log(Jmaxmass) ~ log(Nmass)*Plim_status, data=Acibothfits2, na.action=na.exclude) 

plot(Jmaxmass ~ Nmass, data=Acibothfits2, log = "xy", pch=c(15,6,17,19)[Continent], cex=1.1,
     col=alpha(P_status_set[Acibothfits$Plim_status], 0.3), cex.lab=1.6, cex.axis=1.4, bty="l", 
     xaxs="i",  yaxs="i",# = list(type = "log"), yaxis = list(type = "log"),
     ylab=expression(italic(J)[max_mass]~(nmol~ g^{-1}~s^{-1})),
     xlab=expression("Leaf"~N[mass]~(mg~ g^{-1})),
     xlim=Nma_range, ylim=Jmass_range)  

logmod.Jmaxmass_N2pl <- lm(log(Jmaxmass)~log(Nmass), data=Acinew.Plim) 
logypredJmaxmassNpl <- predict(logmod.Jmaxmass_N2pl)
regr.Jmaxmass_N2pl <- lm(log(Jmaxmass)~log(Nmass), data=Acinew.Plim2, na.action=na.exclude) 

regr.Jmaxmass_Nnlim <- lm(log(Jmaxmass)~log(Nmass), data=Acinew.Pnonlim2, na.action=na.exclude) 
# diagnostics plots
#plot(regr.Jm.Pa.nonlim) #this produces 4 plots; can do par(mfrow = c(2,2)) prior to the plot cmd
summary(logmod.Jmaxmass_N2pl)
# this one non-Plim
lines(exp(logypredJmaxmassN2)~exp(logmod.Jmaxmass_N2$model$`log(Nmass)`), col="black", lty=5, lwd=3)

# now the P-nonlim confidence interval
preds6 <- predict.lm(regr.Jmaxmass_Nnlim, newdata = data.frame(Nmass=newx2), interval = "confidence")
predna6 <- predict.lm(regr.Jmaxmass_Nnlim, newdata = data.frame(Nmass=newx2))
# add transparent fill for the 95% CI
polygon(c(rev(newx2), newx2), c(rev(exp(preds6[ ,3])), exp(preds6[ ,2])), col = greytrans, border = NA)

#the Plim line for Jmax
lines(exp(logypredJmaxmassNpl)~exp(logmod.Jmaxmass_N2pl$model$`log(Nmass)`), col="orchid", lty=5, lwd=3)
summary(regr.Jmaxmass_N2pl)

# now predicts + conf interval
preds5 <- predict.lm(regr.Jmaxmass_N2pl, newdata = data.frame(Nmass=newx), interval = "confidence")
predna5 <- predict.lm(regr.Jmaxmass_N2pl, newdata = data.frame(Nmass=newx))
# add transparent fill for the 95% CI
polygon(c(rev(newx), newx), c(rev(exp(preds5[ ,3])), exp(preds5[ ,2])), col = plumtrans1, border = NA)

#overlay lines again
lines(exp(logypredJmaxmassN2)~exp(logmod.Jmaxmass_N2$model$`log(Nmass)`), col="black", lty=5, lwd=3)
lines(exp(logypredJmaxmassNpl)~exp(logmod.Jmaxmass_N2pl$model$`log(Nmass)`), col="orchid", lty=5, lwd=3)

# legend not needed
mtext("(b)", 3, line=0.5, adj=1, cex=1.4)

# BOTTOM PANELS
# Now for the SECOND FIGURE that plots photosynthesis as a function of P
Pm_range<-c(0.05, 4.5)   # modify axis range here, if needed c(3e+06, 6e+08)

logmod.Vcmaxmass_Pm2 <- lm(log(Vcmaxmass)~log(Pmass), data=Acibothfits2) #removed , na.action=na.exclude
logypredVcmaxmass_Pm2 <- predict(logmod.Vcmaxmass_Pm2)
Regr.Vcmaxmass_Pm2 <- lm(log(Vcmaxmass)~log(Pmass), data=Acibothfits2, na.action=na.exclude)

plot(Vcmaxmass ~ Pmass, data=Acibothfits2, log = "xy", pch=c(15,25,17,19)[Continent], cex=1.1,
     col=alpha("darkorchid4", 0.4), cex.lab=1.6, cex.axis=1.4, bty="l", 
     xaxs="i",  yaxs="i",# = list(type = "log"), yaxis = list(type = "log"),
     ylab=expression(italic(V)[cmax_mass]~(nmol~ g^{-1}~s^{-1})),
     xlab=expression("Leaf"~P[mass]~(mg~ g^{-1})),
     xlim=Pm_range, ylim=Vcmass_range)  

# this one for Vcmax
# predicts + interval for this line
newx7 <- seq(min(Acibothfits2$Pmass), max(Acibothfits2$Pmass), length.out=100)
preds9 <- predict.lm(Regr.Vcmaxmass_Pm2, newdata = data.frame(Pmass=newx7), interval = "confidence")
# add transparent fill for the 95% CI
polygon(c(rev(newx7), newx7), c(rev(exp(preds9[ ,3])), exp(preds9[ ,2])), col = violtrans, border = NA)
# for Vcmax
lines(exp(logypredVcmaxmass_Pm2)~exp(logmod.Vcmaxmass_Pm2$model$`log(Pmass)`), col="darkorchid4", lty=5, lwd=3)

#legend(x="bottomright", pch=19, legend=c(levels(Acinew.Plim$Continent)), col=Contin_set,
#       cex=1.2, bty="n", lwd=1.5)
mtext("(c)", 3, line=0.5, adj=1, cex=1.4)

# third panel this one for Jmax
logmod.Jmaxmass_Pm2 <- lm(log(Jmaxmass)~log(Pmass), data=Acibothfits2)   #, na.action=na.exclude)
logypredJmaxmass_Pm2 <- predict(logmod.Jmaxmass_Pm2)
Regr.Jmaxmass_Pm2 <- lm(log(Jmaxmass)~log(Pmass), data=Acibothfits2, na.action=na.exclude)

plot(Jmaxmass ~ Pmass, data=Acibothfits2, log = "xy", pch=c(15,25,17,19)[Continent], cex=1.1,
     col=alpha("darkorchid4", 0.4), cex.lab=1.6, cex.axis=1.4, bty="l", 
     xaxs="i",  yaxs="i",# = list(type = "log"), yaxis = list(type = "log"),
     ylab=expression(italic(J)[max_mass]~(nmol~ g^{-1}~s^{-1})),
     xlab=expression("Leaf"~P[mass]~(mg~ g^{-1})),
     xlim=Pm_range, ylim=Jmass_range)  

# Same problem as earlier with NaN values
preds8 <- predict.lm(Regr.Jmaxmass_Pm2, newdata = data.frame(Pmass=newx7), interval = "confidence")
# add transparent fill for the 95% CI
polygon(c(rev(newx7), newx7), c(rev(exp(preds8[ ,3])), exp(preds8[ ,2])), col = violtrans, border = NA)
lines(exp(logypredJmaxmass_Pm2)~exp(logmod.Jmaxmass_Pm2$model$`log(Pmass)`), col="darkorchid4", lty=5, lwd=3)

mtext("(d)", 3, line=0.5, adj=1, cex=1.4)

summary(Regr.Jmaxmass_Pm2)

par(op1)

dev.off()

# FIGURE

# Regression analysis for Table 1 manuscript
regr.Amass.Nm.all <- lm(log(Photomass) ~ log(Nmass)*Plim_status, data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.Nm.all<- lm(log(Vcmaxmass)~log(Nmass)*Plim_status, data=Acibothfits2, na.action=na.exclude)
# significant log(Nmass):Plim_status means slope difference, P = 0.00261
regr.Jmaxmass.Nm.all <- lm(log(Jmaxmass) ~ log(Nmass)*Plim_status, data=Acibothfits2, na.action=na.exclude) 
# significant log(Nmass):Plim_status means slope difference, P = 0.00261

## Fill in a set of missing regressions
regr.Vcmaxmass.Nm2.Plim<- lm(log(Vcmaxmass)~log(Nmass), data=Acinew.Plim, na.action=na.exclude)
regr.Amass.Pm.Pnonlim <- lm(log(Photomass) ~ log(Pmass), data=Acinew.Pnonlim, na.action=na.exclude) 
regr.Amass.Pm.Plim <- lm(log(Photomass) ~ log(Pmass), data=Acinew.Plim, na.action=na.exclude) 
regr.Vcmaxmass.Pm.Pnonlim <- lm(log(Vcmaxmass)~log(Pmass), data=Acinew.Pnonlim, na.action=na.exclude) 
regr.Vcmaxmass.Pm.Plim <- lm(log(Vcmaxmass)~log(Pmass), data=Acinew.Plim, na.action=na.exclude) 

regr.Vc.Pa.nonlim <- lm(log(Vcmax) ~ log(Parea), data=Acinew.Pnonlim, na.action=na.exclude) 
regr.Vc.Pa.lim <- lm(log(Vcmax) ~ log(Parea), data=Acinew.Plim, na.action=na.exclude)
# no more regressions
## MAKE TABLE
# see p. 173 "extracting output from ..." Rmanual
TableD1 <- matrix(data=NA, nrow=((25)), ncol=7) 
colnames(TableD1) <- c("Dependent", "Independent","d.f.","Intercept", "Slope", "R2", "P-value")
TableD1 <- data.frame(TableD1)
# row, column
TableD1[1,1:2]  <- names(regr.Amass.Nm.Pnonlim$model[1:2])
TableD1[1,3]  <- regr.Amass.Nm.Pnonlim$df
TableD1[1,4:5]  <- coef(regr.Amass.Nm.Pnonlim)
TableD1[1,6]  <- summary(regr.Jm.Pa.nonlim)$r.squared
TableD1[1,7]  <-anova(regr.Amass.Nm.Pnonlim)$'Pr(>F)'[1]

TableD1[2,1:2]  <- names(regr.Amass.Nm.Plim$model[1:2])
TableD1[2,3]  <- regr.Amass.Nm.Plim$df
TableD1[2,4:5]  <- coef(regr.Amass.Nm.Plim)
TableD1[2,6]  <- summary(regr.Jm.Pa.lim)$r.squared
TableD1[2,7]  <-anova(regr.Amass.Nm.Plim)$'Pr(>F)'[1]

TableD1[3,1:2]  <- names(regr.Amass.Pm.Pnonlim$model[1:2])
TableD1[3,3]  <- regr.Amass.Pm.Pnonlim$df
TableD1[3,4:5]  <- coef(regr.Amass.Pm.Pnonlim)
TableD1[3,6]  <- summary(regr.Amass.Pm.Pnonlim)$r.squared
TableD1[3,7]  <-anova(regr.Amass.Pm.Pnonlim)$'Pr(>F)'[1]

TableD1[4,1:2]  <- names(regr.Amass.Pm.Plim$model[1:2])
TableD1[4,3]  <- regr.Amass.Pm.Plim$df
TableD1[4,4:5]  <- coef(regr.Amass.Pm.Plim)
TableD1[4,6]  <- summary(regr.Amass.Pm.Plim)$r.squared
TableD1[4,7]  <-anova(regr.Amass.Pm.Plim)$'Pr(>F)'[1]

TableD1[5,1:2]  <- names(regr.Vcmaxmass.Nm.Pnonlim$model[1:2])
TableD1[5,3]  <- regr.Vcmaxmass.Nm.Pnonlim$df
TableD1[5,4:5]  <- coef(regr.Vcmaxmass.Nm.Pnonlim)
TableD1[5,6]  <- summary(regr.Vcmaxmass.Nm.Pnonlim)$r.squared
TableD1[5,7]  <-anova(regr.Vcmaxmass.Nm.Pnonlim)$'Pr(>F)'[1]

TableD1[6,1:2]  <- names(regr.Vcmaxmass.Nm2.Plim$model[1:2])
TableD1[6,3]  <- regr.Vcmaxmass.Nm2.Plim$df
TableD1[6,4:5]  <- coef(regr.Vcmaxmass.Nm2.Plim)
TableD1[6,6]  <- summary(regr.Vcmaxmass.Nm2.Plim)$r.squared
TableD1[6,7]  <-anova(regr.Vcmaxmass.Nm2.Plim)$'Pr(>F)'[1]

#
TableD1[7,1:2]  <- names(regr.Vcmaxmass.Pm.Pnonlim$model[1:2])
TableD1[7,3]  <- regr.Vcmaxmass.Pm.Pnonlim$df
TableD1[7,4:5]  <- coef(regr.Vcmaxmass.Pm.Pnonlim)
TableD1[7,6]  <- summary(regr.Vcmaxmass.Pm.Pnonlim)$r.squared
TableD1[7,7]  <-anova(regr.Vcmaxmass.Pm.Pnonlim)$'Pr(>F)'[1]

TableD1[8,1:2]  <- names(regr.Vcmaxmass.Pm.Plim$model[1:2])
TableD1[8,3]  <- regr.Vcmaxmass.Pm.Plim$df
TableD1[8,4:5]  <- coef(regr.Vcmaxmass.Pm.Plim)
TableD1[8,6]  <- summary(regr.Vcmaxmass.Pm.Plim)$r.squared
TableD1[8,7]  <-anova(regr.Vcmaxmass.Pm.Plim)$'Pr(>F)'[1]

TableD1[9,1:2]  <- names(regr.Vc.Na.nonlim$model[1:2])
TableD1[9,3]  <- regr.Vc.Na.nonlim$df
TableD1[9,4:5]  <- coef(regr.Vc.Na.nonlim)
TableD1[9,6]  <- summary(regr.Vc.Na.nonlim)$r.squared
TableD1[9,7]  <-anova(regr.Vc.Na.nonlim)$'Pr(>F)'[1]

TableD1[10,1:2]  <- names(regr.Vc.Na.Plim$model[1:2])
TableD1[10,3]  <- regr.Vc.Na.Plim$df
TableD1[10,4:5]  <- coef(regr.Vc.Na.Plim)
TableD1[10,6]  <- summary(regr.Vc.Na.Plim)$r.squared
TableD1[10,7]  <-anova(regr.Vc.Na.Plim)$'Pr(>F)'[1]

TableD1[11,1:2]  <- names(regr.Vc.Pa.nonlim$model[1:2])
TableD1[11,3]  <- regr.Vc.Pa.nonlim$df
TableD1[11,4:5]  <- coef(regr.Vc.Pa.nonlim)
TableD1[11,6]  <- summary(regr.Vc.Pa.nonlim)$r.squared
TableD1[11,7]  <-anova(regr.Vc.Pa.nonlim)$'Pr(>F)'[1]

TableD1[12,1:2]  <- names(regr.Vc.Pa.lim$model[1:2])
TableD1[12,3]  <- regr.Vc.Pa.lim$df
TableD1[12,4:5]  <- coef(regr.Vc.Pa.lim)
TableD1[12,6]  <- summary(regr.Vc.Pa.lim)$r.squared
TableD1[12,7]  <-anova(regr.Vc.Pa.lim)$'Pr(>F)'[1]

TableD1[13,1:2]  <- names(regr.Jmmass.Pm.nonlim$model[1:2])
TableD1[13,3]  <- regr.Jmmass.Pm.nonlim$df
TableD1[13,4:5]  <- coef(regr.Jmmass.Pm.nonlim)
TableD1[13,6]  <- summary(regr.Jmmass.Pm.nonlim)$r.squared
TableD1[13,7]  <-anova(regr.Jmmass.Pm.nonlim)$'Pr(>F)'[1]

TableD1[14,1:2]  <- names(regr.Jmmass.Pa2.lim$model[1:2])
TableD1[14,3]  <- regr.Jmmass.Pa2.lim$df
TableD1[14,4:5]  <- coef(regr.Jmmass.Pa2.lim)
TableD1[14,6]  <- summary(regr.Jmmass.Pa2.lim)$r.squared
TableD1[14,7]  <-anova(regr.Jmmass.Pa2.lim)$'Pr(>F)'[1]

TableD1[15,1:2]  <- names(regr.An.Na.nonlim$model[1:2])
TableD1[15,3]  <- regr.An.Na.nonlim$df
TableD1[15,4:5]  <- coef(regr.An.Na.nonlim)
TableD1[15,6]  <- summary(regr.An.Na.nonlim)$r.squared
TableD1[15,7]  <-anova(regr.An.Na.nonlim)$'Pr(>F)'[1]

TableD1[16,1:2]  <- names(regr.An.Na.lim$model[1:2])
TableD1[16,3]  <- regr.An.Na.lim$df
TableD1[16,4:5]  <- coef(regr.An.Na.lim)
TableD1[16,6]  <- summary(regr.An.Na.lim)$r.squared
TableD1[16,7]  <-anova(regr.An.Na.lim)$'Pr(>F)'[1]

TableD1[17,1:2]  <- names(regr.Jm.Na.nonlim$model[1:2])
TableD1[17,3]  <- regr.Jm.Na.nonlim$df
TableD1[17,4:5]  <- coef(regr.Jm.Na.nonlim)
TableD1[17,6]  <- summary(regr.Jm.Na.nonlim)$r.squared
TableD1[17,7]  <-anova(regr.Jm.Na.nonlim)$'Pr(>F)'[1]

TableD1[18,1:2]  <- names(regr.Jm.Na.lim$model[1:2])
TableD1[18,3]  <- regr.Jm.Na.lim$df
TableD1[18,4:5]  <- coef(regr.Jm.Na.lim)
TableD1[18,6]  <- summary(regr.Jm.Na.lim)$r.squared
TableD1[18,7]  <-anova(regr.Jm.Na.lim)$'Pr(>F)'[1]

#TableD1[13,1:2]  <- names(regr.Jm.Na.nonlim$model[1:2])
#TableD1[13,3]  <- regr.Jm.Na.nonlim$df
#TableD1[13,4:5]  <- coef(regr.Jm.Na.nonlim)
#TableD1[13,6]  <- summary(regr.Jm.Na.nonlim)$r.squared
#TableD1[13,7]  <-anova(regr.Jm.Na.nonlim)$'Pr(>F)'[1]

TableD1[19,1:2]  <- names(regr.Jm.Pa.nonlim$model[1:2])
TableD1[19,3]  <- regr.Jm.Pa.nonlim$df
TableD1[19,4:5]  <- coef(regr.Jm.Pa.nonlim)
TableD1[19,6]  <- summary(regr.Jm.Pa.nonlim)$r.squared
TableD1[19,7]  <-anova(regr.Jm.Pa.nonlim)$'Pr(>F)'[1]

TableD1[20,1:2]  <- names(regr.Jm.Pa.lim$model[1:2])
TableD1[20,3]  <- regr.Jm.Pa.lim$df
TableD1[20,4:5]  <- coef(regr.Jm.Pa.lim)
TableD1[20,6]  <- summary(regr.Jm.Pa.lim)$r.squared
TableD1[20,7]  <-anova(regr.Jm.Pa.lim)$'Pr(>F)'[1]

TableD1[21,1:2]  <- names(regr.Jmmass.Nm.nonlim$model[1:2])
TableD1[21,3]  <- regr.Jmmass.Nm.nonlim$df
TableD1[21,4:5]  <- coef(regr.Jmmass.Nm.nonlim)
TableD1[21,6]  <- summary(regr.Jmmass.Nm.nonlim)$r.squared
TableD1[21,7]  <-anova(regr.Jmmass.Nm.nonlim)$'Pr(>F)'[1]

TableD1[22,1:2]  <- names(regr.Jmmass.Nm.lim$model[1:2])
TableD1[22,3]  <- regr.Jmmass.Nm.lim$df
TableD1[22,4:5]  <- coef(regr.Jmmass.Nm.lim)
TableD1[22,6]  <- summary(regr.Jmmass.Nm.lim)$r.squared
TableD1[22,7]  <-anova(regr.Jmmass.Nm.lim)$'Pr(>F)'[1]

TableD1[23,1:2]  <- names(reg_Jm_Vc_all.Pnonlim$model[1:2])
TableD1[23,3]  <- reg_Jm_Vc_all.Pnonlim$df
TableD1[23,4:5]  <- coef(reg_Jm_Vc_all.Pnonlim)
TableD1[23,6]  <- summary(reg_Jm_Vc_all.Pnonlim)$r.squared
TableD1[23,7]  <-anova(reg_Jm_Vc_all.Pnonlim)$'Pr(>F)'[1]

TableD1[24,1:2]  <- names(reg_Jm_Vc_all.Plim$model[1:2])
TableD1[24,3]  <- reg_Jm_Vc_all.Plim$df
TableD1[24,4:5]  <- coef(reg_Jm_Vc_all.Plim)
TableD1[24,6]  <- summary(reg_Jm_Vc_all.Plim)$r.squared
TableD1[24,7]  <-anova(reg_Jm_Vc_all.Plim)$'Pr(>F)'[1]

TableD1[25,1:2]  <- names(reg_Jm_Vc_all.allP$model[1:2])
TableD1[25,3]  <- reg_Jm_Vc_all.allP$df
TableD1[25,4:5]  <- coef(reg_Jm_Vc_all.allP)
TableD1[25,6]  <- summary(reg_Jm_Vc_all.allP)$r.squared
TableD1[25,7]  <-anova(reg_Jm_Vc_all.allP)$'Pr(>F)'[1]

TableD1 # give it the once over
write.csv(TableD1, "TableD1_latest_Regr_coefficients.csv")

# NOW THE REGRESSIONS AND 2ND TABLE FOR ALL
# Do regressions on the whole dataset but excluding Northern Hemisphere stuff
# NOTE !!!! this should use Acibothfits2 !!!
# remove obs where there isn't both N and P = use data=Acibothfits2
regr.Amass.Nm.all <- lm(log(Photomass) ~ log(Nmass), data=Acibothfits2, na.action=na.exclude)
regr.Amass.Pm.all <- lm(log(Photomass) ~ log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Amass.Pmlma.all <- lm(log(Photomass) ~ log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Amass.PmNmlma.all <- lm(log(Photomass) ~ log(Nmass) + log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 

regr.Vcmaxmass.Nm.all <- lm(log(Vcmaxmass)~log(Nmass), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.Pm.all <- lm(log(Vcmaxmass)~log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.Pmlma.all <- lm(log(Vcmaxmass) ~ log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.Nmlma.all <- lm(log(Vcmaxmass) ~ log(Nmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.PmNm.all <- lm(log(Vcmaxmass) ~ log(Nmass) + log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmaxmass.PmNmlma.all <- lm(log(Vcmaxmass) ~ log(Nmass) + log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 

regr.Vc.Na.all <- lm(log(Vcmax) ~ log(Narea), data=Acibothfits2, na.action=na.exclude) 
regr.Vc.Pa.all <- lm(log(Vcmax) ~ log(Parea), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.Nm.all <- lm(log(Jmaxmass)~log(Nmass), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.Pm.all <- lm(log(Jmaxmass) ~ log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.Pmlma.all <- lm(log(Jmaxmass) ~ log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.Nmlma.all <- lm(log(Jmaxmass) ~ log(Nmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.PmNm.all <- lm(log(Jmaxmass) ~ log(Nmass) + log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Jmmass.PmNmlma.all <- lm(log(Jmaxmass) ~ log(Nmass) + log(Pmass)+log(LMA), data=Acibothfits2, na.action=na.exclude) 

regr.An.Na.all <- lm(log(Photo) ~ log(Narea), data=Acibothfits2, na.action=na.exclude) 
regr.Jm.Na.all <- lm(log(Jmax) ~ log(Narea), data=Acibothfits2, na.action=na.exclude) 
regr.Jm.Pa.all <- lm(log(Jmax) ~ log(Parea), data=Acibothfits2, na.action=na.exclude) 
# model with both N and P
regr.Jmmass.NPcomp.all <- lm(log(Jmaxmass)~log(Nmass)+log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Vcmass.NPcomp.all <- lm(log(Vcmaxmass)~log(Nmass)+log(Pmass), data=Acibothfits2, na.action=na.exclude) 

regr.Nmass.Pmass.all <- lm(log(Nmass)~log(Pmass), data=Acibothfits2, na.action=na.exclude) 
regr.Pmass.LMA.all <- lm(log(Pmass)~log(LMA), data=Acibothfits2, na.action=na.exclude) 
regr.Nmass.LMA.all <- lm(log(Pmass)~log(LMA), data=Acibothfits2, na.action=na.exclude) 

regr.NmPm.all <- lm(log(Pmass)~log(Nmass), data=Acibothfits2, na.action=na.exclude) 

## MAKE TABLE OF REGR RESULTS
# see p. 173 "extracting output from ..." Rmanual
TableD2 <- matrix(data=NA, nrow=((19)), ncol=7) 
colnames(TableD2) <- c("Dependent", "Independent","d.f.","Intercept", "Slope", "R2", "P-value")
TableD2 <- data.frame(TableD2)
# row, column
TableD2[1,1:2]  <- names(regr.Amass.Nm.all$model[1:2])
TableD2[1,3]  <- regr.Amass.Nm.all$df
TableD2[1,4:5]  <- coef(regr.Amass.Nm.all)
TableD2[1,6]  <- summary(regr.Amass.Nm.all)$r.squared
TableD2[1,7]  <-anova(regr.Amass.Nm.all)$'Pr(>F)'[1]

TableD2[2,1:2]  <- names(regr.Amass.Pm.all$model[1:2])
TableD2[2,3]  <- regr.Amass.Pm.all$df
TableD2[2,4:5]  <- coef(regr.Amass.Pm.all)
TableD2[2,6]  <- summary(regr.Amass.Pm.all)$r.squared
TableD2[2,7]  <-anova(regr.Amass.Pm.all)$'Pr(>F)'[1]

TableD2[3,1:2]  <- names(regr.Vcmaxmass.Nm.all$model[1:2])
TableD2[3,3]  <- regr.Vcmaxmass.Nm.all$df
TableD2[3,4:5]  <- coef(regr.Vcmaxmass.Nm.all)
TableD2[3,6]  <- summary(regr.Vcmaxmass.Nm.all)$r.squared
TableD2[3,7]  <-anova(regr.Vcmaxmass.Nm.all)$'Pr(>F)'[1]

TableD2[4,1:2]  <- names(regr.Vcmaxmass.Pm.all$model[1:2])
TableD2[4,3]  <- regr.Vcmaxmass.Pm.all$df
TableD2[4,4:5]  <- coef(regr.Vcmaxmass.Pm.all)
TableD2[4,6]  <- summary(regr.Vcmaxmass.Pm.all)$r.squared
TableD2[4,7]  <-anova(regr.Vcmaxmass.Pm.all)$'Pr(>F)'[1]

TableD2[5,1:2]  <- names(regr.Vc.Na.all$model[1:2])
TableD2[5,3]  <- regr.Vc.Na.all$df
TableD2[5,4:5]  <- coef(regr.Vc.Na.all)
TableD2[5,6]  <- summary(regr.Vc.Na.all)$r.squared
TableD2[5,7]  <-anova(regr.Vc.Na.all)$'Pr(>F)'[1]

TableD2[6,1:2]  <- names(regr.Vc.Pa.all$model[1:2])
TableD2[6,3]  <- regr.Vc.Pa.all$df
TableD2[6,4:5]  <- coef(regr.Vc.Pa.all)
TableD2[6,6]  <- summary(regr.Vc.Pa.all)$r.squared
TableD2[6,7]  <-anova(regr.Vc.Pa.all)$'Pr(>F)'[1]

TableD2[7,1:2]  <- names(regr.Jmmass.Nm.all$model[1:2])
TableD2[7,3]  <- regr.Jmmass.Nm.all$df
TableD2[7,4:5]  <- coef(regr.Jmmass.Nm.all)
TableD2[7,6]  <- summary(regr.Jmmass.Nm.all)$r.squared
TableD2[7,7]  <-anova(regr.Jmmass.Nm.all)$'Pr(>F)'[1]

TableD2[8,1:2]  <- names(regr.Jmmass.Pm.all$model[1:2])
TableD2[8,3]  <- regr.Jmmass.Pm.all$df
TableD2[8,4:5]  <- coef(regr.Jmmass.Pm.all)
TableD2[8,6]  <- summary(regr.Jmmass.Pm.all)$r.squared
TableD2[8,7]  <-anova(regr.Jmmass.Pm.all)$'Pr(>F)'[1]

TableD2[9,1:2]  <- names(regr.An.Na.all$model[1:2])
TableD2[9,3]  <- regr.An.Na.all$df
TableD2[9,4:5]  <- coef(regr.An.Na.all)
TableD2[9,6]  <- summary(regr.An.Na.all)$r.squared
TableD2[9,7]  <-anova(regr.An.Na.all)$'Pr(>F)'[1]

TableD2[10,1:2]  <- names(regr.Jm.Na.all$model[1:2])
TableD2[10,3]  <- regr.Jm.Na.all$df
TableD2[10,4:5]  <- coef(regr.Jm.Na.all)
TableD2[10,6]  <- summary(regr.Jm.Na.all)$r.squared
TableD2[10,7]  <-anova(regr.Jm.Na.all)$'Pr(>F)'[1]

TableD2[11,1:2]  <- names(regr.Jm.Pa.all$model[1:2])
TableD2[11,3]  <- regr.Jm.Pa.all$df
TableD2[11,4:5]  <- coef(regr.Jm.Pa.all)
TableD2[11,6]  <- summary(regr.Jm.Pa.all)$r.squared
TableD2[11,7]  <-anova(regr.Jm.Pa.all)$'Pr(>F)'[1]

TableD2[12,1:2]  <- names(regr.Vcmass.NPcomp.all$model[1:2])
TableD2[12,3]  <- regr.Vcmass.NPcomp.all$df
TableD2[12,4:5]  <- coef(regr.Vcmass.NPcomp.all)  # this will give the first coeff 
# but note there are two coeffs in this model
TableD2[12,6]  <- summary(regr.Vcmass.NPcomp.all)$r.squared
TableD2[12,7]  <-anova(regr.Vcmass.NPcomp.all)$'Pr(>F)'[1]

TableD2[13,1:2]  <- names(regr.Jmmass.NPcomp.all$model[1:2])
TableD2[13,3]  <- regr.Jmmass.NPcomp.all$df
TableD2[13,4:5]  <- coef(regr.Jmmass.NPcomp.all)  # this will give the first coeff 
# but note there are two coeffs in this model
TableD2[13,6]  <- summary(regr.Jmmass.NPcomp.all)$r.squared
TableD2[13,7]  <-anova(regr.Jmmass.NPcomp.all)$'Pr(>F)'[1]

TableD2[14,1:2]  <- names(regr.Nmass.Pmass.all$model[1:2])
TableD2[14,3]  <- regr.Nmass.Pmass.all$df
TableD2[14,4:5]  <- coef(regr.Nmass.Pmass.all)
TableD2[14,6]  <- summary(regr.Nmass.Pmass.all)$r.squared
TableD2[14,7]  <-anova(regr.Nmass.Pmass.all)$'Pr(>F)'[1]

TableD2[15,1:2]  <- names(regr.Pmass.LMA.all$model[1:2])
TableD2[15,3]  <- regr.Pmass.LMA.all$df
TableD2[15,4:5]  <- coef(regr.Pmass.LMA.all)
TableD2[15,6]  <- summary(regr.Pmass.LMA.all)$r.squared
TableD2[15,7]  <-anova(regr.Pmass.LMA.all)$'Pr(>F)'[1]

TableD2[16,1:2]  <- names(regr.Nmass.LMA.all$model[1:2])
TableD2[16,3]  <- regr.Nmass.LMA.all$df
TableD2[16,4:5]  <- coef(regr.Nmass.LMA.all)
TableD2[16,6]  <- summary(regr.Nmass.LMA.all)$r.squared
TableD2[16,7]  <-anova(regr.Nmass.LMA.all)$'Pr(>F)'[1]

TableD2[17,1:2]  <- names(regr.Jmmass.Nmlma.all$model[1:2])
TableD2[17,3]  <- regr.Jmmass.Nmlma.all$df
TableD2[17,4:5]  <- coef(regr.Jmmass.Nmlma.all)
# there are 2 coefs
TableD2[17,6]  <- summary(regr.Jmmass.Nmlma.all)$r.squared
TableD2[17,7]  <-anova(regr.Jmmass.Nmlma.all)$'Pr(>F)'[1]

TableD2[18,1:2]  <- names(regr.Jmmass.Pmlma.all$model[1:2])
TableD2[18,3]  <- regr.Jmmass.Pmlma.all$df
TableD2[18,4:5]  <- coef(regr.Jmmass.Pmlma.all)
# there are 2 coefs
TableD2[18,6]  <- summary(regr.Jmmass.Pmlma.all)$r.squared
TableD2[18,7]  <-anova(regr.Jmmass.Pmlma.all)$'Pr(>F)'[1]

TableD2[19,1:2]  <- names(regr.Jmmass.PmNmlma.all$model[1:2])
TableD2[19,3]  <- regr.Jmmass.PmNmlma.all$df
TableD2[19,4:5]  <- coef(regr.Jmmass.PmNmlma.all)
# there are 3 coefs
TableD2[19,6]  <- summary(regr.Jmmass.PmNmlma.all)$r.squared
TableD2[19,7]  <-anova(regr.Jmmass.PmNmlma.all)$'Pr(>F)'[1]

TableD2 # give it the once over
write.csv(TableD2, "TableD2_Regr_latest_coefficients_fulldata.csv")

### END

