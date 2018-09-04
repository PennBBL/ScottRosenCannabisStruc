################################################################
# This script will be used to produce figure 3 for:
# Scott et al. "Cannabis Use in Youth is Associated with Limited Alterations in Brain Structure"
# linking the relationship of cannabis consumption to structural
# imaging phenotypes
# This script in particular will run bootstrap analyses and prepare the confidence intervals
# from these analyses when comparing TBV, TBGM, and TBWM across the pairwise group comparisons
################################################################

################################################################
## Load library(s)
################################################################
library('psych')
library('ggplot2')
library('mgcv')
library('foreach')
library('doParallel')
library('simpleboot')

################################################################
## Load the data
################################################################
all.data <- readRDS('mjAnovaData.RDS')
 
################################################################
## Regress out covariates from data
################################################################
orig <- all.data
base.model <- paste("s(ageAtScan1)+sex+averageManualRating+factor(race2)+overall_psychopathology_ar_4factor")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.mod <- gam(tmp.formula, data=all.data)
  index <- which(complete.cases(all.data[,v]))
  all.data[index,name.val] <- NA
  all.data[index,name.val] <- scale(residuals(tmp.mod))
}

################################################################
## Prepare a parallel backend
################################################################
cl <- makeCluster(8)
registerDoParallel(cl)

################################################################
##Now find bootstrapped mean differences
################################################################
output.ci <- foreach(q=1:length(vars.of.interest), .combine='rbind', .packages='simpleboot') %dopar% {
    # First isolate our values
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
    occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    # Now calculate our bootstrapped differences
    n.v.o <- two.boot(non.vals, occ.vals, mean, 1000, student=T, M=50)
    n.v.f <- two.boot(non.vals, fre.vals, mean, 1000, student=T, M=50)
    f.v.o <- two.boot(fre.vals, occ.vals, mean, 1000, student=T, M=50)
    # Now return our confidence intervals
    out.row <- c(name.val,boot.ci(n.v.o)$basic[4:5],boot.ci(n.v.f)$basic[4:5],boot.ci(f.v.o)$basic[4:5], mean(non.vals),mean(occ.vals),mean(fre.vals))
    out.row
}
################################################################
## Remove our parallel backend
################################################################
stopCluster(cl)
colnames(output.ci) <- c("ROI","nvoMin","nvoMax","nvfMin","nvfMax","fvoMin","fvoMax","nonMean","occMean","freMean")

################################################################
## Calculate a standard error from our confidence intervals
################################################################
seDenom <- 2*1.96
nvoSE <- (as.numeric(output.ci[,3]) - as.numeric(output.ci[,2]))/seDenom
nvfSE <- (as.numeric(output.ci[,5]) - as.numeric(output.ci[,4]))/seDenom
fvoSE <- (as.numeric(output.ci[,7]) - as.numeric(output.ci[,6]))/seDenom

################################################################
## Now calculate a z score
################################################################
nvoZ <- (as.numeric(output.ci[,8]) - as.numeric(output.ci[,9])) / nvoSE
nvfZ <- (as.numeric(output.ci[,8]) - as.numeric(output.ci[,10])) / nvfSE
fvoZ <- (as.numeric(output.ci[,9]) - as.numeric(output.ci[,10])) / fvoSE

################################################################
## Now obtain a p value
################################################################
nvoP <- pnorm(-abs(nvoZ))*2
nvfP <- pnorm(-abs(nvfZ))*2
fvoP <- pnorm(-abs(fvoZ))*2

################################################################
## Now apply fdr correction to our p-values
## FDR correction will be applied to all lobular
## and regional values, within the respective specification
## and within each metric of interest i.e. GMD, CT, and Vol
################################################################
lowerLim <- c(1,140,238,361,373,385)
upperLim <- c(139,237,355,372,384,396)
inputPVal <- cbind(nvoP,nvfP,fvoP)
outputFDRPVal <- matrix(NA, ncol=3, nrow=length(nvoP))
for(colVal in 1:3){
    for(limVal in 1:length(lowerLim)){
    outputFDRPVal[lowerLim[limVal]:upperLim[limVal],colVal] <- p.adjust(inputPVal[lowerLim[limVal]:upperLim[limVal],colVal], method='fdr')
  }
}

################################################################
## Now produce histograms for group pairwise differences
## for TBV, TBGM, and TBWM
################################################################
pdf('user.minus.non.pdf')
for(q in c(358,359,360)){
  ## First obtain the BS mean differences
  v <- vars.of.interest[q]
  name.val <- names(all.data)[v]
  tmp.dat <- all.data[,c(name.val, 'marcat')]
  non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
  occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
  n.v.o <- two.boot(non.vals, occ.vals, mean, 1000, student=T, M=50)
  ci.vals <- round(boot.ci(n.v.o)$student[4:5], digits=3)
  ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
  ## Now plot our histogram for these values
  out.plot <- ggplot() +
    geom_histogram(aes(n.v.o$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(out.plot)
}
dev.off()

pdf('freq.minus.non.pdf')
for(q in c(358,359,360)){
    ## First obtain the BS mean differences
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    n.v.f <- two.boot(non.vals, fre.vals, mean, 1000, student=T, M=50)
    ci.vals <- round(boot.ci(n.v.f)$student[4:5], digits=3)
    ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
    ## Now plot our histogram for these values
    out.plot <- ggplot() +
    geom_histogram(aes(n.v.f$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    #xlab(name.val) +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
    print(out.plot)
}
dev.off()

pdf('freq.minus.user.pdf')
for(q in c(358,359,360)){
    ## First obtain the BS mean differences
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    f.v.o <- two.boot(fre.vals, occ.vals, mean, 1000, student=T, M=50)
    ci.vals <- round(boot.ci(f.v.o)$student[4:5], digits=3)
    ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
    ## Now plot our histogram for these values
    out.plot <- ggplot() +
    geom_histogram(aes(f.v.o$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
    print(out.plot)
}
dev.off()
