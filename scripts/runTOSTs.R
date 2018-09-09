################################################################
# This script will be used to run the TOST analyses for:
# Scott et al. "Cannabis Use in Youth is Associated with Limited Alterations in Brain Structure"
# linking the relationship of cannabis consumption to structural
# imaging phenotypes
# This script in particular will run the TOST proceure from the equivalence library
# probing for equivalence between the group comparisons after controlling
# for covariates of interest
################################################################

################################################################
## Load library(s)
################################################################
library('psych')
library('ggplot2')
library('mgcv')
library('equivalence')

################################################################
## Load the data
################################################################
all.data <- readRDS('mjAnovaData.RDS')

################################################################
## Regress out covariates from data
################################################################
equiv.data <- all.data
base.model <- paste("s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
    name.val <- names(all.data)[v]
    tmp.formula <- as.formula(paste(name.val, "~", base.model))
    tmp.col <- rep(NA, 1504)
    tmp.mod <- gam(tmp.formula, data=equiv.data)
    index <- which(complete.cases(all.data[,v]))
    equiv.data[index,name.val] <- NA
    equiv.data[index,name.val] <- scale(residuals(tmp.mod))
}

################################################################
## Run the TOST procedure
################################################################
output.tost.vals <- NULL
for(v in vars.of.interest){
    name.val <- names(all.data)[v]
    ## Run TOST with +/-.5 effect size boundary
    test.val.one <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, var.equal=F)
    test.val.two <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Freq User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, var.equal=F)
    test.val.three <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Freq User"),v], paired=F, var.equal=F)
    print(v)
    ## Run TOST with +/-.3 effect size boundary
    test.val.four <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, epsilon = .6, var.equal=F)
    test.val.five <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Freq User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, epsilon = .6, var.equal=F)
    test.val.six <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Freq User"),v], paired=F, epsilon = .6, var.equal=F)
    ## Now prepare our output values
    output.row <- c(name.val, test.val.one$tost.p.value, test.val.two$tost.p.value, test.val.three$tost.p.value, test.val.four$tost.p.value, test.val.five$tost.p.value, test.val.six$tost.p.value)
    output.tost.vals <- rbind(output.tost.vals, output.row)
}

################################################################
## Now apply fdr correction to our p-values
## FDR correction will be applied to all lobular
## and regional values, within the respective specificity
## and within each metric of interest i.e. GMD, CT, and Vol
################################################################
lowerLim <- c(1,140,238,361,373,385)
upperLim <- c(139,237,355,372,384,396)
inputPVal <- cbind(output.tost.vals[,2],output.tost.vals[,3],output.tost.vals[,4],output.tost.vals[,5],output.tost.vals[,6],output.tost.vals[,7])
outputFDRPVal <- matrix(NA, ncol=6, nrow=dim(output.tost.vals)[1])
for(colVal in 1:6){
    for(limVal in 1:length(lowerLim)){
        outputFDRPVal[lowerLim[limVal]:upperLim[limVal],colVal] <- p.adjust(inputPVal[lowerLim[limVal]:upperLim[limVal],colVal], method='fdr')
    }
}

################################################################
## Now return ROI's that are not FDR sig
################################################################
output.tost.vals <- cbind(output.tost.vals, outputFDRPVal)
all.out <- output.tost.vals[which(outputFDRPVal[,5] > .05 | outputFDRPVal[,6] >.05),]
write.csv(all.out, "outputTOSTNonSig.csv", quote=F,row.names=F)
