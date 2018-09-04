################################################################
# This script will be used to produce figure 2 for:
# Scott et al. "Cannabis Use in Youth is Associated with Limited Alterations in Brain Structure"
# linking the relationship of cannabis consumption to structural
# imaging phenotypes
# This script in particular will prepare volumetric density plots
# for a limbic and basal ganlia regions, as well as associated F stats
################################################################

################################################################
## Load library(s)
################################################################
library('psych')
library('ggplot2')
library('mgcv')

################################################################
## Load the data
################################################################
all.data <- readRDS('mjAnovaData.RDS')

################################################################
## Create function(s)
################################################################
## Cretae a function which will return a density plot for the input variables
writeDensityPlot <- function(dataIn, var.of.interest=NULL, covariates="s(ageAtScan1)+sex+averageManualRating+factor(race2)+overall_psychopathology_ar_4factor",paraMetricValue=NULL, nonParametricValue=NULL, linMod=FALSE){
    ## First check we have a variable of interest
    if (missing(var.of.interest)) { stop("Missing requireed input variable (DV)")}
    ## Regress out covairates
    tmp.form <- as.formula(paste(var.of.interest, covariates, sep="~"))
    tmp.mod <- mgcv::gam(tmp.form, data=dataIn, na.action=na.exclude)
    if(linMod==TRUE){
        tmp.mod <- lm(tmp.form, data=dataIn, na.action=na.exclude)
    }
    
    ## Now prepare our residuals
    all.data$tmpvals <- scale(as.numeric(residuals(tmp.mod)))
    
    ## Now create our density plot
    out.hist <- ggplot(dataIn, aes(tmpvals)) +
    geom_density(data=subset(dataIn,marcat=='MJ Non-User'), fill="#009E73",alpha=.4) +
    geom_density(data=subset(dataIn,marcat=='MJ Occ User'), fill="#9ad0f3", alpha=.4) +
    geom_density(data=subset(dataIn,marcat=='MJ Freq User'), fill="#D55E00",alpha=.4) +
    coord_cartesian(xlim=c(-5, 5), ylim=c(0,.9)) +
    theme_bw() +
    theme(text = element_text(size=30))
    
    ## Now create our text strings to add to the density plot
    if(!identical(paraMetricValue, NULL)){
        f.string <- paste("F-statistic = ",round(as.numeric(paraMetricValue), digits=3), " ")
        out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = f.string, vjust=1.5, hjust=1, parse = F, size=10)
    }
    if(!identical(nonParametricValue, NULL)){
        chi.string <- paste("Kruskal-Wallis H = ",round(as.numeric(nonParametricValue), digits=3), " ")
        out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = chi.string, vjust=3.5, hjust=1, parse=F, size=10)
    }
    # Now return the object
    return(out.hist)
}

## Now make a function which will return a legend color scale
returnColorScale <- function(){
    colfunc1 <- colorRampPalette(c("#009E73"))
    colfunc2 <- colorRampPalette(c("#9ad0f3"))
    colfunc3 <- colorRampPalette(c("#D55E00"))
    
    ## Now return a plot with these colors
    plot(rep(1, 99), col = c(colfunc1(33), colfunc2(33), colfunc3(33)), pch = 15, cex = 10)
}

################################################################
## Now identify our variables of interest
## This will loop through all of the subcortical ROI's
## as well as TBM, TBGM, and TBWM
################################################################
roi.of.interest <- c(1550:1552,109,110,111,112,114,115,127,128,129,130,131,132, 121, 122)

################################################################
## Ensure we do not have any NA's for our models
################################################################
all.data.plot <- all.data[complete.cases(all.data[,roi.of.interest]),]

################################################################
## Now loop through and create our histogram and values
################################################################
output.values.nl <- NULL
for(s in roi.of.interest){
    # Obtain an F value
    mod.one <- gam(all.data.plot[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + marcat, data=all.data.plot)
    aov.mod.one <- anova.gam(mod.one)
    
    # Now obtain the Kruskall-Wallis K value
    all.data.plot$tmpvals <- as.numeric(scale(gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)$residuals))
    mod.two <- kruskal.test(tmpvals ~marcat, data = all.data.plot)
    
    ## Now prepare the output
    output.row <- c(names(all.data.plot)[s], round(aov.mod.one$pTerms.table['marcat',c('F')], digits=2), round(mod.two$statistic,digits=2),round(aov.mod.one$pTerms.table['marcat',c('p-value')], digits=2),round(mod.two$p.value, digits=2))
    output.values.nl <- rbind(output.values.nl, output.row)
}

################################################################
## Now produce the density plots
################################################################
pdf('figure2VolumeDensity.pdf')#,units='in',res=300,width=10,height=8)
for(i in 1:dim(output.values.nl)[1]){
    ## Prepare the ROI's density plot
    tmp.plot <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[i,2], nonParametricValue=output.values.nl[i,3], var.of.interest=output.values.nl[i,1])

    ## Now prepare x axis title
    tmp.string <- gsub(x=output.values.nl[i,1], pattern="mprage_jlf_vol_", replacement="")
    tmp.string <- gsub(x=tmp.string, pattern="L_", replacement="Left ")
    tmp.string <- gsub(x=tmp.string, pattern="R_", replacement="Right ")
    tmp.string <- gsub(x=tmp.string, pattern="_", replacement=" ")
    if(i ==1){tmp.string <- "Total Brain Volume"}
    if(i ==2){tmp.string <- "Gray Matter Volume"}
    if(i ==3){tmp.string <- "White Matter Volume"}
    
    ## Now remove finalize and print the plot
    tmp.plot <- tmp.plot + theme(axis.title.y = element_text(color="white")) + xlab(tmp.string)
    print(tmp.plot)
}
dev.off()
