source("globals.R")
source("GWAA.R")

load(working.data.fname(8))

# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])
rm(genoBim)

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))

# Combine typed and imputed
GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)
head(GWAScomb)

tail(GWAScomb)


# Subset for CETP SNPs
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP)[,c("SNP","p.value","Neg_logP","chr","position","type","gene")]
print(CETP)

write.csv(CETP, CETP.fname, row.names=FALSE) # save for future use

save(genotype, clinical, pcs, imputed, target, rules, phenoSub, support, genes,
     impCETP, impCETPgeno, imputeOut, GWASout, GWAScomb, CETP, file=working.data.fname(9))


load(working.data.fname(9))
########################
# Visualization and QC

########################
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c("black"), col.imputed=c("blue"), col.text="black",
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE,
                           bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
                           Lstringent.adjustment=10000) {
  
  bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
  Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
  xscale <- 1000000
  
  manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
  
  #sort the data by chromosome and then location
  manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
  manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
  
  ##Finding the maximum position for each chromosome
  max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
  max.pos2 <- c(0, cumsum(max.pos))                  
  
  #Add spacing between chromosomes
  max.pos2 <- max.pos2 + c(0:21) * xscale * 10
  
  #defining the positions of each snp in the plot
  manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
  
  # alternate coloring of chromosomes
  manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
  
  # draw the chromosome label roughly in the middle of each chromosome band
  text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
  
  # Plot the data
  plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
       pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab=NA,
       ylab="Negative Log P-value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
  #Add x-label so that it is close to axis
  mtext(side = 1, "Chromosome", line = 1.25)
  
  points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"],
         pch=20, cex=.4, col = col.imputed)
  
  points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
         pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
  
  axis(2)
  abline(h=0)
  
  SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
  
  #Add legend
  legend("topright",c("Bonferroni corrected threshold (p = 5E-8)", "Candidate threshold (p = 5E-6)"),
         border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
         lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.6)
  
  #Add chromosome number
  text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=.8)
  
  #Add bonferroni line
  abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
  
  #Add "less stringent" line
  abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
  
  #Plotting detected genes
  #Were any genes detected?
  if (length(SigNifSNPs)>0){
    
    sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
    
    points(manhat.ord$pos[sig.snps]/xscale,
           manhat.ord$Neg_logP[sig.snps],
           pch=20,col=col.detected, bg=col.detected,cex=0.5)
    
    text(manhat.ord$pos[sig.snps]/xscale,
         manhat.ord$Neg_logP[sig.snps],
         as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.5)
  }
}

# Create Manhattan Plot
GWAS_Manhattan(GWAScomb)


##################### 
# Quantile-quantile plots and the λ-statistic

# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")] # remove all extra factors, leave only phenotype

GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname)

GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
par(mfrow=c(1,2))
par(mfrow=c(1,1))
lambdaAdj <- estlambda(GWASout$t.value^2,plot=TRUE,method="median")
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2,plot=TRUE,method="median")

cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))

# Calculate standardized lambda
lambdaAdj_1000<-1+(lambdaAdj$estimate-1)/nrow(phenoSub)*1000
lambdaUnadj_1000<-1+(lambdaUnadj$estimate-1)/nrow(phenoSub)*1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))


########################
# Heatmap
library(LDheatmap)
library(rtracklayer)

# Add "rs247617" to CETP
CETP <- rbind.fill(GWASout[GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- cbind(genotype[,colnames(genotype) %in% CETP$SNP], impCETPgeno)     # CETP subsets from typed and imputed SNPs

# Subset SNPs for only certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
CETP <- CETP[CETP$SNP %in% colnames(subgen),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames(subgen),CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs

ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
library(ggplot2)
library(grid)

plot.new()
llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

grid.draw(ggplotGrob({
  qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position),
        asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) + 
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(0.75)), legend.position = "none", 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_color_manual(values = c("red", "black"))
}))


################
# Regional Association
# Create regional association plot
# Create data.frame of most significant SNP only
library(postgwas)
snps<-data.frame(SNP=c("rs1532625", "rs247617"))

# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWAScomb, c(p.value="P", chr="CHR", position="BP"))


# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19

myconfig <- biomartConfigs$hsapiens
myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL" 
myconfig$hsapiens$snp$host <- "grch37.ensembl.org"
myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"

# Run regionalplot using HAPMAP data (pop = CEU)
regionalplot(snps, GWAScomb, biomartConfigs$hsapiens, window.size = 400000, draw.snpname = data.frame(
  snps = c("rs1532625", "rs247617"), 
  text = c("rs1532625", "rs247617"),
  angle = c(20, 160),
  length = c(1, 1), 
  cex = c(0.8)
),
ld.options = list(
  gts.source = 2, 
  max.snps.per.window = 2000, 
  rsquare.min = 0.8, 
  show.rsquare.text = FALSE
),
out.format = list(file = NULL, panels.per.page = 2))

