install.packages("wtest")
setwd('/Users/emulciber/MitoClub/fish-strong-purifying-selection/scripts/w-test/')
library(wtest)

geno<- read.table(file="genotype_sample.txt", header=TRUE)
phen<- unlist(read.table(file="phenotype_sample.txt", header=FALSE, row.names=1))

#  Calculate main effect or pairwise interaction test: using default h and f parameters:

#input y as a numeric vector, data as N by P matrix.
w.out1<- wtest(data=geno, y=phen, w.order = 1)

# Step 1.  Estimate h and f parameters: Using bootstrapped samples under null hypothesis to estimate h and f. The option in the function hf.calculation includes n.marker and n.sample. We suggest using all the samples in the data to perform the estimation, which is the default option. When sample size is smaller than 1,000, B should be increased to 400, or 1000, and n.marker is suggested to use 100 or above.

#Estimate h and f for main effect
hf1<-hf(data = geno, B=100, w.order = 1, n.marker=2, n.sample=2775)

#Estimate h and f for interaction effect
hf2<-hf(data = geno, B=100, w.order = 2, n.marker=2, n.sample=2775)

#Step 2A. Calculate pairwise interactions exhaustively using estimated h and f.

w.out2<-wtest(data = geno, y = phen, w.order=2, input.pval= NULL, hf1=hf1, hf2 = hf2)
w.pairs<- w.out2$results

#Step 2B. Calculate pairwise interactions using selected main effect markers:
  
#Input.pval = 0.1:  only input SNPs which main effect p-values are smaller than 0.1 
#Output.pval = 0.0005:  output only the pairs that p-values are smaller than 0.0005.
w.out3<-wtest(data = geno, y = phen, w.order=2, hf1=hf1, hf2 = hf2, input.pval= 0.1, output.pval =0.0005)
w.pairs<- w.out3$results

#output pairs:
write.table(w.pairs, file= "w.pairs.txt", quote=F, col.names=T, row.names=F, sep= "\t")

#Step 3A. Distribution plot:  Diagnostic checking on estimated probability distributions with hypothesis probability distributions. The function will output two distribution curves in the same plot, for each k (number of categories formed by the biomarker), or h and f. 
w.diagnosis(data=geno, w.order = 2, n.rep = 100, hf2=hf2, main=NULL, xlab=NULL, ylab=NULL)

#Step 3B. qq plot:  Diagnostic checking on whether there is inflation of type I error
w.qqplot(geno, w.order = 2, hf2= hf2)
abline(0,1)

