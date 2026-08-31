library(qtl2) 
# 1. Load cross 
cross <- read_cross2("cross.yaml") 
print(summary(cross)) 

##Can remove strains with this code
##cross <- cross[rownames(cross$pheno) != "strain", ]

# 2. Calculate genotype probabilities 
prob <- calc_genoprob(cross, error_prob = 0.002, quiet = FALSE) 

k <- calc_kinship(pr, type = "loco")

# 3. QTL scan across all phenotypes 
out <- scan1(prob, cross$pheno, kinship = k) 

# 4. Permutation test for significance threshold (n_perm takes a long time, reduce to 100 for testing) 
set.seed(1)
perm <- scan1perm(pr, cross$pheno, kinship = k, n_perm = 1000) 
thr_val <- summary(perm, alpha = 0.05)

# 5. Find significant peaks 
peaks <- find_peaks(out, cross$gmap, threshold = thr_val, peakdrop = 1)
colnames(peaks) <- c("Phenotype Index", "Phenotype", "Chromosome", "cM Position", "LOD")
peaks$Threshold <- thr_val[1, as.character(peaks$Phenotype)]
peaks$Margin <- peaks$LOD - peaks$Threshold
print(peaks)

# 6. Plotting
plot_scan1(out, cross$gmap[as.character(1:19)],
           lodcolumn = "phenotype",
           main      = "phenotype",
           ylab      = "LOD")
abline(h = thr_val[, "phenotype"], col = "red", lty = 2)
