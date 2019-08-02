##Load packages
args = commandArgs(trailingOnly=T)

START_TIME <- Sys.time()
INPUT_COL <- as.numeric(args[1])
OUTPUT_DIR_NAME <- args[2]
LOG_FILE <- args[3]
SIGNIFICANT_CHROM <- args[4]

print(INPUT_COL)
print(OUTPUT_DIR_NAME)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(compositions))
suppressPackageStartupMessages(library(readr))

#read in the imputed genotypes data
g_data = t(read.table("impute_genotype.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE,stringsAsFactors=FALSE))

#Phenotype data
pheno <- read.table("../OTU_Data/novo_filtered.txt", header=TRUE, sep="\t", row.names = 1, check.names=FALSE,stringsAsFactors=FALSE,comment.char="@")

pheno$Taxonomy = NULL
pheno <- t(pheno)
samples <- rownames(pheno)
pheno <- as_tibble(pheno)
pheno$samples <- samples

print(colnames(pheno)[INPUT_COL])

# Filter the data
#replace all the . in the colnames with _ so they match the genotype data
pheno$samples <- gsub('\\.','_',pheno$samples)
p2 <- as_tibble(cbind(nms = names(pheno), t(pheno)))

#omit samples that have no genus
pheno <- pheno[,colnames(pheno)!="None"]

dim(g_data)
dim(pheno)

# Arrange Phenotype and Genotype Data
#find the common column and row names
com_id = intersect(colnames(g_data),pheno$samples)
otu_ids <- pheno$samples
taxes <- pheno[nrow(pheno),]

#Arrange the matrices accordingly
g_data = g_data[,com_id]
pheno = pheno[pheno$samples %in% com_id,]

## Pick a Phonetype to Map
#Clr transformation
phenotype_clr = clr(pheno)

#plot it to visualize the normalization
phenotype <- as.numeric(unlist(phenotype_clr[,INPUT_COL]))
ggplot(data=pheno,aes(x=phenotype)) + geom_histogram()

ggsave(paste(OUTPUT_DIR_NAME,"/hist.pdf",sep=""))
NORMALNESS_P_VALUE <- shapiro.test(phenotype)$p.value

# Perform Genome Wide QTL Scan
#get the LOD score for each location in the genome. Each location is a row in the g_data
lod_RIL = apply(g_data, 1, function(x){
  dat = cbind(1,suppressWarnings(as.factor(x)),phenotype)
  dat = dat[complete.cases(dat), ]
  fit1 = .lm.fit(dat[,1:2], dat[,3], tol=1e-12)
  fit0 = .lm.fit(as.matrix(dat[,1]),dat[,3], tol=1e-12)
  rss0 = sum(fit0$residuals^2)
  rss1 = sum(fit1$residuals^2)
  (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
})

#Now we will have to calculate the empirical threshold for the LOD scores. We would randomly shuffle the phenotype values and repeat the procedure n number of time (takes a couple of minute). We then take the 95th quantile of that distribution.
#n is the number of permutation to perform
n = 100
perm = sapply(1:n, function(i){
  if(i%%10==0){
    message(paste0("Bootstrapping ",i))
  }
  pheno_sample = sample(phenotype)
  lod_perm = apply(g_data, 1, function(x){
    dat = cbind(1,suppressWarnings(as.factor(x)),pheno_sample)
    dat = dat[complete.cases(dat), ]
    fit1 = .lm.fit(dat[,1:2], dat[,3], tol=1e-12)
    fit0 = .lm.fit(as.matrix(dat[,1]),dat[,3], tol=1e-12)
    rss0 = sum(fit0$residuals^2)
    rss1 = sum(fit1$residuals^2)
    (dim(dat)[1]/2)*(log10(rss0)-log10(rss1))
  })
  lod_perm[which(lod_perm == max(lod_perm,na.rm=TRUE))]
})
threshold = quantile(as.numeric(unlist(perm)), 0.95, na.rm = TRUE)

#Analyze and Visualize QTLs
#Let's pull out the most significant marker and look at it's effect size. The genotype column should again be edited when swtiching between different genotype coding. The purpose of this is make sure that the most significant marker is actually significant.
lod_RIL_small <- lod_RIL[grepl(SIGNIFICANT_CHROM,unlist(names(lod_RIL)))]
min_marker = names(which(lod_RIL_small == max(lod_RIL_small,na.rm=TRUE)))
effect_data = data.frame(accession = com_id, genotype=as.factor(g_data[min_marker,]),pheno=phenotype)
effect_data = effect_data[complete.cases(effect_data), ]
effect_data$genotype = as.factor(effect_data$genotype)

#Find the p value for the two different genotypes
types <- levels(effect_data$genotype)
g1 <- as.numeric(unlist(pheno[effect_data$genotype==types[1],]))
g2 <- as.numeric(unlist(pheno[effect_data$genotype==types[2],]))
test_result <- t.test(g1,g2)
P_VALUE <- test_result$p.value

#Plot box plots of the two genotypes
ggplot(effect_data, aes(x=genotype,y=pheno)) +
geom_boxplot(width=0.2, outlier.shape=NA, notch = TRUE) +
geom_jitter(width=0.1, size=0.5) +
ggtitle(paste0(min_marker))
ggsave(paste(OUTPUT_DIR_NAME,"/boxplot.pdf",sep=""))

# Create Plot of all the QTL data
# First format the data for plotting
# Some preprocessing of the data
marker = as.data.frame(str_split_fixed(rownames(g_data), "_", 3))
marker$V1 = as.data.frame(str_split_fixed(marker$V1, "-", 2))$V1
marker$V1 = factor(marker$V1, levels=c(paste0("chr",1:12),"scaffold"))
marker$V2 = 1:dim(marker)[1]
marker = marker[,-3]
plotdata = cbind(marker, lod_RIL)
colnames(plotdata) = c("chrom","pos","lod")
head(plotdata)

#Count the number of loci that are over threshhold
OVER_THRESH_COUNT = sum(plotdata$lod > threshold)

# Actually make the QTL Pplot
ggplot(plotdata, aes(pos, lod, group=chrom,color=chrom)) +
geom_point(alpha=0.7) +
geom_line(alpha=0.7) +
facet_wrap(~ chrom, ncol = 13, scales = "free_x", strip.position="bottom") +
geom_hline(yintercept = threshold, linetype = "dotted", color='red') +
theme_minimal() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(0,"line"),
axis.text.x=element_blank(),strip.placement = "outside",legend.position = "none")+
labs(x = "Chromosome", y = "LOD", color = "", linetype = "")
ggsave(paste(OUTPUT_DIR_NAME,"/qtl.pdf",sep=""))

## Find candidate genes around the QTLs
#Now the most important part is look into the candidate genes. We would load up the gene annotation of the new assembly and look at 2.5 Mb on either side of the marker and look at the description of the genes.
plotdata_list <- plotdata
plotdata_list$marker <- rownames(plotdata)

arranged <- arrange(plotdata_list, desc(lod))
head(arranged)

annotation = read.table("genes_description.txt", header=TRUE, sep="\t", quote="\"", row.names = NULL, check.names=FALSE, stringsAsFactors = FALSE)

candidates <- data.frame(matrix(ncol=0,nrow=0))
for (i in 1:10) {
  c <- arranged$chrom[i]
  if (c != SIGNIFICANT_CHROM){
    next
  }
  m <- arranged$marker[i]
  m <- as.numeric(gsub('\\w+\\d?_','',m))
  l <- arranged$lod[i]
  if (l <= threshold){
    next
  }
  print(l)
  if( c != 'scaffold'){
    toAdd <- annotation %>% 
    filter(chrom==c & 
             start>(m-2500000) & 
             end<=(m+2500000)) %>% 
      select(gene, description)
    toAdd$lod <- l
    candidates <- rbind(candidates, toAdd)
  }
}
print(dim(candidates)) 

head(candidates)
aggr <- aggregate(data.frame(count = candidates$gene), list(value = candidates$gene), length)
colnames(aggr) <- c("gene", "count")
head(aggr)

max(aggr$count)

joined <- left_join(candidates, aggr, by = "gene")
head(joined)

unique_candidates <- data.frame(matrix(ncol=0,nrow=0))
temp <- list(matrix(ncol=0, nrow=0))
for (i in 1:nrow(joined)) {
  gene <-joined$gene[i]
  if (gene %nin% temp){
    unique_candidates <- rbind(unique_candidates, joined[i, ])
    temp <- rbind(temp, gene)
  }
  else {
  }
}
dim(unique_candidates)

write.table(unique_candidates, paste(OUTPUT_DIR_NAME,"/CandidateGenes.txt",sep=''), row.names=F, col.names=F, sep="\t", quote=FALSE) 

END_TIME <- Sys.time()
output_line <- paste(INPUT_COL,NORMALNESS_P_VALUE,P_VALUE,OVER_THRESH_COUNT,OUTPUT_DIR_NAME,START_TIME,END_TIME,as.numeric(END_TIME - START_TIME),sep=',')

print(output_line)
write(output_line,file=LOG_FILE,append=TRUE)
