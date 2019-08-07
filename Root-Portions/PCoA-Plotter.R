library(readr)
library(tibble)
library(ggplot2)
library(ape)

a = read_tsv('RootPortions/pcoa_weighted_unifrac.txt',skip=9,col_names = FALSE)
a <- a[,1:6]

a$Group <- c('Top','Fine','Top','Fine','Top','Fine')
colnames(a) <- c('sample','PC1','PC2','PC3','PC4','PC5','Root Portion')

p <- ggplot(data = a, aes(x = PC1,y=PC2, color=`Root Portion`)) + geom_point() + 
  xlab('PC1 79.5%') + 
  ylab('PC2 11.8%') + 
  ggtitle('PCoA Fine and Top Roots, Weighted Unifraq')
p
#ggsave('RootPortions/PCoA-Beta-Diversity-Weighted.png')


b = read_tsv('RootPortions/filtered_top_fine.txt')
bm <- as.matrix(b[,2:ncol(b)])
pbm <- prcomp(t(bm))
s <- summary(pbm)
s$importance[2,1]
tib <- as_tibble(pbm$x)
tib$`Root Portion` <- c('Top','Fine','Top','Fine','Top','Fine')
ggplot(data <- tib, aes(x=PC1,y=PC2,color=`Root Portion`)) + geom_point() +
  xlab(paste('PC1',(s$importance[2,1]*100),'%')) + 
  ylab(paste('PC2',(s$importance[2,2]*100),'%')) + 
  ggtitle('PCA Fine and Top Roots Abundance')
ggsave('RootPortions/PCA-Raw-Abundance.png')
