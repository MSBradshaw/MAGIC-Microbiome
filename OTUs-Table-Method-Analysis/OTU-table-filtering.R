library(readr)
library(tibble)

qiime <- read_tsv('Documents/qiime-taxonomy.tsv')

d <- qiime[!grepl('denovo',unlist(qiime[,1])),]

d$taxonomy <- gsub('D_.*D_','',d$taxonomy)
d <- d[!grepl('uncultured',d$taxonomy),]
d <- d[!grepl('Ambiguous',d$taxonomy),]
d <- d[!grepl('metagenome',d$taxonomy),]
#the total number of reads mapped to a useful OTU using Qiime
qsum <- sum(d[,2:(ncol(d)-1)])

novo <- read_tsv('Documents/novo_filtered.txt')
nsum <- sum(novo[,2:(ncol(novo)-1)])


ninja0 <- read_tsv('Documents/ninja_otus.tsv')
ninja <- ninja0
ninja[,1] <- gsub('D_.*D_','',unlist(ninja[,1]))
ninja <- ninja[!grepl('uncultured',unlist(ninja[,1])),]
ninja <- ninja[!grepl('Ambiguous',unlist(ninja[,1])),]
ninja <- ninja[!grepl('metagenome',unlist(ninja[,1])),]
dim(ninja)
ninjasum <- sum(ninja[,2:ncol(ninja)])


results_tib <- tibble(
  Method=c('Qiime','Ninja','NovoGene'),
  Initial_Number_of_OTUs=c(nrow(qiime),nrow(ninja0),nrow(novo)),
  Filtered_Number_of_OTUs=c(nrow(d),nrow(ninja),nrow(novo)),
  Filtered_Total_Number_of_Reads=c(qsum,ninjasum,nsum)
)
results_tib
