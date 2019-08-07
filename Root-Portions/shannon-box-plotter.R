library(readr)
library(ggplot2)

d <- read_tsv('alpha_shannon_results.txt',skip=1,col_names = FALSE)
d$`Root Portion` <- gsub('\\d','',d$X1)
colnames(d) <- c('Sample','Shannon Index','Root Portion')

p <- ggplot(data=d,aes(y=`Shannon Index`,x=`Root Portion`)) + geom_boxplot() + geom_point() + 
  ggtitle('Shannon Indices of Top and Fine Root')
p
ggsave('shannon-box-plots.png')

#calculate the variance of the various groups
var(d[d$`Root Portion` == 'Top',]$`Shannon Index`)
#0.03758973
var(d[d$`Root Portion` != 'Top',]$`Shannon Index`)
#3.856358