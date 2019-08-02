There are four files for quantitative trait loci (QTL) mapping located in this dirrectory.
The external files refered to in these scripts can all be found in '../Data/'

- qtl-mapping-clr.r
R script containing the QTL mapping code for the data normalized using the CLR normalization methods and fit with a linear model. 
Additionally this is the only file which contains code for extracting gene candidates.

- qtl-mapping-ilr.r
R script containing QTL code for ILR normalized data fit to a linear model.

- qtl-mapping-zi.r
R script containing QTL mapping code for CLR and ILR normalization methods with a zero inflated linear model.
You can switch between the CLR and IRL normalization by commenting or uncommenting the irl normalization function call on line 60.

- qtl-mapping-GLM.Rmd
R markdown file containing more or less the same code as the other script but uses the gerneralized linear model (GLM) and rarefied data.
