"# GeneQC" 
# GeneQC

Tutorial:

1.	Install GeneQC package from GitHub.  This process only needs to happen the first time

if(!require(“devtools”)) install.packages(“devtools”)
devtools::install_github(“apmcdermaid/GeneQC/GeneQC”)

2.	Load GeneQC library
library(GeneQC)

3.	Read in GeneQC Feature Extraction output file

dat <- read.table(“<GeneQC Feature Extraction output source.txt>”, header = TRUE, sep = “\t”) 

4.	Perform GeneQC modeling

dat_mod <- GeneQC(dat)

5.	Save GeneQC results file as CSV

write.csv(dat_mod, file = “GeneQC_results.csv”)
