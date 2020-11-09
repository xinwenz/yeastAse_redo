rm(list=ls())
source("~/cloud/yeastAse_redo/Programs/cis_prmtr_func.R")
load("~/cloud/yeastAse_redo/Z2_out/exph.RData")


# depth
exph_Ct_esum <- colSums(exph[,-1])

# save this dataset and functions 
save(list=ls(),file = '~/cloud/yeastAse_redo/Z3_out/exph_functions.RData')

# random subset pick
sink('~/cloud/yeastAse_redo/Z3_out/exph_random_subset.txt')
i <- 1 
tmplt <- (1:20)[-c(7,19)]
for(siz in 1:18) {
    for(t in 1:25) {
        choRep <- sample(tmplt,size=siz,replace = F)
        input <- c(choRep,choRep + 20)
        cat(c(formatC(i, width = 3, format = "d", flag = "0"),input,'\n'))
        i <- i + 1 
    }
}
sink()