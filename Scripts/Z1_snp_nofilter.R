mydisc <- function(chr1,pos1,chr2,pos2) {
  if(chr1 != chr2) {
    return(Inf)
  }else{
    return(abs(as.numeric(pos1)-as.numeric(pos2)))
  }
}

tf2groupName <- function(x) {
  ans <- vector(length=length(x))
  k <- 1 
  ans[1] <- 1 
  for(i in 2:length(x)) {
    if(x[i] == TRUE) {
      k <- k+1
      ans[i] <- k
    }else{
      ans[i] <- ans[i-1]
    }
  }
  return(ans)
}


setwd('~/cloud/yeastAse_redo/001_snpPos')

yr_pos <- read.table(header=T,file='~/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/yr5B.posmap')
ry_pos <- read.table(header=T,file='~/mnt/nnp/1_1_0_cali_bowtie/0_5Bsnp/ryB5.posmap')

names(yr_pos)[1:6] <- c('ypsChrom','ypsPosit','ypsN','rmChrom','rmPosit','rmN')
names(ry_pos)[1:6] <- c('rmChrom','rmPosit','rmN','ypsChrom','ypsPosit','ypsN')


# remove mummer unmatch snps   
ry_pos_exchg <- ry_pos[,c(4,5,6,1,2,3,7)]
res <- merge.data.frame(yr_pos,ry_pos_exchg,by.x = c(1,2,4,5),by.y=c(1,2,4,5),all=T,sort=F)
res_good <- res[which(res$drct.x == res$drct.y),] # 86351


# combine snps in one segment into one group 
yps_block <- res_good[,1:4] %>% arrange(ypsChrom,ypsPosit)
yps_block_tf <- vector(length=nrow(yps_block))
yps_block_tf[1] <- TRUE  
for(i in 2:nrow(yps_block)) {
  yps_block_tf[i] <- 
    mydisc(yps_block[i-1,"ypsChrom"], yps_block[i-1,"ypsPosit"],
           yps_block[i,"ypsChrom"], yps_block[i,"ypsPosit"]) > 700000
} 

yps_rm_86351_group <- cbind(yps_block,yps_block_tf,gN= tf2groupName(yps_block_tf))

write.table(unique(yps_rm_86351_group[,c(1,2)]),file="~/cloud/yeastAse_redo/001_snpPos/Z1_out/yps128_5_snpls_nofilter",row.names = F,quote=F,col.names = F)

write.table(unique(yps_rm_86351_group[,c(3,4)]),file="~/cloud/yeastAse_redo/001_snpPos/Z1_out/rm11_B_snpls_nofilter",row.names = F,quote=F,col.names = F)

#write.table(unique(yps_rm_86351_group[,c(1,2,6)]),file="~/cloud/yeastAse_redo/001_snpPos/Z1_out/yps128_5_snpls_nofilter_group",row.names = F,quote=F,col.names = F)

#write.table(unique(yps_rm_86351_group[,c(3,4,6)]),file="~/cloud/yeastAse_redo/001_snpPos/Z1_out/rm11_B_snpls_nofilter_group",row.names = F,quote=F,col.names = F)


# combine yps same row. ############## position match of rm and yps/
combn <- function(df,strn) {
    nr <- nrow(df)
    mpool <- df[1,]
    ans  <- data.frame()
    if(strn=='yps'){pos <- c(1,2)}else{pos <- c(3,4)}
    for (i in 2:nr) {
        #print(pos)
        #print(df[i,pos])
        #print(mpool)
        if(all(df[i,pos] == mpool[1,pos])){
            #print("true the same")
            mpool <- rbind(mpool,df[i,]) 
        }else{
            # deal with mpool
            ans <- rbind(ans,mpool[ceiling(nrow(mpool)/2),])
            # re put mpool
            mpool <- df[i,]
        }
    }
    ans <- rbind(ans,mpool[ceiling(nrow(mpool)/2),])
    return(ans)
}  


library(dplyr)
res4 <- res_good %>% arrange(ypsChrom,ypsPosit)
res5 <- combn(res4,"yps")
res6 <- res5 %>% arrange(rmChrom,rmPosit)
res7 <- combn(res6,"rm")

rm(list=grep("res[456]",ls(),value=T))


posmap52 <- res7[,c(1,2,3,4,5,6,8,9,10)]
names(posmap52)[9] <- "drct"
save(posmap52,file = "~/cloud/yeastAse_redo/001_snpPos/Z1_out/posmap52.RData")
