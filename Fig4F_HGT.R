#####包括抗性基因的及其下游500bp序列的滑窗基因组距离###
mls_mge_contigs1=mls_mge_contigs#[!mls_mge_contigs$acc %in% arg_acc2$acc[arg_acc2$label1=="plasmid"],]

mls_mge_gene_500bp=data.frame(data.frame(matrix(NA,ncol = length(mls_mge_contigs1$acc),nrow = 3301)))
colnames(mls_mge_gene_500bp)=mls_mge_contigs1$acc

#col_num=c()
for(i in 1:length(mls_mge_contigs1$acc)){
  seq=fasta_matrix$mls_contigs2.fasta[,mls_mge_contigs1$acc[i]]
  k=mls_mge_contigs1$MSRD_send_arg[i]
  if(k==max(mls_mge_contigs1[i,c("MSRD_sstart_arg","MSRD_send_arg","MEFA_sstart_arg","MEFA_send_arg")])){
    if((mls_mge_contigs1$MSRD_acc_len[i]-k)>500){
      k1=mls_mge_contigs1[i,"MEFA_sstart_arg"]
      mls_mge_gene_500bp[,i]=seq[k1:(k+500)]##gene+downstream 500bp
      #col_num[i]=length(seq[k1:(k+500)])
    }
  }else if(k==min(mls_mge_contigs1[i,c("MSRD_sstart_arg","MSRD_send_arg","MEFA_sstart_arg","MEFA_send_arg")])){
    if(k>500){
      print(i)
      k1=mls_mge_contigs1[i,"MEFA_sstart_arg"]
      ####这里替换成互补序列
      mls_mge_gene_500bp[,i]=seq[k1:(k-500)]
      mls_mge_gene_500bp[,i]=unlist(strsplit(complement(DNAString(paste(mls_mge_gene_500bp[,i],collapse = ""))) %>% as.character(),""))
      #col_num[i]=length(seq[k1:(k-500)])
    }
  }
}
table(col_num)

mls_mge_gene_500bp=mls_mge_gene_500bp[,!is.na(mls_mge_gene_500bp[1,])]

#apply(mls_mge_gene_500bp,2,function(x){length(which(is.na(x)==T))})
##滑窗计算基因组之间的距离
dist1=combn(names(mls_mge_gene_500bp),m = 2) %>% t() %>% as.data.frame()

dist1_res=data.frame(matrix(NA,0,3))
colnames(dist1_res)=c("contig_pair","similarity_scores","site")

for(i in 1:nrow(dist1)){
  print(i)
  seq1 <- DNAString(paste(mls_mge_gene_500bp[,dist1$V1[i]],collapse = ""))
  seq2 <- DNAString(paste(mls_mge_gene_500bp[,dist1$V2[i]],collapse = ""))
  window_size <- 50 
  similarity_scores <- numeric()
  for (j in 1:(length(seq1) - window_size + 1)) {
    window_seq1 <- unlist(strsplit(subseq(seq1, start = j, width = window_size) %>% as.character(),""))
    window_seq2 <- unlist(strsplit(subseq(seq2, start = j, width = window_size) %>% as.character(),""))
    
    similarity_scores[j] <- length(which(window_seq1 == window_seq2)) / window_size
  }
  dist1_res1=matrix(NA,nrow = length(similarity_scores),ncol = 3)
  colnames(dist1_res1)=c("contig_pair","similarity_scores","site")
  dist1_res1=as.data.frame(dist1_res1)
  dist1_res1$similarity_scores=similarity_scores
  dist1_res1$contig_pair=paste(dist1$V1[i],dist1$V2[i],sep = "+")
  dist1_res1$site=1:(length(seq1) - window_size + 1)
  
  dist1_res=rbind(dist1_res,dist1_res1)
}
save(dist1_res,file = "dist1_res.RData")
###
# 计算均值的置信区间

median_func <- function(data, indices) {
  median(data[indices])
}

mean_ci <- function(data) {#, conf.level = 0.95
  # #Q1 <- quantile(data, 0.25)
  # #Q3 <- quantile(data, 0.75)
  # #IQR <- Q3 - Q1
  # #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  # n <- length(data)
  # mean <- mean(data)
  # stderr <- sd(data) / sqrt(n)
  # error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  # return(mean - error_margin)
  
  # 使用引导法计算置信区间
  results <- boot(data, statistic = median_func, R = 1000)
  ci <- boot.ci(results, type = "perc")
  return(ci$percent[4])
}

mean_ci1 <- function(data) {#, conf.level = 0.95
  # #Q1 <- quantile(data, 0.25)
  # #Q3 <- quantile(data, 0.75)
  # #IQR <- Q3 - Q1
  # #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  # n <- length(data)
  # mean <- mean(data)
  # stderr <- sd(data) / sqrt(n)
  # error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  # return(mean - error_margin)
  
  # 使用引导法计算置信区间
  results <- boot(data, statistic = median_func, R = 1000)
  ci <- boot.ci(results, type = "perc")
  return(ci$percent[5])
}

data <- lapply(unique(dist1_res$site),function(x){
  value=dist1_res$similarity_scores[dist1_res$site==x]
  filtered_value=value
  mean( filtered_value)
  }) %>% unlist

data=data.frame(site=unique(dist1_res$site),mean=data)

tmp=lapply(unique(data$site),function(x){print(x);mean_ci(data = dist1_res$similarity_scores[dist1_res$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_lower=tmp)
data=left_join(data,tmp)

tmp=lapply(unique(data$site),function(x){print(x);mean_ci1(data = dist1_res$similarity_scores[dist1_res$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_upper=tmp)
data=left_join(data,tmp)

p <- ggplot() +
  geom_line(data = dist1_res,aes(x = site, y = similarity_scores,color=contig_pair),alpha=0.05) +
  geom_line(data = data, aes(x = site, y = mean),color="blue")+
  ylim(0,1)+xlab("Genomic position")+ylab("Sequence similarity")+
  ggtitle("TnpA-mefA-msrD-carrying contigs")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1217, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1337, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 2800, linetype = "dashed", color = "grey50") +
  #scale_color_simpsons()+
  scale_color_manual(values = rep("black",780))+
  geom_ribbon(data = data,aes(x = site, y = mean,ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    legend.position = "none",
    strip.text = element_text(size = 10))
p

##############################################################
##############without MGE#####################################
mls_other_contigs1=mls_other_contigs#[!mls_other_contigs$acc %in% arg_acc2$acc[arg_acc2$label1=="plasmid"],]
#mls_other_contigs1=mls_other_contigs1[sample(1:nrow(mls_other_contigs1),59),]

mls_other_gene_500bp=data.frame(data.frame(matrix(NA,ncol = length(mls_other_contigs1$acc),nrow = 3301)))
colnames(mls_other_gene_500bp)=mls_other_contigs1$acc

#col_num=c()
for(i in 1:length(mls_other_contigs1$acc)){
  if(i %in% which(col_num==3301)){
    seq=fasta_matrix$mls_contigs1.fasta[,mls_other_contigs1$acc[i]]
    k=mls_other_contigs1$MSRD_send_arg[i]
    if(k==max(mls_other_contigs1[i,c("MSRD_sstart_arg","MSRD_send_arg","MEFA_sstart_arg","MEFA_send_arg")])){
      if((mls_other_contigs1$MSRD_acc_len[i]-k)>500){
        k1=mls_other_contigs1[i,"MEFA_sstart_arg"]
        mls_other_gene_500bp[,i]=seq[k1:(k+500)]##gene+downstream 500bp
        #col_num[i]=length(seq[k1:(k+500)])
      }
    }else if(k==min(mls_other_contigs1[i,c("MSRD_sstart_arg","MSRD_send_arg","MEFA_sstart_arg","MEFA_send_arg")])){
      if(k>500){
        print(i)
        k1=mls_other_contigs1[i,"MEFA_sstart_arg"]
        ####这里替换成互补序列
        mls_other_gene_500bp[,i]=seq[k1:(k-500)]
        mls_other_gene_500bp[,i]=unlist(strsplit(complement(DNAString(paste(mls_other_gene_500bp[,i],collapse = ""))) %>% as.character(),""))
        #col_num[i]=length(seq[k1:(k-500)])
      }
    }
  }
  
}
table(col_num)

mls_other_gene_500bp=mls_other_gene_500bp[,!is.na(mls_other_gene_500bp[1,])]

dist2=combn(names(mls_other_gene_500bp),m = 2) %>% t() %>% as.data.frame()
#dist2=dist2[sample(1:nrow(dist2),500),]

dist2_res=data.frame(matrix(NA,0,3))
colnames(dist2_res)=c("contig_pair","similarity_scores","site")

for(i in 1:nrow(dist2)){
  print(i)
  seq1 <- DNAString(paste(mls_other_gene_500bp[,dist2$V1[i]],collapse = ""))
  seq2 <- DNAString(paste(mls_other_gene_500bp[,dist2$V2[i]],collapse = ""))
  window_size <- 50 
  similarity_scores <- numeric()
  for (j in 1:(length(seq1) - window_size + 1)) {
    window_seq1 <- unlist(strsplit(subseq(seq1, start = j, width = window_size) %>% as.character(),""))
    window_seq2 <- unlist(strsplit(subseq(seq2, start = j, width = window_size) %>% as.character(),""))
    
    similarity_scores[j] <- length(which(window_seq1 == window_seq2)) / window_size
  }
  dist2_res1=matrix(NA,nrow = length(similarity_scores),ncol = 3)
  colnames(dist2_res1)=c("contig_pair","similarity_scores","site")
  dist2_res1=as.data.frame(dist2_res1)
  dist2_res1$similarity_scores=similarity_scores
  dist2_res1$contig_pair=paste(dist2$V1[i],dist2$V2[i],sep = "+")
  dist2_res1$site=1:(length(seq1) - window_size + 1)
  
  dist2_res=rbind(dist2_res,dist2_res1)###储存
}
save(dist2_res,file = "dist2_res.RData")
###
# 计算均值的置信区间
mean_ci <- function(data, conf.level = 0.95) {
  #Q1 <- quantile(data, 0.25)
  #Q3 <- quantile(data, 0.75)
  #IQR <- Q3 - Q1
  #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  n <- length(data)
  mean <- mean(data)
  stderr <- sd(data) / sqrt(n)
  error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  return(mean - error_margin)
}
mean_ci1 <- function(data, conf.level = 0.95) {
  #Q1 <- quantile(data, 0.25)
  #Q3 <- quantile(data, 0.75)
  #IQR <- Q3 - Q1
  #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  n <- length(data)
  mean <- mean(data)
  stderr <- sd(data) / sqrt(n)
  error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  return(mean + error_margin)
}

data <- lapply(unique(dist2_res$site),function(x){
  value=dist2_res$similarity_scores[dist2_res$site==x]
  filtered_value=value
  mean( filtered_value)
}) %>% unlist

data=data.frame(site=unique(dist2_res$site),mean=data)

tmp=lapply(unique(data$site),function(x){mean_ci(data = dist2_res$similarity_scores[dist2_res$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_lower=tmp)
data=left_join(data,tmp)

tmp=lapply(unique(data$site),function(x){mean_ci1(data = dist2_res$similarity_scores[dist2_res$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_upper=tmp)
data=left_join(data,tmp)

p <- ggplot() +
  geom_line(data = dist2_res,aes(x = site, y = similarity_scores,color=contig_pair),alpha=0.05) +
  geom_line(data = data, aes(x = site, y = mean),color="blue")+
  ylim(0,1)+xlab("Genomic position")+ylab("Sequence similarity")+
  ggtitle("mefA-msrD-carrying contigs")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1217, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1337, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 2800, linetype = "dashed", color = "grey50") +
  #scale_color_simpsons()+
  scale_color_manual(values = rep("black",30381))+
  geom_ribbon(data = data,aes(x = site, y = mean,ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    legend.position = "none",
    strip.text = element_text(size = 10))
p

#############################################################
#####包括抗性基因的及其下游500bp序列的滑窗基因组距离###
amin_mge_contigs1=amin_mge_contigs#[!amin_mge_contigs$acc %in% arg_acc2$acc[arg_acc2$label1=="plasmid"],]

amin_mge_gene_500bp=data.frame(data.frame(matrix(NA,ncol = length(amin_mge_contigs1$acc),nrow = 1940)))
colnames(amin_mge_gene_500bp)=amin_mge_contigs1$acc

#col_num=c()
for(i in 1:length(amin_mge_contigs1$acc)){
  if(i %in% which(col_num==1940)){
    seq=fasta_matrix$amin_contigs1.fa[,amin_mge_contigs1$acc[i]]
    k=amin_mge_contigs1$AAC6_send_arg[i]
    if(k==max(amin_mge_contigs1[i,c("AAC6_sstart_arg","AAC6_send_arg","APH2_sstart_arg","APH2_send_arg")])){
      if((amin_mge_contigs1$AAC6_acc_len[i]-k)>500){
        k1=amin_mge_contigs1[i,"AAC6_sstart_arg"]
        amin_mge_gene_500bp[,i]=seq[k1:(k+500)]##gene+downstream 500bp
        #col_num[i]=length(seq[k1:(k+500)])
      }
    }else if(k==min(amin_mge_contigs1[i,c("AAC6_sstart_arg","AAC6_send_arg","APH2_sstart_arg","APH2_send_arg")])){
      if(k>500){
        print(i)
        k1=amin_mge_contigs1[i,"AAC6_sstart_arg"]
        ####这里替换成互补序列
        amin_mge_gene_500bp[,i]=seq[k1:(k-500)]
        amin_mge_gene_500bp[,i]=unlist(strsplit(complement(DNAString(paste(amin_mge_gene_500bp[,i],collapse = ""))) %>% as.character(),""))
        #col_num[i]=length(seq[k1:(k-500)])
      }
    }
  }
  
}

table(col_num)

amin_mge_gene_500bp=amin_mge_gene_500bp[,!is.na(amin_mge_gene_500bp[1,])]

#apply(amin_mge_gene_500bp,2,function(x){length(which(is.na(x)==T))})
##滑窗计算基因组之间的距离
dist1=combn(names(amin_mge_gene_500bp),m = 2) %>% t() %>% as.data.frame()

dist1_res_amin=data.frame(matrix(NA,0,3))
colnames(dist1_res_amin)=c("contig_pair","similarity_scores","site")

for(i in 1:nrow(dist1)){
  print(i)
  seq1 <- DNAString(paste(amin_mge_gene_500bp[,dist1$V1[i]],collapse = ""))
  seq2 <- DNAString(paste(amin_mge_gene_500bp[,dist1$V2[i]],collapse = ""))
  window_size <- 50 
  similarity_scores <- numeric()
  for (j in 1:(length(seq1) - window_size + 1)) {
    window_seq1 <- unlist(strsplit(subseq(seq1, start = j, width = window_size) %>% as.character(),""))
    window_seq2 <- unlist(strsplit(subseq(seq2, start = j, width = window_size) %>% as.character(),""))
    
    similarity_scores[j] <- length(which(window_seq1 == window_seq2)) / window_size
  }
  dist1_res_amin1=matrix(NA,nrow = length(similarity_scores),ncol = 3)
  colnames(dist1_res_amin1)=c("contig_pair","similarity_scores","site")
  dist1_res_amin1=as.data.frame(dist1_res_amin1)
  dist1_res_amin1$similarity_scores=similarity_scores
  dist1_res_amin1$contig_pair=paste(dist1$V1[i],dist1$V2[i],sep = "+")
  dist1_res_amin1$site=1:(length(seq1) - window_size + 1)
  
  dist1_res_amin=rbind(dist1_res_amin,dist1_res_amin1)###储存
}
save(dist1_res_amin,file = "dist1_res_amin.RData")

dist1_res_amin$contig1=lapply(dist1_res_amin$contig_pair,function(x){strsplit(x,split="\\+") %>% unlist() %>% .[1]}) %>% unlist
dist1_res_amin$contig2=lapply(dist1_res_amin$contig_pair,function(x){strsplit(x,split="\\+") %>% unlist() %>% .[2]}) %>% unlist

dist1_res_amin$contig11=lapply(dist1_res_amin$contig1,function(x){arg_acc2$taxa[arg_acc2$acc==x] %>% unique()}) %>% unlist
dist1_res_amin$contig22=lapply(dist1_res_amin$contig2,function(x){arg_acc2$taxa[arg_acc2$acc==x] %>% unique()}) %>% unlist

dist1_res_amin1=dist1_res_amin[!(dist1_res_amin$contig11 == dist1_res_amin$contig22),]

p <- ggplot() +
  geom_line(data = dist1_res_amin1,aes(x = site, y = similarity_scores,color=contig_pair),alpha=0.05) +
  geom_line(data = data, aes(x = site, y = mean),color="blue")+
  ylim(0,1)+xlab("Genomic position")+ylab("Sequence similarity")+
  ggtitle("TnpA-aac(6')-aph(2'')-carrying contigs")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1439, linetype = "dashed", color = "grey50") +
  #scale_color_simpsons()+
  scale_color_manual(values = rep("black",68))+
  geom_ribbon(data = data,aes(x = site, y = mean,ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    legend.position = "none",
    strip.text = element_text(size = 10))
p

###
# 计算均值的置信区间
mean_ci <- function(data) {
  conf.level = 0.95
  n <- length(data)
  mean <- mean(data)
  stderr <- sd(data) / sqrt(n)
  error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  return(mean - error_margin)
  
}

mean_ci1 <- function(data) {
  conf.level = 0.95
  n <- length(data)
  mean <- mean(data)
  stderr <- sd(data) / sqrt(n)
  error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  return(mean + error_margin)
}

data <- lapply(unique(dist1_res_amin1$site),function(x){
  value=dist1_res_amin1$similarity_scores[dist1_res_amin1$site==x]
  #Q1 <- quantile(value, 0.25)
  #Q3 <- quantile(value, 0.75)
  #IQR <- Q3 - Q1
  filtered_value=value#[value >= (Q1 - 1.5 * IQR) & value <= (Q3 + 1.5 * IQR)]
  mean( filtered_value)
}) %>% unlist

data=data.frame(site=unique(dist1_res_amin1$site),mean=data)

tmp=lapply(unique(data$site),function(x){print(x);
  mean_ci(data = dist1_res_amin1$similarity_scores[dist1_res_amin1$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_lower=tmp)
data=left_join(data,tmp)

tmp=lapply(unique(data$site),function(x){
  print(x)
  mean_ci1(data = dist1_res_amin1$similarity_scores[dist1_res_amin1$site==x])}) %>% unlist

tmp=data.frame(site=unique(data$site),y_upper=tmp)
data=left_join(data,tmp)

p <- ggplot(data, aes(x = site, y = mean)) +
  geom_line() +ylim(0,1)+xlab("Genomic position")+ylab("Mean sequence similarity")+
  ggtitle("TnpA-aac(6')-aph(2'')-carrying contigs")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1439, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    strip.text = element_text(size = 10))
p

##############################################################
##############without MGE#####################################
amin_other_contigs1=amin_other_contigs#[!amin_other_contigs$acc %in% arg_acc2$acc[arg_acc2$label1=="plasmid"],]

amin_other_gene_500bp=data.frame(data.frame(matrix(NA,ncol = length(amin_other_contigs1$acc),nrow = 1940)))
colnames(amin_other_gene_500bp)=amin_other_contigs1$acc

#col_num=c()
for(i in 1:length(amin_other_contigs1$acc)){
  if(i %in% which(col_num==1940)){
    seq=fasta_matrix$amin_contigs2.fa[,amin_other_contigs1$acc[i]]
    k=amin_other_contigs1$AAC6_send_arg[i]
    if(k==max(amin_other_contigs1[i,c("AAC6_sstart_arg","AAC6_send_arg","APH2_sstart_arg","APH2_send_arg")])){
      if((amin_other_contigs1$AAC6_acc_len[i]-k)>500){
        k1=amin_other_contigs1[i,"AAC6_sstart_arg"]
        amin_other_gene_500bp[,i]=seq[k1:(k+500)]##gene+downstream 500bp
        #col_num[i]=length(seq[k1:(k+500)])
      }
    }else if(k==min(amin_other_contigs1[i,c("AAC6_sstart_arg","AAC6_send_arg","APH2_sstart_arg","APH2_send_arg")])){
      if(k>500){
        print(i)
        k1=amin_other_contigs1[i,"AAC6_sstart_arg"]
        ####这里替换成互补序列
        amin_other_gene_500bp[,i]=seq[k1:(k-500)]
        amin_other_gene_500bp[,i]=unlist(strsplit(complement(DNAString(paste(amin_other_gene_500bp[,i],collapse = ""))) %>% as.character(),""))
        #col_num[i]=length(seq[k1:(k-500)])
      }
    }
  }
  
}

#table(col_num)

amin_other_gene_500bp=amin_other_gene_500bp[,!is.na(amin_other_gene_500bp[1,])]

dist2=combn(names(amin_other_gene_500bp),m = 2) %>% t() %>% as.data.frame()
#dist2=dist2[sample(1:nrow(dist2),500),]

dist2_res_amin=data.frame(matrix(NA,0,3))
colnames(dist2_res_amin)=c("contig_pair","similarity_scores","site")

for(i in 1:nrow(dist2)){
  print(i)
  seq1 <- DNAString(paste(amin_other_gene_500bp[,dist2$V1[i]],collapse = ""))
  seq2 <- DNAString(paste(amin_other_gene_500bp[,dist2$V2[i]],collapse = ""))
  window_size <- 50 ###窗口不能太小
  similarity_scores <- numeric()
  for (j in 1:(length(seq1) - window_size + 1)) {
    window_seq1 <- unlist(strsplit(subseq(seq1, start = j, width = window_size) %>% as.character(),""))
    window_seq2 <- unlist(strsplit(subseq(seq2, start = j, width = window_size) %>% as.character(),""))
    
    similarity_scores[j] <- length(which(window_seq1 == window_seq2)) / window_size
  }
  dist2_res_amin1=matrix(NA,nrow = length(similarity_scores),ncol = 3)
  colnames(dist2_res_amin1)=c("contig_pair","similarity_scores","site")
  dist2_res_amin1=as.data.frame(dist2_res_amin1)
  dist2_res_amin1$similarity_scores=similarity_scores
  dist2_res_amin1$contig_pair=paste(dist2$V1[i],dist2$V2[i],sep = "+")
  dist2_res_amin1$site=1:(length(seq1) - window_size + 1)
  
  dist2_res_amin=rbind(dist2_res_amin,dist2_res_amin1)###储存
}
save(dist2_res_amin,file = "dist2_res_amin.RData")

dist2_res_amin$contig1=lapply(dist2_res_amin$contig_pair,function(x){strsplit(x,split="\\+") %>% unlist() %>% .[1]})
dist2_res_amin$contig2=lapply(dist2_res_amin$contig_pair,function(x){strsplit(x,split="\\+") %>% unlist() %>% .[2]}) 

dist2_res_amin1=dist2_res_amin[!(dist2_res_amin$contig1 %in% t2$acc[t2$taxa=="Campylobacter concisus"] & dist2_res_amin$contig2 %in% t2$acc[t2$taxa=="Campylobacter concisus"]),]

p <- ggplot() +
  geom_line(data = dist2_res_amin1,aes(x = site, y = similarity_scores,group=contig_pair),alpha=0.05) +
  geom_line(data = data, aes(x = site, y = mean),color="blue")+
  ylim(0,1)+xlab("Genomic position")+ylab("Sequence similarity")+
  ggtitle("aac(6')-aph(2'')-carrying contigs\n(without MGE upstreaming)")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1439, linetype = "dashed", color = "grey50") +
  #scale_color_simpsons()+
  #scale_color_manual(values = rep("black",21))+
  geom_ribbon(data = data,aes(x = site, y = mean,ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    legend.position = "none",
    strip.text = element_text(size = 10))
p

###
# 计算均值的置信区间
mean_ci <- function(data, conf.level = 0.95) {
  Q1 <- quantile(data, 0.35)
  #Q3 <- quantile(data, 0.75)
  #IQR <- Q3 - Q1
  #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  # n <- length(data)
  # mean <- mean(data)
  # stderr <- sd(data) / sqrt(n)
  # error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  # return(mean - error_margin)
  return(Q1)
  #return(median(data)-sd(data))
}
mean_ci1 <- function(data, conf.level = 0.95) {
  #Q1 <- quantile(data, 0.25)
  Q3 <- quantile(data, 0.65)
  #IQR <- Q3 - Q1
  #data=data[data >= (Q1 - 1.5 * IQR) & data <= (Q3 + 1.5 * IQR)]
  # n <- length(data)
  # mean <- mean(data)
  # stderr <- sd(data) / sqrt(n)
  # error_margin <- qt(1 - (1 - conf.level) / 2, df = n - 1) * stderr
  # return(mean + error_margin)
  return(Q3)
  #return(median(data)+sd(data))
}

data <- lapply(unique(dist2_res_amin1$site),function(x){
  value=dist2_res_amin1$similarity_scores[dist2_res_amin1$site==x]
  #Q1 <- quantile(value, 0.25)
  #Q3 <- quantile(value, 0.75)
  #IQR <- Q3 - Q1
  filtered_value=value#[value >= (Q1 - 1.5 * IQR) & value <= (Q3 + 1.5 * IQR)]
  mean( filtered_value)
}) %>% unlist

data=data.frame(site=unique(dist2_res_amin1$site),mean=data)

tmp=lapply(unique(data$site),function(x){
  print(x);mean_ci(data = dist2_res_amin1$similarity_scores[dist2_res_amin1$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_lower=tmp)
data=left_join(data,tmp)

tmp=lapply(unique(data$site),function(x){print(x);mean_ci1(data = dist2_res_amin1$similarity_scores[dist2_res_amin1$site==x])}) %>% unlist
tmp=data.frame(site=unique(data$site),y_upper=tmp)
data=left_join(data,tmp)

p <- ggplot(data, aes(x = site, y = median)) +
  geom_line() +ylim(0,1)+xlab("Genomic position")+ylab("Median sequence similarity")+
  ggtitle("aac(6')-aph(2'')-carrying contigs\n(without MGE upstreaming)")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1439, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2)+
  theme_bw()+
  theme(
    axis.title=element_text(size = 10),
    axis.text.x= element_text(size = 10,angle = 0),
    axis.text.y= element_text(size = 10),
    strip.text = element_text(size = 10))
p









