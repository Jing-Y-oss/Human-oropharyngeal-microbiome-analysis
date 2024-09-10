library(circlize)

View(species_ch_res)

species_ch_res$p.adj=p.adjust(species_ch_res$p,method = "fdr")
species_ch_res$p_spouse.adj=p.adjust(species_ch_res$p_spouse,method = "fdr")
species_ch_res$p_silbling.adj=p.adjust(species_ch_res$p_silbling,method = "fdr")
species_ch_res$p_grantchild.adj=p.adjust(species_ch_res$p_grantchild,method = "fdr")
species_ch_res$p_parentchild.adj=p.adjust(species_ch_res$p_parentchild,method = "fdr")

sig_taxa=Reduce(union,list(species_ch_res$species[species_ch_res$p.adj<0.05],
                               species_ch_res$species[species_ch_res$p_grantchild<0.05],
                               species_ch_res$species[species_ch_res$p_parentchild<0.05],
                               species_ch_res$species[species_ch_res$p_spouse<0.05],
                               species_ch_res$species[species_ch_res$p_silbling<0.05])) %>% na.omit() %>% as.character() %>% sort
sig_taxa[sig_taxa=="Caudoviricetes sp."]="Siphoviridae sp."
sig_taxa[sig_taxa=="Cutibacterium modestum"]="Propionibacterium sp. oral taxon 193"
tmp=species_rela[sig_taxa,]
sig_taxa=which(apply(tmp,1,function(x){length(which(x>0))})>889) %>% rownames(tmp)[.]#在85%的样本中不为0
sig_taxa=setdiff(sig_taxa,"uncultured phage")

sig_taxa_phylum=read.csv("sig_taxa.csv",header = T)
sig_taxa_phylum$species=apply(sig_taxa_phylum,1,function(x){strsplit(x[1],split = ":") %>% unlist %>% .[1] %>%  gsub(" $","",.)})
sig_taxa_phylum$phylum=apply(sig_taxa_phylum,1,function(x){strsplit(x[1],split = ":") %>% unlist %>% .[2] %>%  gsub(" ","",.)})

sig_taxa_phylum=sig_taxa_phylum[sig_taxa_phylum$species %in% sig_taxa,]
sig_taxa_phylum=sig_taxa_phylum[order(sig_taxa_phylum$phylum),]


###get tree

# library(rotl)
sig_taxa_phylum$species[sig_taxa_phylum$species=="Propionibacterium sp. oral taxon 193"]="Cutibacterium modestum"
sig_taxa_phylum$species[sig_taxa_phylum$species=="Siphoviridae sp."]="Caudoviricetes sp."
# sig_taxa_phylum$species[sig_taxa_phylum$species=="Siphoviridae sp."]="Siphoviridae"
# taxa<-tnrs_match_names(names= sig_taxa_phylum$species)
# 
# tree <- tol_induced_subtree(ott_ids = ott_id(taxa))
# 
species_ch_res$species[species_ch_res$species=="Siphoviridae sp."]="Caudoviricetes sp."
species_ch_res$species[species_ch_res$species=="Propionibacterium sp. oral taxon 193"]="Cutibacterium modestum"
rownames(species_ch_res)=species_ch_res$species %>% lapply(.,function(x){gsub(" ","_",x)}) %>% unlist

#######
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(reshape2)
library(dplyr)
library(tidytree)
library(phyloseq)

sig_taxa_phylum$species=lapply(sig_taxa_phylum$species,function(x){gsub(" ","_",x)}) %>% unlist
rownames(sig_taxa_phylum)=sig_taxa_phylum$species
sig_taxa_phylum=sig_taxa_phylum[tree$tip.label,]
sig_taxa_phylum$id=1:42


###phylum name
install.packages("taxizedb")
library(taxizedb)
ids <- name2taxid(sig_taxa_phylum$species, out_type="summary")
all(ids$name==sig_taxa_phylum$species)
sig_taxa_phylum$phylum_correct=NA
for(i in 1:42){
  sig_taxa_phylum$phylum_correct[i]= classification(ids$id)[i] %>% as.data.frame() %>% .[which(.[,2]=="phylum"),1]
}


tree <- read.tree("phyloT_generated_tree_1712760152_newick.txt")

circ <- ggtree(tree, layout = "circular",branch.length = "none")+
  geom_hilight(data=sig_taxa_phylum,aes(node=id,fill=phylum_correct), alpha=.6,align="both")+
  scale_fill_jama( name="Phylum")+geom_tiplab(offset=7.5,size=5)
circ

df <- data.frame(All=species_ch_res[tree$tip.label,"cor"])
df$All[which(species_ch_res[tree$tip.label,"p.adj"]>0.05)]=NA
df$All=as.numeric(df$All)
rownames(df) <- tree$tip.label

circ=circ+new_scale_fill()

p1 <- gheatmap(circ, df, offset=0.2,width=0.1,
               colnames_angle=95, colnames_offset_y = .25,color = "black")+#name="Correlation between\ncohabitants"
  scale_fill_distiller(type = "div",direction = 1,name="Correlation between\ncohabitants",na.value = "white")
p1

df1 <- data.frame(`Grantparent-child`=species_ch_res[tree$tip.label,"cor_grantchild"])
df1$Grantparent.child[which(species_ch_res[tree$tip.label,"p_grantchild"]>0.05)]=NA
df1$Grantparent.child=as.numeric(df1$Grantparent.child)
rownames(df1) <- tree$tip.label

p1 <-gheatmap(p1, df1, offset=1.6, width=0.1,
              colnames_angle=90, colnames_offset_y = .25,color="#EFC000FF")+
  scale_fill_distiller(type = "div",direction = 1,name="Correlation between\ncohabitants",na.value = "white")
p1

df2 <- data.frame(`Parent-child`=species_ch_res[tree$tip.label,"cor_parentchild"])
df2$Parent.child[which(species_ch_res[tree$tip.label,"p_parentchild"]>0.05)]=NA
df2$Parent.child=as.numeric(df2$Parent.child)
rownames(df2) <- tree$tip.label


p2 <-gheatmap(p1, df2, offset=3, width=0.1,
              colnames_angle=90, colnames_offset_y = .25,color="#868686FF")+
  scale_fill_distiller(type = "div",direction = 1,name="Correlation between\ncohabitants",na.value = "white")
p2


df3 <- data.frame(Spouse=species_ch_res[tree$tip.label,"cor_spouse"])
df3$Spouse[which(species_ch_res[tree$tip.label,"p_spouse"]>0.05)]=NA
df3$Spouse=as.numeric(df3$Spouse)
rownames(df3) <- tree$tip.label


p3 <-gheatmap(p2, df3, offset=4.4, width=0.1,
              colnames_angle=90, colnames_offset_y = .25,color="#CD534CFF")+
  scale_fill_distiller(type = "div",direction = 1,name="Correlation between\ncohabitants",na.value = "white")
p3


df4 <- data.frame(Silbling=species_ch_res[tree$tip.label,"cor_silbling"])
df4$Silbling[which(species_ch_res[tree$tip.label,"p_silbling"]>0.05)]=NA
df4$Silbling=as.numeric(df4$Silbling)
rownames(df4) <- tree$tip.label


p4 <-gheatmap(p3, df4, offset=5.8, width=0.1,
              colnames_angle=90, colnames_offset_y = .25,color="#8F7700FF")+
  scale_fill_distiller(type = "div",direction = 1,name="Correlation between\ncohabitants",na.value = "white")
p4

##########version 2.0
tree$tip.label=lapply(tree$tip.label,function(x){gsub("_"," ",x)}) %>% unlist
species_ch_res$species=lapply(species_ch_res$species,function(x){gsub("_"," ",x)}) %>% unlist


circ <- ggtree(tree, layout = "circular",branch.length = "none")+
  geom_hilight(node=111, fill="steelblue", alpha=0.5) +
  geom_hilight(node=90, fill="darkgreen", alpha=0.5) +
  geom_hilight(node=70, fill="gray", alpha=0.5) +
  geom_hilight(node=53, fill="pink", alpha=0.5) +
  geom_hilight(node=62, fill="beige", alpha=0.5) +
  geom_hilight(node=84, fill="yellow", alpha=0.5)+
  geom_hilight(node=44, fill="skyblue", alpha=0.5)+
  #geom_hilight(data=sig_taxa_phylum,aes(node=id,fill=phylum_correct),alpha=0.6,align = "both")+
  #scale_fill_jama(name="Phylum")+
  geom_tiplab(offset=10.5,size=5)
circ <- circ+new_scale_fill()

df <- data.frame(All=species_ch_res[tree$tip.label,"cor"])
df$sig="*"
df$sig[which(species_ch_res[tree$tip.label,"p.adj"]>0.05)]=""
df$All=as.numeric(df$All)
df$spe <- tree$tip.label

df1 <- data.frame(`Grantparent-child`=species_ch_res[tree$tip.label,"cor_grantchild"])
df1$sig="*"
df1$sig[which(species_ch_res[tree$tip.label,"p_grantchild"]>0.05)]=""
df1$Grantparent.child=as.numeric(df1$Grantparent.child)
df1$spe <- tree$tip.label

df2 <- data.frame(`Parent-child`=species_ch_res[tree$tip.label,"cor_parentchild"])
df2$sig="*"
df2$sig[which(species_ch_res[tree$tip.label,"p_parentchild"]>0.05)]=""
df2$Parent.child=as.numeric(df2$Parent.child)
df2$spe <- tree$tip.label

df3 <- data.frame(Spouse=species_ch_res[tree$tip.label,"cor_spouse"])
df3$sig="*"
df3$sig[which(species_ch_res[tree$tip.label,"p_spouse"]>0.05)]=""
df3$Spouse=as.numeric(df3$Spouse)
df3$spe <- tree$tip.label

df4 <- data.frame(Silbling=species_ch_res[tree$tip.label,"cor_silbling"])
df4$sig="*"
df4$sig[which(species_ch_res[tree$tip.label,"p_silbling"]>0.05)]=""
df4$Silbling=as.numeric(df4$Silbling)
df4$spe <- tree$tip.label

circ+
  geom_fruit(data=df,geom = geom_tile,aes(y=spe,fill=All),color="black",
             offset=0.15,width=1.2,lwd = 0.4)+
  geom_fruit(data = df,geom = geom_text,aes(y=spe,x=0,label=sig))+
  geom_fruit(data=df1,geom = geom_tile,aes(y=spe,fill=Grantparent.child),color="#EFC000FF",
             offset=0.15,width=1.2,lwd = 0.7)+
  geom_fruit(data = df1,geom = geom_text,aes(y=spe,x=0,label=sig))+
  geom_fruit(data=df2,geom = geom_tile,aes(y=spe,fill=Parent.child),color="#868686FF",
             offset=0.15,width=1.2,lwd = 0.7)+
  geom_fruit(data = df2,geom = geom_text,aes(y=spe,x=0,label=sig))+
  geom_fruit(data=df3,geom = geom_tile,aes(y=spe,fill=Spouse),color="#CD534CFF",
             offset=0.15,width=1.2,lwd = 0.7)+
  geom_fruit(data = df3,geom = geom_text,aes(y=spe,x=0,label=sig))+
  geom_fruit(data=df4,geom = geom_tile,aes(y=spe,fill=Silbling),color="#8F7700FF",
             offset=0.15,width=1.2,lwd = 0.7)+
  geom_fruit(data = df4,geom = geom_text,aes(y=spe,x=0,label=sig))+
  scale_fill_steps2(low = "grey", mid = "white", high = "brown",  midpoint = 0)








