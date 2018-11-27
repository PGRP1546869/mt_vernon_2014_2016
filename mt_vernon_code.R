library(MASS)
library(devtools)
library(digest)
library(reshape)
library(extrafont)
library(ggpubr)
library(viridis)
library(agricolae)
library(ncf)
library(tidyverse)
library(maSigPro)

#morphometric figures

imagej <- read.table('mt_vernon_morpho_imagej_all_info.txt', header=T,sep="\t")

#keep relevant columns
imagej <- dplyr::select(imagej, circ, ar, round, solidity, rootstock, irrigation, block, year)

#now run the script to include which ones are significant

total_var_mat_p=matrix(,1,11)
colnames(total_var_mat_p)=c("rootstock_var", "irrigation_var","rootstock_irrigation_var", "year_var","block_var", "residual_var","rootstock_p", "irrigation_p","rootstock_irrigation_p", "year_p","block_p")

for (i in which(colnames(imagej)=="circ"):which(colnames(imagej)=="solidity")) {
  imagej_lm=lm(imagej[,i] ~ rootstock + irrigation+rootstock*irrigation+year+block, data = imagej)
  imagej_sums=anova(imagej_lm)$Sum
  total_var=sum(imagej_sums)
  test_total=(imagej_sums/total_var)*100
  imagej_p=anova(imagej_lm)$Pr[1:5]
  var_p=append(test_total,imagej_p)
  total_var_mat_p=rbind(total_var_mat_p,var_p)
  rownames(total_var_mat_p)[i+1]=paste(colnames(imagej)[i])
}

total_var_mat_p=total_var_mat_p[-1,]

total_var_mat_p <- as.data.frame(total_var_mat_p)
total_var_mat_p <- cbind(rownames(total_var_mat_p),total_var_mat_p)
colnames(total_var_mat_p)[1]="phenotype"
write.table(total_var_mat_p, "mt_vernon_imagej_var_p.txt", sep="\t", row.names=F, quote=F, col.names=T)

#then modify the table so it is easy to plot
total_var_mat_p <- read.table('mt_vernon_imagej_var_p.txt', header=T,sep="\t")
#remove residuals since I am not going to visualize it
total_var_mat_no_res_edit <- dplyr::select(total_var_mat_p,-residual_var)
total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("phenotype", "rootstock_var")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="var"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("phenotype", "irrigation_var")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="var"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("phenotype", "rootstock_irrigation_var")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="var"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("phenotype", "year_var")]
colnames(total_var_mat_no_res_edit_year)[2]="var"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("phenotype", "block_var")]
colnames(total_var_mat_no_res_edit_block)[2]="var"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")

total_var_mat_no_res_edit_all <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block)

#then do the same thing with p-values
total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("phenotype", "rootstock_p")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="p-value"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("phenotype", "irrigation_p")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="p-value"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("phenotype", "rootstock_irrigation_p")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="p-value"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("phenotype", "year_p")]
colnames(total_var_mat_no_res_edit_year)[2]="p-value"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("phenotype", "block_p")]
colnames(total_var_mat_no_res_edit_block)[2]="p-value"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")

total_var_mat_no_res_edit_all_p <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block)

var_p=cbind(total_var_mat_no_res_edit_all,total_var_mat_no_res_edit_all_p)
table(var_p[,1]==var_p[,4])
table(var_p[,3]==var_p[,6])
var_p <- var_p[,-4]
var_p <- var_p[,-5]
colnames(var_p)[4]="p_value"

#only plot significant p_values
var_p_sig=dplyr::filter(var_p, p_value<0.05)


write.table(var_p_sig, "mt_vernon_imagej_var_p_sig.txt", sep="\t", row.names=F, quote=F, col.names=T)


#there's no significant effect for block so I need an empty column there 
var_p_sig <- rbind(var_p_sig, c("solidity", NA, "block", NA))
var_p_sig[,"var"] <- as.numeric(as.character(var_p_sig[,"var"]))
var_p_sig[,"p_value"] <- as.numeric(as.character(var_p_sig[,"p_value"]))

var_p_sig[,"phenotype"] <- gsub("circ", "3_circularity", var_p_sig[,"phenotype"])
var_p_sig[,"phenotype"] <- gsub("solidity", "1_solidity", var_p_sig[,"phenotype"])
var_p_sig[,"phenotype"] <- gsub("round", "2_round", var_p_sig[,"phenotype"])
var_p_sig[,"phenotype"] <- gsub("\\<ar\\>", "4_aspect ratio", var_p_sig[,"phenotype"])


pdf("imagej_variance_heatmap_black_viridis_test.pdf",family="Arial", width=6, height=2.5)
t <- ggplot(data=var_p_sig, aes(x=trait, y=phenotype, fill=var))
t + geom_tile(color="black", size=0.2)+theme_minimal()+scale_fill_gradientn(colors=c("white", "black"),limits=c(0, 7))+labs(x = "Factor", y="Phenotype") + theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"),legend.position = "bottom")+ scale_fill_viridis(option="inferno", name="% variance explained", direction = -1,na.value =NA)+scale_x_discrete(position = "top")
dev.off()


#associated boxplots for rootstock and irrigation circularity plots 

rootstock_palette <- c("gray","#1b9e77", "#7570b3",  "#e6ab02")

circ_rootstock <- dplyr::select(imagej, circ, rootstock)
circ_rootstock$rootstock <- factor(circ_rootstock$rootstock,
                                   levels = c('OWN','1103P', '3309C', 'SO4'),ordered = TRUE)
circ_rootstock_plot <- ggplot(circ_rootstock,aes(x = rootstock, y = circ, fill = rootstock)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=rootstock_palette)+theme_bw()+theme(legend.position = "none")+labs(x="rootstock", y="Circularity") 

circ_irrigation <- dplyr::select(imagej, circ, irrigation)
circ_irrigation$irrigation <- factor(circ_irrigation$irrigation,
                                   levels = c('none','RDI', 'full'),ordered = TRUE)
circ_irrigation_plot <- ggplot(circ_irrigation,aes(x = irrigation, y = circ, fill = irrigation)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+theme_bw()+theme(legend.position = "none")+labs(x="irrigation", y="Circularity")+ scale_fill_brewer(direction=1)


dev.off()
require(cowplot)
pdf("imagej_figure.pdf", width=6, height=4,family="Arial")
plot_grid(circ_irrigation_plot,circ_rootstock_plot, ncol = 2,labels="AUTO")
dev.off()


#persistent homology figures
all_morpho <- read.table('mt_vernon_morpho_avg_two_years_all_info.txt', header=T,sep="\t")
#run PCA on all morphometric data together
#save ids which I can't include when doing PCA
all_morpho_id <- dplyr::select(all_morpho,row:block)
all_morpho <- dplyr::select(all_morpho, -(row:block))

#perform PCA, do not scale but do center
all_morpho_pca <- prcomp(all_morpho, center=T)
eigenvalues <- (all_morpho_pca$sdev)^2
pcs_explained <- eigenvalues/sum(eigenvalues)*100
pcs_explained <- cbind(pcs_explained, "sum")
colnames(pcs_explained)[2]="sum"
pcs_explained_prop <- pcs_explained
pcs_explained_prop <- cbind(1:nrow(pcs_explained_prop),pcs_explained_prop[,1])
colnames(pcs_explained_prop) <- c("PC", "PC_variance_explained")
pcs_explained_prop <- as.data.frame(pcs_explained_prop)
write.table(pcs_explained_prop, "mt_vernon_pcs_explained_prop.txt", sep="\t", row.names=F, quote=F, col.names=T)
pcs_explained[,"sum"]<-cumsum(pcs_explained[,"pcs_explained"])
#get total proportion of variance explained by first 20 PCs
pcs_explained <- pcs_explained[1:20, "pcs_explained"]
sum(as.numeric(pcs_explained))
length(pcs_explained)
#20 PCs explain 68.1329% of the variance

#save pcs 
all_morpho_pcs <- all_morpho_pca$x
all_morpho_pcs <- cbind(all_morpho_id,all_morpho_pcs)

write.table(all_morpho_pcs, "mt_vernon_morpho_merged_pcs.txt", sep="\t", row.names=F, quote=F, col.names=T)

#check how much variance is explained by rootstock, irrigation and rootstock by irrigation after we control for block and year. 

morpho_pcs <- read.table('mt_vernon_morpho_merged_pcs.txt', header=T,sep="\t")
morpho_pcs <- mutate(morpho_pcs, column_row =paste(column, row, sep="_"))
morpho_pcs <-mutate(morpho_pcs, rootstock_irr =paste(rootstock, irrigation, sep="_"))
morpho_pcs <- select(morpho_pcs, row:block, column_row:rootstock_irr, PC1:PC566)

#make sure all variables are factors
morpho_pcs[,"year"] <- as.factor(morpho_pcs[,"year"])
total_var_mat=matrix(,1,6)
colnames(total_var_mat)=c("rootstock_var", "irrigation_var","rootstock_irrigation_var", "year_var", "block_var", "residual_var")

for (i in which(colnames(morpho_pcs)=="PC1"):which(colnames(morpho_pcs)=="PC20")) {
  pc_lm=lm(morpho_pcs[,i] ~ rootstock + irrigation+rootstock*irrigation+year+block, data = morpho_pcs)
  pc_sums=anova(pc_lm)$Sum
  total_var=sum(pc_sums)
  test_total=(pc_sums/total_var)*100
  total_var_mat=rbind(total_var_mat,test_total)
  rownames(total_var_mat)[i-8]=paste("PC", i-9, sep="_")
}

total_var_mat=total_var_mat[-1,]

total_var_mat <- as.data.frame(total_var_mat)
#now to plot it what is the max variance explained by any of these factors, ignore residuals
total_var_mat_no_res <- dplyr::select(total_var_mat, -residual_var)

max(total_var_mat_no_res)
min(total_var_mat_no_res)

#now run the script to include which ones are significant

total_var_mat_p=matrix(,1,11)
colnames(total_var_mat_p)=c("rootstock_var", "irrigation_var","rootstock_irrigation_var", "year_var","block_var", "residual_var","rootstock_p", "irrigation_p","rootstock_irrigation_p", "year_p","block_p")

for (i in which(colnames(morpho_pcs)=="PC1"):which(colnames(morpho_pcs)=="PC20")) {
  pc_lm=lm(morpho_pcs[,i] ~ rootstock + irrigation+rootstock*irrigation+year+block, data = morpho_pcs)
  pc_sums=anova(pc_lm)$Sum
  total_var=sum(pc_sums)
  test_total=(pc_sums/total_var)*100
  pc_p=anova(pc_lm)$Pr[1:5]
  var_p=append(test_total,pc_p)
  total_var_mat_p=rbind(total_var_mat_p,var_p)
  rownames(total_var_mat_p)[i-8]=paste("PC", i-9, sep="_")
}

total_var_mat_p=total_var_mat_p[-1,]

total_var_mat_p <- as.data.frame(total_var_mat_p)
rownames(total_var_mat_p) <- 1:nrow(total_var_mat_p)
total_var_mat_p <- cbind(rownames(total_var_mat_p),total_var_mat_p)
colnames(total_var_mat_p)[1]="PC"
write.table(total_var_mat_p, "mt_vernon_pcs_var_p.txt", sep="\t", row.names=F, quote=F, col.names=T)

#then modify the table so it is easy to plot
total_var_mat_p <- read.table('mt_vernon_pcs_var_p.txt', header=T,sep="\t")
#remove residuals since I am not going to visualize it

total_var_mat_no_res_edit <- dplyr::select(total_var_mat_p,-residual_var)
total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("PC", "rootstock_var")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="var"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("PC", "irrigation_var")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="var"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("PC", "rootstock_irrigation_var")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="var"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("PC", "year_var")]
colnames(total_var_mat_no_res_edit_year)[2]="var"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("PC", "block_var")]
colnames(total_var_mat_no_res_edit_block)[2]="var"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")

total_var_mat_no_res_edit_all <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block)

#then do the same thing with p-values
total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("PC", "rootstock_p")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="p-value"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("PC", "irrigation_p")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="p-value"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("PC", "rootstock_irrigation_p")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="p-value"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("PC", "year_p")]
colnames(total_var_mat_no_res_edit_year)[2]="p-value"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("PC", "block_p")]
colnames(total_var_mat_no_res_edit_block)[2]="p-value"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")

total_var_mat_no_res_edit_all_p <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block)

var_p=cbind(total_var_mat_no_res_edit_all,total_var_mat_no_res_edit_all_p)
table(var_p[,1]==var_p[,4])
table(var_p[,3]==var_p[,6])
var_p <- var_p[,-4]
var_p <- var_p[,-5]
colnames(var_p)[4]="p_value"

#create a palette 
morpho_palette <- c(rgb(73,0,146,maxColorValue=255),rgb(146,0,0,maxColorValue=255),rgb(0,146,146,maxColorValue=255),rgb(255,182,119,maxColorValue=255), "lightgray")

#only plot significant p_values
var_p_sig=dplyr::filter(var_p, p_value<0.05)
#arrange by PC
var_p_sig[,"PC"] <- as.numeric(as.character(var_p_sig[,"PC"]))
var_p_sig <- arrange(var_p_sig, PC)
#add in variance explained by PCs
var_explained=read.table('mt_vernon_pcs_explained_prop.txt', header=T,sep="\t")
var_p_sig_info <- left_join(var_p_sig, var_explained)
write.table(var_p_sig_info, "mt_vernon_morpho_var_p_sig.txt", sep="\t", row.names=F, quote=F, col.names=T)


pdf("2morpho_variance_heatmap_black_viridis.pdf",family="Arial", width=6, height=8)
t <- ggplot(data=var_p_sig, aes(x=trait, y=as.numeric(as.character(PC)), fill=var))
t + geom_tile(color="black", size=0.2)+theme_minimal()+scale_fill_gradientn(colors=c("white", "black"),limits=c(0, 8))+labs(x = "Trait", y="PC (% variance explained)") + theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"), legend.position = "bottom")+  scale_y_continuous(trans = "reverse",breaks = seq(1, 20, 1), minor_breaks = NULL, labels=paste("PC",var_explained[1:20,"PC"], " (", round(as.numeric(var_explained[1:20,"PC_variance_explained"]), digits=2), "%)", sep=""))+scale_x_discrete(position = "top")+expand_limits(x = 0, y = 19)+ scale_fill_viridis(option="inferno", name="% variance explained", direction = -1)
dev.off()


#all ionomics data with outliers removed
#test if there is a significant impact on each factor for ionomics 
lm_ionomics_all <- read.table('mt_vernon_ionomics_2014_2016_all_data_with_vineyard_info.txt', header=T,sep="\t")
#make sure all variables are factors
lm_ionomics_all[,"year"] <- as.factor(lm_ionomics_all[,"year"])
total_var_mat=matrix(,1,9)
colnames(total_var_mat)=c("rootstock_var", "irrigation_var", "leaf_var", "rootstock_irrigation_var","rootstock_leaf_var", "irrigation_leaf_var","year_var", "block_var", "residual_var")

for (i in which(colnames(lm_ionomics_all)=="Na"):which(colnames(lm_ionomics_all)=="Cd")) {
  pc_lm=lm(lm_ionomics_all[,i] ~ rootstock + irrigation+leaf+rootstock*irrigation+rootstock*leaf+irrigation*leaf+year+block, data = lm_ionomics_all)
  pc_sums=anova(pc_lm)$Sum
  total_var=sum(pc_sums)
  test_total=(pc_sums/total_var)*100
  total_var_mat=rbind(total_var_mat,test_total)
  rownames(total_var_mat)[i-7]=colnames(lm_ionomics_all)[i]
}

total_var_mat=total_var_mat[-1,]
total_var_mat <- as.data.frame(total_var_mat)
#now to plot it what is the max variance explained by any of these factors, ignore residuals
total_var_mat_no_res <- dplyr::select(total_var_mat, -residual_var)

#now run the script to include which ones are significant

total_var_mat_p=matrix(,1,17)
colnames(total_var_mat_p)=c("rootstock_var", "irrigation_var", "leaf_var", "rootstock_irrigation_var","rootstock_leaf_var", "irrigation_leaf_var","year_var", "block_var", "residual_var","rootstock_p", "irrigation_p", "leaf_p", "rootstock_irrigation_p","rootstock_leaf_p", "irrigation_leaf_p","year_p", "block_p")

for (i in which(colnames(lm_ionomics_all)=="Na"):which(colnames(lm_ionomics_all)=="Cd")) {
  pc_lm=lm(lm_ionomics_all[,i] ~ rootstock + irrigation+leaf+rootstock*irrigation+rootstock*leaf+irrigation*leaf+year+block, data = lm_ionomics_all)
  pc_sums=anova(pc_lm)$Sum
  total_var=sum(pc_sums)
  test_total=(pc_sums/total_var)*100
  pc_p=anova(pc_lm)$Pr[1:8]
  var_p=append(test_total,pc_p)
  total_var_mat_p=rbind(total_var_mat_p,var_p)
  rownames(total_var_mat_p)[i-7]=colnames(lm_ionomics_all)[i]
}

total_var_mat_p=total_var_mat_p[-1,]
total_var_mat_p <- as.data.frame(total_var_mat_p)

total_var_mat_p <- cbind(rownames(total_var_mat_p),total_var_mat_p)
colnames(total_var_mat_p)[1]="element"
rownames(total_var_mat_p)=NULL

write.table(total_var_mat_p, "mt_vernon_ionomics_var_p.txt", sep="\t", row.names=T, quote=F, col.names=T)

#then modify the table so it is easy to plot
total_var_mat_p <- read.table('mt_vernon_ionomics_var_p.txt', header=T,sep="\t")
#remove residuals since I am not going to visualize it
total_var_mat_no_res_edit <- dplyr::select(total_var_mat_p,-residual_var)

total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("element", "rootstock_var")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="var"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("element", "irrigation_var")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="var"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("element", "rootstock_irrigation_var")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="var"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("element", "year_var")]
colnames(total_var_mat_no_res_edit_year)[2]="var"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("element", "block_var")]
colnames(total_var_mat_no_res_edit_block)[2]="var"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")
total_var_mat_no_res_edit_leaf <-  total_var_mat_no_res_edit[,c("element", "leaf_var")]
colnames(total_var_mat_no_res_edit_leaf)[2]="var"
total_var_mat_no_res_edit_leaf <- mutate(total_var_mat_no_res_edit_leaf, trait="leaf")
total_var_mat_no_res_edit_rootstock_leaf <-  total_var_mat_no_res_edit[,c("element", "rootstock_leaf_var")]
colnames(total_var_mat_no_res_edit_rootstock_leaf)[2]="var"
total_var_mat_no_res_edit_rootstock_leaf <- mutate(total_var_mat_no_res_edit_rootstock_leaf, trait="rootstock_leaf")
total_var_mat_no_res_edit_irrigation_leaf <-  total_var_mat_no_res_edit[,c("element", "irrigation_leaf_var")]
colnames(total_var_mat_no_res_edit_irrigation_leaf)[2]="var"
total_var_mat_no_res_edit_irrigation_leaf <- mutate(total_var_mat_no_res_edit_irrigation_leaf, trait="irrigation_leaf")


total_var_mat_no_res_edit_all <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block,total_var_mat_no_res_edit_leaf,total_var_mat_no_res_edit_rootstock_leaf,total_var_mat_no_res_edit_irrigation_leaf)

#then do the same thing with p-values
total_var_mat_no_res_edit_rootstock <- total_var_mat_no_res_edit[,c("element", "rootstock_p")]
colnames(total_var_mat_no_res_edit_rootstock)[2]="p-value"
total_var_mat_no_res_edit_rootstock <- mutate(total_var_mat_no_res_edit_rootstock, trait="rootstock")
total_var_mat_no_res_edit_irrigation<- total_var_mat_no_res_edit[,c("element", "irrigation_p")]
colnames(total_var_mat_no_res_edit_irrigation)[2]="p-value"
total_var_mat_no_res_edit_irrigation <- mutate(total_var_mat_no_res_edit_irrigation, trait="irrigation")
total_var_mat_no_res_edit_rootstock_irrigation<- total_var_mat_no_res_edit[,c("element", "rootstock_irrigation_p")]
colnames(total_var_mat_no_res_edit_rootstock_irrigation)[2]="p-value"
total_var_mat_no_res_edit_rootstock_irrigation <- mutate(total_var_mat_no_res_edit_rootstock_irrigation, trait="rootstock irrigation")
total_var_mat_no_res_edit_year<- total_var_mat_no_res_edit[,c("element", "year_p")]
colnames(total_var_mat_no_res_edit_year)[2]="p-value"
total_var_mat_no_res_edit_year <- mutate(total_var_mat_no_res_edit_year, trait="year")
total_var_mat_no_res_edit_block<- total_var_mat_no_res_edit[,c("element", "block_p")]
colnames(total_var_mat_no_res_edit_block)[2]="p-value"
total_var_mat_no_res_edit_block <- mutate(total_var_mat_no_res_edit_block, trait="block")
total_var_mat_no_res_edit_leaf<- total_var_mat_no_res_edit[,c("element", "leaf_p")]
colnames(total_var_mat_no_res_edit_leaf)[2]="p-value"
total_var_mat_no_res_edit_leaf <- mutate(total_var_mat_no_res_edit_leaf, trait="leaf")
total_var_mat_no_res_edit_rootstock_leaf<- total_var_mat_no_res_edit[,c("element", "rootstock_leaf_p")]
colnames(total_var_mat_no_res_edit_rootstock_leaf)[2]="p-value"
total_var_mat_no_res_edit_rootstock_leaf <- mutate(total_var_mat_no_res_edit_rootstock_leaf, trait="rootstock_leaf")
total_var_mat_no_res_edit_irrigation_leaf<- total_var_mat_no_res_edit[,c("element", "irrigation_leaf_p")]
colnames(total_var_mat_no_res_edit_irrigation_leaf)[2]="p-value"
total_var_mat_no_res_edit_irrigation_leaf <- mutate(total_var_mat_no_res_edit_irrigation_leaf, trait="irrigation_leaf")

total_var_mat_no_res_edit_all_p <- rbind(total_var_mat_no_res_edit_rootstock,total_var_mat_no_res_edit_irrigation,total_var_mat_no_res_edit_rootstock_irrigation,total_var_mat_no_res_edit_year,total_var_mat_no_res_edit_block,total_var_mat_no_res_edit_leaf,total_var_mat_no_res_edit_rootstock_leaf,total_var_mat_no_res_edit_irrigation_leaf)

var_p=cbind(total_var_mat_no_res_edit_all,total_var_mat_no_res_edit_all_p)
table(var_p[,1]==var_p[,4])
table(var_p[,3]==var_p[,6])
var_p <- var_p[,-4]
var_p <- var_p[,-5]
colnames(var_p)[4]="p_value"

#elements need a number for y-axis when plotting 
var_p <- arrange(var_p, element)
var_p <- mutate(var_p, element_number=rep(1:17,each=8))
#only plot significant p_values
var_p_sig=dplyr::filter(var_p, p_value<0.05)
write.table(var_p_sig[,c("element", "trait", "var", "p_value")], "20180309_mt_vernon_ionomics_var_p_sig.txt", sep="\t", row.names=F, quote=F, col.names=T)

pdf("ionomic_variance_heatmap_black_viridis.pdf",family="Arial", width=6, height=8)
t <- ggplot(data=var_p_sig, aes(x=trait, y=as.numeric(as.character(element_number)), fill=var))
t + geom_tile(color="black", size=0.2)+theme_minimal()+scale_fill_gradientn(colors=c("white", "black"),limits=c(0, 61))+labs(x = "Trait", y="Element") + theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"),legend.position = "bottom")+scale_y_continuous(trans="reverse",breaks = seq(1, 17, 1), minor_breaks=NULL, labels=unique(var_p_sig[,"element"]))+scale_x_discrete(position = "top")+ scale_fill_viridis(option="inferno", name="% variance explained", direction = -1)
dev.off()


#boxplots for supplement
ionomics_all_plot=read.table('mt_vernon_ionomics_2014_2016_all_data_with_vineyard_info.txt', header=T,sep="\t")
ionomics_all_plot <- ionomics_all_plot %>% mutate(col_id=paste(row, column, plant, leaf, sep="_"))
ionomics_all_plot <- dplyr::select(ionomics_all_plot, col_id, year, irrigation, rootstock, leaf, Na:Cd)
ionomics_all_plot <- melt(ionomics_all_plot,id = c('col_id','year', 'irrigation', 'rootstock', 'leaf'), variable.name = 'el')

#rootstock plot

ionomics_all_plot$rootstock <- factor(ionomics_all_plot$rootstock ,
                                      levels = c('OWN','1103P', '3309C', 'SO4'),ordered = TRUE)

#reorder elements for Figure. 
ionomics_all_plot$el <- factor(ionomics_all_plot$variable, 
                               levels = c("Al","Ca",
                                          "Cd","Co","Cu","Fe","K",
                                          "Mg","Mn","Mo","Na","Ni","P","Rb", "S", "Sr", "Zn") 
)

#re-order colours since OWN is first now
rootstock_palette <- c("gray","#1b9e77", "#7570b3",  "#e6ab02")

plot1 <- ggplot(ionomics_all_plot,aes(x = interaction(rootstock,year), y = value, colour=rootstock, fill = rootstock)) + geom_boxplot() +  geom_vline(xintercept=c(4.5)) +
  facet_wrap(~el, scales = "free_y",ncol = 3) + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=rootstock_palette)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_color_manual(values=c("black", "black", "black", "black"))+theme(legend.position = "none")

#leaf plot
leaf_palette <- c(rgb(145,102,189,maxColorValue=255),rgb(178,76,76,maxColorValue=255),rgb(0,146,146,maxColorValue=255),rgb(255,182,119,maxColorValue=255))

#Make line drawings
ionomics_all_plot <- mutate(ionomics_all_plot, Numpos = 1)

ionomics_all_plot$Numpos[ionomics_all_plot$leaf == 'Y'] <- 2
ionomics_all_plot$Numpos[ionomics_all_plot$leaf == 'Z'] <- 3

ionomics_all_plot <- group_by(ionomics_all_plot, col_id,year,Numpos, variable)

#reorder elements for Figure. 
ionomics_all_plot$el <- factor(ionomics_all_plot$variable, 
                               levels = c("Al","Ca",
                                          "Cd","Co","Cu","Fe","K",
                                          "Mg","Mn","Mo","Na","Ni","P","Rb", "S", "Sr", "Zn") 
)


plot2 <- ggplot(ionomics_all_plot,aes(x = interaction(leaf,year), y = value, colour=leaf, fill = leaf)) + geom_boxplot() +  geom_vline(xintercept=c(3.5)) +
  facet_wrap(~el, scales = "free_y", ncol=3) + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=leaf_palette)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_color_manual(values=c("black", "black", "black"))+theme(legend.position = "none")



#rootstock by irrigation plot 

#Make line drawings
ionomics_all_ri <- mutate(ionomics_all_plot, root_irr = paste(rootstock, irrigation, sep="_"))

ionomics_all_plot$root_irr <- factor(ionomics_all_ri$root_irr ,
                                     levels = c('OWN_none','OWN_half', 'OWN_full','1103P_none', '1103P_half','1103P_full', '3309C_none','3309C_half','3309C_full', 'SO4_none','SO4_half','SO4_full'),ordered = TRUE)

#reorder elements for Figure. 
ionomics_all_plot$el <- factor(ionomics_all_plot$variable, 
                               levels = c("Al","Ca",
                                          "Cd","Co","Cu","Fe","K",
                                          "Mg","Mn","Mo","Na","Ni","P","Rb", "S", "Sr", "Zn") 
)


plot3 <- ggplot(ionomics_all_plot,aes(x = interaction(root_irr,year), y = value, colour=rootstock, fill = irrigation)) + geom_boxplot() +  geom_vline(xintercept=c(12.5)) +
  facet_wrap(~el, scales = "free_y", ncol=3) + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_brewer(direction=-1)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_colour_manual(values=rootstock_palette)+theme(legend.position = "none")

dev.off()
require(cowplot)
pdf("ionomics_boxplots_supplement1.pdf", width=10, height=12,family="Arial")
plot_grid(plot1, plot2, ncol = 2,labels=NA)
dev.off()
pdf("ionomics_boxplots_supplement2.pdf", width=10, height=12,family="Arial")
plot_grid(plot3, ncol = 1,labels=NA)
dev.off()

#next step is to make the "boxplots of interest" for the main body of the paper, include just 2016 data 

ionomics_all_plot$rootstock <- factor(ionomics_all_plot$rootstock ,
                                      levels = c('OWN','1103P', '3309C', 'SO4'),ordered = TRUE)


ionomics_all_plot$root_irr <- factor(ionomics_all_plot$root_irr ,
                                     levels = c('OWN_none','OWN_half', 'OWN_full','1103P_none', '1103P_half','1103P_full', '3309C_none','3309C_half','3309C_full', 'SO4_none','SO4_half','SO4_full'),ordered = TRUE)

ionomics_k <- filter(ionomics_all_plot, year=="2016", variable=="K")
ionomics_k_plot <- ggplot(ionomics_k,aes(x = leaf, y = value, fill = leaf)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=leaf_palette)+theme_bw()+theme(legend.position = "none")+labs(x="leaf position", y="K (PPM)") 

ionomics_ca <- filter(ionomics_all_plot, year=="2016", variable=="Ca")
ionomics_ca_plot <- ggplot(ionomics_ca,aes(x = leaf, y = value, fill = leaf)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=leaf_palette)+theme_bw()+theme(legend.position = "none")+labs(x="leaf position", y="Ca (PPM)") 

ionomics_mo <- filter(ionomics_all_plot, year=="2016", variable=="Mo")
ionomics_mo_plot <- ggplot(ionomics_mo,aes(x = root_irr, y = value, colour=rootstock, fill = irrigation)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_brewer(direction=-1)+theme_bw()+ scale_colour_manual(values=rootstock_palette)+theme(legend.position = "none")+labs(x="rootstock * irrigation", y="Mo (PPM)") 

ionomics_ni <- filter(ionomics_all_plot, year=="2016", variable=="Ni")
ionomics_ni_plot <- ggplot(ionomics_ni,aes(x = rootstock, y = value, fill = rootstock)) + geom_boxplot() + theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank())+ scale_fill_manual(values=rootstock_palette)+theme_bw()+theme(legend.position = "none")+labs(x="rootstock", y="Ni (PPM)") 


dev.off()
require(cowplot)
pdf("ionomics_boxplots_main_4.pdf", width=10, height=10,family="Arial")
plot_grid(ionomics_ca_plot, ionomics_k_plot,ionomics_ni_plot,ionomics_mo_plot, ncol = 2,labels="AUTO")
dev.off()

#RNA-seq analysis 

exp <- read.table("mt_vernon_rna_counts.txt")
sample <- read.table("mt_vernon_rna_design.txt")
design<-make.design.matrix(sample)
rna_fit<- p.vector(exp, design,  Q = 0.05, MT.adjust = "BH", min.obs = 6, counts=T)
NBt <- T.fit(rna_fit)
#compare to control get list of genes that are significantly different for every rootstock in comparison to own-rooted
control <- get.siggenes(NBt, vars="groups")
control$summary
write.table(control$summary, "gene_expression_results.txt", sep="\t", row.names=F, quote=F, col.names=T)

suma2Venn(control$summary[, c(1:4)])

#I want the lists of the ones that only differ for 1103P, only differ for 3309C and only differ for SO4.
so4_only <- control$summary[!(control$summary[,4] %in% control$summary[,3]),4]
so4_only <- so4_only[!so4_only %in% control$summary[,2]]
so4_only <- so4_only[!so4_only %in% control$summary[,1]]

c3309_only <- control$summary[!(control$summary[,2] %in% control$summary[,3]),2]
c3309_only <- c3309_only[!c3309_only %in% control$summary[,4]]
c3309_only <- c3309_only[!c3309_only %in% control$summary[,1]]

p1103_only <- control$summary[!(control$summary[,3] %in% control$summary[,2]),3]
p1103_only <- p1103_only[!p1103_only %in% control$summary[,4]]
p1103_only <- p1103_only[!p1103_only %in% control$summary[,1]]

own_only <- control$summary[!(control$summary[,1] %in% control$summary[,2]),1]
own_only <- own_only[!own_only %in% control$summary[,3]]
own_only <- own_only[!own_only %in% control$summary[,4]]

write.table(so4_only, "so4_no_overlap.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(c3309_only, "3309c_no_overlap.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(p1103_only, "1103p_no_overlap.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(own_only, "own_no_overlap.txt", sep="\t", row.names=F, quote=F, col.names=F)

#now get ones that are shared amongst all rootstocks but not overlapping with own-rooted (5)
overlap_rootstocks <- control$summary[(control$summary[,3] %in% control$summary[,2]),3]
overlap_rootstocks <- overlap_rootstocks[overlap_rootstocks %in% control$summary[,4]]
overlap_rootstocks <- overlap_rootstocks[!overlap_rootstocks %in% control$summary[,1]]
write.table(overlap_rootstocks, "overlap_rootstocks_no_ownrooted.txt", sep="\t", row.names=F, quote=F, col.names=F)

#get only the ones that are shared amongst all (121)

shared_all <- control$summary[(control$summary[,1] %in% control$summary[,2]),1]
shared_all <- shared_all[shared_all %in% control$summary[,3]]
shared_all <- shared_all[shared_all %in% control$summary[,4]]
write.table(shared_all, "shared_all_overlap.txt", sep="\t", row.names=F, quote=F, col.names=F)
