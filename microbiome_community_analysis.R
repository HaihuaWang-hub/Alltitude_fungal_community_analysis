#call the library
call_library <- function() {
  library(DESeq2)
  library("dplyr")
  library("ggplot2")
  library("calibrate")
  library("RColorBrewer")
  library("glmpca")
  library("gplots")
  library("amap")
  library("ggplot2")
  library("BiocParallel")
  library(ggpubr)
  library("ggpubr")
  library(pheatmap)
  library(limma)
  library(corrplot)
  library(vegan)
  library(tidyverse)
  library(Rtsne)
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(ggrepel)
  library(ggalt)
  library(clusterProfiler)
  library(DOSE)
  library(topGO)
  library(enrichplot)
  library(stringr)
  library(ggnewscale)
  library(ggupset)
  library(ggridges)
  library(europepmc)
  library(Mfuzz)
  library(gdata)
  library(GGally)
  library("gridExtra")
}
call_library()


rm(list=ls())
#1 get the data
work_dir="F:/PostDoc_dataset/16S_XN/workflow_results/799F_1193R-39samples/workflow_results/OtuTaxon_summary/tax_summary_r"
#abundance_out_full = read_tsv(paste(work_dir,"otu_taxon_otu.percent.full.tsv", sep = "/"))
abundance_otu = read_tsv(paste(work_dir,"otu_taxon_otu.percent.tsv", sep = "/"), row.names=1)
abundance_otu = data.frame(row.names = abundance_otu$`OTU ID`, abundance_otu[,2:40])



rm(list=ls())
#setwd("F:/PostDoc_dataset/bioinfomatic/performance evaluation")
setwd("F:/PostDoc_dataset/16S_XN/workflow_results")
list.files(getwd())


ggplot(abundance_otu, aes(x=sample, y=100 * value, fill = Taxonomy)) +
  #数据输入：样本、物种、丰度
  geom_col(position = 'stack', width = 0.6) + # stack：堆叠图
  scale_y_continuous(expand=c(0, 0))+# 调整y轴属性，使柱子与X轴坐标接触
  scale_fill_manual(values =  rev(c('#FF0000', 
                                    '#FF88C2', '#FF00FF', '#9999FF', '#33FFFF',
                                    '#33FF33', '#D1BBFF', '#770077', '#EE7700', 
                                    '#CCEEFF', '#0000AA'))) + #手动修改颜色
  labs(x = 'Samples', y = '相对分度\n Relative Abundance(%)') + #设置X轴和Y轴的信息
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) + #设置主题背景，根据自己的需求定制
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))

p

#2.3.5 #PCA analysis

library(ggplot2)
normalized_counts.table = read.table("normalized_countdata.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
head(normalized_counts.table)
normalized_counts.table <- t(normalized_counts.table)
rownames <-rownames(normalized_counts.table)
normalized_counts.table <- data.frame(normalized_counts.table,rownames)


ord <- prcomp(t(abundance_otu))
summary(ord)
dt <- ord$x
df <- data.frame(dt,sample_name=colnames(abundance_otu))
head(df)
summ <- summary(ord)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

pdf("PCA_normized.pdf")
(p1 <- ggplot(df,aes(df$PC1,df$PC2,color=sample_name))+stat_ellipse(aes(fill=sample_name),
                                                                 type="norm",geom="polygon",alpha=0.2,color=NA)+guides(fill=F)+geom_point()+labs(x=xlab,y=ylab,color=""))
dev.off()



#functional stability
#################################################
rm(list=ls())
setwd("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/4_taxa_abundance")

data = read.csv("OTU_table.csv",header = TRUE, row.names = 1)

df_abundance=as.data.frame(lapply(df, function(x) x / sum(x)))
row.names(df_abundance) <- row.names(df)
colSums(df_abundance)

library(genefilter)
df_abundance$SD_score = rowSds(as.matrix(df_abundance[,1:18]))
df_abundance$mean = rowMeans(df_abundance[1:18])



df_abundance$variation_degree=abs(df_abundance[,1] - df_abundance[,c("mean")])/df_abundance$SD_score
AVD=sum(df_abundance$variation_degree)/(nrow(df_abundance)*ncol(df_abundance))




#function to calculate the community or function stability
fun_stability <- function(df_ori) {
  library(genefilter)
  df=as.data.frame(lapply(df_ori, function(x) x / sum(x)))
  df_value = as.data.frame(matrix(nrow = nrow(df), ncol = 2))
  df_variation_degree = as.data.frame(matrix(nrow = nrow(df), ncol = ncol(df)))
     colnames(df_variation_degree) <- colnames(df)
  df_AVD = as.data.frame(matrix(nrow = 1, ncol = ncol(df)))
     colnames(df_AVD) <- colnames(df)
  df_value$SD_score = rowSds(as.matrix(df))
  df_value$Mean = rowMeans(df)
  for (i in 1:ncol(df)) {
    df_variation_degree[,i] = abs(df[,i] - df_value[,c("Mean")])/df_value$SD_score
    df_variation_degree[,i][which(df_variation_degree[,i] == "NaN" )] <- 1
    df_AVD[,i] = sum(df_variation_degree[,i])/((nrow(df)*ncol(df)))
  }
  data.frame(df_AVD)
}

stability <- fun_stability(data)
root1500_stable <- fun_stability(data[,1:3])
root1900_stable <- fun_stability(data[,4:6])
root2300_stable <- fun_stability(data[,7:9])
soil1500_stable <- fun_stability(data[,10:12])
soil1900_stable <- fun_stability(data[,13:15])
soil2300_stable <- fun_stability(data[,16:18])
stability = cbind(root1500_stable,root1900_stable,root2300_stable,soil1500_stable,soil1900_stable,soil2300_stable)

write.table(stability, "D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/8_stability/stability.txt")

rm(list = ls())



library(corrplot)
stability = read.table("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/8_stability/stability.txt") 
  row.names(stability) = "stability"
env = read.table("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/4_taxa_abundance/envronmental_variables.txt", header = T)
df <- cbind(t(stability),env)
df <- df[,-2:-3]
df <- df[10:18,]
df <- df[1:9,]
matrix <- cor(df,method="pearson",use="everything")
corrplot(corr=matrix)

ggplot(df, aes(x = AP, y = stability)) + 
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') + 
  geom_point(aes(colour = altitude), size = 4) +
  labs(y = 'phosphatase activity', x = 'AP') + 
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12), 
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'), 
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  scale_colour_continuous(high = 'navy', low = 'salmon')

fit <- lm(df$stability~df$AP)
summary(fit)
#############################




#upset:https://zhuanlan.zhihu.com/p/370210775
#upset  and venn
############################################
#
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
rm(list=ls())
setwd("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/4_taxa_abundance")

df = read.csv("OTU_table.csv",header = TRUE, row.names = 1)
df$root1500=rowMeans(as.matrix.data.frame(df[,1:3]))
df$root1900=rowMeans(as.matrix.data.frame(df[,4:6]))
df$root2300=rowMeans(as.matrix.data.frame(df[,7:9]))
df$soil1500=rowMeans(as.matrix.data.frame(df[,10:12]))
df$soil1900=rowMeans(as.matrix.data.frame(df[,13:15]))
df$soil2300=rowMeans(as.matrix.data.frame(df[,16:18]))
Data<-df[,19:24]
df$root_mean_log <- log2(rowMeans(as.matrix.data.frame(df[,1:9])))
df$soil_mean_log <- log2(rowMeans(as.matrix.data.frame(df[,10:18])))

# 将数据表中>0的数值替换为1，数值>0则OTU或各级数据在分组中有出现，替换为1是用于可视化的数据准备
Data[Data>0]=1
head(Data)
tail(Data)
colSums(Data)
Data_1 <- cbind(Data, df[,25:26])
## draw upset diagram
# 为方便代码兼容R3.6和R4.0版本，再进行一次数据处理

# 目的是将'Data4Pic'的行名转换为'numeric'的行号
row.names(Data) <- 1:nrow(Data)

# 随后进行数据可视化
library(UpSetR)
# 集合图的基本图形绘制
pdf(file="Upset.pdf", width=12, height=7, pointsize=8)
upset(Data_1, sets = c("root1500","root1900","root2300","soil1500","soil1900","soil2300"),order.by = "freq",keep.order = TRUE,
      queries = list(
        list(
          query = intersects, 
          params = list("root1900"), 
          active = T,
          color = "blue", active = F),
        list(
          query = intersects, 
          params = list("root2300"), 
          active = T,
          color = "red", active = F),
        list(
          query = intersects, 
          params = list("soil1500"), 
          active = T,
          color = "orange", active = F),
        list(
          query = intersects, 
          params = list("soil1900"), 
          active = T,
          color = "chartreuse3", active = F),
        list(
          query = intersects, 
          params = list("soil2300"), 
          active = T,
          color = "brown", active = F)
      ),
     boxplot.summary = c("root_mean_log","soil_mean_log"),
     sets.bar.color = c("purple","blue","red","orange","chartreuse3","brown")
     )
      
dev.off()     
      
      
     



combine <- Data
combine$ID <- row.names(combine)

  venn.diagram(
    x = list(
      combine %>% filter(root1500 > 0)  %>% select(ID) %>% unlist() , 
      combine %>% filter(root1900 > 0)  %>% select(ID) %>% unlist(),
      combine %>% filter(root2300 > 0)  %>% select(ID) %>% unlist()
    ),
    category.names = c( "root1500" , "root1900","root2300"),
    filename = "venn_root.tiff",
    output = TRUE ,
    imagetype="tiff" ,
    
    height = 4800 , 
    width = 4800 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col=c("#440154ff", '#21908dff', '#fde725ff'),
    fill = c(alpha("#440154ff",0.6), alpha('#21908dff',0.6), alpha('#fde725ff',0.6)),
    cex = 5,
    
    fontfamily = "sans",
    cat.cex = 4,
    cat.default.pos = "outer",
    # cat.pos = c(227, 327, 135,135),
    #cat.dist = c(0.065, 0.055, 0.085, 0.085),
    cat.fontfamily = "sans",
    #cat.col = c("#440154ff", '#21908dff', '#fde725ff','#fde725ff'),
    #rotation = 1
  )
  
  combine$root = rowSums(combine[,1:3])
  combine$soil = rowSums(combine[,4:6])
  venn.diagram(
    x = list(
      combine %>% filter(root > 0)  %>% select(ID) %>% unlist() , 
      combine %>% filter(soil > 0)  %>% select(ID) %>% unlist()
    ),
    category.names = c( "root" , "soil1"),
    filename = "venn_root_soil.tiff",
    output = TRUE ,
    imagetype="tiff" ,
    
    height = 4800 , 
    width = 4800 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col=c("#440154ff", '#21908dff'),
    fill = c(alpha("#440154ff",0.6), alpha('#21908dff',0.6)),
    cex = 5,
    
    fontfamily = "sans",
    cat.cex = 4,
    cat.default.pos = "outer",
    # cat.pos = c(227, 327, 135,135),
    #cat.dist = c(0.065, 0.055, 0.085, 0.085),
    cat.fontfamily = "sans",
    #cat.col = c("#440154ff", '#21908dff', '#fde725ff','#fde725ff'),
    #rotation = 1
  )
##################################




################################
##linear analysis plotting
################################
library(ggplot2)
setwd("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/1_environmental_anova/123")

df = read.csv("envronmental variables.csv",header = TRUE,sep = ",")

xx = ggplot(df, aes(x = phosphatase, y = AP)) + 
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') + 
  geom_point(aes(colour = altitude), size = 4) +
  labs(y = 'phosphatase activity', x = 'AP') + 
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12), 
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'), 
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx

#linear fit
fit <- lm(df$phosphatase~df$AP)
fit <- lm(df$sucrase~df$SOM)
fit <- lm(df$urase~df$NH4)
fit <- lm(df$urase~df$NO3)
summary(fit)
###########










# community diversity
###############################################################

setwd("D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/4_taxa_abundance")

# alpha_diversity
library(vegan)
library(picante)

#read the table
otu = read.csv("OTU_table.csv", header=T, row.names = 1, comment.char = "", )
otu <- t(otu)

##evaluate Richness ָ
richness <- rowSums(otu > 0)
#or
richness <- estimateR(otu)[1, ]

#Shannon
shannon_index <- diversity(otu, index = 'shannon', base = exp(1))

#Shannon diversity
shannon_diversity <- exp(1)^shannon_index

#pielou evenness
pielou <- shannon_index / log(richness, exp(1))

##Simpson
#Gini-Simpson
gini_simpson_index <- diversity(otu, index = 'simpson')

#Simpson ָindex
simpson_index <- 1 - gini_simpson_index

#Invsimpson 
invsimpson_index <- 1 / gini_simpson_index
#?invsimpson_index
invsimpson_index <- diversity(otu, index = 'invsimpson')

#Simpson diversity
simpson_diversity <- 1 / (1 - gini_simpson_index)


#Simpson  equitability 
equitability <- 1 / (richness * (1 - gini_simpson_index))

##Chao1 & ACE
#Chao1 ָ
chao1 <- estimateR(otu)[2, ]
#ACE ָ
ace <- estimateR(otu)[4, ]

##goods_coverage ָ
goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)

##??ϵ?????ԣ??????????ȣ?????ָ?????????ļ???
#????ʱ?????и??????޸?????PD_whole_tree??????????һ???ģ??????޸????ļ?????????
tree <- read.tree('otu_tree.tree')
pd_whole_tree <- pd(otu, tree, include.root = FALSE)



library(picante)

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson ָ??
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}



alpha_all <- alpha(otu, base = 2)

write.csv(alpha_all, 'alpha.csv', quote = FALSE)
#####################################################################





#PERMANOVA
#https://cloud.tencent.com/developer/article/1667540
########################################################################

library(vegan)

otu = read.csv("OTU_table.csv", header=T, row.names = 1, comment.char = "", )
otu <- t(otu)
View(otu)

metadata <- read.csv("metadata.csv")
View(metadata)

#1.root~altitude
otu1 <- otu[1:9,]
metadata1 <- metadata[1:9,1:2]
adonis(otu1~ Group,data = metadata1,permutations = 9999,method="bray")
#2.soil~altitude
otu1 <- otu[10:18,]
metadata1 <- metadata[10:18,1:2]
adonis(otu1~ Group,data = metadata1,permutations = 9999,method="bray")
#3.root~soil
otu1 <- otu
metadata1 <- metadata[,c(1,3)]
adonis(otu1~ Group_2,data = metadata1,permutations = 9999,method="bray")



bray_dis <- vegdist(abc, method = 'bray')
adonis_result_bray <- adonis(bray_dis~fungal*iron, group_data, permutations = 999)
#####################################################################################







#ggClusterNet vennͼ
# install ggClusterNet
devtools::install_github("taowenmicro/ggClusterNet")

library(ggplot2)
library(ggrepel)
library(ggClusterNet)
library(phyloseq)
library(dplyr)
# ????????#-----????????#-------
data(ps)














#correlations analysis
#https://zhuanlan.zhihu.com/p/114008006
#https://blog.csdn.net/lalaxumelala/article/details/86084040
#https://zhuanlan.zhihu.com/p/183880858
#chang color:https://www.thinbug.com/q/30743983
#install.packages("corrplot")
library(corrplot)
otu = read.csv("OTU_table.csv", header=T, row.names = 1, comment.char = "", )
otu <- t(otu)
View(otu)

env = read.table("envronmental_variables.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
View(env)

#1.soil~env
otu1 <- otu[10:18,]
env1 <- env[1:9,2:25]
matrix <- cor(otu1,env1,method="pearson",use="everything")
dim(matrix)
write.csv(matrix,file="soil_env_correlation.csv",quote=F,row.names=T)

#matrix <- cor(data,method="spearman",use="everything")
matrix <- cor(data1,data2,method="spearman",use="complete.obs")
dim(matrix)
matrix[1:5,]
#correlation matrix
write.csv(matrix,file="data/spearman_correlation.csv",quote=F,row.names=T)

#correlation plot
pdf(file="corr_plot.pdf",width=3000,height=3000)
#plotting
corrplot(matrix, type = "lower", order = "original", tl.col = "black", tl.srt = 45,cl.lim = c(0, 1))
#plotting
corrplot(matrix,title = "Correlation Plot",method = "color",outline = T, addgrid.col = "white", order="hclust", mar = c(0,0,1,0), addrect = 4, rect.col = "grey", rect.lwd = 1,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 3)
dev.off()
corrplot(corr=matrix)
corrplot(matrix, type = "lower")


#2.root~env
otu1 <- otu[1:9,]
env1 <- env[10:18,2:25]
matrix <- cor(otu1,env1,method="pearson",use="everything")
dim(matrix)
write.csv(matrix,file="root_env_correlation.csv",quote=F,row.names=T)

#3.root~soil
otu1 <- otu[1:9,]
otu2 <- otu[10:18,]
matrix <- cor(otu1,otu1,method="pearson",use="everything")
dim(matrix)
write.csv(matrix,file="root_soil_correlation.csv",quote=F,row.names=T)



## Five: OTU network data exaction
otu = read.table("OTU_table_average.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
otu <- t(otu)
#change rownames
row.names(otu) <- c("root1500", "root1500","root1500","root1900","root1900","root1900","root2300","root2300","root2300","soil1500","soil1500","soil1500","soil1900","soil1900","soil1900","soil2300","soil2300","soil2300")
#1.??????ת??Ϊ??Ԫ??
# ??Ϊ??Ԫ??
#????1???Ѿ?????????
library(Matrix)
smm = Matrix(otu)
smm
{re = summary(otu, i = c(1,1:18),j = c(1:231,1))
  sm = as.data.frame(re)
  sm
  
  re = summary(smm)
  sm = as.data.frame(re)
  sm
}

#????2?????????????Է?????
#mtrx2cols = function(m1,m2,val1,val2){
lt = lower.tri(m1)  #??ȡ?°???ΪTRUE???ϰ???ΪFALSE??ͬά?Ⱦ?????
res = data.frame(row = row(m1,as.factor = T)[lt],  #???ؾ???m1???°???Ԫ?ض?Ӧ??????
                 col = col(m1,as.factor = T)[lt],  #???ؾ???m1???°??Ƕ?Ӧ??????
                 val1 = m1[lt], val2= m2[lt]) #?????��λ?ȡ?????°??ǵ?Ԫ??
names(res)[3:4] = c(val1,val2) #?Ժ?��??????????֧?ֶ????????ϲ?
return(res)
}  #??ȡ?°???ΪTRUE???ϰ???ΪFALSE??ͬά?Ⱦ??? ???????Ծ???ʱʹ?ã???
#res = mtrx2cols(dst1,dst2,'Eco','logy')

#????3 ????ʹ??
#https://www.jianshu.com/p/a32d12a62819
#install.packages("reshape2")
library(reshape2)
melt(otu)
otu = melt(otu)

#2.filter links
otu_filter <- filter(otu, abs(value) > 0 )
otu_filter_1 <- otu_filter[,c(2,1)]
write.csv(otu_filter_1,file="otu_filter.csv",quote=F,row.names=T)








#???ģD:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/5_network???igraph)
library(psych)
library(dplyr)
library(tibble)
library(Hmisc)
setwd("F:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/5_network"orrelati)on analysis
otu = read.table("OTU_table_network.txt", header=T, row.names = 1, sep = "\t", comment.char = "")

#otu <- t(otu)
#View(otu)
#root_otu <- otu[1:9,]
otu_root <- otu[,109:182]
otu_soil <- otu[,1:108]
otu_1500 <- otu[1:3,]
otu_1900 <- otu[4:6,]
otu_2300 <- otu[7:9,]

occor_otu = corr.test(otu,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor_otu_root = corr.test(otu_root,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor_otu_soil = corr.test(otu_soil,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor_otu_1500 = corr.test(otu_1500,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor_otu_1900 = corr.test(otu_1900,use="pairwise",method="pearson",adjust="fdr",alpha=.05)
occor_otu_2300 = corr.test(otu_2300,use="pairwise",method="pearson",adjust="fdr",alpha=.05)


#4.2 extract R and P values
#
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame( row = rownames(cormat)[row(cormat)[ut]],
              column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}

otu_output=flattenCorrMatrix(occor_otu$r, occor_otu$p)
otu_root_output=flattenCorrMatrix(occor_otu_root$r, occor_otu_root$p)
otu_soil_output=flattenCorrMatrix(occor_otu_soil$r, occor_otu_soil$p)
otu_1500_output=flattenCorrMatrix(occor_otu_1500$r, occor_otu_1500$p)
otu_1900_output=flattenCorrMatrix(occor_otu_1900$r, occor_otu_1900$p)
otu_2300_output=flattenCorrMatrix(occor_otu_2300$r, occor_otu_2300$p)

otu_output <- unique(otu_output)
otu_root_output <- unique(otu_root_output)
otu_soil_output <- unique(otu_soil_output)
otu_1500_output <- unique(otu_1500_output)
otu_1900_output <- unique(otu_1900_output)
otu_2300_output <- unique(otu_2300_output)

#4.3 filter relations
otu_output_filter <- filter(otu_output, abs(cor) >= 0.6 & p <= 0.05)
otu_root_output_filter <- filter(otu_root_output, abs(cor) >= 0.6 & p <= 0.05)
otu_soil_output_filter <- filter(otu_soil_output, abs(cor) >= 0.6 & p <= 0.05)
otu_1500_output_filter <- filter(otu_1500_output, abs(cor) >= 0.6 & p <= 0.05)
otu_1900_output_filter <- filter(otu_1900_output, abs(cor) >= 0.6 & p <= 0.05)
otu_2300_output_filter <- filter(otu_2300_output, abs(cor) >= 0.6 & p <= 0.05)

write.csv(otu_output_filter, 'otu_o,t"output_filter.csv"csv(otu_root_output_filter, 'otu_roo"otu_roo_filter.csv"csv(otu_)
write.soil_output_filter, 'otu_soil_output_filter.csv')
write.csv(otu_1500_output_filter, 'otu_1500_output_filter.csv')
write.csv(otu_1900_output_filter, 'otu_1900_output_filter.csv')
write.csv(otu_2300_output_filter, 'otu_2300_output_filter.csv')

#CYTOSCAPE visualization






# set the working directory
setwd("F:/Dropox/luohuanhuanhuan/altitude_fungal_????/????/2_correlation_2")
#??????
?correlation analysis//zhuanlan.zhihu.com/p/114008006
#https://blog.csdn.net/lalaxumelala/article/details/86084040
#https://zhuanlan.zhihu.com/p/183880858
#chang color:https://www.thinbug.com/q/30743983
#install.packages("corrplot")
library(corrplot)
#otu = read.table("genus_abundance_root.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
otu_root = read.csv("otu_table.Genus.absolute_root.csv",header = TRUE,sep = ",", row.names = 1)
otu_soil = read.csv("otu_table.Genus.absolute_soil.csv",header = TRUE,sep = ",", row.names = 1)
otu_root <- t(otu_root)
otu_soil <- t(otu_soil)
otu_soil <- otu_soil[,1:47]
otu_root <- otu_root[,1:28]

#otu <- t(otu)
View(otu)

#env = read.table("envronmental_variables.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
env = read.csv("envronmental variables.csv",header = TRUE,sep = ",", row.names = 1)

env_soil <- env[1:9,]
env_root <- env[10:18,]
View(env)


matrix_soil_env <- cor(otu_soil,env_soil,method="pearson",use="everything")
matrix_root_env <- cor(otu_root,env_root,method="pearson",use="everything")
matrix_root_root <- cor(otu_root,otu_root,method="pearson",use="everything")
matrix_soil_soil <- cor(otu_soil,otu_soil,method="pearson",use="everything")
matrix_soil_root <- cor(otu_soil,otu_root,method="pearson",use="everything")

dim(matrix)

write.csv(matrix_soil_env,file="matrix_soil_env_correlation.csv",quote=F,row.names=T)
write.csv(matrix_root_env,file="matrix_root_env_correlation.csv",quote=F,row.names=T)
write.csv(matrix_root_root,file="matrix_root_root_correlation.csv",quote=F,row.names=T)
write.csv(matrix_soil_soil,file="matrix_soil_soil_correlation.csv",quote=F,row.names=T)
write.csv(matrix_soil_root,file="matrix_soil_root_correlation.csv",quote=F,row.names=T)




#correlation plot
#pdf(file = "soil.pdf",width=600,height=600)
png(file="Fig.S1_env_root_corr.png",width=8000,height=8000)
png(file="Fig.S2_env_soil_corr.png",width=8000,height=8000)
png(file="Fig.S4_soil_soil_corr.png",width=8000,height=8000)
png(file="Fig.S5_root_root_corr.png",width=10000,height=10000)
png(file="Fig.S6_soil_root_corr.png",width=10000,height=10000)

#?????٣????ߵ?ͼ
#corrplot(matrix, type = "lower", order = "original", tl.col = "black", tl.srt = 45,cl.lim = c(0, 1))
#???ݶ࣬??????ͼ
#corrplot(matrix,title = "Correlation Plot",method = "color",outline = T, addgrid.col = "white", order="hclust", mar = c(0,0,1,0), addrect = 4, rect.col = "grey", rect.lwd = 1,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 3)

corrplot(matrix_root_env,tl.cex = 12, cl.cex = 16, order = "original",cl.pos = "b",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
corrplot(matrix_soil_env,tl.cex = 10, cl.cex = 10, order = "original",cl.pos = "b",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))

corrplot(matrix_soil_soil,tl.cex = 10, cl.cex = 13, order = "original",cl.pos = "r",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
corrplot(matrix_soil_root,tl.cex = 12, cl.cex = 13, order = "original",cl.pos = "b",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))
corrplot(matrix_root_root,type = "upper",tl.cex = 12, cl.cex = 15, order = "original",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))

dev.off()


#corrplot(matrix, type = "full",tl.cex = 1, cl.cex = 3)
help("corrplot")


#Ⱥ??????F:/DropbD:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/ation_2')

library(ggplot2)
otu_root = read.csv("abundance_Genus.absolute_root.csv",header = TRUE,sep = ",")
otu_soil = read.csv("abundance_Genus.absolute_soil.csv",header = TRUE,sep = ",")


#????ת??
otu_root
otu_root$richness=log2(otu_root$richness)
otu_root$richness=log10(otu_root$richness)
otu_soil$richness=log2(otu_soil$richness)
otu_soil$richness=log10(otu_soil$richness)
#????abundance
otu_soil$richness <- lapply(otu_soil$richness, function(x) x / 1000)

#set order
otu_root$Taxonomy <- factor(otu_root$Taxonomy, levels= c("Unspecified_Agaricomycetes","Unspecified_Helotiales","Phialophora","Piloderma","Unspecified_Ascomycota","Pseudogymnoascus","Lactarius","Russula","Unspecified_Embryophyta","Helvella","Unspecified_Cupressales","Taxus","Unspecified_Pezizales","Tricholoma","Unspecified_Pleosporales","Cortinarius","Chlorencoelia","Termitomyces","Craterocolla","Unspecified_Hydnodontaceae","Cladophialophora","Unspecified_Pezizaceae","Didymella","Clavulina","Unspecified_Saccharomycetales","Penicillium","Unspecified_Boletaceae","Metarhizium"),ordered = T)
otu_soil$Taxonomy <- factor(otu_soil$Taxonomy, levels= c("Unspecified_Agaricomycetes","Piloderma","Unspecified_Helotiales","Phialophora","Lactarius","Termitomyces","Unspecified_Ascomycota","Pseudogymnoascus","Tricholoma","Russula","Cortinarius","Unclassified_Mortierellales","Eimeriidae","Archaeorhizomyces","Unspecified_Pleosporales","Unspecified_Hydnodontaceae","Penicillium","Inocephalus","Unspecified_Pezizales","Clavulina","Unspecified_Eurotiales","Unspecified_Sordariales","Chaetomium","Cladophialophora","Macrolepiota","Metarhizium","Craterocolla","Cryptococcus","Dipodascopsis","Didymella","Geminibasidium","Unspecified_Tremellales","Unspecified_Mytilinidiales","Cladosporium","Unclassified_Tremellales","Helvella","Unspecified_Tremellaceae","Umbelopsis","Unspecified_Pezizaceae","Unspecified_Saccharomycetales","Anomoporia","Sarcosphaera","Glaciozyma","Unspecified_Capnodiales","Multiclavula","Peziza","Pseudoplatyophyra"),ordered = T)
otu_root$sample <- factor(otu_root$sample, levels= c("root23003","root23002","root23001","root19003","root19002","root19001","root15003","root15002","root15001"),ordered = T)
otu_soil$sample <- factor(otu_soil$sample, levels= c("soil23003","soil23002","soil23001","soil19003","soil19002","soil19001","soil15003","soil15002","soil15001"),ordered = T)


mm <- ggplot(otu_soil,aes(x=sample,y=Taxonomy),order(na.last = TURE)) + 
  geom_point(alpha = 1.75, aes(color=richness/1000, size=richness/1000)) + 
  scale_color_gradient(low = "white",high = "red",space = "Lab",na.value = "white") +
  labs(x = "Sample", y = "Taxonomy", fill = "Richness") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) 

mm



p <- ggplot(otu_root,aes(x=sample,y=Taxonomy),order(na.last = TURE))
p + geom_point(aes(color=richness/1000, size=richness/1000)) +
  scale_color_continuous(low="white",high="red") 

p <- ggplot(otu_soil,aes(x=sample,y=Taxonomy),order(na.last = TURE))
p + geom_point(aes(color=richness, size=richness)) +
  scale_color_gradient(low = "white",high = "red",space = "Lab",na.value = "white")









###mantel test
#https://cloud.tencent.com/developer/article/1667543
#https://www.jianshu.com/p/6de6debbcd67
#https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
#https://jkzorz.github.io/2019/07/08/mantel-test.html
library(vegan)
setwd('D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/6_differential_analysis')

#1. read the table
df = read.table("soil+root+env.txt", header=T, row.names = 1, sep = "\t", comment.char = "")

#2. calculate the distance
abund_root <- df[,137:201]
abund_soil <- df[,23:136]
altitude <- df[,1]
enzyme <- df[,13:16]
nutrient_0 <- df[,5:12]
env <- df[,2:4]
colonize <- df[,18:21]
soil_N <- nutrient_0[,c(6:7)]
soil_P <- nutrient_0[,3:4]

#ѭ??ɸѡ????
nutrient <- nutrient_0[,-c(1,3,4,8)]
nutrient <- nutrient_0[,3:4]

#2.1 taxa abundance distance (Bray-curtis)
dist.abund_soil <- vegdist(abund_soil, method = 'bray')
dist.abund_root <- vegdist(abund_root, method = 'bray')
dist.colonize <- vegdist(colonize, method = 'bray')
dist.enzyme <- vegdist(enzyme, method = 'bray')

#2.2 env distance
#single
dist.altitude <- dist(altitude, method = 'euclidean')

#multi
scale.env <- scale(env, center = TRUE, scale = TRUE)
dist.env <- dist(scale.env, method = 'euclidean')

scale.soil_N <- scale(soil_N, center = TRUE, scale = TRUE)
dist.soil_N <- dist(scale.soil_N, method = 'euclidean')

scale.soil_P <- scale(soil_P, center = TRUE, scale = TRUE)
dist.soil_P <- dist(scale.soil_P, method = 'euclidean')

#2.3 mantel test
mantel <- mantel(dist.soil_N,dist.abund_soil,method = 'pearson',permutations = 9999, na.rm = TRUE)
mantel

abund_soil_altitude <- mantel(dist.abund_soil,dist.altitude,method = 'pearson',permutations = 9999, na.rm = TRUE)
abund_soil_altitude
abund_soil_env <- mantel(dist.abund_soil,dist.env,method = 'pearson',permutations = 9999, na.rm = TRUE)
abund_soil_env
abund_soil_nutrient <- mantel(dist.abund_soil,dist.nutrient,method = 'pearson',permutations = 9999, na.rm = TRUE)
abund_soil_nutrient

##3. plotting
library(ggplot2)

#3.0 mentel test heatmap
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
install.packages("ComplexHeatmap")
mantel_test_result = read.table("mantal_test_result.txt", header=T, row.names = 1, sep = "\t", comment.char = "")
mantel_test_result <- as.matrix.data.frame(mantel_test_result)

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)#?????ÿ?????ɫ
my_palette <- colorRampPalette(c("green", "white", "red"))(n = 1000) #?????Լ?????ɫ??
hM <- format(round(mantel_test_result, 2))#?????ݱ???2λС??

#plot
pheatmap(mantel_test_result,display_numbers = TRUE,order = "original",type = "upper",cluster_rows = FALSE,cluster_cols = FALSE)


M<-cor(mtcars)
corrplot.mixed(mantel_test_result)
corrplot(mantel_test_result,tl.cex = 3, cl.cex = 3, order = "original",cl.pos = "b",cl.lim=c(-1,1), col=colorRampPalette(c("blue","white","red"))(200))


heatmap.2(mantel_test_result,cluster_rows = FALSE,cluster_cols = FALSE,
          trace="none",#????ʾtrace
          col=coul,#?޸???ͼ??ɫ
          density.info = "none",#ͼ??ȡ??density
          key.xlab ='mantel_test',
          key.title = "",
          cexRow = 1,cexCol = 1,#?޸ĺ???????????
          Rowv = F,Colv = F, #ȥ??????
          margins = c(6, 6),
          cellnote = hM,notecol='black'#????????ϵ????ֵ???޸???????ɫ
)


heatmap.2(mantel_test_result,cluster_rows = FALSE,cluster_cols = FALSE,
          trace="none",#????ʾtrace
          col=my_palette,#ʹ???Լ?????ɫ??
          density.info = "none",#ͼ??ȡ??density
          key.xlab ='mantel_test',
          key.title = "",
          cexRow = 1,cexCol = 1,#?޸ĺ???????????
          Rowv = F,Colv = F, #ȥ??????
          margins = c(6, 6),
          cellnote = hM,notecol='black'#????????ϵ????ֵ???޸???????ɫ
)

help(heatmap.2)


#3.1 #ĳ???????¶ȵ??????ԣ??????¶ȣ????????ַ??ȣ???ɫ??ʾ??????γ??

xx = ggplot(df, aes(x = ST, y = Soil_Piloderma)) + 
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') + 
  geom_point(aes(colour = SMC), size = 4) +
  labs(y = 'Soil_Piloderma', x = 'Temperature (C)') + 
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12), 
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'), 
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx
#????ͼ?е????Իع?
fit <- lm(df$ST~df$Soil_Piloderma)
summary(fit)

#3.2 #????ͼչʾ??????��???????ԣ??????¶ȣ?????Ϊ???????ֵķ??ȣ???ɫ??ʾ??????γ??
library(reshape2)

otus <- df[ ,1:11]
otus_melt <- melt(otus, id = c('Station', 'Salinity', 'Temperature', 'Oxygen', 'Nitrate', 'Latitude', 'Longitude'))

xx <- ggplot(otus_melt, aes(x = Temperature, y = value)) + 
  facet_wrap(.~variable, scales = 'free_y') +
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') + 
  geom_point(aes(colour = Latitude), size = 4) +
  labs(y = 'Relative Abundance (%)', x = 'Temperature (C)') + 
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12), 
         axis.text.y = element_text(face = 'bold', size = 10, colour = 'black'), 
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         legend.text = element_text(size = 10, face = 'bold', colour = 'black'), 
         legend.position = 'top', strip.background = element_rect(fill = 'grey90', colour = 'black'),
         strip.text = element_text(size = 9, face = 'bold')) +
  scale_colour_continuous(high = 'navy', low = 'salmon')

xx

##3.3???????ȼ?????????
#?????Ļ??õľ??????ȣ?ת??Ϊ???ݿ???һһ??Ӧ??��
altitude <- as.vector(dist.altitude)
abund.root <- as.vector(dist.abund_root)
abund.soil <- as.vector(dist.abund_soil)
envs <- as.vector(dist.env)
colonizes <- as.vector(dist.colonize)
enzymes <- as.vector(dist.enzyme)
soil.N <- as.vector(dist.soil_N)
soil.P <- as.vector(dist.soil_P)


mat <- data.frame(altitude, abund.root, abund.soil, envs, colonizes, enzymes, soil.N, soil.P)

#???????ַ??ȵľ??????????¶?ָ???ľ???֮??????????ɢ??ͼ????????֪???????????أ?ͬʱ??ɫ??ʾ??????????????
mm <- ggplot(mat, aes(y = abund.root, x = envs)) + 
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = altitude)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  labs(x = "Difference in Environmental Variables", y = "Bray-Curtis Dissimilarity", fill = "Altitude gradient (m)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue")

mm


fit <- lm(mat$abund.soil~mat$soil.N)
fit <- lm(mat$abund.soil~mat$soil.P)
fit <- lm(mat$abund.soil~mat$envs)
fit <- lm(mat$abund.root~mat$soil.N)
fit <- lm(mat$abund.root~mat$soil.P)
fit <- lm(mat$abund.root~mat$envs)

summary(fit)



#???????ַ??ȵľ???????????????????֮??????????ɢ??ͼ????????֪????????????
mm <- ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 3, alpha = 0.5) + 
  labs(x = "Physical separation (km)", y = "Bray-Curtis Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))

mm













###SEM-PLS-PM·??????
# Xiaofei Gao, Huihuang Chen, Lynn Govaert,  Wenping Wang, Jun Yang. (2019). Responses of zooplankton body size and community trophic  structure to temperature change in a subtropical reservoir. Ecology and Evolution, 22(9), 12544-12555
# https://mp.weixin.qq.com/s?__biz=MzUzMjA4Njc1MA==&mid=2247493909&idx=1&sn=17b110d253921325e3ac26b657f92744&chksm=faba03a4cdcd8ab24de5d4955311fafa450646ad4adbc2627b805887de0e1d70f86fde0dab13&mpshare=1&scene=1&srcid=1122WglaP0z9AQRdJFF6RpXD&sharer_sharetime=1606010323196&sharer_shareid=5e26caa3a2ed75e9d7c29e213ad5207a&key=d073af33ae0d2ce990573d5b62c77e680802ae6f972fad4efb7b4905fbee337039c91d7fa78dbb6316e547c20a88ed317392af0628b788f99b6daa96d904aaaa87bb1e10bfa562ae3be3955fe967bb4b92ce5a58064edaa68c5dca52d5a68839a28c1926808e0b5a58b77ff883b504ca7f8e29d5a2bd6d3534dcc9163edcb48b&ascene=1&uin=MjY2MjY3OTUzMw%3D%3D&devicetype=Windows+7+x64&version=6300002f&lang=en&exportkey=AxLDj5RAely3brnBkDzr6ik%3D&pass_ticket=gpcotgpJ%2BOzYT5eQ0CNOlrCRbreZ96rDruHpJDAMuvrP%2F54qtxa4PUTfl4vj7oaH&wx_header=0

# install.packages("ggplot2")
# install.packages("plspm")
# install.packages("vegan")
# install.packages("ggrepel")
install.packages("tester")
install.packages("turner")
install.packages("diagram")
install.packages("shape")
install.packages("amap")

# R install plspm
# install "devtools"
install.packages("devtools") 
library(devtools)
# install "plspm"
install_github("gastonstat/plspm")
#download package
# wget https://cran.r-project.org/src/contrib/Archive/plspm/plspm_0.4.9.tar.gz
#mannual instal plspm
# INSTALL plspm_0.4.9.tar.gz

library(ggplot2)
library(plspm)
library(vegan)
library(ggrepel)


setwd('D:/Dropbox Folder/Dropbox/luohuanhuanhuan/altitude_fungal_analysis/R_analysis/6_differential_analysis')

#1. read the table
df = read.table("soil+root+env.txt", header=T, row.names = 1, sep = "\t", comment.char = "")

# 1. draw the full path
altitude = c(0,0,0,0,0,0,0,0)
env = c(1,0,0,0,0,0,0,0)
abund_root = c(0,1,0,0,0,0,0,0)
abund_soil = c(1,0,0,0,0,0,0,0)
enzyme = c(1,1,1,0,0,0,0,0)
colonize = c(1,1,1,1,1,0,0,0)
soil_N = c(1,1,1,1,1,1,0,0)
soil_P = c(1,1,1,1,1,1,1,0)
path_mat = rbind(abund_root, abund_soil, altitude, enzyme, env, colonize, soil_N,soil_P)
innerplot(path_mat)


# 2.???????????ӣ???��??????????VIF??<10??????20??
#ȥ??block??ģ?飩?ڲ????ӹ?????
DATA=df
spe.1 <- rda(DATA ~ SMC+ST+pH+TK+TN+TP+AP+AK+NO3.N+NH4.N+SOM, data = DATA)
vif.cca(spe.1)
spe.1 <- rda(DATA ~ Sucrase+Urase+Dehydrogenase+Phosphatase+Catalase, data = DATA)
vif.cca(spe.1)
spe.1 <- rda(DATA ~ Total_colonization_rate+ECMF_colonization_rate+DSE_colonization_rate+Microsclerotia_colonization_rate+Soil_mycelium_density, data = DATA)
vif.cca(spe.1)




abund_root <- df[,137:201]
abund_soil <- df[,23:136]
altitude <- df[,1]
enzyme <- df[,13:17]
nutrient_0 <- df[,5:12]
env <- df[,2:4]
colonize <- df[,18:21]
soil_N <- df[,c(6,10:11)]
soil_P <- df[,7:8]



#3.4 ·??????
#3.4.1 total
#????ÿ??ģ???ı?��??��???????ݴ??????ݱ??е?????????????????VIF<10
blocks=list(137:201, 23:136, 1, 13:16,  2:4, 18:21,c(10:11),7:8)
modes = c("A","A","A","A","A","A","A","A")
sat_pls = plspm(df, path_mat, blocks, modes=modes)
summary(sat_pls)

##Outer Model??????Loading??????0.7?????ݽ???????ȥ??ÿ??ģ????LoadingֵС??0.7?ı?��??ֱ?????б?��Loading > 0.7??????????·??????ģ??
blocks=list(2, 3:5, 6:7, 8, 9, 11, 13)
modes = c("A","A","A","A","A","A","A")
sat_pls = plspm(DATA, path_mat, blocks, modes=modes)
summary(sat_pls)

#Loading >0.7??????Loadingֵ??Ϊ??Loadingֵ????????????·??????ģ??
blocks=list(1, c(4,20), c(7,9,10), 11, 12, c(13,14,15,17,18), 19)
modes = c("A","A","A","A","A","A","A")
sat_pls = plspm(DATA, path_mat, blocks, modes=modes)
summary(sat_pls)

#3.4.2 soil
altitude = c(0,0,0,0,0,0,0,0)
env = c(1,0,0,0,0,0,0,0)
soil_N = c(0,1,0,0,0,0,0,0)
soil_P = c(0,1,0,0,0,0,0,0)
SOM = c(0,1,0,0,0,0,0,0)
abund_soil = c(0,1,1,1,1,0,0,0)
colonize = c(0,1,1,1,1,1,0,0)
enzyme = c(0,1,1,1,1,1,1,0)

path_mat = rbind(altitude, env, soil_N, soil_P, SOM, abund_soil, colonize, enzyme)
innerplot(path_mat)

blocks=list(1, 2:4, 10,7:8, 12, 23:136,   19:20, c(13,15:16))
modes = c("A","A","A","A","A","A","A","A")
sat_pls = plspm(df, path_mat, blocks, modes=modes)
summary(sat_pls)




#3.4.3 root
altitude = c(0,0,0,0,0,0,0,0)
env = c(1,0,0,0,0,0,0,0)
soil_N = c(0,1,0,0,0,0,0,0)
soil_P = c(0,1,0,0,0,0,0,0)
SOM = c(0,1,0,0,0,0,0,0)
abund_root = c(0,1,1,1,1,0,0,0)
colonize = c(0,1,1,1,1,1,0,0)
enzyme = c(0,1,1,1,1,1,1,0)

path_mat = rbind(altitude, env, soil_N, soil_P, SOM, abund_root, colonize, enzyme)
innerplot(path_mat)

blocks=list(1, 2:4, c(10:11),7:8, 12,  137:201,   19:20, c(13,15:16))
modes = c("A","A","A","A","A","A","A","A")
sat_pls = plspm(df, path_mat, blocks, modes=modes)
summary(sat_pls)



#MAKE PLOT
























rm(list=ls())

DATA=read.csv("AMOS.csv",row.names=1,header=T)
head(DATA)

# 1. 设置路径???
colonization = c(0,0,0,0,0,0,0,0)
Envir = c(0,0,0,0,0,0,0,0)
nutrient = c(1,0,0,0,0,0,0,0)
enzyme = c(1,1,1,0,0,0,0,0)
soil_div = c(1,1,1,1,0,0,0,0)
root_div = c(1,1,1,1,1,0,0,0)
soil_taxa = c(1,1,1,1,1,1,0,0)
root_taxa = c(1,1,1,1,1,1,1,0)

path_mat = rbind(colonization, Envir, nutrient, enzyme, soil_div, root_div, soil_taxa, root_taxa)
innerplot(path_mat)




#3.4 路径分析
##设置每个模块的变量（括号中数据代表数据表中的列数），膨胀因子VIF<10
blocks=list(1:4, 5:7, 8:17, 18:22, 23:31, 32:40, 41:95, 96:122)
modes = c("A","A","A","A","A","A","A","A")
sat_pls = plspm(DATA, path_mat, blocks, modes=modes)
summary = summary(sat_pls)


blocks=list(3:4, 5:6, c(8,10,11,13,14), 21:22, c(26:29,31), 36:37, c(42,47,48,53,55,58,61,63,64,66,68,69,72,75,78:82,85,87,89,91:93,95),c(99,104,107,110,118))
modes = c("A","A","A","A","A","A","A","A")
sat_pls = plspm(DATA, path_mat, blocks, modes=modes)
summary = summary(sat_pls)











m = read.table("D:/Desktop/ZHAO RUIHUA(1).txt", header = F, fill = T)

after <-m[,1:6]
for(i in 7:27) {
  tmp <- as.matrix(strsplit(as.character(m[,i]), '')) %>% str_replace("A", "1") %>% str_replace("C", "2") %>% str_replace("G", "3") %>% str_replace("T", "4") 
  after <- cbind(after, as.data.frame(tmp))
}
after[is.na(after)] = 0



tmp <- as.data.frame(strsplit(as.character(m[,i]), '')) %>% str_replace("A", "1") %>% str_replace("C", "2") %>% str_replace("G", "3") %>% str_replace("T", "4") 


m %>% separate(V10, sep = "", into = c("V28", "V29"), remove = FALSE)

m %>% separate(V10, into = c("V8", "V9"), sep = "")


cbind(m, stringr::str_split_fixed(as.character(m$V8), " ", ncols))




df <- data.frame(ID=11:13, FOO=c('a|b','b|c','x|y'))


after <-m[,1:6]
cbind(after, strsplit(as.character(m[,i]), ''))










#安装 plspm 包，由于plspm没有上CRAN，而在github上有，因此需通过github安装
install.packages('devtools') #安装github包
devtools::install_github('gastonstat/plspm') #如果没有安装
#加载 plspm 包和vegan包
library(plspm)
library(vegan)
data(varechem) #查看数据
#指定潜变量及关系
latent_r <- list(
  Ab = c('N', 'P', 'K'), 
  Rare = c('Ca', 'Mg', 'S'), 
  Ma= c('Fe', 'Mn', 'Zn',"Mo"), 
  SH = c('Baresoil', "Humdepth"),
  Soil = 'pH'
)
latent_r
