
#Feb 10 2019

# modify 20200821-----------------------------------
setwd("E:/Project/2019__proliferation/3.GSEA_explaination")

preribosome <- read.table("./GO_PRERIBOSOME.tsv",header = T,sep = "\t")
head(preribosome)

fib <- read.table("./FIB_CFSE_entrez.gct",header = T,sep = "\t")
head(fib)


esc <- read.table("./ESC_CFSE_entrez.gct",header = T,sep = "\t")
head(esc)


# FIGURE 3B FIB heatmap for GO_PRERIBOSOME---------------------
preribosome <- fib[match(preribosome$PROBE,fib$NAME),]
head(preribosome)



# entrez id to human symbol id---------------
# require("biomaRt")
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# 
# annote <- getBM(
#   mart=mart,
#   attributes=c("entrezgene","external_gene_name"),
#   filter="entrezgene",
#   values=preribosome$NAME,
#   uniqueRows=T)
# head(annote)
# annote[annote$external_gene_name=="FBLL1",]
# 
# saveRDS(annote,"./annote.rds")
###########################
annote <- readRDS("./annote.rds")

rownames(preribosome) <- annote$external_gene_name[match(preribosome$NAME,annote$entrezgene)]
preribosome <- preribosome[,c(3:10)]
library(pheatmap)
pheatmap(preribosome)


# order by signal2noise-------

signal2noise <- function(data) {
  a <- data[1:4]
  b <- data[5:8]
  res <- (mean(a)-mean(b))/(var(a)+var(b))
  
  return(res)
}

preribosome$signal2noise <- apply(preribosome,1,signal2noise)


preribosome <- preribosome[order(preribosome$signal2noise,decreasing = T),]


preribosome <- preribosome[,1:8]

preribosome_row <- t(scale(t(log(preribosome,2))))
preribosome_row[is.nan(preribosome_row)] <- 0




pdf("./FIB preribo heatmap 20200821.pdf", width = 7, height = 14)

pheatmap(preribosome_row,clustering_method = "complete",
         kmeans_k = NA,
         clustering_distance_rows = "euclidean",
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         breaks = NA, 
         cellwidth = 28, cellheight = 12,
         scale = "none",
         fontsize_col=15)

dev.off()




# FIGURE 3S B ESC heatmap for GO_PRERIBOSOME-----------------------

preribosome <- esc[match(preribosome$PROBE,esc$NAME),]

head(preribosome)


###########################
annote <- readRDS("./annote.rds")

rownames(preribosome) <- annote$external_gene_name[match(preribosome$NAME,annote$entrezgene)]
preribosome <- preribosome[,c(3:10)]
library(pheatmap)
pheatmap(preribosome)


# order by signal2noise--------------------------

preribosome$signal2noise <- apply(preribosome,1,signal2noise)


preribosome <- preribosome[order(preribosome$signal2noise,decreasing = T),]


preribosome <- preribosome[,1:8]

preribosome_row <- t(scale(t(log(preribosome,2))))
preribosome_row[is.nan(preribosome_row)] <- 0




pdf("./ESC preribo heatmap 20200821.pdf", width = 7, height = 14)

pheatmap(preribosome_row,clustering_method = "complete",
         kmeans_k = NA,
         clustering_distance_rows = "euclidean",
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         breaks = NA, 
         cellwidth = 28, cellheight = 12,
         scale = "none",
         fontsize_col=15)

dev.off()



# FIGURE 3S C myc-----------------------



setwd("E:/Project/2019__proliferation/14.new 20200925/correlation with proliferation for ESC and FIB/")
# 20200229 spearman correlation with CFSE


# read in fib data
fib_data <- read.table("./fib_all_data_20190517.tsv",header = T, sep='\t')
str(fib_data)

head(fib_data)
dim(fib_data)

esc_data <- read.table("./esc_all_data_20190517.tsv",header = T, sep='\t')
str(esc_data)

head(esc_data)
dim(esc_data)


# normalize by sequencing depth
colSums(fib_data)

colSums(esc_data)

fib_data_norm <- sweep(fib_data,2,FUN = '/',STATS = colSums(fib_data))*10^6
colSums(fib_data_norm)


esc_data_norm <- sweep(esc_data,2,FUN = '/',STATS = colSums(esc_data))*10^6
colSums(esc_data_norm)


all(rownames(fib_data_norm)==rownames(esc_data_norm))
# transfer to human entrez id


# require("biomaRt")
# mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

#View(listAttributes(mart1))

# human / mouse
# human_gene <- getLDS(attributes=c("ensembl_gene_id"),
#                      filters="ensembl_gene_id", values=rownames(fib_data_norm), mart=mart2,
#                      attributesL=c("entrezgene_id","external_gene_name"), martL=mart1)
# 
# 
# saveRDS(human_gene,'./mouse_human_id.rds')


human_gene <- readRDS("./mouse_human_id.rds")


fib_data_norm$human_entrez <- human_gene$NCBI.gene.ID[match(rownames(fib_data_norm),human_gene$Gene.stable.ID)]

dim(fib_data_norm)
fib_data_norm <- na.omit(fib_data_norm)

length(unique(fib_data_norm$human_entrez))

# get mean of same entrez
library(data.table)

fib_data_norm <- as.data.table(fib_data_norm)

fib_data_norm <- fib_data_norm[, lapply(.SD,mean), by=human_entrez]

dim(fib_data_norm)


esc_data_norm$human_entrez <- human_gene$NCBI.gene.ID[match(rownames(esc_data_norm),human_gene$Gene.stable.ID)]

dim(esc_data_norm)
esc_data_norm <- na.omit(esc_data_norm)

length(unique(esc_data_norm$human_entrez))

# get mean of same entrez
library(data.table)

esc_data_norm <- as.data.table(esc_data_norm)

esc_data_norm <- esc_data_norm[, lapply(.SD,mean), by=human_entrez]

dim(esc_data_norm)


# mytheme-----------------
ggplot_theme_pdf <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=10),
  axis.title.y = element_text(color="black", size=10),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=0.5),
  #刻度线
  axis.ticks = element_line(size=0.5),
  axis.ticks.length=unit(-0.1,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=10,color='black',margin = margin(t=0.2, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=10,color='black',margin = margin(r=0.2, unit = "cm")),
  #图例
  legend.title = element_text(colour="black", size=10),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=10)
  # remove legend
  ,legend.position="none"
)



# MYC 4609---------


MYC_fib <- fib_data_norm[fib_data_norm$human_entrez=='4609',][c(2:5,8:11)]
log2fc_MYC_fib <- log2(MYC_fib[5:8]/MYC_fib[1:4])

log2fc_MYC_fib <- as.data.frame(t(log2fc_MYC_fib))
colnames(log2fc_MYC_fib) <- 'log2FC'


MYC_esc <- esc_data_norm[esc_data_norm$human_entrez=='4609',][c(2:5,8:11)]
log2fc_MYC_esc <- log2(MYC_esc[5:8]/MYC_esc[1:4])

log2fc_MYC_esc <- as.data.frame(t(log2fc_MYC_esc))
colnames(log2fc_MYC_esc) <- 'log2FC'


all_myc <- rbind(log2fc_MYC_fib,log2fc_MYC_esc)
all_myc$condition <- substr(rownames(all_myc),1,3)

ggplot(data=all_myc,aes(x=condition,y=log2FC,color=condition)) + 
  geom_boxplot(size=0.5,outlier.shape = NA) +
  geom_point(size=3) +
  #geom_jitter(shape=16, position=position_jitter(0.2),size=0.5) +
  scale_color_manual(values =c("#7C285C","#2171B5")) +
  labs(x='',y='Log2FoldChange(Fast/Slow)') +
  ggplot_theme_pdf


ggsave2('E:/Project/2019__proliferation/14.new 20200925/correlation with proliferation for ESC and FIB/20201001 MYC expression.pdf', width = 60, height = 70, units = "mm")













