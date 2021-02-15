esc <- read.table("./data/ESC_TMRE_h5.tsv",header = T,sep="\t")
esc <- esc[,c(1,3,5,7,2,4,6,8)]
colSums(esc)
esc <- sweep(esc,2,colSums(esc),`/`)*10^6
head(esc)

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./esc_TMRE_col_8.csv", header=T,sep = ',')
head(coldata)

esc_all <- round(esc)
nrow(esc_all) #48177
head(esc_all)
dds <- DESeqDataSetFromMatrix(esc_all,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","batch")) 


# ggplot2-------
esc_pca <- plotPCA( dds_rLogTransformed ,intgroup=c( "condition","batch"),returnData=T)
colnames(esc_pca) <- c("PC1: 64% variance","PC2: 27% variance","group","TMRE","replicate","name")
# + geom_text(aes(label=name),vjust=2)
# use ggplot2 to plot it
esc_pca$replicate <- as.factor(esc_pca$replicate)
esc_pca$TMRE <- c(rep("high",4),rep("low",4))

p <- ggplot(esc_pca,aes(`PC1: 64% variance`,`PC2: 27% variance`)) +
  geom_point(aes(shape=TMRE,color=replicate),size=3) +
  labs(#title = "ESCs TMRE data before removing batch effect", 
    x="PC1: 64% variance", y="PC2: 27% variance") +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # set legend point size
  guides(colour = guide_legend(override.aes = list(size=5))) +
  guides(shape = guide_legend(override.aes = list(size=5))) +
  guides(size = 'none') + 
  theme(
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
    #,legend.position="none"
  )


p

library(cowplot)
ggsave2('E:/Project/2019__proliferation/figures/20200929 figures code and table/Figure6/tmre_esc_pca.pdf', width = 80, height = 60, units = "mm")

