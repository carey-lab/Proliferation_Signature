
# 20190605 CFSE doubling time with measured doubling time

library(flowCore)

# read in all files under a folder------------------
filenames <- list.files("temp", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
res <- lapply(ldf, summary)
names(res) <- substr(filenames, 6, 30)
###########################################################

# fib data-------------------------------------
library(openxlsx)
fib_data <- read.xlsx("./FIB_Sorting_fcs/DS.xlsx")
colnames(fib_data)

fib_data <- fib_data[,c(2,41)]
head(fib_data)

calculate_doubling_time <- function(data,name) {
  data <-  data[grep(name,data$filename),]
  data$time <- seq(24,24*5,24)
  result <- lm(log2(data$median_FITC_A)~data$time)
  return(-1/result$coefficients[2])
}

all_results <- c()
all_names <- c("H1","H2","H3","H4",
               "HA","HB","HC","HD",
               "M1","M2","M3","M4",
               "MA","MB","MC","MD",
               "L1","L2","L3","L4",
               "LA","LB","LC","LD")

for (i in all_names) {
  temp <- calculate_doubling_time(fib_data,paste0(i,".fcs"))
  names(temp) <- i
  all_results <- c(all_results,temp)
}
all_results

all_results <- as.data.frame(all_results)

all_results$conditions <- c("HL0","HL0","HL1","HL1",
                        "HN0","HN0","HN1","HN1",
                        "ML0","ML0","ML1","ML1",
                        "MN0","MN0","MN1","MN1",
                        "LL0","LL0","LL1","LL1",
                        "LN0","LN0","LN1","LN1")

all_results$TMRE <- substr(all_results$conditions, start = 1, stop = 1)
all_results$O2 <- substr(all_results$conditions, start = 2, stop = 2)
all_results$VitC <- substr(all_results$conditions, start = 3, stop = 3)


all_results$conditions <- factor(all_results$conditions,levels = unique(all_results$conditions))
all_results$TMRE <- factor(all_results$TMRE,levels = unique(all_results$TMRE))

all_results$VitC <- factor(all_results$VitC)



library(ggplot2)
library(ggbeeswarm)
fib_all_results <- ggplot(data = all_results, aes(y = all_results, x = conditions)) +

  geom_beeswarm(aes(col=TMRE, shape=O2, alpha=VitC),size=3) +

  scale_alpha_discrete(range = c(0.5,1)) +

  scale_colour_manual(values = c("#c51b8a", "#31a354", "#3182bd")) +

  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +

  labs(x="",y="Doubling time") +

  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         alpha = guide_legend(order = 3)) +

  theme(
    plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
    axis.title.x = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=18),
    #delete background
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank(),

    #加上坐标轴
    axis.line = element_line(colour = "black",size=0.5),
    #刻度线
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(-0.2,"cm"),

    #x轴标签
    axis.text.x = element_text(size=14,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=14,color='black',margin = margin(r=0.3, unit = "cm")),

    #边框
    panel.border = element_rect(colour = "black", fill=NA, size=1),

    #图例
    #legend.title = element_text(colour="black", size=14),
    #legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour="black", size=10)
    # remove legend
    ,legend.justification=c(0.5,0),
    legend.background = element_rect(fill="white",
                                     size=0.5, linetype="solid",
                                     colour ="white")
    ,legend.position="bottom"
    #,legend.position=c(0.8,0.1)
  )
fib_all_results
ggsave(fib_all_results, filename = "fib_treatment.pdf", width = 8, height = 6)


# esc data--------------------------------------
library(openxlsx)
esc_data <- read.xlsx("./ESC_Sorting_fcs/DS.xlsx")
colnames(esc_data)

esc_data <- esc_data[,c(2,42)]
head(esc_data)

calculate_doubling_time <- function(data,name) {
  data <-  data[c(2,grep(name,data$filename)),]
  print(data)
  data$time <- c(0,48,72)
  result <- lm(log2(data$median_FITC_A)~data$time)
  return(-1/result$coefficients[2])
}

calculate_doubling_time(esc_data,paste0("P6 1",".fcs"))

grep(paste0("P6 1",".fcs"),esc_data$filename)

all_results <- c()
all_names <- c("P6 1","P6 2","P6 AA","P6 AA2",
               "P6 LOW 1","P6 LOW 2","P6 LOW AA1","P6 LOW AA2",
               "P7 1","P7 2","P7 AA","P7 AA 2",
               "P7 LOW 1","P7 LOW 2","P7 LOW AA","P7 LOW AA2",
               "P8 1","P8 2","P8 AA","P8 AA2",
               "P8 LOW 1","P8 LOW 2","P8 LOW AA","P8 LOW AA2")

for (i in all_names) {
  temp <- calculate_doubling_time(esc_data,paste0(i,".fcs"))
  names(temp) <- i
  all_results <- c(all_results,temp)
}
all_results

all_results <- as.data.frame(all_results)

all_results$conditions <- c("LN0","LN0","LN1","LN1",
                            "LL0","LL0","LL1","LL1",
                            "MN0","MN0","MN1","MN1",
                            "ML0","ML0","ML1","ML1",
                            "HN0","HN0","HN1","HN1",
                            "HL0","HL0","HL1","HL1")

all_results <- all_results[order(all_results$conditions),]
all_results <- all_results[c(1:8,17:24,9:16),]



all_results$TMRE <- substr(all_results$conditions, start = 1, stop = 1)
all_results$O2 <- substr(all_results$conditions, start = 2, stop = 2)
all_results$VitC <- substr(all_results$conditions, start = 3, stop = 3)


all_results$conditions <- factor(all_results$conditions,levels = unique(all_results$conditions))
all_results$TMRE <- factor(all_results$TMRE,levels = unique(all_results$TMRE))

all_results$VitC <- factor(all_results$VitC,levels = unique(all_results$VitC))



library(ggplot2)
library(ggbeeswarm)
esc_all_results <- ggplot(data = all_results, aes(y = all_results, x = conditions)) +

  geom_beeswarm(aes(col=TMRE, shape=O2, alpha=VitC),size=3) +

  scale_alpha_discrete(range = c(0.5,1)) +

  scale_colour_manual(values = c("#c51b8a", "#31a354", "#3182bd")) +

  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +

  labs(x="",y="Doubling time") +

  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         alpha = guide_legend(order = 3)) +

  theme(
    plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
    axis.title.x = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=18),
    #delete background
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
    panel.grid.minor = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
    panel.background = element_blank(),

    #加上坐标轴
    axis.line = element_line(colour = "black",size=0.5),
    #刻度线
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(-0.2,"cm"),

    #x轴标签
    axis.text.x = element_text(size=14,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=14,color='black',margin = margin(r=0.3, unit = "cm")),

    #边框
    panel.border = element_rect(colour = "black", fill=NA, size=1),

    #图例
    #legend.title = element_text(colour="black", size=14),
    #legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour="black", size=10)
    # remove legend
    ,legend.justification=c(0.5,0),
    legend.background = element_rect(fill="white",
                                     size=0.5, linetype="solid",
                                     colour ="white")
    ,legend.position="bottom"
    #,legend.position=c(0.8,0.1)
  )
esc_all_results
ggsave(esc_all_results, filename = "esc_treatment.pdf", width = 8, height = 6)


# ANNOVA to show the difference------------------------------

# fib  treat TMRE as continuous value
annova_fib <- read.table("./fib_dt.tsv",sep = "\t",header = T,check.names = F)
annova_fib$VitC <- as.factor(annova_fib$VitC)

annova_fib$TMRE_n <- c(rep(3,8),rep(2,8),rep(1,8))

res_annova_fib <- aov(doubling_time ~ TMRE_n, data = annova_fib)

summary(res_annova_fib)
writeLines(unlist(summary(res_annova_fib)),"./fib_annova.tsv",sep = "\t")

TukeyHSD(res_annova_fib)

# esc
annova_esc <- read.table("./esc_dt.tsv",sep = "\t",header = T,check.names = F)
annova_esc$VitC <- as.factor(annova_esc$VitC)

res_annova_esc <- aov(doubling_time ~ TMRE+O2+VitC, data = annova_esc)

summary(res_annova_esc)

TukeyHSD(res_annova_esc)


boxplot(annova_esc$doubling_time~annova_esc$TMRE)

points(unique(annova_esc$TMRE), aggregate(annova_esc$doubling_time,by=list(annova_esc$TMRE),FUN=mean)$x, col = "orange", pch = 18)

annova_esc$TMRE <- factor(annova_esc$TMRE,levels = c("L","M","H"))
library(ggpubr)
p <- ggboxplot(annova_esc, x = "TMRE", y = "doubling_time",
               color = "TMRE", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "TMRE")
p + stat_compare_means(comparisons = list(c("L","M")))





# 1. Homogeneity of variances------------
plot(res_annova_fib, 1)
# 2. Normality
plot(res_annova_fib, 2)

# ANOVA test with no assumption of equal variances-----
res_annova_fib_welch <- oneway.test(doubling_time ~ TMRE+O2+VitC, data = annova_fib)

summary(res_annova_fib_welch)

res_annova_fib <- aov(doubling_time ~ VitC, data = annova_fib)

summary(res_annova_fib)

# Pairewise t-test-------
pairwise.t.test(rep(1,10), rep(2,10),
                p.adjust.method = "BH",pool.sd = F)




# figure S5-------------------------------------------





setwd("E:/Project/2019__proliferation/14.new 20200925/Vocalno plot")

# apply deseq2################----------------------------------


# CFSE------
library(data.table)
fib <- read.table("./fib gene id counts.csv",sep=",",header=T)
head(fib)
fib <- na.omit(fib)
colSums(fib)

esc <- read.table("./ESC gene id counts.csv",sep=",",header=T)
head(esc)
esc <- na.omit(esc)
colSums(esc)

rownames(fib)==rownames(esc)


# # get human symbol name
# library(biomaRt)
#
# mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",host = "asia.ensembl.org", ensemblRedirect = FALSE)
# # get human symbol for mouse ensembl gene id
#
# human_symbol<- getBM(
#   mart=mart,
#   attributes=c("ensembl_gene_id", "gene_biotype", "description","external_gene_name","entrezgene"),
#   filter="ensembl_gene_id",
#   values=rownames(fib),
#   uniqueRows=T)
#
# head(human_symbol)
# saveRDS(human_symbol,"./human_symbol.rds")




# fibroblasts-----------------------------------------------------


library("DESeq2")

coldata<-read.table(file="./fib_coldata.tsv", header=T,sep = "\t")

head(coldata)

fib <- round(fib)

dds <- DESeqDataSetFromMatrix(fib,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(dds$id);head(dds$condition)

# EXTRACT NORMALIZED COUNTS
# dds <- estimateSizeFactors(dds)
# hos3_counts <- counts(dds, normalized=TRUE)
# colSums(hos3_counts)


# wt vs. deltahos3
res_fib = results(dds, contrast=c("condition", "fast", "slow"))

res_fib = res_fib[order(res_fib$padj),]

head(res_fib)
summary(res_fib)

table(res_fib$padj<0.05 & abs(res_fib$log2FoldChange) > 1)


# add gene symbol and gene description
human_symbol <- readRDS('./human_symbol.rds')
res_fib$gene_symbol <- human_symbol$external_gene_name[match(rownames(res_fib),human_symbol$ensembl_gene_id)]
res_fib$gene_description <- human_symbol$description[match(rownames(res_fib),human_symbol$ensembl_gene_id)]
res_fib$entrez <- human_symbol$entrezgene[match(rownames(res_fib),human_symbol$ensembl_gene_id)]

head(res_fib)


res_fib <- res_fib[,c(7,8,9,1:6)]

write.table(res_fib,file="./fib_CFSE_deseq2.tsv",sep = "\t",col.names = T,row.names = T,quote = F)


#------------------------------------------------------

ggplot_theme <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14),
  #delete background
  panel.grid.major = element_line(colour="grey92", size = (0.5)),
  panel.grid.minor = element_line(size = (0.5), colour="grey92"),
  # panel.grid.major =element_blank(),
  # panel.grid.minor = element_blank(),
  panel.background = element_blank(),

  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  #加上坐标轴
  axis.line = element_line(colour = "black",size=0.5),
  #刻度线
  axis.ticks = element_line(size=0.5),
  #axis.ticks.length=unit(-0.2,"cm"),
  axis.ticks.length.x = unit(-0.2,"cm"),
  axis.ticks.length.y = unit(-0.2,"cm"),
  #x轴标签
  axis.text.x = element_text(size=14,color='black',margin = margin(t=0.3, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=14,color='black',margin = margin(r=0.3, unit = "cm")),
  #图例
  legend.title = element_blank(),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=10)
  # remove legend
  ,legend.position="none"
)

#---------------------------------
volcano_plot <- function(res) {
  library("ggplot2")

  diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

  dim(diff_gene_deseq2)

  res$Significant<-'Not change'

  res$Significant[res$padj< 0.05 & res$log2FoldChange>1]<-'Up'

  res$Significant[res$padj < 0.05 & res$log2FoldChange<(-1)]<-'Down'

  head(res)

  res <- as.data.frame(res)


  res$log2FoldChange[res$log2FoldChange > 5] <- 5
  res$log2FoldChange[res$log2FoldChange < -5] <- -5
  res$padj[! is.na(res$padj) & -log10(res$padj) > 5] <- 10^(-5)

  ggplot(data=res,aes(x=log2FoldChange,y=-log10(padj)))+

    geom_point(size=0.5,aes(color=Significant))+

    #geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+

    coord_cartesian(xlim=c(-5,5),ylim=c(0,5))+

    scale_color_manual(values=c("grey","grey","grey"),guide=guide_legend(title="Significant"))+

    labs(x="Log2FoldChange",y="-Log10(adjusted p-value)")+

    guides(colour = guide_legend(override.aes = list(size=4))) +

    ggplot_theme
}





volcano_plot(res_fib)


ggsave2('E:/Project/2019__proliferation/14.new 20200925/Vocalno plot/20201001 vocalno fib cfse.pdf', width = 80, height = 80, units = "mm")


###############################
############check gene sets#####################################
res <-  res_fib
head(res)


# add more proteasome and ribosome gene set
proteasome_core_complex <- c("26441","26440","73902","19177","19173","19172","16912","26445","19175","19171","26443","16913","26442","19170","26446","73677","26444","19167","19166")

no_col <- max(count.fields("./ribosome&protesome_mouse_entrez.tsv", sep = "\t"))
ribi <- read.table("./ribosome&protesome_mouse_entrez.tsv",sep="\t",fill=TRUE,col.names=1:no_col)
head(ribi$X1)



RIBOSOME_BIOGENESIS <- as.character(ribi[1,][-c(1,2)])
PRERIBOSOME <- as.character(ribi[2,][-c(1,2)])
PROTEASOME_COMPLEX <- as.character(ribi[3,][-c(1,2)])
BIOCARTA_PROTEASOME_PATHWAY <- as.character(ribi[4,][-c(1,2)])
CYTOSOLIC_RIBOSOME <- as.character(ribi[5,][-c(1,2)])



library("ggplot2")

diff_gene_deseq2 <- subset(res_fib, padj < 0.05 & abs(log2FoldChange) > 1)

dim(diff_gene_deseq2)

res$Significant<-'Not change'

res$Significant[res$padj< 0.05 & res$log2FoldChange>1]<-'Up'

res$Significant[res$padj < 0.05 & res$log2FoldChange<(-1)]<-'Down'

res$Significant[res$entrez %in% proteasome_core_complex]<-'v_proteasome_core_complex'

res$Significant[res$entrez %in% RIBOSOME_BIOGENESIS]<-'w_RIBOSOME_BIOGENESIS'

res$Significant[res$entrez %in% PRERIBOSOME]<-'x_PRERIBOSOME'

res$Significant[res$entrez %in% PROTEASOME_COMPLEX]<-'y_PROTEASOME_COMPLEX'

res$Significant[res$entrez %in% BIOCARTA_PROTEASOME_PATHWAY]<-'y_BIOCARTA_PROTEASOME_PATHWAY'

res$Significant[res$entrez %in% CYTOSOLIC_RIBOSOME]<-'y_CYTOSOLIC_RIBOSOME'



head(res)

res <- as.data.frame(res)


ggplot(data=res,aes(x=log2FoldChange,y=-log10(padj)))+

  geom_point(size=0.5,aes(color=Significant))+




  geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+

  coord_cartesian(xlim=c(-7,7),ylim=c(0,30))+

  scale_color_manual(values=c("blue","grey","red","green"),guide=guide_legend(title="Significant"))+

  labs(x="Log2FoldChange",y="-Log10(padj)")+

  theme_bw(base_size=15) +

  geom_point(data= res[res$Significant=="y_CYTOSOLIC_RIBOSOME",],size=2,color="green")

#########################################################


# esc------------------------------------


library("DESeq2")

coldata<-read.table(file="./esc_coldata.tsv", header=T,sep = "\t")

head(coldata)

esc <- round(esc)

dds <- DESeqDataSetFromMatrix(esc,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(dds$id);head(dds$condition)

# EXTRACT NORMALIZED COUNTS
# dds <- estimateSizeFactors(dds)
# hos3_counts <- counts(dds, normalized=TRUE)
# colSums(hos3_counts)


# wt vs. deltahos3
res = results(dds, contrast=c("condition", "fast", "slow"))

res = res[order(res$padj),]

head(res)
summary(res)

table(res$padj<0.05 & abs(res$log2FoldChange) > 1)
# add gene symbol and gene description

res$gene_symbol <- human_symbol$external_gene_name[match(rownames(res),human_symbol$ensembl_gene_id)]
res$gene_description <- human_symbol$description[match(rownames(res),human_symbol$ensembl_gene_id)]
res <- res[,c(7,8,1:6)]

write.table(res,file="./esc_CFSE_de2.tsv",sep = "\t",col.names = T,row.names = T,quote = F)

volcano_plot(res)

ggsave2('E:/Project/2019__proliferation/14.new 20200925/Vocalno plot/20201001 vocalno esc cfse.pdf', width = 80, height = 80, units = "mm")
# TMRE -------------------------------

fib_TMRE <- read.table("./FIB_TMRE_h5.tsv",sep="\t",header=T)
head(fib_TMRE)
fib <- na.omit(fib_TMRE)
colSums(fib_TMRE)

fib_TMRE <- fib_TMRE[,-3]

coldata<-read.table(file="./fib_TMRE_coldata.tsv", header=T,sep = "\t")

head(coldata)

fib_TMRE <- round(fib_TMRE)

dds <- DESeqDataSetFromMatrix(fib_TMRE,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(dds$id);head(dds$condition)



# wt vs. deltahos3
res = results(dds, contrast=c("condition", "high", "low"))

res = res[order(res$padj),]

head(res)
summary(res)

table(res$padj<0.05 & abs(res$log2FoldChange) > 1)
# add gene symbol and gene dfib_TMREription

res$gene_symbol <- human_symbol$external_gene_name[match(rownames(res),human_symbol$ensembl_gene_id)]
res$gene_description <- human_symbol$description[match(rownames(res),human_symbol$ensembl_gene_id)]
res <- res[,c(7,8,1:6)]

write.table(res,file="./fib_TMRE_de2.tsv",sep = "\t",col.names = T,row.names = T,quote = F)

volcano_plot(res)


ggsave2('E:/Project/2019__proliferation/14.new 20200925/Vocalno plot/20201001 vocalno fib tmre.pdf', width = 80, height = 80, units = "mm")
# esc

esc_TMRE <- read.table("./esc_TMRE_h5.tsv",sep="\t",header=T)
head(esc_TMRE)
esc <- na.omit(esc_TMRE)
colSums(esc_TMRE)

esc_TMRE <- esc_TMRE[,-c(1:2)]

coldata<-read.table(file="./esc_TMRE_coldata.tsv", header=T,sep = "\t")

head(coldata)

esc_TMRE <- round(esc_TMRE)

dds <- DESeqDataSetFromMatrix(esc_TMRE,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(dds$id);head(dds$condition)



# wt vs. deltahos3
res = results(dds, contrast=c("condition", "high", "low"))

res = res[order(res$padj),]

head(res)
summary(res)

table(res$padj<0.05 & abs(res$log2FoldChange) > 1)
# add gene symbol and gene desc_TMREription

res$gene_symbol <- human_symbol$external_gene_name[match(rownames(res),human_symbol$ensembl_gene_id)]
res$gene_description <- human_symbol$description[match(rownames(res),human_symbol$ensembl_gene_id)]
res <- res[,c(7,8,1:6)]

write.table(res,file="./esc_TMRE_de2.tsv",sep = "\t",col.names = T,row.names = T,quote = F)

volcano_plot(res)


ggsave2('E:/Project/2019__proliferation/14.new 20200925/Vocalno plot/20201001 vocalno esc tmre.pdf', width = 80, height = 80, units = "mm")
#-------------------20190425-------------------------------------

fib_de <- read.table("./fib_de2.tsv", sep = "\t", header = T)
fib_de <- fib_de[! is.na(fib_de$log2FoldChange),]
head(fib_de)
tail(fib_de)

esc_de <- read.table("./esc_de2.tsv", sep = "\t", header = T)
esc_de <- esc_de[! is.na(esc_de$log2FoldChange),]
head(esc_de)
nrow(fib_de)


merge <- merge(fib_de,esc_de,by="gene_symbol")
nrow(merge)
head(merge)




fib_de_sign <- fib_de[abs(fib_de$log2FoldChange) >=1 & fib_de$padj<0.05,] #3166
fib_de_sign <- na.omit(fib_de_sign)
esc_de_sign <- esc_de[abs(esc_de$log2FoldChange) >=1 & esc_de$padj<0.05,] #2892
esc_de_sign <- na.omit(esc_de_sign)

library(tidyverse)
fib_de_sign <- as.tibble(fib_de_sign)
fib_de_sign <- arrange(fib_de_sign, desc(abs(log2FoldChange)) )

g <- filter(fib_de_sign,fib_de_sign$log2FoldChange>0)


write.table(fib_de_sign$entrez,"./fib_de_sign_entrez.tsv",sep = "\t",row.names = F,col.names = F,quote = F)


overlap_sign <- merge(fib_de_sign,esc_de_sign,by="gene_symbol")

overlap_sign <- na.omit(overlap_sign)

overlap_fib <- fib_de[match(esc_de$gene_symbol,fib_de$gene_symbol) ,]
overlap_fib <- na.omit(overlap_fib)

overlap_esc <- esc_de[match(fib_de$gene_symbol,esc_de$gene_symbol) ,]
overlap_esc <- na.omit(overlap_esc)


setdiff(overlap_fib$gene_symbol,overlap_esc$gene_symbol)


merge <- merge(fib_de,esc_de,by="gene_symbol")
head(merge)

temp <- merge[merge$log2FoldChange.x*merge$log2FoldChange.y>0,]
head(temp)

temp[temp$log2FoldChange.x > 0.5 & temp$log2FoldChange.y > 0.5 ,][,c(4,11)]




























