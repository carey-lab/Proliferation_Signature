
# FIGURE 4A  Scatter plot of all sorted genes----------

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






# calculate spearman correlation with CFSE-------
head(fib_data_norm)
fib_data_norm <- as.data.frame(fib_data_norm)
esc_data_norm <- as.data.frame(esc_data_norm)


CFSE <- c(1,2,3)

cor_fib_1 <- cor(t(fib_data_norm[,c(2,6,8)]),CFSE,method = 'spearman')
cor_fib_2 <- cor(t(fib_data_norm[,c(3,7,9)]),CFSE,method = 'spearman')

cor_fib_1[is.na(cor_fib_1),] <- 0 
cor_fib_2[is.na(cor_fib_2),] <- 0 

fib_data_norm$cor_fib_1 <- cor_fib_1
fib_data_norm$cor_fib_2 <- cor_fib_2
fib_data_norm$mean_cor_fib <- rowMeans(fib_data_norm[,c(12:13)])


cor_esc_1 <- cor(t(esc_data_norm[,c(2,6,8)]),CFSE,method = 'spearman')
cor_esc_2 <- cor(t(esc_data_norm[,c(3,7,9)]),CFSE,method = 'spearman')

cor_esc_1[is.na(cor_esc_1),] <- 0 
cor_esc_2[is.na(cor_esc_2),] <- 0 

esc_data_norm$cor_esc_1 <- cor_esc_1
esc_data_norm$cor_esc_2 <- cor_esc_2
esc_data_norm$mean_cor_esc <- rowMeans(esc_data_norm[,c(12:13)])



proliferation_signature <- readRDS('./proliferation_signature_combined_human.rds')


plot_data <- cor_table[,c(1,2,9)]

# sample proliferatino siganture
sample_data <- cbind(fib_data_norm[,c(1:3,6:9)],esc_data_norm[,c(2:3,6:9)])
# sample_data <- sample_data[sample_data$human_entrez %in% proliferation_signature$proliferation_signature,]

rownames(sample_data) <- sample_data$human_entrez
sample_data <- sample_data[,-1]
# ssGSEA get proliferatino siganture score

library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(GSEABase)



gene_set <- list(as.character(na.omit(proliferation_signature$proliferation_signature)))

gsva_matrix<- gsva(as.matrix(sample_data),
                   gene_set, # gene set 用list格式
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)


head(plot_data)

gsva_matrix <- rbind(gsva_matrix[c(1,3,5)],gsva_matrix[c(2,4,6)],gsva_matrix[c(7,9,11)],gsva_matrix[c(8,10,12)])

ps_mean <- mean(cor(t(gsva_matrix),CFSE,method = 'spearman'))


ps_mean <- c('prolife','prolife',ps_mean)

plot_data <- rbind(ps_mean,plot_data)

plot_data$mean_cor_fib_esc <- as.numeric(plot_data$mean_cor_fib_esc)

plot_data <- plot_data[order(plot_data$mean_cor_fib_esc,decreasing = T),]

plot_data$rank <- 1:nrow(plot_data)

plot(plot_data$rank,plot_data$mean_cor_fib_esc)

grep('PCNA',plot_data$human_symbol)
grep('MKI67',plot_data$human_symbol)
grep('prolife',plot_data$human_symbol)


plot_data$group <- 'Other genes'
plot_data$group[plot_data$human_symbol=='prolife'] <- 'Proliferation siganture'
plot_data$group[plot_data$human_symbol=='PCNA'] <- 'PCNA'
plot_data$group[plot_data$human_symbol=='MKI67'] <- 'KI67'

plot_data$size <- 0.5
plot_data$size[plot_data$human_symbol=='prolife'] <- 4
plot_data$size[plot_data$human_symbol=='PCNA'] <- 4
plot_data$size[plot_data$human_symbol=='MKI67'] <- 4


cor_table[cor_table$human_symbol=='PCNA',]
cor_table[cor_table$human_symbol=='MKI67',]


fib_data_norm[fib_data_norm$human_entrez=='5111',]
esc_data_norm




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


ggplot_theme <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14),
  #delete background
  # panel.grid.major = element_line(colour="grey", size = (1.5)),
  # panel.grid.minor = element_line(size = (0.2), colour="grey"),
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
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
  #,legend.position="none"
)
#--------------------------------

library(grid)
ggplot(plot_data, aes(x = rank, y = mean_cor_fib_esc,
                      col=group)) +
  geom_point(size=plot_data$size) +
  labs(x='Rank',y='Correlation with proliferation rate') +
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_hline(yintercept=seq(-1,1,0.5),linetype="dashed",colour="grey") +
  geom_vline(xintercept=c(0,5000,10000,15000),linetype="dashed",colour="grey") +
  scale_x_continuous(breaks= seq(0,20000,by=1000),
                     labels = c(0,rep("",4),
                                5000,rep("",4),
                                10000,rep("",4),
                                15000,rep("",4),
                                20000)
  ) +
  #annotation_ticks(sides = 'b', scale = 'identity',ticks_per_base=5) +
  ggplot_theme

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 gene rank PCNA KI67 .pdf', width = 140, height = 90, units = "mm")


# To validate the proliferation signature use 2019PNAS yeast data-------------------------
setwd('E:/坚果云同步private/literature-article-papers/cell_proliferation/2019 02 22/input data')

#----preprocess data----------------------------------------
load('./1.2009 Dataset S1 SuppFullArchive.RData')


frmeDataCharles <- as.data.frame(frmeDataCharles)
dim(frmeDataCharles)
dim(frmeDataGresham)
head(frmeDataCharles)

merge_1_2 <- merge(frmeDataCharles,frmeDataGresham,by="row.names",all.x=TRUE)
dim(merge_1_2) #5795   63



rownames(merge_1_2) <- merge_1_2$Row.names
merge_1_2 <- merge_1_2[,-1]

ugrr_data <- read.table("./3.ugrr_data.txt",sep = '\t',header = T)
dim(ugrr_data)
head(ugrr_data)

rownames(ugrr_data) <- ugrr_data$YORF
ugrr_data <- ugrr_data[,-1]

yeast_107_data <- merge(merge_1_2,ugrr_data,by="row.names",all.x=TRUE)

dim(yeast_107_data)

yeast_107_data <- na.omit(yeast_107_data)

saveRDS(yeast_107_data,"./yeast_107_data.rds")

#------load processed data----------------------------------------------------

setwd('E:/坚果云同步private/literature-article-papers/cell_proliferation/2019 02 22/input data')

yeast_107_data <- readRDS("./yeast_107_data.rds")
yeast_107_data[1:5,1:5]

rownames(yeast_107_data) <- yeast_107_data$Row.names
yeast_107_data <- yeast_107_data[,-1]




proliferation_signature <- readRDS('./proliferation_signature_combined_human_celegans_yeast.rds')


# use GSVA ssGSEA to calculate proliferation signature-----


# ssGSEA
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(GSEABase)


yeast_107_data[1:5,1:5]


gene_set <- list(as.character(na.omit(proliferation_signature$yeast_ensembl_id)))

gsva_matrix<- gsva(as.matrix(yeast_107_data),
                   gene_set, # gene set 用list格式
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)

# gsva_matrix<- gsva(as.matrix(yeast_107_data),
#                    gene_set, # gene set 用list格式
#                    method='gsva',
#                    kcdf='Gaussian')


plot(gsva_matrix[1,],growth_rate)




# create growth rate------------------------
colnames(frmeDataCharles)
colnames(frmeDataGresham)

colnames(ugrr_data)

frmeDataCharles_growth_rate <- substr(colnames(frmeDataCharles),11,14)

frmeDataCharles_growth_rate[frmeDataCharles_growth_rate=='0.1.'] <- '0.1'
frmeDataCharles_growth_rate[frmeDataCharles_growth_rate=='0.2H'] <- '0.2'


frmeDataGresham_growth_rate <- gsub(".*?\\..*?\\.",'',colnames(frmeDataGresham))

frmeDataGresham_growth_rate[2] <- '0.12'


ugrr_data_growth_rate <- gsub(".*?\\.\\..*?\\.",'',colnames(ugrr_data))


growth_rate <- as.numeric(c(frmeDataCharles_growth_rate,frmeDataGresham_growth_rate,ugrr_data_growth_rate))

saveRDS(growth_rate,'./growth_rate.rds')
# ----------20200309-----------------


growth_rate <- readRDS('./growth_rate.rds')

# ssGSEA--
proliferation_score <- as.data.frame(gsva_matrix[1,])
colnames(proliferation_score)[1] <- 'Proliferation_signature'
proliferation_score$growth_rate <- growth_rate


plot(proliferation_score)
abline(lm(growth_rate~proliferation_score,data = proliferation_score))
text(-30,0.28,paste0('r = ',round(cor(proliferation_score)[1,2],2)))



cor(proliferation_score[1:38,])

cor(proliferation_score[63:107,])

cor(proliferation_score[-c(1:38,63:107),])



cor.test(proliferation_score[-c(1:38,63:107),]$Proliferation_signature,proliferation_score[-c(1:38,63:107),]$growth_rate)
# p-value = 8.901e-07

cor.test(proliferation_score[63:107,]$Proliferation_signature,proliferation_score[63:107,]$growth_rate)
# p-value = 1.338e-08


plot(proliferation_score[63:107,],main='r=0.73')
plot(proliferation_score[1:38,], main='r=0.47')

plot(proliferation_score[-c(1:38,63:107),],main='r=0.82')


Data_Charles_cor <- cor(proliferation_score[1:38,])[2]
Data_Gresham_cor <- cor(proliferation_score[-c(1:38,63:107),])[2]
Data_ugrr_cor <- cor(proliferation_score[63:107,])[2]


#--------ggplot_theme-----------------------------
ggplot_theme <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
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
  #图例
  #legend.title = element_text(colour="black", size=14),
  legend.title = element_blank(),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=10)
  # remove legend
  #,legend.position="none"
)
#--------ggplot_theme_pdf-------------------------

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
  #,legend.position="none"
)



# FIGURE 4B Data_Gresham_cor  --------------------------
ggplot(proliferation_score[-c(1:38,63:107),],
       aes(x = Proliferation_signature, y = growth_rate)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='Proliferation signature score') +
  annotate("text", x=0.05, y=0.22, label= paste0("R = ",round(Data_Gresham_cor,2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification1_Data_Gresham_ssgesa.pdf', width = 60, height = 50, units = "mm")



# FIGURE 4C Data_ugrr_cor ------------
ggplot(proliferation_score[63:107,],
       aes(x = Proliferation_signature, y = growth_rate)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='Proliferation signature score') +
  annotate("text", x=-0.2, y=0.26, label= paste0("R = ",round(Data_ugrr_cor,2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification2_Data_ugrr_ssgsea.pdf', width = 60, height = 50, units = "mm")



# directly sum-----------------------------------

yeast_107_data_proliferation <- yeast_107_data[rownames(yeast_107_data) %in% proliferation_signature$yeast_ensembl_id,]

proliferation_score <- colSums(yeast_107_data_proliferation)

proliferation_score <- as.data.frame(proliferation_score)


proliferation_score$growth_rate <- growth_rate


plot(proliferation_score)
abline(lm(growth_rate~proliferation_score,data = proliferation_score))
text(-30,0.28,paste0('r = ',round(cor(proliferation_score)[1,2],2)))

colnames(proliferation_score)[1] <- 'Proliferation_signature'

cor.test(proliferation_score[-c(1:38,63:107),]$Proliferation_signature,proliferation_score[-c(1:38,63:107),]$growth_rate)

cor.test(proliferation_score[63:107,]$Proliferation_signature,proliferation_score[63:107,]$growth_rate)



plot(proliferation_score[63:107,],main='r=0.66')

plot(proliferation_score[-c(1:38,63:107),],main='r=0.78')


#FIGURE 4S A-----------------------------------

Data_Charles_cor <- cor(proliferation_score[1:38,])[2]
Data_Gresham_cor <- cor(proliferation_score[-c(1:38,63:107),])[2]
Data_ugrr_cor <- cor(proliferation_score[63:107,])[2]

cor.test(proliferation_score[1:38,]$Proliferation_signature,proliferation_score[1:38,]$growth_rate)




# Data_Gresham_cor
ggplot(proliferation_score[-c(1:38,63:107),],
       aes(x = Proliferation_signature, y = growth_rate)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='Proliferation signature score') +
  annotate("text", x=-10, y=0.22, label= paste0("R = ",round(Data_Gresham_cor,2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification1_Data_Gresham_sum.pdf', width = 60, height = 50, units = "mm")




# #FIGURE 4S B--Data_ugrr_cor---------

ggplot(proliferation_score[63:107,],
       aes(x = Proliferation_signature, y = growth_rate)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='Proliferation signature score') +
  annotate("text", x=-80, y=0.26, label= paste0("R = ",round(Data_ugrr_cor,2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification2_Data_ugrr_sum.pdf', width = 60, height = 50, units = "mm")


# #FIGURE 4S C-----2006 GB -------


setwd("E:/坚果云同步private/literature-article-papers/cell_proliferation/2006GenomeBiologydata/")

# rm(list=ls())

load('./E-MEXP-593.eSet.r')



phenoData <- study@phenoData@data

assayData <- study@assayData$exprs

colnames(phenoData)


colnames(assayData)


dim(assayData)

study@annotation


# read in processed data-------------------------------------

GB_processed_data <- read.table('./rescale.txt',sep = '\t',header = T)
dim(GB_processed_data)
str(GB_processed_data)


head(GB_processed_data)

GB_processed_data <- GB_processed_data[,c(1,seq(2,37,3))]

GB_processed_data <- GB_processed_data[-1,]

# read in annotation data------------------------

options(stringsAsFactors = F)
GB_annotation_data <- read.table('./A-AFFY-27.adf.txt',sep = '\t',header = T,skip = 16)

head(GB_annotation_data$Composite.Element.Database.Entry.sgd.,100)


GB_processed_data$yeast_sgd <- GB_annotation_data$Composite.Element.Database.Entry.sgd.[match(GB_processed_data$Scan.REF,GB_annotation_data$Composite.Element.Name)]



GB_processed_data <- GB_processed_data[nchar(GB_processed_data$yeast_sgd)>1,]
head(GB_processed_data)


rownames(GB_processed_data) 

# 
# GB_processed_data[GB_processed_data$yeast_sgd=='YJL116C',]

length(unique(GB_processed_data$yeast_sgd)) == length(GB_processed_data$yeast_sgd) # FALSE


library(data.table)

GB_processed_data_table <- GB_processed_data[,-1]

GB_processed_data_table[1:5,1:5]


GB_processed_data_table[,1:12] <- as.data.frame(lapply(GB_processed_data_table[,1:12],function(x) {gsub(',','.',x)}))

GB_processed_data_table[,1:12] <- apply(GB_processed_data_table[,1:12],2,as.numeric)
GB_processed_data_table <- as.data.table(GB_processed_data_table)

GB_processed_data_table <- GB_processed_data_table[,lapply(.SD, mean),by='yeast_sgd']

GB_processed_data_table <- as.data.frame(GB_processed_data_table)

colSums(GB_processed_data_table[,-1])


# read in yeast proliferation genes---------------------------

#genes <- readRDS('./YEAST_proliferation_genes.rds')

genes <- readRDS('./proliferation_signature_combined_human_celegans_yeast.rds')




proliferation_signature <- GB_processed_data_table[GB_processed_data_table$yeast_sgd %in% genes$yeast_ensembl_id,] 


rownames(proliferation_signature) <- proliferation_signature$yeast_sgd

proliferation_signature <- proliferation_signature[,-c(1)]



# directly sum up------------
growth_rate_data <- data.frame(proliferation_signature=colSums(proliferation_signature))

growth_rate_data$growth_rate <- c(0.02,0.05,rep(0.1,3),rep(0.2,3),0.25,rep(0.33,3))


# 20 out of 28 genes -----------------------------
cor(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature) # 0.9794513   # 0.5474304


cor.test(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature)


plot(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature)

ggplot(growth_rate_data,
       aes(y = growth_rate, x = proliferation_signature)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='growth_rate',x='Proliferation signature score') +
  annotate("text", y=0.3, x=270000, label= paste0('R = ',round(cor(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature),2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification1_Data_2006GB_sum.pdf', width = 60, height = 50, units = "mm")


# use ssgsea -------
GB_processed_data_table_ssgsea <- GB_processed_data_table
rownames(GB_processed_data_table_ssgsea) <- GB_processed_data_table_ssgsea$yeast_sgd
GB_processed_data_table_ssgsea <- GB_processed_data_table_ssgsea[,-1]

# ssGSEA
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(GSEABase)


gene_set <- list(as.character(na.omit(genes$yeast_ensembl_id)))

gsva_matrix<- gsva(as.matrix(GB_processed_data_table_ssgsea),
                   gene_set, # gene set 用list格式
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)


growth_rate_data <- as.data.frame(t(gsva_matrix))
colnames(growth_rate_data) <- 'proliferation_signature'

growth_rate_data$growth_rate <- c(0.02,0.05,rep(0.1,3),rep(0.2,3),0.25,rep(0.33,3))

cor(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature)   # 0.7661911
# #FIGURE 4D----------------------------------------------------------------------------------------------------------

library(ggplot2)
ggplot(growth_rate_data,
       aes(y = growth_rate, x = proliferation_signature)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='Proliferation signature score') +
  annotate("text", y=0.3, x=3.3, label= paste0('R = ',round(cor(growth_rate_data$growth_rate,growth_rate_data$proliferation_signature),2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification1_Data_2006GB_ssgesa.pdf', width = 60, height = 50, units = "mm")


# PCNA ----------------------
PCNA <- as.data.frame(t(GB_processed_data_table[GB_processed_data_table$yeast_sgd=='YBR088C',][-1]))
colnames(PCNA) <- 'PCNA'
PCNA$growth_rate <- c(0.02,0.05,rep(0.1,3),rep(0.2,3),0.25,rep(0.33,3))



cor.test(PCNA$PCNA,PCNA$growth_rate)


ggplot(PCNA,
       aes(y = growth_rate, x = PCNA)) +
  geom_point(size=1.5,shape=21,fill='red',col='red') +
  geom_smooth(method='lm', formula= y~x, se=0.95) +   # gam formula = y ~ s(x, bs = "cs")
  #geom_line(size=1) +
  #geom_errorbar(aes(ymin = overlap_cor_gene - se, ymax = overlap_cor_gene + se),width = 0.2) +
  labs(y='Growth rate',x='PCNA') +
  annotate("text", y=0.35, x=1200, label= paste0('R = ',round(cor(PCNA$PCNA,PCNA$growth_rate),2)),size=4,col='black') +
  #ylim(0,80) + 
  ggplot_theme_pdf

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 yeast_verification1_PCNA_2006GB.pdf', width = 60, height = 50, units = "mm")
















