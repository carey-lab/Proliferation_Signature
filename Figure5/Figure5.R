
#######################################################

# terminal and preterminal----------------------

setwd("E:/Project/2019__proliferation/Single_cell_data_analysis/C.elegant_single_cell/20190914_new_data/")

#-ggplot_theme_pdf-------
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


expression_data <- read.table("./celegant_id_singlecell_20201027_max_as_terminal.tsv",header = T,sep = "\t")


expression_data[1:5,1:4]

ncol(expression_data)

rownames(expression_data) <- expression_data$gene.short.id
expression_data <- expression_data[,-c(1:2)]
# proliferation_signature_combined----------------------------
proliferation_signature <- readRDS('./proliferation_signature_combined_human_celegans.rds')


# ssGSEA
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(GSEABase)


gene_set <- list(as.character(na.omit(proliferation_signature$c.elegant_id)))



gsva_matrix<- gsva(as.matrix(expression_data),
                   gene_set, # gene set 用list格式
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)



proliferation_signature <- as.data.frame(t(gsva_matrix))
colnames(proliferation_signature) <- "proliferation_signature"
proliferation_signature$condition <- "na"

proliferation_signature$condition[c(1:140)] <- "Terminal cell type"
proliferation_signature$condition[c(141:654)] <- "Preterminal cell lineage"

head(proliferation_signature)

library(cowplot)
library(RColorBrewer)
brewer.pal(11,"RdGy")
display.brewer.all(type="div")

Reds <- brewer.pal(n = 9, name = "Reds")[7]
Blues <- brewer.pal(n = 9, name = "Blues")[7]

proliferation_signature$condition <- as.factor(proliferation_signature$condition)


# ggplot is better than the ggboxplot when need to add point to the boxplot-----------
ggplot(data=proliferation_signature,aes(x=condition,y=proliferation_signature,color=condition)) + 
  geom_boxplot(size=0.5,outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=0.5) +
  scale_color_manual(values =c(Reds,Blues)) +
  labs(x='',y='Proliferation signature score') +
  ggplot_theme_pdf


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201027 2 C.elegans terminal boxplot max(terminal) ssGSEA.pdf', width = 50, height = 40, units = "mm")

# dev.off()


a=t.test(proliferation_signature~condition,data=proliferation_signature)
a$p.value

# 4.899673e-41

#####################################################
setwd('E:/Project/2019__proliferation/11.celegant_scRNAseq')




########################################################
# load data----
library(VisCello.celegans)
cello()


# ggplot theme-----

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
  legend.title = element_text(colour="black", size=14),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=10)
  # remove legend
  #,legend.position="none"
)



ggplot_theme_pdf <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7),
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
  axis.text.x = element_text(size=8,color='black',margin = margin(t=0.2, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=8,color='black',margin = margin(r=0.2, unit = "cm")),
  #图例
  legend.title = element_text(colour="black", size=7),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=7)
  # remove legend
  ,legend.position="none"
)

ggplot_theme_pdf_temp <- theme(
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
  axis.text.x = element_text(size=10,color='black',margin = margin(t=0.2, unit = "cm")
                             ,angle = 45
                             ,hjust = 1),
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
#----------------------------------------------------------------------------------
# get normalized expression data -------
assayData_eset <- assayData(eset)

assayData_eset_norm_exprs <- assayData_eset$norm_exprs


# check the raw expression data and normalized expression data -----
class(assayData_eset_norm_exprs)
sum(assayData_eset$exprs[,1]) # 1630
sum(assayData_eset$exprs[,2]) # 2319
sum(assayData_eset$exprs[,3]) # 2319


assayData_eset_norm_exprs
sum(assayData_eset_norm_exprs[,1])
sum(assayData_eset_norm_exprs[,2])
sum(assayData_eset_norm_exprs[,3])


# get metadata -----
metadata <- r_data$cmeta$df
head(metadata)
nrow(metadata) # 89701

metadata <- metadata[ metadata$to.filter=="FALSE",]

plot_value <- pvals$proj


# set up plot color ---------
library(ggplot2)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())

display.brewer.all()
brewer.pal()

set1 <- brewer.pal(5,'Set1')

color_set1 <- rgb(colorRamp(c(set1[2],set1[2],set1[3],set1[3],set1[5],set1[1]))( (0:12)/12 )/256)


plot(plot_value$V1,plot_value$V2,col=color_set1)




# get filtered data-----------------

dim(assayData_eset_norm_exprs)

assayData_eset_norm_exprs_filter <- assayData_eset_norm_exprs[,r_data$cmeta$df$to.filter=="FALSE"]
dim(assayData_eset_norm_exprs_filter)  # 20222 89701

all(colnames(assayData_eset_norm_exprs)==rownames(r_data$cmeta$df)) # true

dim(assayData_eset_norm_exprs_filter)

assayData_eset_norm_exprs_filter[1:5,1:5]



assayData_eset_norm_exprs_filter_depth <- colSums(assayData_eset_norm_exprs_filter)

plot_value$depth <- assayData_eset_norm_exprs_filter_depth


# human id to celegans id------------------

# require("biomaRt")
# mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# mart2 = useMart("ensembl", dataset="celegans_gene_ensembl")
# 
# # get human symbol for celegans transcript id
# # listAttributes
# human_entrez <- getLDS(attributes=c("external_gene_name"),
#                        filters="external_gene_name",
#                        values=rownames(assayData_eset_norm_exprs_filter), mart=mart2,
#                        attributesL=c("ensembl_gene_id","entrezgene_id"), martL=mart1                         ,uniqueRows = F)
# 
# head(human_entrez)
# nrow(human_entrez)
# 
# saveRDS(human_entrez,"./celegant_human_entrez.rds")


human_entrez <- readRDS("./celegant_human_entrez.rds")
# proliferation_signature <- readRDS('./proliferation_signature_combined_human.rds')

# proliferation_signature$c.elegant_id <- human_entrez$Gene.name[match(proliferation_signature$proliferation_signature,human_entrez$NCBI.gene.ID)]
# 
# 
# saveRDS(proliferation_signature,'./proliferation_signature_combined_human_celegans.rds')

proliferation_signature <- readRDS('./proliferation_signature_combined_human_celegans.rds')

# pick celegans genes that have corresponding human genes  ----------------------------------------------------------

class(assayData_eset_norm_exprs_filter)

assayData_eset_norm_exprs_filter_pick <- 
  human_entrez$NCBI.gene.ID[match(rownames(assayData_eset_norm_exprs_filter), human_entrez$Gene.name)]

assayData_eset_norm_exprs_filter_ma <- assayData_eset_norm_exprs_filter[! is.na(assayData_eset_norm_exprs_filter_pick),]

dim(assayData_eset_norm_exprs_filter_ma)

r_data$cmeta$df$embryo.time.bin

# table of clegans gene and human entrez----
assayData_names <- as.data.frame(rownames(assayData_eset_norm_exprs_filter_ma))

head(assayData_names)

assayData_names$entrez <- human_entrez$NCBI.gene.ID[match(assayData_names$`rownames(assayData_eset_norm_exprs_filter_ma)`,human_entrez$Gene.name)]

all(! is.na(assayData_names$entrez))

assayData_eset_norm_exprs_filter_ma[1:10,1:10]



####################################################################  this is main figure #############################
proliferation_signature <- readRDS('./proliferation_signature_combined_human_celegans_yeast.rds')

require("biomaRt")
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "http://asia.ensembl.org/")
listmart<-listDatasets(mart)
View(listAttributes(mart))

annotLookup <- getBM(
  mart=mart,
  attributes=c("entrezgene_id","ensembl_gene_id", "external_gene_name"),
  filter="entrezgene_id",
  values=proliferation_signature$proliferation_signature,
  uniqueRows=TRUE)

proliferation_signature$human_symbol <- annotLookup$external_gene_name[match(proliferation_signature$proliferation_signature,annotLookup$entrezgene_id)]

# search by hands in NCBI for NA genes


proliferation_signature[is.na(proliferation_signature$human_symbol),]$human_symbol <- c('PSMA6P4','RPLP0P6','PSMC1P4','LOC652826','LOC100132108')


saveRDS(proliferation_signature,'./proliferation_signature_combined_human_celegans_yeast2.rds')

write.table(proliferation_signature,'./proliferation_signature_combined_human_celegans_yeast2.tsv',sep = '\t')

# figure G -------------------------
# boxplot  Boxplot of C.elegant (left panel) --------------------

proliferation_signature_expr <- assayData_eset_norm_exprs_filter_ma[assayData_names$entrez %in% proliferation_signature$proliferation_signature,]

dim(proliferation_signature_expr)

plot_value$proliferation_signature_expr <- colSums(as.array(proliferation_signature_expr))



plot_value_1 <- plot_value[!plot_value$cell.type %in% c('Germline','Intestine','M_cell'),]



ggplot(data=plot_value_1,aes(x=embryo.time.bin,y=proliferation_signature_expr,color=embryo.time.bin)) + 
  geom_boxplot(size=0.5,outlier.size = 0.3
               #,outlier.shape = NA
  ) +
  #geom_jitter(shape=16, position=position_jitter(0.2),size=0.5) +
  #scale_color_manual(values =c(Reds,Blues)) +
  labs(x="Estimated embryo time (minutes)",y="Proliferation signature score") +
  ggplot_theme_pdf_temp


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans bin embryo time boxplot.pdf', width = 55, height = 45, units = "mm")


# boxplot  Boxplot of human (right panel) ------------------

# EMTAB3929 ##########################################################
proliferation_signature <- readRDS('./proliferation_signature_combined_human_celegans_yeast2.rds')


setwd('E:/坚果云同步private/literature-article-papers/single cell development/data_analysis/EMTAB3929')

# Genes were filtered, keeping 15,633 out of 26,178 genes that were expressed in at least 5 out of 1,919 sequenced cells (RPKM ≥ 10) and for which cells with expression came from at least two different embryos. 
rpkmHuman <- read.table("./E-MTAB-3929.processed/rpkm.txt",sep='\t',row.names = 1,header = T)

dim(rpkmHuman)

rpkm_filter <- function(x) {
  flag1 <- length(x[x>=10])
  flag2 <- length(unique(substr(colnames(rpkmHuman)[x>=10],1,2)))
  
  if (flag1>=5 & flag2>=2) return(TRUE) else return(FALSE)
}

index <- apply(rpkmHuman, 1, rpkm_filter)
table(index)

rpkmHuman <- rpkmHuman[index,]

proliferation_score <- colSums(log10(rpkmHuman[rownames(rpkmHuman) %in% proliferation_signature$human_symbol,]+1))

label <- names(proliferation_score)

get_label <- function(x) {
  if (grepl("early",x) | grepl("late",x)) {
    return(substr(x,1,8))
  } else {
    return(substr(x,1,2))
  }
}

label <- as.character(sapply(label,get_label))

proliferation_score <- as.data.frame(proliferation_score)

colnames(proliferation_score)[1] <- 'proliferation_score'

proliferation_score$embryo_time <- label

proliferation_score$embryo_time[proliferation_score$embryo_time=='E4.late.'] <- 'E4'#'E4 late'

proliferation_score$embryo_time[proliferation_score$embryo_time=='E5.early'] <- 'E5'#'E5 early'


# proliferation_score$embryo_time[proliferation_score$embryo_time=='E4.late.'] <- 'E4 late'
# 
# proliferation_score$embryo_time[proliferation_score$embryo_time=='E5.early'] <- 'E5 early'

unique(proliferation_score$embryo_time)

proliferation_score$embryo_time <- factor(proliferation_score$embryo_time,levels = c('E3','E4','E5','E6','E7')
                                          # levels = c('E3,'E4','E4 late','E5 early','E5','E6','E7')
)




ggplot(data=proliferation_score,aes(x=embryo_time,y=proliferation_score,color=embryo_time)) + 
  geom_boxplot(size=0.5,outlier.size = 0.3
               #,outlier.shape = NA
  ) +
  #geom_jitter(shape=16, position=position_jitter(0.2),size=0.5) +
  #scale_color_manual(values =c(Reds,Blues)) +
  labs(x="Embryonic days",y="Proliferation signature score") +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 human bin embryo time boxplot.pdf', width = 55, height = 45, units = "mm")





# 20200318 combine proliferation siganture and proliferation_signature_expr-----------
# load('./proliferation_signature.Rdata')
# proliferation_signature_1 <- proliferation_signature
# load('./proliferation_signature_overlap_gene1.Rdata')
# proliferation_signature_2 <- proliferation_signature
# 
# proliferation_signature <- rbind(proliferation_signature_1,proliferation_signature_2)
# 
# proliferation_signature <- proliferation_signature[order(as.numeric(proliferation_signature$proliferation_signature)),]
# 
# proliferation_signature <- proliferation_signature[!duplicated(proliferation_signature$proliferation_signature),]
# 
# saveRDS(proliferation_signature,'./proliferation_signature_combined.rds')


# check correaltion -----

cor(plot_value$embryo.time,plot_value$proliferation_signature_expr,method='spearman')  # -0.5333214

cor(plot_value$embryo.time,plot_value$proliferation_signature_expr,method='pearson') # -0.4703898


# check special cell types-----------
unique(plot_value$cell.type)

Germline_plot <- plot_value[plot_value$cell.type=='Germline',]

plot(Germline_plot$embryo.time,Germline_plot$proliferation_signature_expr,data=Germline_plot)

M_cell_plot <- plot_value[plot_value$cell.type=='M_cell',]

plot(M_cell_plot$embryo.time,M_cell_plot$proliferation_signature_expr,data=M_cell_plot)


Intestine_plot <- plot_value[plot_value$cell.type=='Intestine',]

plot(Intestine_plot$embryo.time,Intestine_plot$proliferation_signature_expr,data=Intestine_plot)


colnames(plot_value)
head(plot_value)





#  figure E Boxplot of proliferation signature for all cell types with embryo time > 650 ---------------------

all_bigger_than_650 <- plot_value[plot_value$embryo.time>650,][,c(12,19)]
all_bigger_than_650$cell.type <- 'All cell'



for (cell in unique(plot_value$cell.type)) {
  temp <- plot_value[plot_value$embryo.time>650 & plot_value$cell.type==cell,][,c(12,19,8)]
  all_bigger_than_650 <- rbind(all_bigger_than_650,temp)
}


library(data.table)
all_bigger_than_650 <- as.data.table(all_bigger_than_650)
all_bigger_than_650_order <- all_bigger_than_650[,lapply(.SD,median),by=cell.type]
all_bigger_than_650_order <- all_bigger_than_650_order$cell.type[order(all_bigger_than_650_order$proliferation_signature_expr,decreasing = T)]
all_bigger_than_650_order <- c(all_bigger_than_650_order[9],all_bigger_than_650_order[-9])

all_bigger_than_650 <- as.data.frame(all_bigger_than_650)
all_bigger_than_650$cell.type <- factor(all_bigger_than_650$cell.type,levels = all_bigger_than_650_order)


# library("ggpubr")
# ggboxplot(all_bigger_than_650, x = 'cell.type', y = "proliferation_signature_expr", 
#           #color = "Generation_number",
#           #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#           #order = c("ctrl", "trt1", "trt2"),
#           #ylab = "Weight", 
#           #title='Proliferation signature of all cell type with embryo time > 650',
#           xlab = "",
#           ylab = "Proliferation signature score"
# ) + rotate_x_text(angle = 45
#                   #, hjust = NULL, vjust = NULL
#                   )


ggplot_theme_pdf_temp2 <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7),
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
  axis.text.x = element_text(size=5,color='black',margin = margin(t=0.2, unit = "cm")
                             ,angle = 45
                             ,hjust = 1),
  #y轴标签
  axis.text.y = element_text(size=7,color='black',margin = margin(r=0.2, unit = "cm")),
  #图例
  legend.title = element_text(colour="black", size=10),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=10)
  # remove legend
  ,legend.position="none"
)


ggplot(data=all_bigger_than_650,aes(x=cell.type,y=proliferation_signature_expr)) + 
  geom_boxplot(size=0.3,outlier.size = 0.3
               #,outlier.shape = NA
  ) +
  #geom_jitter(shape=16, position=position_jitter(0.2),size=0.5) +
  #scale_color_manual(values =c(Reds,Blues)) +
  labs(y="Proliferation signature score",x='') +
  ggplot_theme_pdf_temp2

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans cell type boxplot.pdf', width = 120, height = 60, units = "mm")


boxplot(proliferation_signature_expr~cell.type,data=all_bigger_than_650,las=3)

table(plot_value[plot_value$cell.type=='M_cell',]$cell.subtype)

#  figure F --line plot for Germline, M_Cell and intestine-------------

# Germline-----------------------
library(beeswarm)
Germline_plot <- plot_value[ plot_value$cell.type=='Germline',][,c(12,13,19,8)]

beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Germline_plot, main='Germline',pch = 16, col = rainbow(12))
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Germline_plot, probs=0.5, add = TRUE)


boxplot(proliferation_signature_expr~embryo.time.bin,data=Germline_plot,main='Germline')

Germline_plot$embryo.time.bin <- as.character(Germline_plot$embryo.time.bin)

Germline_plot$embryo.time.bin[Germline_plot$embryo.time < 50] <- '< 50' 
Germline_plot$embryo.time.bin[Germline_plot$embryo.time >= 50 & Germline_plot$embryo.time < 75] <- '50-75' 
Germline_plot$embryo.time.bin[Germline_plot$embryo.time >= 75 & Germline_plot$embryo.time < 100] <- '75-100'
Germline_plot$embryo.time.bin[Germline_plot$embryo.time >= 100 & Germline_plot$embryo.time < 500] <- '100-500'
Germline_plot$embryo.time.bin[Germline_plot$embryo.time >= 500 ] <- '> 500'

Germline_plot$embryo.time.bin <- factor(Germline_plot$embryo.time.bin, levels = c('< 50','50-75','75-100','100-500','> 500'))

library(beeswarm)
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Germline_plot,
         # ylim=c(0,6),
         #cex = 2,
         log = F, pch = 16, col = rainbow(12),
         main = 'Germline',
         xlab="embryo time bin",
         cex.axis=1,
         cex.lab=1.5,
         cex.main=2)
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Germline_plot, probs=0.5, add = TRUE)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

library(Rmisc)
Germline_plot_se <- summarySE(Germline_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))




ggplot_theme_pdf_temp <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=8),
  axis.title.y = element_text(color="black", size=8),
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
  axis.text.x = element_text(size=7,color='black',margin = margin(t=0.2, unit = "cm")
                             ,angle = 45
                             ,hjust = 1),
  #y轴标签
  axis.text.y = element_text(size=7,color='black',margin = margin(r=0.2, unit = "cm")),
  #图例
  legend.title = element_text(colour="black", size=7),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=7)
  # remove legend
  ,legend.position="none"
)

ggplot(Germline_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans germline boxplot.pdf', width = 45, height = 45, units = "mm")

# M_Cell----

M_cell_plot <- plot_value[ plot_value$cell.type=='M_cell',][,c(12,13,19,8)]


boxplot(proliferation_signature_expr~embryo.time.bin,data=M_cell_plot,main='M_cell')
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = M_cell_plot, main='M_cell',pch = 16, col = rainbow(12))
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = M_cell_plot, probs=0.5, add = TRUE)


M_cell_plot$embryo.time.bin <- as.character(M_cell_plot$embryo.time.bin)

M_cell_plot$embryo.time.bin[M_cell_plot$embryo.time < 400] <- '< 400' 

M_cell_plot$embryo.time.bin[M_cell_plot$embryo.time >= 400 & M_cell_plot$embryo.time < 450] <- '400-450' 
M_cell_plot$embryo.time.bin[M_cell_plot$embryo.time >= 450 & M_cell_plot$embryo.time < 500] <- '450-500'
M_cell_plot$embryo.time.bin[M_cell_plot$embryo.time >= 500 & M_cell_plot$embryo.time < 700] <- '500-700'
M_cell_plot$embryo.time.bin[M_cell_plot$embryo.time >= 700 ] <- '> 700'

M_cell_plot$embryo.time.bin <- factor(M_cell_plot$embryo.time.bin, levels = c('< 400','400-450','450-500','500-700','> 700'))

library(beeswarm)
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = M_cell_plot,
         # ylim=c(0,6),
         #cex = 2,
         log = F, pch = 16, col = rainbow(12),
         main = 'M_cell',
         xlab="embryo time bin",
         cex.axis=1,
         cex.lab=1.5,
         cex.main=2)
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = M_cell_plot, probs=0.5, add = TRUE)


library(Rmisc)
M_cell_plot_se <- summarySE(M_cell_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))

ggplot(M_cell_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans M_cell boxplot.pdf', width = 45, height = 45, units = "mm")
# Intestine------

Intestine_plot <- plot_value[ plot_value$cell.type=='Intestine',][,c(12,13,19,8)]


boxplot(proliferation_signature_expr~embryo.time.bin,data=Intestine_plot,main='Intestine')
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Intestine_plot, main='Intestine',pch = 16, col = rainbow(12))
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Intestine_plot, probs=0.5, add = TRUE)


Intestine_plot$embryo.time.bin <- as.character(Intestine_plot$embryo.time.bin)

Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time < 200] <- '< 200' 

Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 200 & Intestine_plot$embryo.time < 300] <- '200-300' 
Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 300 & Intestine_plot$embryo.time < 400] <- '300-400'
Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 400 & Intestine_plot$embryo.time < 500] <- '400-500'
Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 500 & Intestine_plot$embryo.time < 600] <- '500-600'
Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 600 & Intestine_plot$embryo.time < 700] <- '600-700'
Intestine_plot$embryo.time.bin[Intestine_plot$embryo.time >= 700 ] <- '> 700'

Intestine_plot$embryo.time.bin <- factor(Intestine_plot$embryo.time.bin, levels = c('< 200','200-300','300-400','400-500','500-600','600-700','> 700'))

library(beeswarm)
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Intestine_plot,
         # ylim=c(0,6),
         #cex = 2,
         log = F, pch = 16, col = rainbow(12),
         main = 'Intestine',
         xlab="embryo time bin",
         cex.axis=1,
         cex.lab=1.5,
         cex.main=2)
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Intestine_plot, probs=0.5, add = TRUE)


library(Rmisc)
Intestine_plot_se <- summarySE(Intestine_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))

ggplot(Intestine_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans intestine boxplot.pdf', width = 45, height = 45, units = "mm")
# decrease example-----

# Body_wall_muscle-------------
sort(table(plot_value$cell.type))


Body_wall_muscle_plot <- plot_value[ plot_value$cell.type=='Body_wall_muscle',][,c(12,13,19,8)]


boxplot(proliferation_signature_expr~embryo.time.bin,data=Body_wall_muscle_plot,main='Body_wall_muscle')
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Body_wall_muscle_plot, main='Body_wall_muscle',pch = 16, col = rainbow(12))
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Body_wall_muscle_plot, probs=0.5, add = TRUE)


Body_wall_muscle_plot$embryo.time.bin <- as.character(Body_wall_muscle_plot$embryo.time.bin)

Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time < 200] <- '< 200' 

Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 200 & Body_wall_muscle_plot$embryo.time < 300] <- '200-300' 
Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 300 & Body_wall_muscle_plot$embryo.time < 400] <- '300-400'
Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 400 & Body_wall_muscle_plot$embryo.time < 500] <- '400-500'
Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 500 & Body_wall_muscle_plot$embryo.time < 600] <- '500-600'
Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 600 & Body_wall_muscle_plot$embryo.time < 700] <- '600-700'
Body_wall_muscle_plot$embryo.time.bin[Body_wall_muscle_plot$embryo.time >= 700 ] <- '> 700'

Body_wall_muscle_plot$embryo.time.bin <- factor(Body_wall_muscle_plot$embryo.time.bin, levels = c('< 200','200-300','300-400','400-500','500-600','600-700','> 700'))

library(beeswarm)
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Body_wall_muscle_plot,
         # ylim=c(0,6),
         #cex = 2,
         log = F, pch = 16, col = rainbow(12),
         main = 'Body_wall_muscle',
         xlab="embryo time bin",
         cex.axis=1,
         cex.lab=1.5,
         cex.main=2)
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Body_wall_muscle_plot, probs=0.5, add = TRUE)


library(Rmisc)
Body_wall_muscle_plot_se <- summarySE(Body_wall_muscle_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))

ggplot(Body_wall_muscle_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans body wall boxplot.pdf', width = 45, height = 45, units = "mm")


# Ciliated_amphid_neuron----------
Ciliated_amphid_neuron_plot <- plot_value[ plot_value$cell.type=='Ciliated_amphid_neuron',][,c(12,13,19,8)]

Ciliated_amphid_neuron_plot$embryo.time.bin <- as.character(Ciliated_amphid_neuron_plot$embryo.time.bin)

Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time < 200] <- '< 200' 

Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 200 & Ciliated_amphid_neuron_plot$embryo.time < 300] <- '200-300' 
Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 300 & Ciliated_amphid_neuron_plot$embryo.time < 400] <- '300-400'
Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 400 & Ciliated_amphid_neuron_plot$embryo.time < 500] <- '400-500'
Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 500 & Ciliated_amphid_neuron_plot$embryo.time < 600] <- '500-600'
Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 600 & Ciliated_amphid_neuron_plot$embryo.time < 700] <- '600-700'
Ciliated_amphid_neuron_plot$embryo.time.bin[Ciliated_amphid_neuron_plot$embryo.time >= 700 ] <- '> 700'

Ciliated_amphid_neuron_plot$embryo.time.bin <- factor(Ciliated_amphid_neuron_plot$embryo.time.bin, levels = c('< 200','200-300','300-400','400-500','500-600','600-700','> 700'))


Ciliated_amphid_neuron_plot_se <- summarySE(Ciliated_amphid_neuron_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))

ggplot(Ciliated_amphid_neuron_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans Ciliated_amphid_neuron boxplot.pdf', width = 45, height = 45, units = "mm")


# Hypodermis-------------
sort(table(plot_value$cell.type))


Hypodermis_plot <- plot_value[ plot_value$cell.type=='Hypodermis',][,c(12,13,19,8)]


boxplot(proliferation_signature_expr~embryo.time.bin,data=Hypodermis_plot,main='Hypodermis')
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Hypodermis_plot, main='Hypodermis',pch = 16, col = rainbow(12))
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Hypodermis_plot, probs=0.5, add = TRUE)


Hypodermis_plot$embryo.time.bin <- as.character(Hypodermis_plot$embryo.time.bin)

Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time < 200] <- '< 200' 

Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 200 & Hypodermis_plot$embryo.time < 300] <- '200-300' 
Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 300 & Hypodermis_plot$embryo.time < 400] <- '300-400'
Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 400 & Hypodermis_plot$embryo.time < 500] <- '400-500'
Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 500 & Hypodermis_plot$embryo.time < 600] <- '500-600'
Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 600 & Hypodermis_plot$embryo.time < 700] <- '600-700'
Hypodermis_plot$embryo.time.bin[Hypodermis_plot$embryo.time >= 700 ] <- '> 700'

Hypodermis_plot$embryo.time.bin <- factor(Hypodermis_plot$embryo.time.bin, levels = c('< 200','200-300','300-400','400-500','500-600','600-700','> 700'))

library(beeswarm)
beeswarm(proliferation_signature_expr ~ embryo.time.bin, data = Hypodermis_plot,
         # ylim=c(0,6),
         #cex = 2,
         log = F, pch = 16, col = rainbow(12),
         main = 'Hypodermis',
         xlab="embryo time bin",
         cex.axis=1,
         cex.lab=1.5,
         cex.main=2)
bxplot(proliferation_signature_expr ~ embryo.time.bin, data = Hypodermis_plot, probs=0.5, add = TRUE)


library(Rmisc)
Hypodermis_plot_se <- summarySE(Hypodermis_plot, measurevar="proliferation_signature_expr", groupvars=c("embryo.time.bin"))

ggplot(Hypodermis_plot_se, aes(x = embryo.time.bin, y = proliferation_signature_expr,group=1)) +
  geom_point(size=1) +
  geom_line(size=0.5) +
  geom_errorbar(aes(ymin = proliferation_signature_expr - se, ymax = proliferation_signature_expr + se),width = 0.3) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  
  ylim(0,100) +
  geom_hline(yintercept=50,linetype="dashed", 
             color = "red", size=0.5) +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans Hypodermis boxplot.pdf', width = 45, height = 45, units = "mm")







# figure C the UMAP PLOT -----proliferation_genes_exps-----------------------------------------------------
library(ggplot2)
head(plot_value)
hist(plot_value$proliferation_signature_expr)



ggplot_theme <- theme(
  plot.title = element_text(lineheight=18, size=20,hjust = 0.5),
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
  axis.ticks.length=unit(-0.2,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=8,color='black',margin = margin(t=0.3, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=8,color='black',margin = margin(r=0.3, unit = "cm")),
  #图例
  legend.title = element_text(colour="black", size=10),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=8)
  # remove legend
  #,legend.position="none"
)

ggplot(data = plot_value,aes(x=V1,y=V2)) +
  geom_point(aes(col=proliferation_signature_expr),size=0.01) +
  #scale_color_continuous(low = "#132B43", high = "#56B1F7",) +
  scale_colour_gradientn(name="Proliferation\nsignature",colours = c(rep(set1[1],3),
                                                                     rep(set1[5],2),
                                                                     rep(set1[3],5),
                                                                     rep(set1[2],5),rep(set1[4],6))) +
  #scale_color_manual(values = color_set1) +
  
  labs(x='UMAP1',y='UMAP2') +
  ggplot_theme

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans UMAP proliferation boxplot.pdf', width = 70, height = 60, units = "mm")



ggplot(data = plot_value,aes(x=V1,y=V2)) +
  geom_point(aes(col=embryo.time.bin),size=0.01) +
  #scale_color_continuous(low = "#132B43", high = "#56B1F7",) +
  # scale_colour_discrete(name="Estimated\nembryo time\n(minutes) ",
  #                        colours = c(rep(set1[1],2),
  #                                    rep(set1[5],1),
  #                                    rep(set1[3],4),
  #                                    rep(set1[2],4),rep(set1[4],6))) +
  
  scale_color_manual(name="Estimated\nembryo time\n(minutes)",
                     values = color_set1) +
  
  labs(x='UMAP1',y='UMAP2') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggplot_theme

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 C.elegans UMAP embryo time boxplot.pdf', width = 200, height = 160, units = "mm")




# FIGURE D Proliferation signature decrease with the embryo time. --------------------------
cor.test(plot_value_1$raw.embryo.time,plot_value_1$proliferation_signature_expr,method='spearman',exact=FALSE)

temp <- data.frame(embryo.time=plot_value$embryo.time,proliferation_signature_expr=plot_value$proliferation_signature_expr,type=plot_value$cell.type)

temp2 <- temp[!temp$type %in% c('Germline','Intestine','M_cell'),]

re=cor.test(temp2$embryo.time,temp2$proliferation_signature_expr,method='spearman',exact=F)
re$p.value  # -0.452308   p-value < 2.2e-16
cor.test(temp2$embryo.time,temp2$proliferation_signature_expr,method='pearson',exact=F)




library(data.table)
temp3 <- as.data.table(temp2[,1:2])

temp3 <- temp3[,lapply(.SD,median),by=embryo.time]

temp3 <- as.data.frame(temp3)
temp3 <- temp3[order(temp3$embryo.time),]

plot(temp3$embryo.time,temp3$proliferation_signature_expr)


re=cor.test(temp3$embryo.time,temp3$proliferation_signature_expr,method='spearman',exact=FALSE)
re$p.value  # -0.7342559  8.672783e-24



ggplot(data=temp3,aes(x=embryo.time,y=proliferation_signature_expr)) +
  geom_point(size=1) +
  geom_smooth(method='loess',se=T,size=0.5) +
  labs(x='Embryo time (minutes)',y='Proliferation signature score') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ggplot_theme_pdf


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201127 2 C.elegans smooth decrease median exclude.pdf', width = 65, height = 50, units = "mm")





# all cells -------
re=cor.test(temp3$embryo.time,temp3$proliferation_signature_expr,method='spearman',exact=FALSE)
re$p.value # -0.6498347  9.315872e-19


ggplot(data=temp3,aes(x=embryo.time,y=proliferation_signature_expr)) +
  geom_point(size=1) +
  geom_smooth(method='loess',se=T,size=0.5) +
  labs(x='Embryo time (minutes)',y='Proliferation signature score') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ggplot_theme_pdf


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201127 2 C.elegans smooth decrease median.pdf', width = 65, height = 50, units = "mm")

# after 650 -----

temp3 <- temp3[temp3$embryo.time>=650,]

re=cor.test(temp3$embryo.time,temp3$proliferation_signature_expr,method='spearman',exact=FALSE)
re$p.value # 0.5   p-value = 0.01512



ggplot(data=temp3,aes(x=embryo.time,y=proliferation_signature_expr)) +
  geom_point(size=1) +
  geom_smooth(method='loess',se=T,size=0.5) +
  labs(x='Embryo time (minutes)',y='Proliferation signature score') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  xlim(650,830) +
  ggplot_theme_pdf


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans smooth decrease after 650.pdf', width = 65, height = 50, units = "mm")


# FIGURE D change divide into terminal and preterminal---------------------------------

pre_terminal_cell_lineage <- read.table('./pre_terminal_cell_lineage.tsv',sep = '\t',quote="",header = T)

pre_terminal_cell_lineage[1:10,1:8]

table(pre_terminal_cell_lineage$UMAP.group)


terminal_cell_type <- read.table('./terminal_cell_type.tsv',sep = '\t',quote="",header = T)

terminal_cell_type[1:10,1:5]


head(plot_value)
head(plot_value[plot_value$lineage=='MSxpappp',])

table(plot_value$batch)
table(plot_value$cell.type)

table(pre_terminal_cell_lineage$IsTerminal)
pre_terminal_cell_lineage_yes <- pre_terminal_cell_lineage$Lineage[pre_terminal_cell_lineage$IsTerminal==0]

pre_terminal_cell_lineage_no <- pre_terminal_cell_lineage$Lineage[pre_terminal_cell_lineage$IsTerminal==1]


plot_value$cell_kinds <- 'Unannotated'

plot_value$cell_kinds[plot_value$lineage %in% pre_terminal_cell_lineage_yes] <- 'Preterminal'

plot_value$cell_kinds[plot_value$lineage %in% pre_terminal_cell_lineage_no] <- 'Terminal'


table(plot_value$cell_kinds)


#---

temp <- data.frame(embryo.time=plot_value$embryo.time,proliferation_signature_expr=plot_value$proliferation_signature_expr,type=plot_value$cell.type,kind=plot_value$cell_kinds)

temp2 <- temp[!temp$type %in% c('Germline','Intestine','M_cell'),]

re=cor.test(temp2$embryo.time,temp2$proliferation_signature_expr,method='spearman',exact=FALSE)
re$p.value  # -0.5789683 

# Ciliated_amphid_neuron----
Ciliated_amphid_neuron <- temp2[temp2$type=='Ciliated_amphid_neuron',]
plot(Ciliated_amphid_neuron$embryo.time,Ciliated_amphid_neuron$proliferation_signature_expr,col=as.factor(Ciliated_amphid_neuron$kind))

#Ciliated_amphid_neuron <- aggregate(proliferation_signature_expr ~ kind + embryo.time, Ciliated_amphid_neuron, mean)
ggplot_theme2 <- theme(
  plot.title = element_text(lineheight=.8, size=20,hjust = 0.5),
  axis.title.x = element_text(color="black", size=8),
  axis.title.y = element_text(color="black", size=8),
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
  axis.text.x = element_text(size=7,color='black',margin = margin(t=0.3, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=7,color='black',margin = margin(r=0.3, unit = "cm")),
  #图例
  legend.title = element_blank(),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank()
  legend.text = element_text(colour="black", size=7)
  # remove legend
  #,legend.position="none"
)


ggplot(data=Ciliated_amphid_neuron,aes(x=embryo.time,y=proliferation_signature_expr,group=kind,col=kind)) +
  geom_point(size=0.1) +
  #geom_smooth(method='loess',se=F,size=0.5) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ylim(0,100) +
  ggplot_theme2



ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans scatter Ciliated_amphid_neuron.pdf', width = 75, height = 45, units = "mm")


# Body_wall_muscle----
Body_wall_muscle <- temp2[temp2$type=='Body_wall_muscle',]
plot(Body_wall_muscle$embryo.time,Body_wall_muscle$proliferation_signature_expr,col=as.factor(Body_wall_muscle$kind))

ggplot(data=Body_wall_muscle,aes(x=embryo.time,y=proliferation_signature_expr,group=kind,col=kind)) +
  geom_point(size=0.1) +
  #geom_smooth(method='loess',se=F,size=0.5) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ylim(0,100) +
  ggplot_theme2



ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans scatter Body_wall_muscle.pdf', width = 75, height = 45, units = "mm")


# Hypodermis---------
Hypodermis <- temp2[temp2$type=='Hypodermis',]
plot(Hypodermis$embryo.time,Hypodermis$proliferation_signature_expr,col=as.factor(Hypodermis$kind))

#Hypodermis <- aggregate(proliferation_signature_expr ~ kind + embryo.time, Hypodermis, mean)


ggplot(data=Hypodermis,aes(x=embryo.time,y=proliferation_signature_expr,group=kind,col=kind)) +
  geom_point(size=0.1) +
  #geom_smooth(method='loess',se=F,size=0.5) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ylim(0,100) +
  ggplot_theme2



ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 2 C.elegans scatter Hypodermis.pdf', width = 75, height = 45, units = "mm")

# library(openxlsx)
# write.xlsx(temp2,'./temp2.xlsx')
# write.xlsx(plot_value,'./plot_value_with_unannotated.xlsx',row.names=T)
# library(data.table)

# all terminal and preterminal

temp <- data.frame(embryo.time=plot_value$embryo.time,proliferation_signature_expr=plot_value$proliferation_signature_expr,type=plot_value$cell.type,kind=plot_value$cell_kinds)

temp2 <- temp[!temp$type %in% c('Germline','Intestine','M_cell'),]
temp3 <- as.data.table(temp2[,c(1,2,4)])


temp3 %>% group_by(kind, embryo.time) %>%
  mutate(mean = mean(proliferation_signature_expr))


temp3 <- aggregate(proliferation_signature_expr ~ kind + embryo.time, temp3, mean)


temp3 <- temp3[order(temp3$embryo.time),]

plot(temp3$embryo.time,temp3$proliferation_signature_expr)
temp2 <- temp2[order(temp2$embryo.time)]

re=cor.test(temp3$embryo.time,temp3$proliferation_signature_expr,method='spearman',exact=FALSE)
re$p.value # 1.386724e-37

ggplot(data=temp3[temp3$kind=='Terminal',],aes(x=embryo.time,y=proliferation_signature_expr,group=kind,col=kind)) +
  geom_point(size=0.5) +
  geom_smooth(method='loess',se=F,size=0.5) +
  labs(x='Estimated embryo time (minutes)',y='Proliferation signature score') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ggplot_theme

library(plyr)
mu <- ddply(temp3, "kind", summarise, grp.mean=mean(proliferation_signature_expr))
ggplot(temp3, aes(x=proliferation_signature_expr, color=kind)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=kind),
             linetype="dashed") +
  ggplot_theme


ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/C.elegans smooth decrease.pdf', width = 60, height = 60, units = "mm")











# FIGURE S4A ------------------------------------
# figure, density plot to show most of the gene have 0 expression-----------------

set.seed(1000)
cell_1000 <- data.frame(
  cell_1000=as.vector(assayData_eset_norm_exprs_filter[,sample(1:89701, 1000, replace=FALSE)])
)

ggplot(cell_1000, aes(x=cell_1000)) + 
  geom_histogram(aes(y=..density..),colour="black", fill="white")+
  #geom_density(alpha=.2, fill="#FF6666")+
  labs(x="Counts(UMI)",y="Density") +
  ggplot_theme_pdf

table(cell_1000==0)
19326178/(895822+19326178)

# 0.9557006

ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201028 1000cellUMI density.pdf', width = 65, height = 65, units = "mm")
# FIGURE S4D and S4E ------------------------------------

setwd('E:/Project/2019__proliferation/11.celegant_scRNAseq/C.elegant cell number vs time')



# read the webplot data--------

webplotdata <- read.table('./webplot.csv',sep=',')


colnames(webplotdata) <- c('time','cellnumber')

webplotdata <- webplotdata[order(webplotdata$time),]

rownames(webplotdata) <- 1:nrow(webplotdata)


plot(webplotdata)

webplotdata <- webplotdata[c(10,15,24,37,47,54,57,61,63,67,72,75,78,81,84),]

webplotdata <- rbind(c(0,1),webplotdata)
webplotdata <- as.data.frame(webplotdata)
webplotdata$log2cellnumber <- log2(webplotdata$cellnumber)

#webplotdata <- webplotdata[c(1,seq(2,16,2)),]



plot(webplotdata$time,webplotdata$log2cellnumber, xlab='embryo time (min)',ylab='Log2 of cell number',type = "b", lty = 1,lwd=2,pch=19,col='purple',cex=1.5)

plot(diff(webplotdata$log2cellnumber)/diff(webplotdata$time),type = "b", lty = 1,lwd=2,pch=19,col='purple',cex=1.5,
     ylab = "", xaxt='n',xlab = "", yaxt='n')




library(ggplot2)


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




ggplot(data=webplotdata,aes(x=time,y=log2cellnumber,group=1)) +
  geom_line() +
  geom_point(size=1.5,col='purple') +
  labs(x='Time (minutes)',y='Log2 of cell number') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ggplot_theme_pdf

library(cowplot)
ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 C.elegans cell number line plot 1.pdf', width = 60, height = 50, units = "mm")



rollmean(x, 3)

plot(diff(data$cellnumber)/diff(data$etime),type = "b", lty = 1,lwd=2,pch=19,col='purple',cex=1.5,
     ylab = "", xaxt='n',xlab = "", yaxt='n')


library(zoo)
diff_webplotdata <- data.frame(matrix(nrow=15))
diff_webplotdata$diff_time <- rollmean(webplotdata$time,2)
diff_webplotdata$diff_log2cellnumber_change <- diff(webplotdata$log2cellnumber)/diff(webplotdata$time)



ggplot(data=diff_webplotdata,aes(x=diff_time,y=diff_log2cellnumber_change,group=1)) +
  geom_line() +
  geom_point(size=1.5,col='purple') +
  labs(x='Time (minutes)',y='Proliferation rate') +
  #annotate("text", x=400, y=50, label= "R = -0.85",size=6,col='black') +
  ggplot_theme_pdf

library(cowplot)
ggsave2('E:/Project/2019__proliferation/figures/C.elegans part figures/20201001 C.elegans cell number line plot2.pdf', width = 60, height = 50, units = "mm")





