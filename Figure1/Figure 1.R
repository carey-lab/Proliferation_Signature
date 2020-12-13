
# 20190605 CFSE doubling time with measured doubling time


# modify 20200820 ----------------

setwd('E:/Project/2019__proliferation/1.Flow_Cytometry')

library(flowCore)


#Load packages.

require(ggplot2)
require(ggridges)

#Load files and uniformize column names.

uM5 <- read.table("AllFCS5uM.tsv", header = TRUE, sep = "\t", row.names = 1)
colnames(uM5) <- c("FITC", "Generation", "CellType")
colnames(uM5)
uM5$log2FITC <- log2(uM5$FITC)
head(uM5)
uM5$CombinedID <- paste(uM5$CellType, uM5$Generation)
uM5$Generation <- factor(uM5$Generation)

uM5_line <- uM5

# use um10 to show the density plot--
#----um10-----------------------------
uM10 <- read.table("AllFCS10nM.tsv", header = TRUE, sep = "\t", row.names = 1)
colnames(uM10) <- c("FITC", "Generation", "CellType", "log2FITC", "CombinedID")

# minus blank
fib_unstained <- as.data.frame(exprs(read.FCS("./raw files/fib20171023_unstained.fcs", which.lines=NULL, transformation="linearize")))
#fib_unstained <- median(fib_unstained$`FITC-A`)

# minus blank
esc_unstained <- as.data.frame(exprs(read.FCS("./raw files/esc05122017_Unstained.fcs", which.lines=NULL, transformation="linearize")))
#esc_unstained <- median(esc_unstained$`FITC-A`)

fib_unstained$Generation <- 'AutoFlo'
esc_unstained$Generation <- 'AutoFlo'
fib_unstained$CellType <- "FIB"
esc_unstained$CellType <- "ESC"
fib_unstained$log2FITC <- log2(fib_unstained$`FITC-A`)
esc_unstained$log2FITC <- log2(esc_unstained$`FITC-A`)
fib_unstained$CombinedID <- paste(fib_unstained$CellType,fib_unstained$Generation)
esc_unstained$CombinedID <- paste(esc_unstained$CellType,esc_unstained$Generation)


esc_unstained <- esc_unstained[,c(4,33:36)]
fib_unstained <- fib_unstained[,c(4,33:36)]

colnames(esc_unstained)[1] <- "FITC"
colnames(fib_unstained)[1] <- "FITC"


esc_10_0 <- as.data.frame(exprs(read.FCS("./esc_CFSE_Generation0 10uM.fcs", which.lines=NULL, transformation="linearize")))

fib_10_0 <- as.data.frame(exprs(read.FCS("./fib_cfse_Generation0 10uM.fcs", which.lines=NULL, transformation="linearize")))


head(fib_10_0)
head(esc_10_0)
fib_10_0$Generation <- 0
esc_10_0$Generation <- 0
fib_10_0$CellType <- "FIB"
esc_10_0$CellType <- "ESC"
fib_10_0$log2FITC <- log2(fib_10_0$`FITC-A`)
esc_10_0$log2FITC <- log2(esc_10_0$`FITC-A`)
fib_10_0$CombinedID <- paste(fib_10_0$CellType,fib_10_0$Generation)
esc_10_0$CombinedID <- paste(esc_10_0$CellType,esc_10_0$Generation)


esc_10_0 <- esc_10_0[,c(5,12:15)]
fib_10_0 <- fib_10_0[,c(5,12:15)]

colnames(esc_10_0)[1] <- "FITC"
colnames(fib_10_0)[1] <- "FITC"

uM10 <- rbind(uM10,esc_10_0)
uM10 <- rbind(uM10,fib_10_0)

uM10 <- rbind(uM10,esc_unstained)
uM10 <- rbind(uM10,fib_unstained)

head(uM10)

table(uM10$Generation)


saveRDS(uM10,"./uM10_2.rds")
#-----------------------------------------------
uM10 <- readRDS('./uM10_2.rds')
uM10 <- uM10[order(uM10$Generation),]

uM10$Generation <- factor(uM10$Generation)

uM10_line <- uM10
#Removal of -Inf values from the log2 column.

uM5 <- uM5[is.finite(uM5$log2FITC), ]
uM10 <- uM10[is.finite(uM10$log2FITC), ]

#Plot generation.
library(RColorBrewer)

display.brewer.all()
display.brewer.pal(9,"Blues")

brewer.pal(9,"Greens")
brewer.pal(9,"Purples")
brewer.pal(9,"Purples")
brewer.pal(9,"Reds")
brewer.pal(9,"Set1")
brewer.pal(9,"Blues")

colorRampPalette(brewer.pal(9,"Blues"))(100)


uM5$CellType <- factor(uM5$CellType, levels = rev(levels(uM5$CellType)))
uM10$CellType <- factor(uM10$CellType, levels = rev(levels(uM10$CellType)))

uM10$CellType <- factor(uM10$CellType,levels = c("ESC","FIB"))


# figure 1b--------------

plot10uM <- ggplot(uM10, aes(x = uM10$log2FITC, y = uM10$Generation, fill = uM10$CellType)) + geom_density_ridges(alpha = 0.7) + scale_fill_manual(values = c("#2171B5", "#7C285C"), name = "Cell type") + ylab("Number of doublings") + xlab("log2(FITC)") + 
  ggtitle("Concentration: 10 nM") + 
  theme_classic() + 
  theme(axis.line = element_blank(), panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(color="black", size=18),
        axis.title.y = element_text(color="black", size=18),
        axis.text = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
        legend.position = "none") 
plot10uM

ggsave(plot10uM, filename = "Distribution_10uM5 20200820.pdf", width = 5, height = 4.7)



# plot the line plot for 5uM and 10uM---------------------------------------
library(data.table)
head(uM5_line)
head(uM10_line)


fib_5 <- uM5_line[,c(1,2,4)][uM5_line$CellType=="FIB",]
fib_5 <- as.data.table(fib_5)
fib_5 <- fib_5[,lapply(.SD, median),by=Generation]
# minus blank
fib_unstained <- as.data.frame(exprs(read.FCS("./raw files/fib20171023_unstained.fcs", which.lines=NULL, transformation="linearize")))
fib_unstained <- median(fib_unstained$`FITC-A`)
fib_5$FITC <- fib_5$FITC - fib_unstained
fib_5$log2FITC <- log2(fib_5$FITC)



esc_5 <- uM5_line[,c(1,2,4)][uM5_line$CellType=="ESC",]
esc_5 <- as.data.table(esc_5)
esc_5 <- esc_5[,lapply(.SD, median),by=Generation]
# minus blank
esc_unstained <- as.data.frame(exprs(read.FCS("./raw files/esc05122017_Unstained.fcs", which.lines=NULL, transformation="linearize")))
esc_unstained <- median(esc_unstained$`FITC-A`)
esc_5$FITC <- esc_5$FITC - esc_unstained
esc_5$log2FITC <- log2(esc_5$FITC)





fib_10 <- uM10_line[,c(1,2,4)][uM10_line$CellType=="FIB",]
fib_10 <- as.data.table(fib_10)
fib_10 <- fib_10[,lapply(.SD, median),by=Generation]
# minus blank

fib_10$FITC <- fib_10$FITC - fib_unstained
fib_10$log2FITC <- log2(fib_10$FITC)



esc_10 <- uM10_line[,c(1,2,4)][uM10_line$CellType=="ESC",]
esc_10 <- as.data.table(esc_10)
esc_10 <- esc_10[,lapply(.SD, median),by=Generation]
# minus blank

esc_10$FITC <- esc_10$FITC - esc_unstained
esc_10$log2FITC <- log2(esc_10$FITC)


fib <- merge(fib_5,fib_10,by="Generation")

esc <- merge(esc_5,esc_10,by="Generation")

fib$time <- seq(0,7*24,24)
esc$time <- seq(0,6*24,24)

colnames(fib)[c(3,5)] <- c("log2FITC_5uM","log2FITC_10uM")


fib$Generation <- as.numeric(fib$Generation)
head(fib)
fib <- fib[,-c(1,2,4)]
colnames(fib) <- c("5uM CFSE","10uM CFSE","time")

library(reshape2)
fib_melt <- melt(fib,id.vars="time")
head(fib_melt)
fib_melt$variable <- factor(fib_melt$variable)


fib_dt <- fib[2:7,]
fib_dt <- lm(`5uM CFSE`~time,data=fib_dt)
summary(fib_dt)
fib_dt <- 1/fib_dt$coefficients[2]  # 18.52635 

fib_dt <- fib[2:7,]
fib_dt <- lm(`10uM CFSE`~time,data=fib_dt)
summary(fib_dt)
fib_dt <- 1/fib_dt$coefficients[2]  # 18.22631 

##############################################
# figure 1c--------------

library(cowplot)
library(ggpubr)

fib_line <- ggplot(data = fib_melt,aes(x=time,y=value)) +
  geom_point(aes(x=time,y=value,shape=variable),fill="#7C285C",col="#7C285C",size=5, alpha = 0.8) + 
  geom_smooth(data=subset(fib_melt, time %in% 24:165),aes(group=variable),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid", colour = c("#7C285C")) +
  
  # scale_colour_manual(values = c("#7C285C", "#377EB8"),
  #                     #aesthetics = c("colour", "fill"),
  #                     labels = c("5μM CFSE (DT=21hrs)", "10μM CFSE (DT=19hrs)")) + # legned labels
  
  scale_shape_manual(values = c(21, 24),
                     #aesthetics = c("colour", "fill"),
                     labels = c("5μM CFSE (DT=19hrs)", "10μM CFSE (DT=18hrs)")) +
  
  scale_x_continuous(breaks = seq(0, 168, by = 24)) +
  
  scale_y_continuous(breaks = seq(4, 16, by = 2)) +
  
  labs(x="Time (hours after initial sort) ",y="log2(median CFSE)") +
  
  theme_bw(base_rect_size=0.5) +
  
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
  axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
  #图例
  #legend.title = element_text(colour="black", size=14),
  legend.title = element_blank(),
  #legend.text = element_blank()
  legend.text = element_text(colour=c("black"), size=16)
  # remove legend
  ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
  legend.background = element_rect(fill="white",
                                   size=0.5, linetype="solid", 
                                   colour ="black")
  #,legend.position="none"
) 
fib_line
ggsave(fib_line, filename = "fib_doubling_time8.pdf", width = 5, height = 4.7)




esc$Generation <- as.numeric(esc$Generation)
head(esc)
esc <- esc[,-c(1,2,4)]
colnames(esc) <- c("5uM CFSE","10uM CFSE","time")

library(reshape2)
esc_melt <- melt(esc,id.vars="time")
head(esc_melt)
esc_melt$variable <- factor(esc_melt$variable)


esc_dt <- esc[2:6,]
esc_dt <- lm(`5uM CFSE`~time,data=esc_dt)
summary(esc_dt)
esc_dt <- 1/esc_dt$coefficients[2]  # 10.01226  

esc_dt <- esc[2:6,]
esc_dt <- lm(`10uM CFSE`~time,data=esc_dt)
summary(esc_dt)
esc_dt <- 1/esc_dt$coefficients[2]  # 10.37069 


#################################################
# figure 1d--------------
esc_line <- ggplot(data = esc_melt,aes(x=time,y=value)) +
  geom_point(aes(x=time,y=value,shape=variable),fill="#2171B5",col="#2171B5",size=5, alpha = 0.8) + 
  geom_smooth(data=subset(esc_melt, time %in% 24:120),aes(group=variable),method='lm',formula=y~x,se=F, fullrange=FALSE, level=0.95, linetype = "solid", colour = "#2171B5") +
  
  # scale_colour_manual(values = c("#7C285C", "#377EB8"),
  #                     aesthetics = c("colour", "fill"),
  #                     labels = c("5μM CFSE (DT=12hrs)", "10μM CFSE (DT=10hrs)")) + # legned labels
  
  
  scale_shape_manual(values = c(21, 24),
                     #aesthetics = c("colour", "fill"),
                     labels = c("5μM CFSE (DT=10hrs)", "10μM CFSE (DT=10hrs)")) +
  
  
  scale_x_continuous(breaks = seq(0, 168, by = 24)) +
  
  scale_y_continuous(breaks = seq(4, 16, by = 2)) +
  
  labs(x="Time",y="log2( median CFSE )") +
  
  theme_bw(base_rect_size=0.5) +
  
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
    axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
    #图例
    #legend.title = element_text(colour="black", size=14),
    legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour=c("black"), size=16)
    # remove legend
    ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
    legend.background = element_rect(fill="white",
                                     size=0.5, linetype="solid", 
                                     colour ="black")
    #,legend.position="none"
  )

esc_line

ggsave(esc_line, filename = "esc_doubling_time8.pdf", width = 5, height = 4.7)

##################################################
# figure 1g--------------


# re check the genes that genes have spearman corelation 1 
# 20200820----------

require(ggplot2)
require(reshape2)
require(gtools)
library(cowplot)

FIB2016 <- read.csv("FIBBehavGroups.tsv", sep = "\t", header = TRUE)
FIB2016$X <- substr(FIB2016$X, 1, 18)
ESC2016 <- read.csv("ESCBehavGroups.tsv", sep = "\t", header = TRUE)
ESC2016$X <- substr(ESC2016$X, 1, 18)



# create a table to check genes------------

FIB2016_check <- FIB2016
FIB2016_check$spearman_1 <- t(cor(c(1,2,3),t(FIB2016_check[,c(4,7,10)]),method = 'spearman'))

FIB2016_check$spearman_2 <- t(cor(c(1,2,3),t(FIB2016_check[,c(13,16,19)]),method = 'spearman'))


FIB2016_check_res <- FIB2016_check[FIB2016_check$spearman_1==1 & FIB2016_check$spearman_2==1,]

FIB2016_check_res <-  na.omit(FIB2016_check_res)

dim(FIB2016_check_res) # 1616   31

colnames(FIB2016_check)




ESC2016_check <- ESC2016
ESC2016_check$spearman_1 <- t(cor(c(1,2,3),t(ESC2016_check[,c(4,7,10)]),method = 'spearman'))

ESC2016_check$spearman_2 <- t(cor(c(1,2,3),t(ESC2016_check[,c(13,16,19)]),method = 'spearman'))


ESC2016_check_res <- ESC2016_check[ESC2016_check$spearman_1==1 & ESC2016_check$spearman_2==1,]

ESC2016_check_res <-  na.omit(ESC2016_check_res)

dim(ESC2016_check_res) # 1616   31

colnames(ESC2016_check)



FIB2016_check_te <- FIB2016_check_res[,c(1,4,7,10,13,16,19)]
ESC2016_check_te <- ESC2016_check_res[,c(1,4,7,10,13,16,19)]

colnames(ESC2016_check_te) <- paste0('ESC',colnames(ESC2016_check_te))
colnames(ESC2016_check_te)[1] <- 'X'

FIB_ESC_2016 <- merge(FIB2016_check_te,ESC2016_check_te,by = 'X')

dim((FIB_ESC_2016))

#-----------------------------------
require("biomaRt")
listMarts()
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# , host = "http://asia.ensembl.org/"

listmart<-listDatasets(mart)
View(listmart)


annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "description", "external_gene_name"),
  filter="ensembl_gene_id",
  values=FIB_ESC_2016$X,
  uniqueRows=TRUE)

#--------------------------------------

head(FIB_ESC_2016)

FIB_ESC_2016$gene_name <- annotLookup$external_gene_name[match(FIB_ESC_2016$X,annotLookup$ensembl_gene_id)]

FIB_ESC_2016$description <- annotLookup$description[match(FIB_ESC_2016$X,annotLookup$ensembl_gene_id)]


ribosome_biogenesis <- read.table('./ribi_geneset.txt',header = T)

head(ribosome_biogenesis)


FIB_ESC_2016$gene_name <- toupper(FIB_ESC_2016$gene_name)


FIB_ESC_2016_save <- FIB_ESC_2016[FIB_ESC_2016$gene_name %in% ribosome_biogenesis$GO_RIBOSOME_BIOGENESIS,]



line_fit <- function(data) {
  x = c(1,2,3)
  res <- lm(data~x)
  res <- res[[1]][2]
  return(res)
}



FIB_ESC_2016_save$slope_fib_1 <- apply(FIB_ESC_2016_save[,c(2:4)],1,line_fit)

FIB_ESC_2016_save$slope_fib_2 <- apply(FIB_ESC_2016_save[,c(5:7)],1,line_fit)


FIB_ESC_2016_save$slope_esc_1 <- apply(FIB_ESC_2016_save[,c(8:10)],1,line_fit)

FIB_ESC_2016_save$slope_esc_2 <- apply(FIB_ESC_2016_save[,c(11:13)],1,line_fit)




head(FIB_ESC_2016_save)

write.table(FIB_ESC_2016_save,'./FIB_ESC_2016_save.csv',sep = ",",col.names = T)



#---------------------------------------------------



Columns2Keep2016S1 <- c("S1_TPM", "S1_TPM_95CI_low", "S1_TPM_95CI_high") #Slow replicate 1
Columns2Keep2016M1 <- c("M1_TPM", "M1_TPM_95CI_low", "M1_TPM_95CI_high") #Med replicate 1
Columns2Keep2016F1 <- c("F1_TPM", "F1_TPM_95CI_low", "F1_TPM_95CI_high") #Fast replicate 1
Columns2Keep2016S2 <- c("S2_TPM", "S2_TPM_95CI_low", "S2_TPM_95CI_high") #Slow replicate 2
Columns2Keep2016M2 <- c("M2_TPM", "M2_TPM_95CI_low", "M2_TPM_95CI_high") #Med replicate 2
Columns2Keep2016F2 <- c("F2_TPM", "F2_TPM_95CI_low", "F2_TPM_95CI_high") #Fast replicate 2

NewColNames <- c("median", "lowerCI", "higherCI")



#RRS1------------------------------------

RRS1FIB16S1 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016S1]
RRS1FIB16M1 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016M1]
RRS1FIB16F1 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016F1]
RRS1FIB16S2 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016S2]
RRS1FIB16M2 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016M2]
RRS1FIB16F2 <- FIB2016[FIB2016$X == "ENSMUSG00000061024", Columns2Keep2016F2]

colnames(RRS1FIB16S1) <- NewColNames
colnames(RRS1FIB16M1) <- NewColNames
colnames(RRS1FIB16F1) <- NewColNames
colnames(RRS1FIB16S2) <- NewColNames
colnames(RRS1FIB16M2) <- NewColNames
colnames(RRS1FIB16F2) <- NewColNames

RRS1ESC16S1 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016S1]
RRS1ESC16M1 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016M1]
RRS1ESC16F1 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016F1]
RRS1ESC16S2 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016S2]
RRS1ESC16M2 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016M2]
RRS1ESC16F2 <- ESC2016[ESC2016$X == "ENSMUSG00000061024", Columns2Keep2016F2]
colnames(RRS1ESC16S1) <- NewColNames
colnames(RRS1ESC16M1) <- NewColNames
colnames(RRS1ESC16F1) <- NewColNames
colnames(RRS1ESC16S2) <- NewColNames
colnames(RRS1ESC16M2) <- NewColNames
colnames(RRS1ESC16F2) <- NewColNames

#Build a dataframe to plot. 
RRS1 <- smartbind(RRS1FIB16S1, RRS1FIB16M1, RRS1FIB16F1, RRS1FIB16S2, RRS1FIB16M2, RRS1FIB16F2, RRS1ESC16S1, RRS1ESC16M1, RRS1ESC16F1, RRS1ESC16S2, RRS1ESC16M2, RRS1ESC16F2)

ProlifIndex <- rep(c("Slow", "Med", "Fast"), 4)
CellType <- rep(c(rep("FIB", 6), rep("ESC", 6)))
Replicate <- rep(c(rep(1, 3), rep(2, 3)), 2)

RRS1$ProlifIndex <- ProlifIndex
RRS1$CellType <- CellType
RRS1$Replicate <- Replicate 
RRS1$all <- paste(CellType, Replicate)



RRS1Plot <- ggplot(RRS1, aes(x = factor(RRS1$ProlifIndex, levels = rev(levels(factor(RRS1$ProlifIndex)))), y = log2(RRS1$median), group = all, color = CellType)) + 
  geom_line(size=1.5) +  #aes(linetype = Replicate),
  geom_point() + 
  geom_errorbar(width = 0.1, aes(ymin = log2(RRS1$lowerCI), ymax = log2(RRS1$higherCI)),size=1.5) + 
  scale_color_manual(values = c("#7C285C", "#2171B5")) + 
  ylab("log2(TPM)") + xlab("") +  ggtitle("RRS1") +
  theme_bw(base_rect_size=0.5) +
  
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
    axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
    #图例
    #legend.title = element_text(colour="black", size=14),
    #legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour=c("black"), size=16)
    # remove legend
    # ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
    # legend.background = element_rect(fill="white",
    #                                  size=0.5, linetype="solid", 
    #                                  colour ="black")
    #,legend.position="none"
  )

RRS1Plot
ggsave(RRS1Plot, filename = "RRS1_2016_20200821.pdf", width = 5, height = 5)


#*Tnfrsf21*-------------------------------

Tnfrsf21FIB16S1 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016S1]
Tnfrsf21FIB16M1 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016M1]
Tnfrsf21FIB16F1 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016F1]
Tnfrsf21FIB16S2 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016S2]
Tnfrsf21FIB16M2 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016M2]
Tnfrsf21FIB16F2 <- FIB2016[FIB2016$X == "ENSMUSG00000023915", Columns2Keep2016F2]

colnames(Tnfrsf21FIB16S1) <- NewColNames
colnames(Tnfrsf21FIB16M1) <- NewColNames
colnames(Tnfrsf21FIB16F1) <- NewColNames
colnames(Tnfrsf21FIB16S2) <- NewColNames
colnames(Tnfrsf21FIB16M2) <- NewColNames
colnames(Tnfrsf21FIB16F2) <- NewColNames

Tnfrsf21ESC16S1 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016S1]
Tnfrsf21ESC16M1 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016M1]
Tnfrsf21ESC16F1 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016F1]
Tnfrsf21ESC16S2 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016S2]
Tnfrsf21ESC16M2 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016M2]
Tnfrsf21ESC16F2 <- ESC2016[ESC2016$X == "ENSMUSG00000023915", Columns2Keep2016F2]
colnames(Tnfrsf21ESC16S1) <- NewColNames
colnames(Tnfrsf21ESC16M1) <- NewColNames
colnames(Tnfrsf21ESC16F1) <- NewColNames
colnames(Tnfrsf21ESC16S2) <- NewColNames
colnames(Tnfrsf21ESC16M2) <- NewColNames
colnames(Tnfrsf21ESC16F2) <- NewColNames

#Build a dataframe to plot. 

Tnfrsf21 <- smartbind(Tnfrsf21FIB16S1, Tnfrsf21FIB16M1, Tnfrsf21FIB16F1, Tnfrsf21FIB16S2, Tnfrsf21FIB16M2, Tnfrsf21FIB16F2, Tnfrsf21ESC16S1, Tnfrsf21ESC16M1, Tnfrsf21ESC16F1, Tnfrsf21ESC16S2, Tnfrsf21ESC16M2, Tnfrsf21ESC16F2)

ProlifIndex <- rep(c("Slow", "Med", "Fast"), 4)
CellType <- rep(c(rep("FIB", 6), rep("ESC", 6)))
Replicate <- rep(c(rep(1, 3), rep(2, 3)), 2)

Tnfrsf21$ProlifIndex <- ProlifIndex
Tnfrsf21$CellType <- CellType
Tnfrsf21$Replicate <- Replicate 
Tnfrsf21$all <- paste(CellType, Replicate)



Tnfrsf21$CellType <- factor(Tnfrsf21$CellType,levels = c("FIB","ESC"))
Tnfrsf21$Replicate <- factor(Tnfrsf21$Replicate)


Tnfrsf21Plot <- ggplot(Tnfrsf21, aes(x = factor(Tnfrsf21$ProlifIndex, levels = rev(levels(factor(Tnfrsf21$ProlifIndex)))), y = log2(Tnfrsf21$median), group = all, color = CellType)) + 
  geom_line(size=1.5) +  #aes(linetype = Replicate),
  geom_point() + 
  geom_errorbar(width = 0.1, aes(ymin = log2(Tnfrsf21$lowerCI), ymax = log2(Tnfrsf21$higherCI)),size=1.5) + 
  scale_color_manual(values = c("#7C285C", "#2171B5")) + 
  ylab("log2(TPM)") + xlab("") +  ggtitle("Tnfrsf21") +
  theme_bw(base_rect_size=0.5) +
  
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
    axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
    #图例
    #legend.title = element_text(colour="black", size=14),
    #legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour=c("black"), size=16)
    # remove legend
    # ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
    # legend.background = element_rect(fill="white",
    #                                  size=0.5, linetype="solid", 
    #                                  colour ="black")
    #,legend.position="none"
  )

Tnfrsf21Plot
ggsave(Tnfrsf21Plot, filename = "Tnfrsf21_2016_20190630.pdf", width = 5, height = 5)




#*Pmaip*----------------------------------------

PmaipFIB16S1 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016S1]
PmaipFIB16M1 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016M1]
PmaipFIB16F1 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016F1]
PmaipFIB16S2 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016S2]
PmaipFIB16M2 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016M2]
PmaipFIB16F2 <- FIB2016[FIB2016$X == "ENSMUSG00000024521", Columns2Keep2016F2]

colnames(PmaipFIB16S1) <- NewColNames
colnames(PmaipFIB16M1) <- NewColNames
colnames(PmaipFIB16F1) <- NewColNames
colnames(PmaipFIB16S2) <- NewColNames
colnames(PmaipFIB16M2) <- NewColNames
colnames(PmaipFIB16F2) <- NewColNames

PmaipESC16S1 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016S1]
PmaipESC16M1 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016M1]
PmaipESC16F1 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016F1]
PmaipESC16S2 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016S2]
PmaipESC16M2 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016M2]
PmaipESC16F2 <- ESC2016[ESC2016$X == "ENSMUSG00000024521", Columns2Keep2016F2]
colnames(PmaipESC16S1) <- NewColNames
colnames(PmaipESC16M1) <- NewColNames
colnames(PmaipESC16F1) <- NewColNames
colnames(PmaipESC16S2) <- NewColNames
colnames(PmaipESC16M2) <- NewColNames
colnames(PmaipESC16F2) <- NewColNames

#Build a dataframe to plot. 

Pmaip <- smartbind(PmaipFIB16S1, PmaipFIB16M1, PmaipFIB16F1, PmaipFIB16S2, PmaipFIB16M2, PmaipFIB16F2, PmaipESC16S1, PmaipESC16M1, PmaipESC16F1, PmaipESC16S2, PmaipESC16M2, PmaipESC16F2)

ProlifIndex <- rep(c("Slow", "Med", "Fast"), 4)
CellType <- rep(c(rep("FIB", 6), rep("ESC", 6)))
Replicate <- rep(c(rep(1, 3), rep(2, 3)), 2)

Pmaip$ProlifIndex <- ProlifIndex
Pmaip$CellType <- CellType
Pmaip$Replicate <- Replicate 
Pmaip$all <- paste(CellType, Replicate)



Pmaip$CellType <- factor(Pmaip$CellType,levels = c("FIB","ESC"))
Pmaip$Replicate <- factor(Pmaip$Replicate)


PmaipPlot <- ggplot(Pmaip, aes(x = factor(Pmaip$ProlifIndex, levels = rev(levels(factor(Pmaip$ProlifIndex)))), y = log2(Pmaip$median), group = all, color = CellType)) + 
  geom_line(size=1.5) +  #aes(linetype = Replicate),
  geom_point() + 
  geom_errorbar(width = 0.1, aes(ymin = log2(Pmaip$lowerCI), ymax = log2(Pmaip$higherCI)),size=1.5) + 
  scale_color_manual(values = c("#7C285C", "#2171B5")) + 
  ylab("log2(TPM)") + xlab("") +  ggtitle("Pmaip") +
  theme_bw(base_rect_size=0.5) +
  
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
    axis.text.x = element_text(size=18,color='black',margin = margin(t=0.3, unit = "cm")),
    #y轴标签
    axis.text.y = element_text(size=18,color='black',margin = margin(r=0.3, unit = "cm")),
    #图例
    #legend.title = element_text(colour="black", size=14),
    #legend.title = element_blank(),
    #legend.text = element_blank()
    legend.text = element_text(colour=c("black"), size=16)
    # remove legend
    # ,legend.justification=c(1,0), legend.position=c(0.95,0.8),
    # legend.background = element_rect(fill="white",
    #                                  size=0.5, linetype="solid", 
    #                                  colour ="black")
    #,legend.position="none"
  )

PmaipPlot
ggsave(PmaipPlot, filename = "Pmaip_2016_20190630.pdf", width = 5, height = 5)


















# figure 1e,f ------------------------------------


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
  axis.text.x = element_text(size=10,color='black',margin = margin(t=0.2, unit = "cm")),
  #y轴标签
  axis.text.y = element_text(size=10,color='black',margin = margin(r=0.2, unit = "cm")),
  #图例
  legend.title =element_blank(),
  #legend.title = element_blank(colour="black", size=10),
  #legend.text = element_blank(),
  legend.text = element_text(colour="black", size=10)
  # remove legend
  #,legend.position="none"
)

setwd("X:/Projects/2019__ProliferationCorrelatedExpression/Manuscript__ProlifCorrExpr/20200929 figures code and table/Figure1/data/figure 1f data brdu")

library(openxlsx)

FIB_Brdu_data <- read.xlsx('./BrdU quantification.xlsx',colNames = TRUE,rowNames = TRUE,sheet = 1)
ESC_Brdu_data <- read.xlsx('./BrdU quantification.xlsx',colNames = TRUE,rowNames = TRUE,sheet = 2)

library(reshape2)

FIB_Brdu_data  <- as.data.frame(t(FIB_Brdu_data))
FIB_Brdu_data$day <- rownames(FIB_Brdu_data)

FIB_Brdu_data <- melt(FIB_Brdu_data,id.var='day')


library(ggplot2)


ggplot(data = FIB_Brdu_data,aes(x=day,y=value,fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#e67cbd", "#7C285C")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_fib20201125.pdf', width = 80, height = 50, units = "mm")



ESC_Brdu_data  <- as.data.frame(t(ESC_Brdu_data))
ESC_Brdu_data$day <- rownames(ESC_Brdu_data)

ESC_Brdu_data <- melt(ESC_Brdu_data,id.var='day')


library(ggplot2)
library(cowplot)

ggplot(data = ESC_Brdu_data,aes(x=day,y=value,fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#87bbe8", "#2171B5")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_ESC20201125.pdf', width = 80, height = 50, units = "mm")

# figure 1 s -----------

FIB_Brdu_data2 <- read.xlsx('./Brdu_data.xlsx',colNames = TRUE,rowNames = F,sheet = 2)
ESC_Brdu_data2 <- read.xlsx('./Brdu_data.xlsx',colNames = TRUE,rowNames = F,sheet = 1)


FIB_Brdu_data2[,c(3:5)] <- FIB_Brdu_data2[,c(3:5)]*100
ESC_Brdu_data2[,c(3:4)] <- ESC_Brdu_data2[,c(3:4)]*100

FIB_Brdu_data2 <- FIB_Brdu_data2[FIB_Brdu_data2$treatment=="Fast cells DMSO" | FIB_Brdu_data2$treatment=="Slow cells DMSO",]
ESC_Brdu_data2 <- ESC_Brdu_data2[ESC_Brdu_data2$treatment=="Fast cells DMSO" | ESC_Brdu_data2$treatment=="Slow cells DMSO",]


# fib Replicate.1
ggplot(data = FIB_Brdu_data2,aes(x=day,y=Replicate.1,fill=treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#e67cbd", "#7C285C"),labels = c("Fast", "Slow")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_fib Replicate.1.20201125.pdf', width = 80, height = 50, units = "mm")

# fib Replicate.2
ggplot(data = FIB_Brdu_data2,aes(x=day,y=Replicate.2,fill=treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#e67cbd", "#7C285C"),labels = c("Fast", "Slow")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_fib Replicate.2.20201125.pdf', width = 80, height = 50, units = "mm")

# fib Replicate.3
ggplot(data = FIB_Brdu_data2,aes(x=day,y=Replicate.3,fill=treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#e67cbd", "#7C285C"),labels = c("Fast", "Slow")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_fib Replicate.3.20201125.pdf', width = 80, height = 50, units = "mm")


# esc Replicate.1
ggplot(data = ESC_Brdu_data2,aes(x=day,y=Replicate.1,fill=treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#87bbe8", "#2171B5"),labels = c("Fast", "Slow")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_esc Replicate.1.20201125.pdf', width = 80, height = 50, units = "mm")


# esc Replicate.2
ggplot(data = ESC_Brdu_data2,aes(x=day,y=Replicate.2,fill=treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#87bbe8", "#2171B5"),labels = c("Fast", "Slow")) +
  labs(x='',y='% of cells in S-phase') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_esc Replicate.2.20201125.pdf', width = 80, height = 50, units = "mm")

# binomal test for brdu fib and esc----

prob <- 0.5^3

binom.test(4, 4, p = prob,
           alternative = c("greater"),
           conf.level = 0.95)
# p-value = 0.0002441

prob <- 0.5^3

binom.test(3, 3, p = prob,
           alternative = c("greater"),
           conf.level = 0.95)
# p-value = 0.001953



# viability plot---------------------

setwd("E:/Project/2019__proliferation/figures/20200929 figures code and table/Figure1")






library(openxlsx)

esc1_viability_data <- read.xlsx('./data/figure 1f data brdu/viability_data.xlsx',colNames = TRUE,rowNames = F,sheet = 1)

esc2_viability_data <- read.xlsx('./data/figure 1f data brdu/viability_data.xlsx',colNames = TRUE,rowNames = F,sheet = 2)

fib1_viability_data <- read.xlsx('./data/figure 1f data brdu/viability_data.xlsx',colNames = TRUE,rowNames = F,sheet = 3)

fib2_viability_data <- read.xlsx('./data/figure 1f data brdu/viability_data.xlsx',colNames = TRUE,rowNames = F,sheet = 4)


esc_viability <- rbind(esc1_viability_data[c(1,5),],
                       esc2_viability_data[c(1,5),])

esc_viability$rep <- c(1,1,2,2)

esc_viability <- esc_viability[,c(1,3,4,5,6)]

library(reshape)


esc_viability <- melt(esc_viability,id.vars=c("celltype","rep"))
esc_viability$rep <- as.factor(esc_viability$rep)

esc_viability$celltype_rep <- paste(esc_viability$celltype,
                                     esc_viability$rep,
                                    sep = ' rep. ')


esc_viability$celltype_rep_day <- paste(esc_viability$celltype_rep,
                                    esc_viability$variable,
                                    sep = ' ')

esc_viability$celltype_rep_day <- factor(esc_viability$celltype_rep_day,levels =esc_viability$celltype_rep_day[c(1,3,2,4,5,7,6,8,9,11,10,12)] )

library(ggplot2)
library(cowplot)



ggplot(data = esc_viability,aes(x=variable,y=value,fill=celltype_rep_day)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c(rep("#e67cbd",2),rep("#7C285C",2),rep("grey",2),rep("red",2),rep("blue",2),rep("green",2))) +
  labs(x='',y='% of viabiliy') +
  ggplot_theme_pdf_temp

ggsave2('E:/Project/2019__proliferation/figures/Figure1/brdu_ESC.pdf', width = 80, height = 50, units = "mm")



