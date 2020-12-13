

# 2019 06 29 to check the ESCs cell pluripotency and differentiation related genes expression------
library(Hmisc)
setwd("E:/Project/2019__proliferation/2.Metabolism/stem cell")

# 2015 cell stem cell--------------------------
pluripotency_gene_stem <- c("Dppa3","Esrrb","Fgf4","Klf4",
                            "Nanog","Oct4","Pou5f1","Sall4",
                            "Sox2","Stat3","Tbx3","Zfp42")

pluripotency_gene_stem <- capitalize(pluripotency_gene_stem)


twoc_like_gene_stem <- c("zscan4c","zscan4d","zscan4f",
                         "tmem92","usp17lc","usp17la",
                         "usp17lb","zscan4b","usp17le",
                         "pramef25","zscan4e",'zfp352')


twoc_like_gene_stem <- capitalize(twoc_like_gene_stem)



differentiation_gene_stem <- c("runx1","bhlhe40","fosl2",
                               "klf7","klf6","zfhx3",
                               "hes1","tgfb1i1","ikzf2",
                               "meis1","zkscan3","hdac7",
                               "gata3","id1","rsl1")

differentiation_gene_stem <- capitalize(differentiation_gene_stem)


#-  stem  ------------------------------

library(org.Mm.eg.db)

keytypes(org.Mm.eg.db)



mouse_ensembl <- select(org.Mm.eg.db,keys=twoc_like_gene_stem,columns = c("SYMBOL","ENSEMBL"),keytype = "SYMBOL")

twoc_like_gene_stem <- mouse_ensembl
twoc_like_gene_stem$ENSEMBL[11] <- "ENSMUSG00000095936"


mouse_ensembl <- select(org.Mm.eg.db,keys=pluripotency_gene_stem,columns = c("SYMBOL","ENSEMBL","ENTREZID"),keytype = "SYMBOL")

pluripotency_gene_stem <- mouse_ensembl
pluripotency_gene_stem$ENSEMBL[6] <- "ENSMUSG00000024406"


mouse_ensembl <- select(org.Mm.eg.db,keys=differentiation_gene_stem,columns = c("SYMBOL","ENSEMBL"),keytype = "SYMBOL")

differentiation_gene_stem <- mouse_ensembl


# read in CFSE ESCs expression file-------
mESC <- read.table("./ESC_all_data_20190517.tsv",header = T,sep = "\t")
head(mESC)
colSums(mESC)

mESC <- sweep(mESC,2,colSums(mESC),`/`)*10^6
dim(mESC)
mESC <- mESC[rowMeans(mESC[,c(1:4,7:10)])>0.001,]

mESC$log2fast2slow <- log2((rowMeans(mESC[,7:10])+1)/(rowMeans(mESC[,1:4])+1))

table(mESC$log2fast2slow==0)
mESC[mESC$log2fast2slow==0,]

table(mESC$log2fast2slow>0)[2]/length(mESC$log2fast2slow)# 0.3282727

table(mESC$log2fast2slow<0)[2]/length(mESC$log2fast2slow) # 0.370684

#--figure 2A stem differentiation----------------------------------------------

mESC_differentiation_gene_stem <- mESC[rownames(mESC) %in% differentiation_gene_stem$ENSEMBL, ]

mESC_differentiation_gene_stem$symbol <- differentiation_gene_stem$SYMBOL[match(rownames(mESC_differentiation_gene_stem),differentiation_gene_stem$ENSEMBL)]

mESC_differentiation_gene_stem$log2fast2slow <- log2(rowMeans(mESC_differentiation_gene_stem[,7:10])/rowMeans(mESC_differentiation_gene_stem[,1:4]))

mESC_differentiation_gene_stem <- mESC_differentiation_gene_stem[order(mESC_differentiation_gene_stem$log2fast2slow),]



library(ggsci)

pdf("./mESC_differentiation_regulators_stem20201021.pdf", width = 7, height = 6.5)

nrow(mESC_differentiation_gene_stem)

dotchart(mESC_differentiation_gene_stem$log2fast2slow, labels =mESC_differentiation_gene_stem$symbol,
         cex = 1.5,pt.cex=2, xlab = "log2(fast/slow)",main = "Genes upregulated in lineage commitment (CFSE)",
         pch=21,col="black",bg = "#6BAED6",font.main=1)
abline(v=0,lty="dashed")

dev.off()



binom_test_data <- mESC_differentiation_gene_stem$log2fast2slow

binom.test(table(binom_test_data>0)[1], length(binom_test_data), p = 0.4647868,
           alternative = c("two.sided"),
           conf.level = 0.95)

# p-value = 0.3132

# figure 2B  -stem pluripotency---------------------------------

library(RColorBrewer)
brewer.pal(n = 9, name = "Blues")

brewer.pal(n = 9, name = "Greens")

brewer.pal(n = 9, name = "Purples")

display.brewer.pal(n = 9, name = "Purples")
rowmedian <- function(x) {
  apply(x,1,median)
}

mESC_pluripotency_regulators_stem <- mESC[rownames(mESC) %in% pluripotency_gene_stem$ENSEMBL, ]

mESC_pluripotency_regulators_stem$symbol <- pluripotency_gene_stem$SYMBOL[match(rownames(mESC_pluripotency_regulators_stem),pluripotency_gene_stem$ENSEMBL)]

mESC_pluripotency_regulators_stem$log2fast2slow <- log2(rowMeans(mESC_pluripotency_regulators_stem[,7:10])/rowMeans(mESC_pluripotency_regulators_stem[,1:4]))


mESC_pluripotency_regulators_stem <- mESC_pluripotency_regulators_stem[order(mESC_pluripotency_regulators_stem$log2fast2slow),]


binom_test_data <- mESC_pluripotency_regulators_stem$log2fast2slow

binom.test(table(binom_test_data>0)[1], length(binom_test_data), p = 0.4647868,
           alternative = c("two.sided"),
           conf.level = 0.95)

# p-value = 0.0002187


library(ggsci)

pdf("./mESC_pluripotency_genes_stem20201021.pdf", width = 7, height = 6)

dotchart(mESC_pluripotency_regulators_stem$log2fast2slow,
         labels =mESC_pluripotency_regulators_stem$symbol,
         cex = 1.5,pt.cex=2,xlab = "log2(fast/slow)",main = "Genes upregulated in naive pluripotent state (CFSE)",
         pch=21,col="black",bg = "#9E9AC8",font.main=1,xlim = c(-0.2,0.05))
abline(v=0,lty="dashed")

dev.off()

#-figure 2C stem 2c-like---------------------------------

mESC_2clike_regulators_stem <- mESC[rownames(mESC) %in% twoc_like_gene_stem$ENSEMBL, ]


mESC_2clike_regulators_stem$symbol <- twoc_like_gene_stem$SYMBOL[match(rownames(mESC_2clike_regulators_stem),twoc_like_gene_stem$ENSEMBL)]

mESC_2clike_regulators_stem$log2fast2slow <- log2(rowMeans(mESC_2clike_regulators_stem[,7:10])/rowMeans(mESC_2clike_regulators_stem[,1:4]))


mESC_2clike_regulators_stem <- mESC_2clike_regulators_stem[order(mESC_2clike_regulators_stem$log2fast2slow),]



binom_test_data <- mESC_2clike_regulators_stem$log2fast2slow

binom.test(table(binom_test_data>0)[1], length(binom_test_data), p = 0.4647868,
           alternative = c("two.sided"),
           conf.level = 0.95)

# p-value = 0.0001016





rowMeans(mESC_2clike_regulators_stem[,c(1:4,7:10)])


library(ggsci)

pdf("./mESC_2c-like_stem20201021.pdf", width = 7, height = 6)

dotchart(mESC_2clike_regulators_stem$log2fast2slow,
         labels =mESC_2clike_regulators_stem$symbol,
         cex = 1.5,pt.cex=2,xlab = "log2(fast/slow)",main = "Genes upregulated in 2C-like state (CFSE)",
         pch=21,col="black",bg = "#fa9fb5",
         xlim = c(-4.5,0),
         font.main=1)
abline(v=0,lty="dashed")

dev.off()


















