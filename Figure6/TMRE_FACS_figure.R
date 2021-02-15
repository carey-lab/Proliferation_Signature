

# 20201219 ------------------------------

library(flowCore)

require(ggplot2)
require(ggridges)



#  esc tmre --------------

esc_tmre_unstained <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY0/05122017_Unstained.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day0 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY0/05122017_Day 0.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day2_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P6 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day2_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P6 2.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day2_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P7 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day2_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P7 2.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day2_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P8 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day2_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY2/08122017_P8 2.fcs", which.lines=NULL, transformation="linearize")))


esc_tmre_day3_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P6 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day3_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P6 2.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day3_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P7 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day3_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P7 2.fcs", which.lines=NULL, transformation="linearize")))

esc_tmre_day3_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P8 1.fcs", which.lines=NULL, transformation="linearize")))
esc_tmre_day3_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/ESC_Sorting_fcs/DAY3/08122017_P8 2.fcs", which.lines=NULL, transformation="linearize")))




esc_tmre_unstained <- data.frame(days="esc_tmre_unstained",PE_A=esc_tmre_unstained[,22])
esc_tmre_day0 <- data.frame(days="esc_tmre_day0",PE_A=esc_tmre_day0[,22])

esc_tmre_day2_low1 <- data.frame(days="esc_tmre_day2_low1",PE_A=esc_tmre_day2_low1[,22])
esc_tmre_day2_low2 <- data.frame(days="esc_tmre_day2_low2",PE_A=esc_tmre_day2_low2[,22])
esc_tmre_day2_med1 <- data.frame(days="esc_tmre_day2_med1",PE_A=esc_tmre_day2_med1[,22])
esc_tmre_day2_med2 <- data.frame(days="esc_tmre_day2_med2",PE_A=esc_tmre_day2_med2[,22])
esc_tmre_day2_high1 <- data.frame(days="esc_tmre_day2_high1",PE_A=esc_tmre_day2_high1[,22])
esc_tmre_day2_high2 <- data.frame(days="esc_tmre_day2_high2",PE_A=esc_tmre_day2_high2[,22])


esc_tmre_day3_low1 <- data.frame(days="esc_tmre_day3_low1",PE_A=esc_tmre_day3_low1[,22])
esc_tmre_day3_low2 <- data.frame(days="esc_tmre_day3_low2",PE_A=esc_tmre_day3_low2[,22])
esc_tmre_day3_med1 <- data.frame(days="esc_tmre_day3_med1",PE_A=esc_tmre_day3_med1[,22])
esc_tmre_day3_med2 <- data.frame(days="esc_tmre_day3_med2",PE_A=esc_tmre_day3_med2[,22])
esc_tmre_day3_high1 <- data.frame(days="esc_tmre_day3_high1",PE_A=esc_tmre_day3_high1[,22])
esc_tmre_day3_high2 <- data.frame(days="esc_tmre_day3_high2",PE_A=esc_tmre_day3_high2[,22])


esc_tmre_all <- rbind(esc_tmre_unstained,
                      esc_tmre_day0,
                      esc_tmre_day2_low1,
                      esc_tmre_day2_low2,
                      esc_tmre_day2_med1,
                      esc_tmre_day2_med2,
                      esc_tmre_day2_high1,
                      esc_tmre_day2_high2,
                      esc_tmre_day3_low1,
                      esc_tmre_day3_low2,
                      esc_tmre_day3_med1,
                      esc_tmre_day3_med2,
                      esc_tmre_day3_high1,
                      esc_tmre_day3_high2)



esc_tmre_all <- esc_tmre_all[esc_tmre_all$PE_A>0,]
esc_tmre_all$log10PE_A <- log10(esc_tmre_all$PE_A+1)

library(ggplot2)
library(ggridges)
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


esc_tmre_all_day2_all <- esc_tmre_all[esc_tmre_all$days %in% c("esc_tmre_unstained"
                                                               ,"esc_tmre_day0",
                                                           "esc_tmre_day2_low1",
                                                           "esc_tmre_day2_low2",
                                                           "esc_tmre_day2_med1",
                                                           "esc_tmre_day2_med2",
                                                           "esc_tmre_day2_high1",
                                                           "esc_tmre_day2_high2"),]

ggplot(esc_tmre_all_day2_all, aes(x=log10PE_A, color=days)) +
  geom_density() +
  xlim(0,6) +
  title("esc_tmre_all_day2_all")

esc_tmre_all_day3_all <- esc_tmre_all[esc_tmre_all$days %in% c("esc_tmre_day0",
                                                               "esc_tmre_day3_low1",
                                                               "esc_tmre_day3_low2",
                                                               "esc_tmre_day3_med1",
                                                               "esc_tmre_day3_med2",
                                                               "esc_tmre_day3_high1",
                                                               "esc_tmre_day3_high2"),]

ggplot(esc_tmre_all_day3_all, aes(x=log10PE_A, color=days)) +
  geom_density() +
  xlim(0,6) +
  title("esc_tmre_all_day3_all")

esc_tmre_all_day3_rep1 <- esc_tmre_all[esc_tmre_all$days %in% c("esc_tmre_day0",
                                                           "esc_tmre_day3_low1",
                                                           "esc_tmre_day3_med1",
                                                           "esc_tmre_day3_high1"),]

ggplot(esc_tmre_all_day3_rep1, aes(x=log10PE_A, color=days)) +
  geom_density() +
  xlim(0,6) +
  title("esc_tmre_all_day3_rep1")

esc_tmre_all_day3_rep2 <- esc_tmre_all[esc_tmre_all$days %in% c("esc_tmre_day0",
                                                                "esc_tmre_day3_low2",
                                                                "esc_tmre_day3_med2",
                                                                "esc_tmre_day3_high2"),]
ggplot(esc_tmre_all_day3_rep2, aes(x=log10PE_A, color=days)) +
  geom_density() +
  xlim(0,6) +
  title("esc_tmre_all_day3_rep2")


# density ------
ggplot(esc_tmre_all_day3_rep2, aes(x=log10PE_A, color=days)) +
  geom_density() +
  title("esc_tmre_all_day2_rep1")

  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")


# geom_density_ridges -----

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





#  fib tmre --------------

fib_tmre_unstained <- as.data.frame(exprs(read.FCS("./data/rawdata/FIB_Sorting_fcs/DAY 0 + SINGLE CONTROLS/20171023_unstained.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day0 <- as.data.frame(exprs(read.FCS("./data/rawdata/FIB_Sorting_fcs/DAY 0 + SINGLE CONTROLS/20171023_day 0.fcs", which.lines=NULL, transformation="linearize")))




fib_tmre_day1_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_LA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day1_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_LB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day1_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_MA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day1_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_MB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day1_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_HA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day1_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 1/20171024_HB.fcs", which.lines=NULL, transformation="linearize")))




fib_tmre_day2_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_LA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day2_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_LB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day2_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_MA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day2_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_MB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day2_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_HA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day2_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 2/20171025_HB.fcs", which.lines=NULL, transformation="linearize")))


fib_tmre_day3_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_LA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day3_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_LB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day3_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_MA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day3_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_MB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day3_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_HA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day3_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 3/20171026_HB.fcs", which.lines=NULL, transformation="linearize")))


fib_tmre_day4_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_LA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day4_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_LB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day4_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_MA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day4_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_MB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day4_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_HA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day4_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 4/20171027_HB.fcs", which.lines=NULL, transformation="linearize")))



fib_tmre_day5_low1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_LA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day5_low2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_LB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day5_med1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_MA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day5_med2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_MB.fcs", which.lines=NULL, transformation="linearize")))

fib_tmre_day5_high1 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_HA.fcs", which.lines=NULL, transformation="linearize")))
fib_tmre_day5_high2 <- as.data.frame(exprs(read.FCS("./data/rawdata/fib_Sorting_fcs/DAY 5/20171028_HB.fcs", which.lines=NULL, transformation="linearize")))











fib_tmre_unstained <- data.frame(days="fib_tmre_unstained",PE_A=fib_tmre_unstained[,22])
fib_tmre_day0 <- data.frame(days="fib_tmre_day0",PE_A=fib_tmre_day0[,22])



fib_tmre_day1_low1 <- data.frame(days="fib_tmre_day1_low1",PE_A=fib_tmre_day1_low1[,22])
fib_tmre_day1_low2 <- data.frame(days="fib_tmre_day1_low2",PE_A=fib_tmre_day1_low2[,22])
fib_tmre_day1_med1 <- data.frame(days="fib_tmre_day1_med1",PE_A=fib_tmre_day1_med1[,22])
fib_tmre_day1_med2 <- data.frame(days="fib_tmre_day1_med2",PE_A=fib_tmre_day1_med2[,22])
fib_tmre_day1_high1 <- data.frame(days="fib_tmre_day1_high1",PE_A=fib_tmre_day1_high1[,22])
fib_tmre_day1_high2 <- data.frame(days="fib_tmre_day1_high2",PE_A=fib_tmre_day1_high2[,22])



fib_tmre_day2_low1 <- data.frame(days="fib_tmre_day2_low1",PE_A=fib_tmre_day2_low1[,22])
fib_tmre_day2_low2 <- data.frame(days="fib_tmre_day2_low2",PE_A=fib_tmre_day2_low2[,22])
fib_tmre_day2_med1 <- data.frame(days="fib_tmre_day2_med1",PE_A=fib_tmre_day2_med1[,22])
fib_tmre_day2_med2 <- data.frame(days="fib_tmre_day2_med2",PE_A=fib_tmre_day2_med2[,22])
fib_tmre_day2_high1 <- data.frame(days="fib_tmre_day2_high1",PE_A=fib_tmre_day2_high1[,22])
fib_tmre_day2_high2 <- data.frame(days="fib_tmre_day2_high2",PE_A=fib_tmre_day2_high2[,22])

fib_tmre_day3_low1 <- data.frame(days="fib_tmre_day3_low1",PE_A=fib_tmre_day3_low1[,22])
fib_tmre_day3_low2 <- data.frame(days="fib_tmre_day3_low2",PE_A=fib_tmre_day3_low2[,22])
fib_tmre_day3_med1 <- data.frame(days="fib_tmre_day3_med1",PE_A=fib_tmre_day3_med1[,22])
fib_tmre_day3_med2 <- data.frame(days="fib_tmre_day3_med2",PE_A=fib_tmre_day3_med2[,22])
fib_tmre_day3_high1 <- data.frame(days="fib_tmre_day3_high1",PE_A=fib_tmre_day3_high1[,22])
fib_tmre_day3_high2 <- data.frame(days="fib_tmre_day3_high2",PE_A=fib_tmre_day3_high2[,22])


fib_tmre_day4_low1 <- data.frame(days="fib_tmre_day4_low1",PE_A=fib_tmre_day4_low1[,22])
fib_tmre_day4_low2 <- data.frame(days="fib_tmre_day4_low2",PE_A=fib_tmre_day4_low2[,22])
fib_tmre_day4_med1 <- data.frame(days="fib_tmre_day4_med1",PE_A=fib_tmre_day4_med1[,22])
fib_tmre_day4_med2 <- data.frame(days="fib_tmre_day4_med2",PE_A=fib_tmre_day4_med2[,22])
fib_tmre_day4_high1 <- data.frame(days="fib_tmre_day4_high1",PE_A=fib_tmre_day4_high1[,22])
fib_tmre_day4_high2 <- data.frame(days="fib_tmre_day4_high2",PE_A=fib_tmre_day4_high2[,22])



fib_tmre_day5_low1 <- data.frame(days="fib_tmre_day5_low1",PE_A=fib_tmre_day5_low1[,22])
fib_tmre_day5_low2 <- data.frame(days="fib_tmre_day5_low2",PE_A=fib_tmre_day5_low2[,22])
fib_tmre_day5_med1 <- data.frame(days="fib_tmre_day5_med1",PE_A=fib_tmre_day5_med1[,22])
fib_tmre_day5_med2 <- data.frame(days="fib_tmre_day5_med2",PE_A=fib_tmre_day5_med2[,22])
fib_tmre_day5_high1 <- data.frame(days="fib_tmre_day5_high1",PE_A=fib_tmre_day5_high1[,22])
fib_tmre_day5_high2 <- data.frame(days="fib_tmre_day5_high2",PE_A=fib_tmre_day5_high2[,22])






fib_tmre_all <- rbind(fib_tmre_unstained,
                      fib_tmre_day0,
                      fib_tmre_day1_low1,
                      fib_tmre_day1_low2,
                      fib_tmre_day1_med1,
                      fib_tmre_day1_med2,
                      fib_tmre_day1_high1,
                      fib_tmre_day1_high2,
                      fib_tmre_day2_low1,
                      fib_tmre_day2_low2,
                      fib_tmre_day2_med1,
                      fib_tmre_day2_med2,
                      fib_tmre_day2_high1,
                      fib_tmre_day2_high2,
                      fib_tmre_day3_low1,
                      fib_tmre_day3_low2,
                      fib_tmre_day3_med1,
                      fib_tmre_day3_med2,
                      fib_tmre_day3_high1,
                      fib_tmre_day3_high2,
                      fib_tmre_day4_low1,
                      fib_tmre_day4_low2,
                      fib_tmre_day4_med1,
                      fib_tmre_day4_med2,
                      fib_tmre_day4_high1,
                      fib_tmre_day4_high2,
                      fib_tmre_day5_low1,
                      fib_tmre_day5_low2,
                      fib_tmre_day5_med1,
                      fib_tmre_day5_med2,
                      fib_tmre_day5_high1,
                      fib_tmre_day5_high2)



fib_tmre_all <- fib_tmre_all[fib_tmre_all$PE_A>0,]
fib_tmre_all$log10PE_A <- log10(fib_tmre_all$PE_A+1)


display.brewer.all()
display.brewer.pal(9,"Blues")
brewer.pal(9,"Greens")
brewer.pal(9,"Purples")
brewer.pal(9,"Purples")
brewer.pal(9,"Reds")
brewer.pal(9,"Set1")
brewer.pal(9,"Blues")

colorRampPalette(brewer.pal(9,"Blues"))(100)


fib_tmre_all_day1_all <- fib_tmre_all[fib_tmre_all$days %in% c(
                                                               "fib_tmre_day0",
                                                               "fib_tmre_day1_low1",
                                                               "fib_tmre_day1_low2",
                                                               "fib_tmre_day1_med1",
                                                               "fib_tmre_day1_med2",
                                                               "fib_tmre_day1_high1",
                                                               "fib_tmre_day1_high2"),]

ggplot(fib_tmre_all_day1_all, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_day1_all")
  

fib_tmre_all_day2_all <- fib_tmre_all[fib_tmre_all$days %in% c("fib_tmre_day0",
                                                               "fib_tmre_day2_low1",
                                                               "fib_tmre_day2_low2",
                                                               "fib_tmre_day2_med1",
                                                               "fib_tmre_day2_med2",
                                                               "fib_tmre_day2_high1",
                                                               "fib_tmre_day2_high2"),]

ggplot(fib_tmre_all_day2_all, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_day2_all")


fib_tmre_all_day3_all <- fib_tmre_all[fib_tmre_all$days %in% c(
                                                               "fib_tmre_day0",
                                                               "fib_tmre_day3_low1",
                                                               "fib_tmre_day3_low2",
                                                               "fib_tmre_day3_med1",
                                                               "fib_tmre_day3_med2",
                                                               "fib_tmre_day3_high1",
                                                               "fib_tmre_day3_high2"),]

ggplot(fib_tmre_all_day3_all, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_day3_all")


fib_tmre_all_day4_all <- fib_tmre_all[fib_tmre_all$days %in% c(
                                                               "fib_tmre_day0",
                                                               "fib_tmre_day4_low1",
                                                               "fib_tmre_day4_low2",
                                                               "fib_tmre_day4_med1",
                                                               "fib_tmre_day4_med2",
                                                               "fib_tmre_day4_high1",
                                                               "fib_tmre_day4_high2"),]

ggplot(fib_tmre_all_day4_all, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_day4_all")

fib_tmre_all_day5_all <- fib_tmre_all[fib_tmre_all$days %in% c(
                                                               "fib_tmre_day0",
                                                               "fib_tmre_day5_low1",
                                                               "fib_tmre_day5_low2",
                                                               "fib_tmre_day5_med1",
                                                               "fib_tmre_day5_med2",
                                                               "fib_tmre_day5_high1",
                                                               "fib_tmre_day5_high2"),]



ggplot(fib_tmre_all_day5_all, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_day5_all")



# all_low1 ------

fib_tmre_all_low1 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_low1",
  "fib_tmre_day2_low1",
  "fib_tmre_day3_low1",
  "fib_tmre_day4_low1",
  "fib_tmre_day5_low1"),]

ggplot(fib_tmre_all_low1, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_low1")


# all_low2 ------

fib_tmre_all_low2 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_low2",
  "fib_tmre_day2_low2",
  "fib_tmre_day3_low2",
  "fib_tmre_day4_low2",
  "fib_tmre_day5_low2"),]

ggplot(fib_tmre_all_low2, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_low2")


# all_high1 ------

fib_tmre_all_high1 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_high1",
  "fib_tmre_day2_high1",
  "fib_tmre_day3_high1",
  "fib_tmre_day4_high1",
  "fib_tmre_day5_high1"),]

ggplot(fib_tmre_all_high1, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_high1")


# all_high2 ------

fib_tmre_all_high2 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_high2",
  "fib_tmre_day2_high2",
  "fib_tmre_day3_high2",
  "fib_tmre_day4_high2",
  "fib_tmre_day5_high2"),]

ggplot(fib_tmre_all_high2, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_high2")



# all_med1 ------

fib_tmre_all_med1 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_med1",
  "fib_tmre_day2_med1",
  "fib_tmre_day3_med1",
  "fib_tmre_day4_med1",
  "fib_tmre_day5_med1"),]

ggplot(fib_tmre_all_med1, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_med1")


# all_med2 ------

fib_tmre_all_med2 <- fib_tmre_all[fib_tmre_all$days %in% c(
  "fib_tmre_unstained",
  "fib_tmre_day0",
  "fib_tmre_day1_med2",
  "fib_tmre_day2_med2",
  "fib_tmre_day3_med2",
  "fib_tmre_day4_med2",
  "fib_tmre_day5_med2"),]

ggplot(fib_tmre_all_med2, aes(x=log10PE_A, color=days)) +
  geom_density() + 
  xlim(0,6) +
  labs(title = "fib_tmre_all_med2")






# geom_density_ridges -----


all_data <- list(
                 fib_tmre_all_low1,
                 fib_tmre_all_low2,
                 fib_tmre_all_med1,
                 fib_tmre_all_med2,
                 fib_tmre_all_high1,
                 fib_tmre_all_high2)


all_data_name <- c(
                   "fib_tmre_all_low1",
                   "fib_tmre_all_low2",
                   "fib_tmre_all_med1",
                   "fib_tmre_all_med2",
                   "fib_tmre_all_high1",
                   "fib_tmre_all_high2")

for (i in 1:6) {
  
  
  plotfibtmre <- ggplot(all_data[[i]], aes(x = all_data[[i]]$log10PE_A, 
                                          y = all_data[[i]]$days, 
                                          fill = all_data[[i]]$days)) + 
    geom_density_ridges(alpha = 0.7) + 
    #scale_fill_manual(values = c("#2171B5"), name = "Cell type") + ylab("Number of doublings") + xlab("log2(FITC)") + 
    #ggtitle("Concentration: 10 nM") + 
    labs(x="TMRE (log10(PE))",y="Days after sorting",title = all_data_name[i]) +
    theme_classic() + 
    theme(axis.line = element_blank(), panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=12),
          axis.text = element_text(size=10,color='black',
                                   margin = margin(t=0.3, unit = "cm")),
          legend.position = "none") 
  
  show(plotfibtmre)
}

