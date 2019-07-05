setwd(getwd())

library(tidyr)
library(lawstat)
library(stringr)
library(ggsignif)
library(dplyr)
library(nlme)

library(ggplot2)
library(ggpubr)
library(multcomp)
library(xlsx)
my_comparisons_snr<- list( c("Raw", "SI-GEVD"), c("CCA", "SI-GEVD") )

transparent_legend =  theme(plot.title = element_text(hjust = 0.5),
                            legend.background = element_rect(fill = "transparent"),
                            legend.key = element_rect(fill = "transparent", 
                                                      color = "transparent") )

data_all <- read.xlsx("relMSE_LR_TRFs_120s.xlsx",sheetIndex = 1)
getwd()
data_all$SNR = as.factor(data_all$SNR)
data_all$Triallength = as.factor(data_all$Triallength)
data_all$Method = as.factor(data_all$Method)
data_all$Method  <- factor(data_all$Method , levels = c("Raw", "CCA", "SI-GEVD"))

# log scale ticks
ggplot(data_all, aes(x=Method, y=relMSE, colour=Method))+  theme_bw(base_size = 20)+ #20) +
  ylab('relative MSE')+ # ylim(0, 150)+
  ggtitle('EEG SNR (dB)')+ 
  geom_boxplot() +
  scale_y_log10()+
  facet_grid(.~SNR)+
  # scale_x_discrete(limits = rev(levels(data_all$SNR)))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  transparent_legend+
  theme(legend.position = "bottom")+
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..),p.adjust.method = "bonferroni",comparisons = my_comparisons_snr, paired = TRUE) +
  annotation_logticks(base = 10, sides = "l", alpha = 0.5)

