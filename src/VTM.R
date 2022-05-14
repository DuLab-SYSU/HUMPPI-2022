library(RColorBrewer)
library(tidyverse)
mytheme <- theme(panel.grid.major = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                                 linetype = "longdash"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 13),
                 legend.title = element_text(size = 16),
                 legend.background = element_blank())

###########################
### Virally-targeted genes
vtgs <- read.table("../result/virally-targeted_genes29.txt",
                   header = TRUE, sep = '\t', quote = "")
IRGs <- read.table("../data/IRGs/all_IRGs.txt", header = FALSE)
vtgs$Immunity <- apply(vtgs[,1,drop = FALSE],1,
                       function(x) ifelse(x %in% IRGs$V1, "IRG", "non-IRG"))

all_vtg_num <- vtgs %>%
  group_by(Abbre, Baltimore, Family) %>%
  dplyr::summarise(Number = n())
virus_limit <- all_vtg_num$Abbre[rev(order(all_vtg_num$Number))]

library(plyr)
vtg_num <- ddply(vtgs, .(Abbre,Baltimore,Family,Immunity), nrow)
vtg_num <- dplyr::rename(vtg_num, Number=V1)

############################
### Virally-targeted modules

vtgs <- read.table("../intermediate/virally-targeted_genes.txt", 
                   header = TRUE, sep = '\t', quote = "")
virus2abbre <- read.table("../result/virus2abbre.txt",
                          header = TRUE, sep = "\t", quote = "")

### Number of virus-targeted modules (VTMs)
library(plyr)
vtms <- vtgs[,2:5]
vtms <- vtms[vtms$Module!=0,]
vtms <- plyr::ddply(vtms, .(Virus, Baltimore, Family, Module), nrow)
vtms <- dplyr::rename(vtms, Overlap_num=V1)

### hypergeometric test
vhf <- vtgs[,2,drop=FALSE]
vhf <- plyr::ddply(vhf, .(Virus), nrow)
vhf <- dplyr::rename(vhf, VTG_num=V1)

vtms2 <- merge(vtms, vhf, by.x = 'Virus', by.y = 'Virus')

modules <- read.table("../intermediate/batch_MCL_out/out.HUMPPI2022.I228", 
                      header = FALSE, sep = " ")
modules$Size <- apply(modules[,1,drop = FALSE], 1,
                      function(x) length(stringr::str_split(x, "\t")[[1]]))
modules <- modules %>%
  filter(Size >= 3) %>%
  dplyr::select(Size) %>%
  mutate(Module = 1:1225)

vtms2 <- merge(vtms2, modules, by = 'Module')


vtms2$P_value <- phyper(vtms2$Overlap_num-1, vtms2$VTG_num, 17402, vtms2$Size,
                       lower.tail = F)
vtms2$ADJ_pvalue <-p.adjust(vtms2[,8], method="BH", n=length(vtms2[,8]))

vtms3 <- merge(vtms2, virus2abbre, by.x = 'Virus', by.y = 'Virus')
vtms3 <- vtms3[vtms3$ADJ_pvalue<=0.05,c(10,1,3,4,2,8,9)]

write.table(vtms3[,c(2,1,3:7)],
            '../result/virally-targeted_modules29.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)


library(plyr)
vtm_num <- ddply(vtms3[,c(1,3,4)], .(Abbre, Baltimore, Family), nrow) #total number of VTC
vtm_num <- dplyr::rename(vtm_num, Number=V1)
virus_limit2 <- vtm_num[,1][rev(order(vtm_num[,4]))]

module_info <- read.table("../intermediate/module_info.txt",
                          header = TRUE, sep = "\t")
vtms4 <- dplyr::select(module_info, Module, Enriched, Immune_enriched) %>%
  right_join(vtms3, by = "Module")

#Number of Immunity and  non-immunity related VTM
vtm_num2 <- ddply(vtms4[,c(4,6,7,3)], .(Abbre, Baltimore, Family, Immune_enriched), nrow) 
vtm_num2 <- dplyr::rename(vtm_num2, Number=V1)

# VTG and VTM number in a common barplot
p_vtg2 <- ggplot(vtg_num, aes(Abbre, Number, fill = factor(Immunity, levels = c("non-IRG", "IRG")))) +
  geom_bar(stat = 'identity', alpha = 0.75) +
  scale_x_discrete(limits = rev(virus_limit)) +
  scale_y_continuous(expand = c(0.02,0.02)) +
  scale_fill_manual(values = c('#999999', '#FF3300')) +
  labs(fill = "Gene Class", y = "Gene Number", x = "") +
  coord_flip() +
  mytheme +
  theme(panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.7, 0.1),
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))

p_vtm2 <- ggplot(vtm_num2, aes(Abbre, Number, fill = Immune_enriched)) +
  geom_bar(stat = 'identity', alpha = 0.75) +
  scale_x_discrete(limits = rev(virus_limit)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = c('#999999', '#FF3300')) +
  labs(y = "Community Number", x = "") +
  coord_flip() +
  mytheme +
  theme(panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.7, 0.1),
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))

p_s3b <- ggpubr::ggarrange(p_vtg2,p_vtm2, nrow=1, labels = c("B", ""))
ggsave(p_s3b,
       filename = "../result/Fig_S3B.pdf",
       width = 10.5,
       height = 7.5)


### coronaviruses targeted communities
covs <- c("MERS-CoV", "SARS-CoV-1", "SARS-CoV-2")
ctms <- filter(vtms4, Abbre %in% covs) %>%
  dplyr::select(Virus, Module, Enriched, Immune_enriched) %>%
  arrange(Virus)

write.table(ctms, 
            '../intermediate/CoVs-modules.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

