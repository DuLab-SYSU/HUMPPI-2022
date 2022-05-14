### Module-module association network

edge_info <- read.table("../intermediate/edges_info.txt", 
                        header = TRUE, sep = "\t")

module_link <- subset(edge_info, Category != "-", select = c(3,4))
module_link$Module_A <- as.numeric(as.character(module_link$Module_A))
module_link$Module_B <- as.numeric(as.character(module_link$Module_B))

module_link2 <- module_link

for(i in 1:nrow(module_link2)){
  if(module_link2[i,1]>module_link2[i,2]){
    module_link2[i,1] = module_link[i,2]
    module_link2[i,2] = module_link[i,1]
  }
}

library(plyr)
module_link2 <- ddply(module_link2, .(Module_A, Module_B), nrow)
module_link2 <- dplyr::rename(module_link2, Number=V1)

intra_link <- subset(module_link2, Module_A == Module_B, select = c(1,3))
intra_link <- dplyr::rename(intra_link, Module = Module_A)

inter_link <- subset(module_link2, Module_A != Module_B)
inter_link <- dplyr::rename(inter_link, nAB=Number)

out_link <- function(id){
  num <- sum(c(inter_link[inter_link$Module_A==id,3], inter_link[inter_link$Module_B==id,3]))
  return(num)
}

in_link <- function(id){
  num <- intra_link[intra_link$Module==id,2]
  return(num)
}

inter_link$nA <- apply(inter_link[,1,drop = FALSE], 1, out_link)
inter_link$nB <- apply(inter_link[,2,drop = FALSE], 1, out_link)
inter_link$na <- apply(inter_link[,1,drop = FALSE], 1, in_link)
inter_link$nb <- apply(inter_link[,2,drop = FALSE], 1, in_link)

# Fisher exact test
nT <- nrow(edge_info) # all the links in HUMPPI-2022

fisherfun <- function(x){
  tmp<- fisher.test(matrix(c(x[3], x[5]-x[3], x[4]-x[3],
                             nT-x[6]-x[7]-x[5]-x[4]+x[3]), 
                           nrow = 2), alternative = "greater")
  return(tmp$p.value)
}

inter_link$P_value <- apply(inter_link, 1, 
                            function(x) fisher.test(matrix(c(x[3], x[5]-x[3], x[4]-x[3],
                                                              nT-x[6]-x[7]-x[5]-x[4]+x[3]), 
                                                            nrow = 2), alternative = "greater")$p.value)

# multiple testing correction
inter_link$ADJ_pvalue <- p.adjust(inter_link[,8], method="BH",
                                     n=length(inter_link[,8]))
inter_link$`-log10(ADJ_pvalue)` <- -log10(inter_link$ADJ_pvalue)

inter_link2 <- inter_link[inter_link$ADJ_pvalue<0.01,]
write.table(inter_link2[,c(1,2,9)], "../result/module_associations.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(inter_link2[,1:2],
            "../intermediate/module_associations.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

### Modules enrichment of IRGs
modules <- read.table("../result/Module_members.txt",
                      header = TRUE, sep = "\t")
module_size <- modules %>%
  filter(Module <= 1225) %>%
  group_by(Module) %>%
  dplyr::summarise(Size = n())

irgs <- read.table('../data/IRGs/all_IRGs.txt', header = FALSE)
colnames(irgs) <- c('Gene')

irgs2 <- modules %>%
  filter(Module <= 1225) %>%
  filter(Gene %in% irgs$Gene) %>%
  dplyr::rename(IRG = Gene)

IRGs_module <- plyr::ddply(irgs2, .(Module), nrow)
IRGs_module <- dplyr::rename(IRGs_module, IRGs_num=V1)

IRGs_module2 <- merge(IRGs_module, module_size, by = 'Module', 
                      all = TRUE)
IRGs_module2[is.na(IRGs_module2)] <- 0
IRGs_module2 <- dplyr::rename(IRGs_module2, Module_size=Size)

IRGs_module2$P_value <- phyper(IRGs_module2$IRGs_num-1, 3724, 7478, 
                                IRGs_module2$Module_size, lower.tail = F)
IRGs_module2$ADJ_pvalue <- p.adjust(IRGs_module2[,4], method = "BH",
                                     n = length(IRGs_module2[,4]))

myFun <- function(IRGs_module2, c1, c2) {
  ifelse(IRGs_module2[c1]<=0.01, 'Enriched', 
         ifelse(IRGs_module2[c2]>1, '2+', '0-1'))
}

IRGs_module2$Enriched <- apply(IRGs_module2, 1, myFun, c1 = 'ADJ_pvalue', c2 = 'IRGs_num')

module_count <- plyr::ddply(IRGs_module2[,6,drop=FALSE], .(Enriched), nrow)
module_count <- dplyr::rename(module_count, module_num=V1)
module_count$module_percent <- paste0(
  round(module_count$module_num/sum(module_count$module_num),4)*100,"%")

p_IRGs <- ggplot(IRGs_module2, aes(Enriched,fill=Enriched)) +
  geom_bar(stat = 'count') +
    scale_fill_manual(values = c('#999999', '#00CCCC', '#006666')) +
  geom_text(data = module_count, size = 4, color = "white",
            aes(x = Enriched, y = module_num - 10, label = module_percent)) +
  scale_y_continuous(expand = c(0.01,0.01), limits = c(0, 620)) +
  labs(x="", y='Module count', fill='IRG number') +
  mytheme +
  theme(legend.position = "none",
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))


### Modules enrichment of immune-related processes
module_imm <- read.table("../intermediate/module_IR_process.txt", 
                         header = TRUE, sep = "\t")

module_imm <- left_join(module_imm, IRGs_module2[,c(1,6)], by = "Module") %>%
  dplyr::select(Module, Enriched, Description:geneID) %>%
  dplyr::rename(IRG_contained = Enriched)

write.table(module_imm, "../result/module_IR_process.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

### Modules infomation
IRGs_module2$Immune_enriched <- apply(IRGs_module2[,1,drop = FALSE], 1,
                                      function(x) ifelse(x %in% module_imm$Module, "Yes", "No"))
centrality <- read.table("../intermediate/Module_centrality.txt",
                         header = TRUE, sep = '\t')

module_info <- left_join(IRGs_module2, centrality, by = "Module") %>%
  dplyr::select(-c(P_value, ADJ_pvalue))
module_info[is.na(module_info)] <- "-"

write.table(module_info,
            '../intermediate/module_info.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

### module centrality
module_centrality <- subset(module_info, Degree != '-')

# Two-sample Wilcoxon tests
p1 <- wilcox.test(as.numeric(as.character(module_centrality[module_centrality$Enriched == 'Enriched',6])),
            as.numeric(as.character(module_centrality[module_centrality$Enriched != 'Enriched',6])))
p1 <- p1$p.value

p2 <- wilcox.test(as.numeric(as.character(module_centrality[module_centrality$Enriched == '2+',6])),
            as.numeric(as.character(module_centrality[module_centrality$Enriched != '0-1',6])))
p2 <- p2$p.value

p3 <- wilcox.test(as.numeric(as.character(module_centrality[module_centrality$Immune_enriched == 'Yes',6])),
            as.numeric(as.character(module_centrality[module_centrality$Immune_enriched == 'No',6])))
p3 <- p3$p.value

# plot
p_degree <- ggplot(module_centrality, aes(Enriched, as.numeric(Degree), fill = Enriched)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#999999', '#00CCCC', '#006666')) +
  labs(x="", y = 'Degree') +
  mytheme +
  theme(legend.position = "none",
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))

p_degree2 <- ggplot(module_centrality, aes(Immune_enriched, as.numeric(Degree), color = Immune_enriched)) +
  geom_boxplot(lwd=1.2) +
  scale_color_manual(values = c('#999999', '#FF3300')) +
  labs(x="", y = 'Degree') +
  mytheme +
  theme(legend.position = "none",
        plot.margin = unit(x=c(6,6,6,6),units = 'mm'))

p_2hij <- ggpubr::ggarrange(p_IRGs, p_degree, p_degree2, nrow=1, labels = c("H", "I", "J"))
ggsave(p_2hij,
       filename = "../result/Fig_2HIJ.pdf",
       width = 12,
       height = 6.5)
