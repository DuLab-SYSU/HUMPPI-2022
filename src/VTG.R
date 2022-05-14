library(ggplot2)
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

###################################
### Human-Virus PPIs and Human PPIN

### Number of host factors for each virus
vtgs <- read.table("../intermediate/virally-targeted_genes.txt", 
                   header = TRUE, sep = '\t', quote = "")
IRGs <- read.table("../data/IRGs/all_IRGs.txt", header = FALSE)
vtgs$Immunity <- apply(vtgs[,1,drop = FALSE],1,
                       function(x) ifelse(x %in% IRGs$V1, "IRG", "non-IRG"))

gene_num <- as.data.frame(table(vtgs$Virus))
colnames(gene_num) <- c("Virus", "Number")

virus <- gene_num[gene_num$Number>=200,1]
virus <- virus[-which(virus == 'Macaca mulatta polyomavirus 1')] # remove Macaca mulatta polyomavirus 1
virus <- virus[-which(virus == 'Murid herpesvirus 4')] # remove Macaca mulatta polyomavirus 1
gene_num <- gene_num[gene_num$Virus %in% virus,]

virus_abbre <- c('DENV-2', 'EBOV', 'EBV', 'HCV', 'HAdV', 'HCMV', 'HHV-1', 'HIV-1', 'HPV-11', 'HPV-16',
                 'HPV-18', 'HPV-31', 'HPV-33', 'HPV-5', 'HPV-6b', 'HPV-8', 'HPV-9', 'HRSV', 'H1N1', 
                 'H3N2', 'H5N1', 'KSHV', 'LCMV', 'MV', 'MERS-CoV', 'SARS-CoV-1', 'SARS-CoV-2', 'VACV', 'ZIKV')

virus_Baltimore <- read.table("../data/virus_description.txt", 
                              header = TRUE, sep = "\t", quote = "")
virus_Baltimore <- dplyr::select(virus_Baltimore, Parent.Viral.Name:Family) %>%
  distinct() %>%
  dplyr::rename(Virus = Parent.Viral.Name)

virus2abbre <- data.frame(Virus = virus, Abbre = virus_abbre)
factor_num <- merge(gene_num, virus2abbre, by.x = "Virus", by.y = "Virus")
virus_limit <- gene_num[,1][rev(order(gene_num[,2]))]

library(plyr)
vtg_num <- ddply(vtgs, .(Virus,Baltimore,Family,Immunity), nrow)
vtg_num <- dplyr::rename(vtg_num, Number=V1)
vtg_num <- merge(virus2abbre, vtg_num, by.x = "Virus", by.y = "Virus")

vtgs2 <- vtgs[vtgs$Virus %in% virus,]
vtgs2 <- merge(vtgs2, virus2abbre,
              by.x = 'Virus', by.y = 'Virus')
write.table(vtgs2[,c(2,1,7,3,4)], '../result/virally-targeted_genes29.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)
  
write.table(virus2abbre, "../result/virus2abbre.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Number of virus targeted genes
gene_virus_num <- vtgs[order(vtgs[,'Gene'],
                             vtgs[,'Virus']),c(1,2)]
gene_virus_num <- gene_virus_num %>% distinct()

gene_virus_num2 <- aggregate(gene_virus_num[2], gene_virus_num[1],
                             FUN = function(X) paste(unique(X), collapse=","))
gene_virus_num2$Virus_num <- apply(gene_virus_num2[,2,drop=FALSE], 1, 
                                    function(x) length(strsplit(x, split = ',')[[1]]))
gene_virus_num2$Virally_targeted <- apply(subset(gene_virus_num2, select = Virus_num),
                                            1, function(x) ifelse(x>quantile(gene_virus_num2[,3], 0.99)[[1]], 
                                                                  'Pan', 'Specific'))
write.table(gene_virus_num2, 
            '../intermediate/viralNum_target_a_gene.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

### Number of host factors for each family
gene_family_num <- vtgs %>%
  group_by(Baltimore, Family) %>%
  dplyr::summarise(Gene_num = n())

p_s3a <- ggplot(gene_family_num, aes(Family, Gene_num)) +
  geom_bar(stat = 'identity', aes(fill = Baltimore)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_discrete(limits = gene_family_num$Family[order(gene_family_num[,3])]) +
  scale_fill_manual(values = c(brewer.pal(9, 'Set1')[1:7])) +
  labs(y = "Gene Number") +
  coord_flip() +
  geom_text(data = gene_family_num, aes(x=Family, y=50, label = Family), 
            size=3,vjust=0.3, hjust=0) +
  mytheme +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.position = c(0.7, 0.2))

ggsave(p_s3a,
       filename = "../result/Fig_S3A.pdf",
       width = 4,
       height = 8)

