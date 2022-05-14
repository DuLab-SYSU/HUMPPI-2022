### Coronaviridae-targeted genes and functions 
## Categories based on the gene overlap
gene_virus_num <- read.table("../intermediate/viralNum_target_a_gene.txt",
                             header = TRUE, sep = "\t", quote = "")
node_info <- read.table('../intermediate/nodes_info.txt', 
                        header = TRUE, sep = '\t')

vtg_info <- merge(gene_virus_num, node_info, by = 'Gene')


vtgs <- read.table("../result/virally-targeted_genes29.txt", 
                   header = TRUE, sep = '\t', quote = "")

vtg_cov <- vtgs %>%
  dplyr::select(Gene, Abbre) %>%
  mutate(Number = 1) %>%
  pivot_wider(names_from = Abbre, values_from = Number)
vtg_cov <- vtg_cov %>% 
  dplyr::select(colnames(vtg_cov)[c(1, 26:28, 2:25, 29, 30)])

vtg_cov$Other_num <- apply(vtg_cov[,5:30], 1, 
                           function(x) sum(x, na.rm = TRUE))

vtg_cov2 <- filter(vtg_info[,c(1,4)], Virally_targeted == "Pan") %>%
  mutate(Virally_targeted = 1) %>%
  dplyr::rename(Pan = Virally_targeted) %>%
  full_join(vtg_cov[,c(1:4,31)], by = "Gene")

vtg_cov2[is.na(vtg_cov2)] <- 0
vtg_cov2$CoV_num <- apply(vtg_cov2[,3:5], 1, 
                          function(x) sum(x, na.rm = TRUE))

vtg_cov2 <- rbind(filter(vtg_cov2, Pan == 1),
                  filter(vtg_cov2, Pan == 0 & CoV_num != 0))


Labels <- c()
for (i in 1:nrow(vtg_cov2)) {
  if (vtg_cov2[i, 2] == 1){
    if (sum(vtg_cov2[i, 3:5]) > 0) {
      Labels <- c(Labels, "Pan")
    } else {
      Labels <- c(Labels, "Pan_Not_CoV")
    }
  } else if (vtg_cov2[i, 6] == 0) {
    if (vtg_cov2[i, 3] == 1) {
      if (sum(vtg_cov2[i, 4:5]) == 2) {
        Labels <- c(Labels, "CoV-pan")
      } else if (vtg_cov2[i, 4] == 1) {
        Labels <- c(Labels, "SARS-CoV-1 & MERS-CoV")
      } else if (vtg_cov2[i, 5] == 1) {
        Labels <- c(Labels, "SARS-CoV-2 & MERS-CoV")
      } else {
        Labels <- c(Labels, "MERS-CoV")
      }
    } else if (vtg_cov2[i, 4] == 1) {
      if (vtg_cov2[i, 5] == 1) {
        Labels <- c(Labels, "SARS-CoV-2 & SARS-CoV-1")
      } else {
        Labels <- c(Labels, "SARS-CoV-1")
      }
    } else {
      Labels <- c(Labels, "SARS-CoV-2")
    }
  } else {
    Labels <- c(Labels, "CoV & Others")
  }
}

vtg_cov2$Label <- Labels
write.table(vtg_cov2, "../intermediate/CTG_category.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

vtg_cov2$Label <- factor(vtg_cov2$Label,
                         levels = rev(c("Pan", "Pan_Not_CoV", 
                                        "CoV & Others", "CoV-pan",
                                        "SARS-CoV-2 & SARS-CoV-1",
                                        "SARS-CoV-2 & MERS-CoV",
                                        "SARS-CoV-1 & MERS-CoV",
                                        "SARS-CoV-2", "SARS-CoV-1", "MERS-CoV")))

cov_num <- vtg_cov2 %>%
  group_by(Label) %>%
  dplyr::summarise(Number = n())

## Functional enrichment
bp <- read.table('../data/MSigDB7.4/c5.bp.v7.4.symbols_GO.ID.txt',
                 header = FALSE, sep = '\t', quote = "")

colnames(bp) <- c("GO_ID", 'Description', 'Gene')
bp2 <- bp[bp$Gene %in% node_info$Gene,]

library(clusterProfiler)
cov_go <- data.frame(Category = character(), Description = character(), GeneRatio = character(), 
                     BgRatio = character(), pvalue = character(), p.adjust = character(), 
                     qvalue = character(), geneID = character())

for (i in unique(vtg_cov2$Label)){
  go <- enricher(vtg_cov2[vtg_cov2$Label == i,1],
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 universe=as.character(node_info$Gene),
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 maxGSSize = 50000,
                 TERM2GENE= bp2[,2:3])
  
  go2 <- go@result
  go2 <- go2[go2$pvalue < 0.01 & go2$p.adjust < 0.01,-1]
  go2 <- go2[,1:7]
  go2$Cluster <- rep(i, nrow(go2))
  go2 <- go2[,c(8,1:7)]
  
  if(nrow(go2)>0){
    reduRow <- c()
    for(n in 1:nrow(go2)){
      for(m in 1:nrow(go2)){
        if(n != m){
          list1 = strsplit(go2$geneID[n], split = '/', fixed = T)[[1]]
          list2 = strsplit(go2$geneID[m], split = '/', fixed = T)[[1]]
          ji = jaccard(list1, list2)
          if(ji >= 0.99){
            if(go2$p.adjust[n] > go2$p.adjust[m]){
              reduRow <- c(reduRow, n)
            } else {
              reduRow <- c(reduRow, m)
            }
          }
        }
      }
    }
    
    if(length(reduRow)>0){
      go3 <- go2[-unique(reduRow),]
      cov_go <- rbind(cov_go, go3)
    }else{
      cov_go <- rbind(cov_go, go2)
    }
  }
  print(i)
}

cov_go <- cov_go %>% 
  mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")"))
write.table(cov_go[,c(1:4,6)], "../result/CTG_category_GO_bp.txt",
            sep = "\t", 
            row.names = FALSE, quote = FALSE)

## Network topological properties
cov_net <- vtg_cov2 %>%
  inner_join(node_centrality, by = "Gene")

cov_num_net <- data.frame(table(cov_net$Label))
colnames(cov_num_net) <- c("Label", "Number")

p_cov_num <- ggplot(cov_num_net, aes(Label, Number)) +
  geom_bar(stat = 'identity', width = .7, fill = "steelblue") +
  scale_y_log10() +
  geom_text(aes(x = Label, y = Number, label = Number, angle = 270), color = "white") +
  labs(x = "", y = "Gene Number") +
  coord_flip() +
  mytheme

cov_random_centrality <- data.frame(Number = numeric(), Label = character(),
                                    Degree = numeric(), Eigenvector_centrality = numeric())
for (v in cov_num_net$Label) {
  random_node <- data.frame(Number = 1:10000, Label = v)
  for (attr in colnames(node_centrality)[2:3]) {
    var <- assign(attr, c())
    
    for (i in 1:10000) {
      set.seed(i)
      node_samp <- node_centrality[sample(1:nrow(node_centrality),
                                          size = cov_num_net[cov_num_net$Label == v, 2], replace = FALSE),]
      var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
    }
    
    random_node[attr] <- var
  }
  
  cov_random_centrality <- rbind(cov_random_centrality, random_node)
}

cov_random_centrality$Label <- factor(cov_random_centrality$Label,
                                      levels = rev(c("Pan", "Pan_Not_CoV", 
                                                     "CoV & Others", "CoV-pan",
                                                     "SARS-CoV-2 & SARS-CoV-1",
                                                     "SARS-CoV-2 & MERS-CoV",
                                                     "SARS-CoV-1 & MERS-CoV",
                                                     "SARS-CoV-2", "SARS-CoV-1", "MERS-CoV")))

cov_random_centrality2 <- cov_random_centrality %>%
  group_by(Label) %>%
  dplyr::summarise(Degree_mean = mean(Degree),
                   Degree_sd = sd(Degree),
                   Eigenvector_centrality_mean = mean(Eigenvector_centrality),
                   Eigenvector_centrality_sd = sd(Eigenvector_centrality))

cov_centrality <- cov_net %>%
  group_by(Label) %>%
  dplyr::summarise(Degree = mean(Degree),
                   Eigenvector_centrality = mean(Eigenvector_centrality)) %>%
  inner_join(cov_random_centrality2, by = "Label")

cov_centrality <- cov_centrality %>%
  mutate(Degree_Z = (Degree - Degree_mean)/Degree_sd,
         Degree_P = 2*pnorm(-abs(Degree_Z)),
         Eigenvector_centrality_Z = (Eigenvector_centrality - 
                                       Eigenvector_centrality_mean)/Eigenvector_centrality_sd,
         Eigenvector_centrality_P = 2*pnorm(-abs(Eigenvector_centrality_Z)))

signif_fun <- function(x){
  s <- c()
  for (i in x) {
    ifelse(i < 0.001, sn <- "***", 
           ifelse(i < 0.01, sn <- "**", 
                  ifelse(i < 0.05, sn <- "*", 
                         sn <- "ns")))
    s <- c(s, sn)
  }
  s
}

cov_centrality <- cov_centrality %>%
  mutate(Degree_Signif = signif_fun(Degree_P),
         Eigenvector_centrality_Signif = signif_fun(Eigenvector_centrality_P))

p_cov_d <- ggplot(cov_random_centrality, aes(Label, Degree)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1) +
  geom_point(data = cov_centrality, aes(Label, Degree), size = 2, 
             color = "red") +
  geom_text(data = cov_centrality, aes(x = Label, y = Degree + 0.5, label = Degree_Signif)) + 
  labs(x="",y="Degree") +
  coord_flip() +
  mytheme

p_cov_e <- ggplot(cov_random_centrality, aes(Label, 1000*Eigenvector_centrality)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1) +
  geom_point(data = cov_centrality, aes(Label, 1000*Eigenvector_centrality), size = 2, 
             color = "red") +
  geom_text(data = cov_centrality, aes(x = Label, y = 1000*Eigenvector_centrality + 0.2, 
                                       label = Eigenvector_centrality_Signif)) + 
  labs(x="",y="Eigenvector centrality (10-3)") +
  coord_flip() +
  mytheme

### Functional domains (SAFE)
safe_anno <- read.table("../intermediate/SAFE/human_gobp-attribute_properties_annotation-highest.txt", 
                        header = TRUE, sep = "\t", quote = "")
safe_anno$Domain.Id <- factor(safe_anno$Domain.Id,
                              labels = letters[1:19])

safe_anno2 <- safe_anno %>%
  dplyr::select(Attribute.name, Domain.Id)

safe_node_property <- read.table("../intermediate/SAFE/human_gobp-node_properties_annotation-highest.txt",
                                 header = TRUE, sep = "\t", quote = "")
domain_num <- safe_node_property %>%
  dplyr::select(Node.label, Domain..predominant., Neighborhood.score..max.1..min.0...predominant., 
                Total.number.of.enriched.domains) %>%
  dplyr::rename(Gene = Node.label,
                Predominant_domain = Domain..predominant.,
                Neighborhood_score = Neighborhood.score..max.1..min.0...predominant.,
                Domain_num = Total.number.of.enriched.domains)
domain_num$Predominant_domain <- factor(domain_num$Predominant_domain,
                                        labels = c(1, letters[1:19]))

safe_neighbor_score <- read.table("../intermediate/SAFE/human_gobp-neighborhood_scores_annotation-highest.txt",
                                  header = TRUE, sep = "\t", quote = "", skip = 3)
safe_neighbor_score$Function_num <- apply(safe_neighbor_score[,-1], 1, 
                                          function(x) sum(x >= 0.323))

cov_domain <- inner_join(vtg_cov2, domain_num, by = "Gene")
cov_fun <- inner_join(vtg_cov2, safe_neighbor_score[,c(1, ncol(safe_neighbor_score))],
                      by = c("Gene" = "name"))

library(ggridges)
library(viridis)
Colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
p_d <- ggplot(cov_domain, aes(Domain_num, Label, fill = ..x..)) +
  geom_density_ridges_gradient() +
  scale_x_continuous(expand = c(0.01, 0), breaks = seq(0, 10, 2)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  #scale_fill_viridis(name = "Domain_num", option = "C") +
  scale_fill_gradientn(colours=Colormap) +
  labs(x = "Domain Number") +
  theme_ridges(font_size = 13, grid = FALSE) +
  mytheme +
  theme(axis.title.y = element_blank(),
        legend.position = "none")

p_fig_3de <- ggpubr::ggarrange(p_cov_num, p_cov_d, p_cov_e, p_d, nrow=1,
                  labels = c("D", "", "", "E"))
ggsave(p_fig_3de,
       filename = "../result/Fig_3DE.pdf",
       width = 18,
       height = 7)

## Comparison between SARS-CoV-2 and other labels
compare_list <- list()
for (l in unique(as.character(cov_domain$Label))){
  if (l != "SARS-CoV-2") {
    compare_list <- c(compare_list, list(c('SARS-CoV-2', l)))
  }
}

ggplot(cov_domain, aes(Label, Domain_num)) + 
  geom_boxplot() +
  stat_compare_means(comparisons = compare_list) +
  coord_flip()

cov_fun2 <- cov_fun %>%
  dplyr::select(Gene, Label) %>%
  left_join(safe_neighbor_score[,1:ncol(safe_neighbor_score)-1], 
            by = c("Gene" = "name")) %>%
  pivot_longer(-(Gene:Label), names_to = "Function", 
               values_to = "Enrichment_score") %>%
  filter(Enrichment_score >= 0.323) %>%
  mutate(Function = stringr::str_replace_all(Function, "\\.", " "))

cov_domain2 <- safe_anno %>%
  dplyr::select(Attribute.name, Domain.Id) %>%
  dplyr::rename(Function = Attribute.name) %>%
  inner_join(cov_fun2, by = "Function") %>%
  dplyr::select(Domain.Id:Label) %>%
  distinct()

cov_domain2 <- aggregate(cov_domain2[2], cov_domain2[c(1,3)],
                         FUN = function(X) paste(unique(X), collapse=","))
cov_domain2$Gene_num <- apply(cov_domain2[,3,drop = FALSE], 1, 
                              function(x) length(stringr::str_split(x, ",")[[1]]))

cov_domain2_2 <- safe_anno %>%
  dplyr::select(Attribute.name, Domain.Id) %>%
  dplyr::rename(Function = Attribute.name) %>%
  inner_join(cov_fun2, by = "Function") %>%
  dplyr::select(Domain.Id:Label) %>%
  distinct()

cov_domain2_2 <- aggregate(cov_domain2_2[1], cov_domain2_2[c(2,3)],
                           FUN = function(X) paste(unique(X), collapse=","))
cov_domain2_2$Domain_num <- apply(cov_domain2_2[,3,drop = FALSE], 1, 
                                  function(x) length(stringr::str_split(x, ",")[[1]]))
cov_domain2_2 <- dplyr::rename(cov_domain2_2,
                               Category = Label)

write.table(cov_domain2_2[,c(2,1,3,4)],
            "../result/CoV-targets_domain_num.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)



