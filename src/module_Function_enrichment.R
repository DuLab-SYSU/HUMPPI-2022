############################################
# Gene ontology over-representation analysis
### GO enrichment
library(clusterProfiler)
node_info <- read.table('../intermediate/nodes_info.txt', 
                        sep = '\t', header = TRUE)
module_list <- read.table('../intermediate/batch_MCL_out/out.HUMPPI2022.I228', 
                           header = FALSE, sep = ' ')
module_list$Module <- rownames(module_list)
module_list <- dplyr::rename(module_list, Gene_list=V1)
module_list <- module_list[1:1225,]

jaccard <- function(x,y){
  jaccard_index <- length(intersect(x, y))/length(union(x, y))
  jaccard_index
}

for (j in c('bp', 'cc', 'mf')){
  annotation <- read.table(paste0('../data/MSigDB7.4/c5.general_', j, '.v7.4.symbols_GO.ID.txt'),
                           header = FALSE, sep = '\t', quote = "")
  colnames(annotation) <- c("GO_ID", 'Description', 'Gene')
  annotation2 <- annotation[annotation$Gene %in% node_info[node_info$Degree!="-",1],]
  
  module_go <- data.frame(Description = character(), GeneRatio = character(), BgRatio = character(), 
                           pvalue = character(), p.adjust = character(), qvalue = character(), 
                           geneID = character(), Module = character())
  
  for (i in 1:nrow(module_list)) {
    gene_list = unlist(strsplit(as.character(module_list[i,1]), split = '\t', fixed = T))
    if (length(intersect(gene_list, unique(annotation2$Gene))) == 0) next
    go <- enricher(gene_list,
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   universe=node_info$Gene,
                   qvalueCutoff = 0.01,
                   minGSSize = 1,
                   maxGSSize = 50000,
                   TERM2GENE= annotation2[,2:3])
    go2 <- go@result
    go2 <- go2[go2$pvalue < 0.01 & go2$p.adjust < 0.01,-1]
    go2 <- go2[,1:7]
    go2$Module <- rep(i, nrow(go2))
    
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
        module_go <- rbind(module_go, go3)
      }else{
        module_go <- rbind(module_go, go2)
      }
    }
    print(paste("Module", i, "finished!"))
  }
  
  module_go <- module_go %>% 
    mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")"))
  write.table(module_go[,c(8,1:7)], paste0("../result/module_GO_", j, ".txt"),
              quote = FALSE, sep = '\t', row.names = FALSE)
  print(paste("###", j, "finished!"))
}
