######################################################
# Immune-related function over-representation analysis

library(clusterProfiler)
node_info <- read.table('../intermediate/nodes_info.txt', 
                        sep = '\t', header = TRUE)
module_list <- read.table('../intermediate/batch_MCL_out/out.HUMPPI2022.I228', 
                          header = FALSE, sep = ' ')
module_list$Module <- rownames(module_list)
module_list <- dplyr::rename(module_list, Gene_list=V1)
module_list <- module_list[1:1225,]

immunity <- read.table('../data/IRGs/Immune-related_process.txt',
                       header = FALSE, sep = '\t', quote = "")
colnames(immunity) <- c('Gene', 'Description')
immunity2 <- immunity[immunity$Gene %in% node_info[node_info$Degree!="-",1],]
immunity2 <- immunity2[,c(2,1)]

module_imm <- data.frame(Description = character(), GeneRatio = character(), BgRatio = character(), 
                         pvalue = character(), p.adjust = character(), qvalue = character(), 
                         geneID = character(), Module = character())

jaccard <- function(x,y){
  jaccard_index <- length(intersect(x, y))/length(union(x, y))
  jaccard_index
}

for (i in 1:nrow(module_list)) {
  gene_list = unlist(strsplit(as.character(module_list[i,1]), split = '\t', fixed = T))
  if (length(intersect(gene_list, unique(immunity2$Gene))) == 0) next
  imm <- enricher(gene_list,
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 universe=as.character(node_info$Gene),
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 maxGSSize = 50000,
                 TERM2GENE= immunity2)
  imm2 <- imm@result
  imm2 <- imm2[imm2$pvalue < 0.01 & imm2$p.adjust < 0.01,-1]
  imm2 <- imm2[,1:7]
  imm2$Module <- rep(i, nrow(imm2))
  
  if(nrow(imm2)>0){
    reduRow <- c()
    for(n in 1:nrow(imm2)){
      for(m in 1:nrow(imm2)){
        if(n != m){
          list1 = strsplit(imm2$geneID[n], split = '/', fixed = T)[[1]]
          list2 = strsplit(imm2$geneID[m], split = '/', fixed = T)[[1]]
          ji = jaccard(list1, list2)
          if(ji >= 0.99){
            if(imm2$p.adjust[n] > imm2$p.adjust[m]){
              reduRow <- c(reduRow, n)
            } else {
              reduRow <- c(reduRow, m)
            }
          }
        }
      }
    }
    
    if(length(reduRow)>0){
      imm3 <- imm2[-unique(reduRow),]
      module_imm <- rbind(module_imm, imm3)
    }else{
      module_imm <- rbind(module_imm, imm2)
    }
  }
  print(paste("Module", i, "finished!"))
}

module_imm <- module_imm %>% 
  mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")"))
write.table(module_imm[,c(8,1:7)], "../intermediate/module_IR_process.txt",
            quote = FALSE, sep = '\t', row.names = FALSE)

