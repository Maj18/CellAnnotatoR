#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType <- function(p2, ann.by.level, clf.data){
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts %>% as.matrix() #(normalized counts, row: cells, col: genes)
  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cells

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T)) #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene)
  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  max.total.GE.per.cell <- sapply(colsum, function(n) nrow(n)/n[1,]) %>% unlist() %>% max() #The maximum of the per-cell total gene expression levels.
  tfidf.per.type.per.gene <- list()
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-(mean.GE.per.type.per.gene[[i]]/rowsum[[i]])+log(ncol(colsum[[i]])/max.total.GE.per.cell*colsum[[i]]) #The point to have max.total.GE.per.cel here is that the gene expression data here is normalized (max may be only 1.7), not read count, so a bit deviate from the orignal tf-idf idea (where here it is the total number of documents that contain a particular word).
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.

  markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
  names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
  names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

  s.pos.v3 <-
    sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
      mapply(function(m, n)
        x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
    sapply(function(n)
      sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.

  return(s.pos.v3) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}
