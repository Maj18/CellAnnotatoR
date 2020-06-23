RrawCountNormalization <- function(p2) {
  cm <- p2$counts # A single matrix; row: cell; col: gene
  cm %<>% as("dgCMatrix") %>% Matrix::t()
  #The normalizaiton before TF-IDF: combine CPM and part of UQ idea.

  # Total count normalization (i.e. TF-step)
  cm@x <- cm@x / rep(Matrix::colSums(cm), diff(cm@p)) #Single counts within a cell are divided by the total read counts in that cell
  cm <- Matrix::t(cm) #now row, cell, col, gene

  # Factors for min-max gene normalization
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95) %>% `+`(1e-10) # Robust alternative to maximum #Here we get the 95% quantile of the above-pre-normalized counts among all cells for each gene, here we will get a vector of norm_factors, each is for a single gene.

  cm@x <- cm@x / rep(max.vals, diff(cm@p)) #Divide the above-pre-normalized counts (different cells) of a gene (all in a column) by the norm_factor for that gene.
  return(cm)
}

#tfidfOnSingleDataframe <- function(){}

#tfidfOnTypeByLevel <- function(){}

#tfidfOnTypeByParentByLevel <- function(){}






#This version prefer smaller clusters
tfIdfOnType10 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$misc$rawCounts
  cm %<>% as("dgCMatrix") %>% Matrix::t()
  #The normalizaiton before TF-IDF: 1.UQ normalization

  norm.factor <- as.data.frame(as.matrix(cm)) %>% sapply(quantile, 0.99) # the 99% is for the genes within each cell.
  #we will get a vector of quantile values, each for a cell
  norm.factor <- norm.factor/median(norm.factor) # find then median of the quantile values, and divide all quantile values by this median, so we get norm_factor
  cm <- Matrix::t(cm)
  cm <- cm/norm.factor # divide the values for the cells of a gene (column) by the norm.factor. Then we get the normalized values

  # Get tf (Nij/sum within cells)
  cm %<>% Matrix::t()
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf) %>% as.matrix()

  #cm%<>% t(.)%>%as.matrix() #(normalized counts, row: cells, col: genes)
  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) tf[intersect(names(n), rownames(tf)),] %>% colMeans)) %>% lapply(function(x) as.data.frame(x))

  #The freq of cells (belonging to one particular type) that have one particular gene expressed
  gene.presence.freq.per.type <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) tf[intersect(names(n), rownames(tf)),] %>% colNon0Count)) %>% lapply(function(x) as.data.frame(x))

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T))
  len.by.level <- ann.per.type.by.level %>% lapply(function(n) unlist(n)%>%length)
  rowsum <- list()
  for (i in names(gene.presence.freq.per.type)){
    rowsum[[i]] <- rep(rowSums(gene.presence.freq.per.type[[i]])/len.by.level[[i]], ncol(gene.presence.freq.per.type[[i]]))%>%matrix(nrow=nrow(gene.presence.freq.per.type[[i]]), ncol=ncol(gene.presence.freq.per.type[[i]]), byrow=F)
  }

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-mean.GE.per.type.per.gene[[i]]*log(1+(1/(p+rowsum[[i]])))
  }

  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                              function(m) head(sort(m, decreasing = T), topgene)%>%sum)) %>% lapply(unlist)

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}












#THIS IS NO GOOD EITHER
#This one is the best function so far (theoretically), works kind of fine as well
tfIdfOnType9 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts
  cm %<>% as("dgCMatrix") %>% Matrix::t()

  # Get tf (Nij/sum within cells)
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf)

  #Normalize the counts further, so that the counts mean for each gene among all cells are at the same scale as (slightly bigger than) the total number of cells/cell types.
  cm  %>% Matrix::t()
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.75) #this is to ensure he counts sum for each gene among all cells are at the same scale as the total number of cells/cell types
  cm@x <- cm@x / rep(max.vals, diff(cm@p))
  cm@x[cm@x>1] <- 1
  cm %<>% Matrix::t()

  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  tf <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(tf)[intersect(names(n), rownames(tf)),] %>% colMeans)) %>% lapply(function(x) as.data.frame(x)) #For tf #divide tf into levels

  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(cm)[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% # For the further normalzed counts #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-mean.GE.per.type.per.gene[[i]]*log(1+ncol(rowsum[[i]])/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                                                function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}








#Here for each clusters, we get score*size, to compare between trees, we should take sum_j(score*size)
tfIdfOnType8 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts
  cm %<>% as("dgCMatrix") %>% Matrix::t()

  # Get tf (Nij/sum within cells)
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf)

  #Normalize the counts further, so that the counts mean for each gene among all cells are at the same scale as (slightly bigger than) the total number of cells/cell types.
  cm  %>% Matrix::t()
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95) %>% `+`(1e-10) #this is to ensure he counts sum for each gene among all cells are at the same scale as the total number of cells/cell types
  cm@x <- cm@x / rep(max.vals, diff(cm@p))
  cm %<>% Matrix::t() %>% as.matrix()

  annotation <- list()
  annotation[["root"]] <- rep("root", length( ann.by.level$annotation[[1]])) %>% setNames(., names(ann.by.level$annotation[[1]]))
  annotation %<>% append(.,ann.by.level$annotation)

  ann.by.parents <- list()
  for (i in 1:((names(annotation) %>% length)-1)){
    ann.by.parents[[names(annotation[i])]] <- split(annotation[[i+1]], annotation[[i]])
  }

  #ann.by.parents %<>% lapply(function(m) lapply(m, function(n) split(n,n)))

  GE.sum.per.type.per.gene.by.parent <- lapply(ann.by.parents, function(m) lapply(m, function(n) split(n,n)%>%lapply(function(x) cm[intersect(names(x), rownames(cm)),] %>% colMeans))) %>% lapply(function(x) lapply(x, function(k) as.data.frame(k))) #*

  len.per.type.by.parent.by.level <- lapply(ann.by.parents, function(m) lapply(m, function(n) split(n,n)%>%lapply(function(x) cm[intersect(names(x), rownames(cm)),] %>% nrow))) %>%
    lapply(function(x) lapply(x, function(k) as.data.frame(k)))

  colsum <- list()
  for (i in names(GE.sum.per.type.per.gene.by.parent)){
    for (j in names(GE.sum.per.type.per.gene.by.parent[[i]])){
      cs <- colSums(GE.sum.per.type.per.gene.by.parent[[i]][[j]])#/len.per.type.by.parent.by.level[[i]][[j]]
      cs <- rep(cs, nrow(GE.sum.per.type.per.gene.by.parent[[i]][[j]]))%>%matrix(nrow=nrow(GE.sum.per.type.per.gene.by.parent[[i]][[j]]), byrow=T)%>%data.frame%>% `rownames<-`(rownames(GE.sum.per.type.per.gene.by.parent[[i]][[j]]))
      colsum[[i]][[j]] <- cs
    }
  }

  rowsum <- list()
  for (i in names(GE.sum.per.type.per.gene.by.parent)){
    for (j in names(GE.sum.per.type.per.gene.by.parent[[i]])){
      cs <- rowSums(GE.sum.per.type.per.gene.by.parent[[i]][[j]])#/sum(len.per.type.by.parent.by.level[[i]][[j]])
      cs <- rep(cs, ncol(GE.sum.per.type.per.gene.by.parent[[i]][[j]]))%>%matrix(.,nrow=nrow(GE.sum.per.type.per.gene.by.parent[[i]][[j]]), byrow=F)%>%data.frame%>% `colnames<-`(colnames(GE.sum.per.type.per.gene.by.parent[[i]][[j]]))
      rowsum[[i]][[j]] <- cs
    }
  }

  ann.by.parents %<>% lapply(function(m) lapply(m, function(n) split(n,n)))
  tf <- as.matrix(tf)
  tf <- lapply(ann.by.parents, function(m) lapply(m, function(n) lapply(n, function(h) tf[intersect(names(h), rownames(tf)),] %>% colMeans))) %>% lapply(function(x) lapply(x, function(y) as.data.frame(y))) #For tf #divide tf into levels

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(GE.sum.per.type.per.gene.by.parent)) { #i: hierarchical level
    for (j in names(GE.sum.per.type.per.gene.by.parent[[i]])){
      tfidf.per.type.per.gene[[i]][[j]]<-tf[[i]][[j]]*log(1+ncol(rowsum[[i]][[j]])/(p+rowsum[[i]][[j]]))
    }
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  tfidf.type.size <- list()
  for (i in names(tfidf.per.type.per.gene)){
    for (j in names(tfidf.per.type.per.gene[[i]])){
      tfidf.type.size[[i]][[j]] <- t(t(tfidf.per.type.per.gene[[i]][[j]]%>%as.data.frame)*(len.per.type.by.parent.by.level[[i]][[j]]%>%as.matrix%>%as.vector)) %>% as.data.frame()
    }
  }

  top.ifidf.sum.per.type.by.level <- tfidf.type.size %>% lapply(function(n) lapply(n,function(m) lapply(m, function(k) head(sort(k, decreasing = T), topgene)%>%sum))) %>% lapply(function(m) lapply(m, function(n) as.data.frame(n)))
  #Take the genes with top tfidf values for each cell type. #Above, the size of the cell type is also taken into consideration. mean*size for each gene and each type
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(m) lapply(m, function(n) topgenesName(n, topgene)))
  #here, the size of cell type is not considered

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}









#In this version, I started from the raw count data, normalize it first to exclude variation in library size
tfIdfOnType7 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$misc$rawCounts
  cm %<>% as("dgCMatrix") %>% Matrix::t()
  #The normalizaiton before TF-IDF: 1.UQ normalization

  norm.factor <- as.data.frame(as.matrix(cm)) %>% sapply(quantile, 0.99) # the 99% is for the genes within each cell.
  #we will get a vector of quantile values, each for a cell
  norm.factor <- norm.factor/median(norm.factor) # find then median of the quantile values, and divide all quantile values by this median, so we get norm_factor
  cm <- Matrix::t(cm)
  cm <- cm/norm.factor # divide the values for the cells of a gene (column) by the norm.factor. Then we get the normalized values

  #Or we can exclude 0 prior to calculating the 75% quantile (There are too many 0 for single cell transcriptomic data)

  # Get tf (Nij/sum within cells)
  cm %<>% Matrix::t()
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf)

  #Normalize the counts further, so that the count sum reflects the spread of the data from the 95% quantile (the smaller the sum is, the more spread the data is)
  #while the sum never exceeds the total number of cells/cell types.
  cm  %>% Matrix::t()
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95) %>% `+`(1e-10)
  cm@x <- cm@x / rep(max.vals, diff(cm@p))
  cm %<>% Matrix::t()

  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  tf <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(tf)[intersect(names(n), rownames(tf)),] %>% colMeans)) %>% lapply(function(x) as.data.frame(x)) #For tf #divide tf into levels

  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(cm)[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% # For the further normalzed counts #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-tf[[i]]*log(1+ncol(rowsum[[i]])/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                                function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}








#This version also prefer bigger clusters, THIS IS NOT TRUE, IF WE TAKE CLUSTER SIZE INTO CONSIDERATION
#This one is the best function so far (theoretically), works kind of fine as well
tfIdfOnType6 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts
  cm %<>% as("dgCMatrix") %>% Matrix::t()

  # Get tf (Nij/sum within cells)
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf)

  #Normalize the counts further, so that the counts mean for each gene among all cells are at the same scale as (slightly bigger than) the total number of cells/cell types.
  cm  %>% Matrix::t()
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95) %>% `+`(1e-10) #this is to ensure he counts sum for each gene among all cells are at the same scale as the total number of cells/cell types
  cm@x <- cm@x / rep(max.vals, diff(cm@p))
  cm %<>% Matrix::t()
  #ann.per.type.by.level <- ann.by.level@annotation %>% lapply(function(n) split(n,n)) 
  ann.per.type.by.level <- ann.by.level %>% lapply(function(n) split(n,n)) %>% lapply(function(m) m[lapply(m, length)>1]) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  tf <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(tf)[intersect(names(n), rownames(tf)),] %>% colMeans)) %>% lapply(function(x) as.data.frame(x)) #For tf #divide tf into levels

  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(cm)[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% # For the further normalzed counts #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-tf[[i]]*log(1+ncol(rowsum[[i]])/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  
  top.tfidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,                                                 function(m) head(sort(m, decreasing = T), topgene)%>%sum) %>% unlist)
  len.per.type.by.level <- ann.per.type.by.level %>% lapply(function(m) sapply(m, length))
  
  top.tfidf.sum.per.type.by.level.size <- list()
  for (i in names(top.tfidf.sum.per.type.by.level)) {
    ji <- top.tfidf.sum.per.type.by.level[[i]]* len.per.type.by.level[[i]]
    top.tfidf.sum.per.type.by.level.size[[i]] <- ji
  }

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(top.tfidf.sum.per.type.by.level.size=top.tfidf.sum.per.type.by.level.size, s.pos.v3=s.pos.v3, top.tfidf.sum.per.type.by.level=top.tfidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}









#' @description This version, the markers are aimed for separating cell types between all cell types at the same hierarchical level, maybe next version, we can try to aim for separating cell types among subtypes under the same cell type.!!!!!
#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @param topgene For each cell, the genes have been sort in a descending order regarding the TF-IDF values, topgene here decides how many top genes will be collected for calculating top5ifidf.sum.per.type.by.level.
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType5 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$misc$rawCounts
  cm %<>% as("dgCMatrix") %>% Matrix::t()
  #The normalizaiton before TF-IDF: 1.UQ normalization

  norm.factor <- as.data.frame(as.matrix(cm)) %>% sapply(quantile, 0.99)
  norm.factor <- norm.factor/median(norm.factor)
  cm <- Matrix::t(cm)
  cm <- cm/norm.factor

  # Get tf (Nij/sum within cells)
  cm %<>% Matrix::t()
  tf <- cm
  tf@x <- tf@x / rep(Matrix::colSums(tf), diff(tf@p))
  tf <- Matrix::t(tf)

  #Normalize the counts further, so that the counts mean for each gene among all cells are close to one (a little bigger)
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95)
  cm@x <- cm@x / rep(max.vals, diff(cm@p))
  cm %<>% Matrix::t()

  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  tf <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(tf)[intersect(names(n), rownames(tf)),] %>% colMeans)) %>% lapply(function(x) as.data.frame(x))

  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) as.matrix(cm)[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  #The freq of cells (belonging to one particular type) that have one particular gene expressed
  #colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T)) #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene)
  #need to add column names here!!!!!!
  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-tf[[i]]*log(ncol(rowsum[[i]])/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                                                                           function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}


















#' @description This version, the markers are aimed for separating cell types between all cell types at the same hierarchical level, maybe next version, we can try to aim for separating cell types among subtypes under the same cell type.!!!!!
#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @param topgene For each cell, the genes have been sort in a descending order regarding the TF-IDF values, topgene here decides how many top genes will be collected for calculating top5ifidf.sum.per.type.by.level.
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType4 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$misc$rawCounts %<>% as("dgCMatrix") %>% Matrix::drop0()
  cm@x <- cm@x/rep(Matrix::colSums(cm), diff(cm@p)) #Here the raw counts are normalized by being devided by the total counts of a gene of all cells. Aren't these values too small? Yes, they are!
  cm%<>% as.matrix() #(normalized counts, row: cells, col: genes)
  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  #The freq of cells (belonging to one particular type) that have one particular gene expressed
  gene.presence.freq.per.type <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colNon0Freq)) %>% lapply(function(x) as.data.frame(x))

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T)) #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene)
  #need to add column names here!!!!!!
  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-(mean.GE.per.type.per.gene[[i]]/colsum[[i]])*log((max(rowsum[[i]]))/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                                                                           function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}






#******GOOD*****
#' @description This version, the markers are aimed for separating cell types between all cell types at the same hierarchical level, maybe next version, we can try to aim for separating cell types among subtypes under the same cell type.!!!!!
#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @param topgene For each cell, the genes have been sort in a descending order regarding the TF-IDF values, topgene here decides how many top genes will be collected for calculating top5ifidf.sum.per.type.by.level.
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType3 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$misc$rawCounts
  cm %<>% as("dgCMatrix") %>% Matrix::t()
  #The normalizaiton before TF-IDF: combine CPM and part of UQ idea.

  # Total count normalization (i.e. TF-step)
  cm@x <- cm@x / rep(Matrix::colSums(cm), diff(cm@p))
  cm <- Matrix::t(cm) #now row, genes, col, cells

  # Factors for min-max gene normalization
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))]
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, 0.95) %>% `+`(1e-10) # Robust alternative to maximum

  cm@x <- cm@x / rep(max.vals, diff(cm@p)) # fast way to do column (row?)-wise normalization for sparse matrices
  cm %<>% as.matrix()
  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cell types?

  #The freq of cells (belonging to one particular type) that have one particular gene expressed
  #*gene.presence.freq.per.type <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colNon0Freq)) %>% lapply(function(x) as.data.frame(x))

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  #colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T)) #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene)
  #need to add column names here!!!!!!
  rowsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))
  #rowsum <- gene.presence.freq.per.type %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  p <- 1e-10 #pseudonumber
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-mean.GE.per.type.per.gene[[i]]*log((max(rowsum[[i]]))/(p+rowsum[[i]]))
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.#rowsum[[i]] is a dataframe.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                                                                           function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}













#' @description This version, the markers are aimed for separating cell types between all cell types at the same hierarchical level, maybe next version, we can try to aim for separating cell types among subtypes under the same cell type.!!!!!
#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @param topgene For each cell, the genes have been sort in a descending order regarding the TF-IDF values, topgene here decides how many top genes will be collected for calculating top5ifidf.sum.per.type.by.level.
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType2 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #library(pivot)
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts %>% as.matrix() #(normalized counts, row: cells, col: genes)
  annotation <- list()
  annotation[["root"]] <- rep("root", length( ann.by.level$annotation[[1]])) %>% setNames(., names(ann.by.level$annotation[[1]]))
  annotation %<>% append(.,ann.by.level$annotation) #NO, this is not alright

  ann.by.parents <- list()
  for (i in 1:((names(annotation) %>% length)-1)){
    ann.by.parents[[names(annotation[i])]] <- split(annotation[[i+1]], annotation[[i]])
  }

  mean.GE.per.type.per.gene <- lapply(ann.by.parents, function(m) lapply(m, function(n) split(n,n)%>%lapply(function(x) cm[intersect(names(x), rownames(cm)),] %>% colMeans))) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) lapply(x, function(k) as.data.frame(k))) #*

  gene.presence.freq.per.type.by.parent <- lapply(ann.by.parents, function(m) lapply(m, function(n) split(n,n)%>%lapply(function(x) cm[intersect(names(x), rownames(cm)),] %>% colNon0Freq))) %>% lapply(function(x) lapply(x, function(k) as.data.frame(k))) #*

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  colsum <- mean.GE.per.type.per.gene %>% lapply(function(m) lapply(m, function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), byrow=T)%>%data.frame%>% `rownames<-`(rownames(n))))  #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene) #*
  #need to add column names here!!!!!!
  colsum %<>% lapply(function(m) lapply(m, function(n) {n[n==0]<-1e-10; n}))
  rowsum <- gene.presence.freq.per.type.by.parent %>% lapply(function(m) lapply(m,function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), byrow=F)%>%data.frame%>%setNames(colnames(n))))
  rowsum %<>% lapply(function(m) lapply(m, function(n) {n[n==0]<-1e-10; n}))  #replace all 0 with 1e-10

  tfidf.per.type.per.gene <- list() #*
  for (i in 1:(names(mean.GE.per.type.per.gene)%>%length)) { #i: hierarchical level
    l <- length(names(mean.GE.per.type.per.gene[[i]]))
    for (j in 1:l){
      tfidf.per.type.per.gene[[names(mean.GE.per.type.per.gene[i])]][[names(mean.GE.per.type.per.gene[[i]][j])]] <- (mean.GE.per.type.per.gene[[i]][[j]]/colsum[[i]][[j]])*log(1+ncol(rowsum[[i]][[j]])/rowsum[[i]][[j]])
    }
  }#tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes. #Here, the point to have 0.1 is avoid the second part becomes 0   #*

  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up. #?????
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                              function(m) sapply(m, function(k) head(sort(k, decreasing = T), topgene)%>%sum)))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) lapply(n, function(l) topgenesName(l, topgene)))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}









#' @description This version, the markers are aimed for separating cell types between all cell types at the same hierarchical level, maybe next version, we can try to aim for separating cell types among subtypes under the same cell type.!!!!!
#' @description Calculate S(pos) (corresponding to "scores.raw" in score.info) on the basis of TF-IDF values on cell types (instead of on cells)
#' @describeIn The TF-IDF function here has been adjusted to fit the dataset.
#' @param p2 the counts object within data variable p2 (p2$counts) is what need here. counts is a sparse matrix of normalized gene expresson, col: cells, row: genes.
#' @param ann.by.level Output from re-annotation (assignCellsByScores()); a list of 3 objects: "annotation" (this is what we need here), "scores", and "annoation.filt"; annotation is list of hierarchial levels (e.g. l1, l2, l3); at each level, there are annotation (cell types at this hierarchical level are assigned to cells) on all the studied cells.
#' @param clf.data Output from getClassificationData(); a list of 3 objects: "cm", "classification.tree", "gene.table" and "marker.list" (this is what we need here in this function); marker.list is a list of cell types, with each type containing 3 objects: "expressed", "not_expressed" and "parent", all three objects each has a vector of markers
#' @param topgene For each cell, the genes have been sort in a descending order regarding the TF-IDF values, topgene here decides how many top genes will be collected for calculating top5ifidf.sum.per.type.by.level.
#' @export s.pos.v3 a list of hierarchial levels, each level is a vector of per-cell total gene expression levels.
tfIdfOnType0 <- function(p2, ann.by.level, clf.data=NULL, topgene=10){
  #Mean gene expression per cell type for each gene:
  cm <- p2$counts %>% as.matrix() #(normalized counts, row: cells, col: genes)
  ann.per.type.by.level <- ann.by.level$annotation %>% lapply(function(n) split(n,n)) # a list of hierarchial levels (e.g. l1, l2, l3; each level has its own cell types and each type has all the cells [names] that belong to it)
  mean.GE.per.type.per.gene <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colMeans)) %>% #GE: gene expression # m: hierarchical level, n: cell type. Here we get the mean gene expression levels of the cells that belong to each cell type for each gene.
    lapply(function(x) as.data.frame(x)) # list of 3 levels, each level has a dataframe, row: genes, col: cells

  #The freq of cells (belonging to one particular type) that have one particular gene expressed
  gene.presence.freq.per.type <- lapply(ann.per.type.by.level, function(m) lapply(m, function(n) cm[intersect(names(n), rownames(cm)),] %>% colNon0Freq)) %>% lapply(function(x) as.data.frame(x))

  ##tf-idf transformation on cell types (the value is mean GE for each type)
  colsum <- mean.GE.per.type.per.gene %>% lapply(function(n) rep(colSums(n), nrow(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=T)) #GE: gene expression; colsum has the same data structure as mean.GE.per.type.per.gene. each matrix object within it has the same nrow and ncol as the corresponding object in mean.GE.per.type.per.gene; all rows have the same values (col sums of the corresponding object in mean.GE.per.type.per.gene)
  #need to add column names here!!!!!!
  rowsum <- gene.presence.freq.per.type %>% lapply(function(n) rep(rowSums(n), ncol(n))%>%matrix(nrow=nrow(n), ncol=ncol(n), byrow=F))

  tfidf.per.type.per.gene <- list()
  for (i in names(mean.GE.per.type.per.gene)) { #i: hierarchical level
    tfidf.per.type.per.gene[[i]]<-(mean.GE.per.type.per.gene[[i]]/colsum[[i]])*log(ncol(colsum[[i]])/rowsum[[i]])
  } #tfidf.per.type.per.gene: a list of hierarchial levels, each level has a dataframe: col, cell types; row, genes.
  if (!is.null(clf.data)) {
    markers.per.type <- clf.data$marker.list %>% lapply(function(n) n$"expressed") #here, we only take positive markers!!!!! markers.per.type: a list of cell types, each type has positive markers.
    names(markers.per.type)[names(markers.per.type) %in% c("0", "1")] <- c("X0", "X1") #change the cell type names, so they are the same as those in tfidf.per.type.per.gene
    names(markers.per.type) <- gsub(" ", ".", names(markers.per.type))

    s.pos.v3 <-
      sapply(tfidf.per.type.per.gene, function(x) #x (a hierarchical level of tfidf.per.type.per.gene): a dataframe, col, cell types; row, genes
        mapply(function(m, n)
          x[m, n], markers.per.type[intersect(colnames(x), names(markers.per.type))], intersect(colnames(x), names(markers.per.type)))) %>% # so far the data: a list of hierarchial levels, each levels contains a dataset: col, genes; row, cells.
      sapply(function(n)
        sapply(n, sum)) # here we calcualte col sum: sum of the expression of all genes within each cell.
  }else{
    s.pos.v3 <- NULL
  }

  #sort the tfidf values for each cell type, and take the top topgene values and sum it up.
  top.ifidf.sum.per.type.by.level <- tfidf.per.type.per.gene %>% lapply(function(n) lapply(n,
                                            function(m) head(sort(m, decreasing = T), topgene)%>%sum))

  #Take the genes with top tfidf values for each cell type.
  top.genes.per.type.by.level <- lapply(tfidf.per.type.per.gene, function(n) topgenesName(n, topgene))

  return(list(s.pos.v3=s.pos.v3, top.ifidf.sum.per.type.by.level=top.ifidf.sum.per.type.by.level, top.genes.per.type.by.level=top.genes.per.type.by.level)) # a list of hierarchial levels, each level is a vector of per-cell GE sums.
}


#' Calculate the freq of non-0 number for each column of a dataframe, row: cells belonging to one type; col: gene
colNon0Freq <- function(df){
  col.non.0.freq<-c()
  for (i in 1:ncol(df)){
    f<-sum(df[,i]!=0)/length(df[,i])
    col.non.0.freq <- c(col.non.0.freq, f)
  }
  return(col.non.0.freq)
}

colNon0Count <- function(df){
  col.non.0.count<-c()
  for (i in 1:ncol(df)){
    f<-sum(df[,i]!=0)
    col.non.0.count <- c(col.non.0.count, f)
  }
  return(col.non.0.count)
}

#' Grab the topgens regarding their tfidf values
topgenesName <- function(df, topgene){
  topgenes <- list()
  for (i in 1:ncol(df)){
    column <- df[, i, drop=F]
    topgenes.name <- row.names(column)[order(column[1], decreasing=TRUE)][1:topgene]
    topgenes.tfidf <- column[order(column[1], decreasing=TRUE),][1:topgene] %>% setNames(.,topgenes.name)
    topgenes[colnames(column)]= list(topgenes.tfidf)
  }
  return(topgenes)
}

