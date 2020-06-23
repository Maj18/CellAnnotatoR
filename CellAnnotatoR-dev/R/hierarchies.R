#' Merge Annotation to Level
#' @param annotation annotation for high-resolution. Cell type names must correspond to nodes in the `classification.tree`
mergeAnnotationToLevel <- function(level, annotation, classification.tree) {
  parent.types <- classificationTreeToDf(classification.tree) %$% Node[PathLen == level] %>% unique()
  if (length(parent.types) == 0) {
    warning("Nothing to merge at level ", level)
    return(annotation)
  }

  if (is.factor(annotation)) {
    annotation <- setNames(as.character(annotation), names(annotation))
  }

  type.map <- unique(annotation) %>% setNames(., .)
  for (pt in parent.types) {
    for (st in getAllSubtypes(pt, classification.tree)) {
      type.map[st] = pt
    }
  }

  return(setNames(type.map[annotation], names(annotation)))
}

#' Merge Annotation By Levels
#' @description merge provided annotation using `classification.tree` to get annotations for different levels of hierarchy
#'
#' @inheritParams classificationTreeToDf
#' @inheritParams mergeAnnotationToLevel
#' @return list of annotations where each element correspond to some hierarchy level
#'
#' @examples
#'   ann_by_level <- mergeAnnotationByLevels(annotation, clf_data$classification.tree)
#'
#' @export
mergeAnnotationByLevels <- function(annotation, classification.tree) {
  max.level <- classificationTreeToDf(classification.tree)$PathLen %>% max()
  anns <- 1:max.level %>% setNames(., paste0("l", .)) %>%
    lapply(mergeAnnotationToLevel, annotation, classification.tree)
  return(anns)
}

simplifyHierarchy <- function(branch) {
  if (!is.list(branch))
    return(branch)

  for (i in 1:length(branch)) {
    while (is.list(branch[[i]]) && length(branch[[i]]) == 1) {
      branch[[i]] <- branch[[i]][[1]]
    }
  }

  if (length(branch) == 1)
    return(branch[[1]])

  return(lapply(branch, simplifyHierarchy))
}

appendHierarchyBranch <- function(branch, parent.name) {
  branch %<>% simplifyHierarchy()
  is.leaf <- (sapply(branch, is.character) & (sapply(branch, length) == 1))

  current.nodes <- c(unlist(branch[is.leaf], use.names=F), names(branch[!is.leaf])) %>%
    setNames(., .) %>% lapply(function(x) list(expressed=c(), not_expressed=c(), parent=parent.name))

  sub.branches <- names(branch)[!is.leaf] %>%
    lapply(function(n) appendHierarchyBranch(branch[[n]], n)) %>% unlist(recursive=F)
  return(c(current.nodes, sub.branches))
}

#' Hierarchy to Classification Tree
#' @description Convert list with cell type hierarchy to classification tree in igraph format
#' @param hierarchy list of cell types, where each element represent a cell type. If cell type has some subtypes it's represented as another list, otherwise it's just a string
#' @inherit createClassificationTree return
#' @examples
#'   hierarchy <- list(
#'     Alveolar=c("AT1 Cell", "AT2 Cell"),
#'     B=c("B Cell", "Ig-producing B cell"),
#'     NK_T=list(`T`=c("Dividing T cells", "T Cell_Cd8b1 high", "Nuocyte"), "NK Cell"),
#'     "Endothelial"
#'   )
#'
#'   clf_tree <- hierarchyToClassificationTree(hierarchy)
#'
#' @export
hierarchyToClassificationTree <- function(hierarchy) {appendHierarchyBranch(hierarchy, "root") %>% createClassificationTree()}

splitClusteringDf <- function(df) {
  if (ncol(df) == 1)
    return(df[[1]])

  return(split(df, df[,1]) %>% lapply(function(x) splitClusteringDf(x[,2:ncol(x)]))) # 这个就是循环套循环的情况
}
#split the dataframe into list layer after layer; a hierarchical list of clusters, the clusters in the last layer has all their cell types

#' Derive Hierarchy
#' @description derive hierarchy from the data using hclust
#' @param feature.matrix matrix where rows represent cells and columns represent either genes or some embedded space (e.g. PCA)
#' @param annotation vector with cell type label per cell
#' @param dist.method method for pairwise distance estimation. Either "cor" for correlation distance or any method supported by `dist` function
#' @param max.depth maximal depth of the hierarchy
#' @return list with cell type hierarchy
#' @examples
#'   hierarchy <- deriveHierarchy(pca_mtx, annotation, max.depth=3)
#'   clf_tree <- hierarchyToClassificationTree(hierarchy)
#'   plotTypeHierarchy(clf_tree)
#'
#' @export
deriveHierarchy <- function(feature.matrix, annotation, dist.method="cor", levels=c(0.8, (1-1e-5))) {
  feature.matrix %<>% as("dgCMatrix") %>% conos:::collapseCellsByType(groups=annotation, min.cell.count=0) # take PCA components sum of each cluster of cells, 100 (100 PCA components) values for each cluster. # The clusters here function like cells
  # Some study take of cluster median for each PCA component, while here we use sum, I think sum is better
  # I don't think I will need this step, because I need the cells to calcualte S(t) scores

  if (dist.method == "cor") {
    c.dist <- (1 - cor(Matrix::t(feature.matrix))) %>% as.dist()
  } else {
    c.dist <- dist(feature.matrix, method=dist.method)
  }

  clusts <- hclust(c.dist)
  #heights <- quantile(clusts$height, seq(0, 1, length.out=max.depth + 1) %>% .[2:(length(.)-1)]) %>% rev()
  #heights <- quantile(clusts$height, c(0.5, 0.6, 0.7, 0.8, 0.9, (1-1e-5))) %>% rev() 
  heights <- quantile(clusts$height, levels) %>% rev() 
  clusts %<>% cutree(h=heights) %>% as.matrix()
  for (i in 1:ncol(clusts)) {
    clusts[,i] %<>% paste0("l", i, "_", .)
  }  # rename the clusters at each layer to e.g. l1_1, l1_2, l1_3 ...

  clusts %<>% cbind(rownames(.)) %>% `colnames<-`(paste0("l", 1:ncol(.))) %>% tibble::as_tibble() #name the columns as e.g. l1, l2, l3...

  return(clusts %>% splitClusteringDf() %>% simplifyHierarchy())
}

# modified from stats:::cuttree  (only has stats::: added)
# At each cutting line, this program assign all the cells to the available clusters at the cutting lines.
# Here the output is like the annotation at different layers.
cuttreeLi <- function (tree, k = NULL, h = NULL) 
{
  if (is.null(n1 <- nrow(tree$merge)) || n1 < 1) 
    stop("invalid 'tree' ('merge' component)")
  n <- n1 + 1
  if (is.null(k) && is.null(h)) 
    stop("either 'k' or 'h' must be specified")
  if (is.null(k)) {
    if (is.unsorted(tree$height)) 
      stop("the 'height' component of 'tree' is not sorted (increasingly)")
    k <- n + 1L - apply(outer(c(tree$height, Inf), h, ">"), 2, which.max)
    #apply(outer(c(tree$height, Inf), h, ">"), 2, which.max), here since the height is in ascending order, we get the first height that is larger than the cutting threshold provided by h.
    #k refers to, if height is in a decending order, the last height that is above the threshold given by h.
    if (getOption("verbose")) 
      message("cutree(): k(h) = ", k, domain = NA)
  } else {
    k <- as.integer(k)
    if (min(k) < 1 || max(k) > n) 
      stop(gettextf("elements of 'k' must be between 1 and %d", 
                    n), domain = NA)
  }
  ans <- .Call(stats:::C_cutree, tree$merge, k) #Get the cluster labels for each cell at the cutting (i.e. all the branches that were cut through)  # The cell order here is the same as that of clusts$labels
  if (length(k) == 1L) {
    ans <- setNames(as.vector(ans), tree$labels)
  }
  else {
    colnames(ans) <- if (!is.null(h)) 
      h
    else k
    rownames(ans) <- tree$labels  #
  } 
  return(ans)
}


#Modified from conos:::collapseCellsByType (only added conos:::)
collapseCellsByTypeLi <- function (cm, groups, min.cell.count = 10) 
{
  groups <- as.factor(groups)
  cl <- factor(groups[match(rownames(cm), names(groups))], 
               levels = levels(groups))
  tc <- conos:::colSumByFactor(cm, cl)
  tc <- tc[-1, , drop = F]
  tc[table(cl) >= min.cell.count, ]
}