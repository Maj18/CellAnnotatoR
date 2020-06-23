#' Get annotation for both parent types (higher hierarchical levels) and the lowest-level clusters
#' @param ann.by.levelM annotation at each of the 4 hierarchical levels (l0, l1, l2, and l3), each level has annotation (with cell types only from the target level) for all cells.
#' @export ann.by.all.typeM: a list, with names: root, intermediate parent cell types and the lowest-level clusters; each type may have different numbers of cells, the lower level one target type is at, the less cells it contains.  "values": annotation (assigned with the subtypes under the target type, for clusters at the lowest level, all cells there belong to the target cluster); ann.by.parentM: comparing to ann.by.all.typeM, ann.by.parentM lacks the annotations for the lowest-level clusters (the annotation there for each cell is the target type itself).
getAnnotationAllType <- function(ann.by.levelM){
  ann.by.levelM$l0 <- rep("root", length(ann.by.levelM$l1)) %>% setNames(names(ann.by.levelM$l1))
  L <- names(ann.by.levelM) %>% length()
  ann.by.parentM <- lapply(0:(L-2), function(i) # 3
    split(ann.by.levelM[[paste0("l", i+1)]], ann.by.levelM[[paste0("l", i)]])) %>% Reduce(c, .) %>% .[sapply(., length) > 5] %>%
        .[sapply(., function(x) length(unique(x)) > 1)]

  bottom.level <- paste0("l",L-1)
  ann.by.clusters <- ann.by.levelM$l0 %>% split(ann.by.levelM[[bottom.level]])
  ann.by.all.typeM <- c(ann.by.parentM, ann.by.clusters)
  return(list(ann.by.parentM=ann.by.parentM, ann.by.all.typeM=ann.by.all.typeM))
}



#' Get DE genes and pre-select markers from them: #The DE genes here are not only for parent types (intermediate nodes), but also for the lowest-level clusters.
#' @param p2
#' @param ann.all.type contains ann.by.parentM and ann.by.all.type, among which ann.by.levelM is what we will use in this function: annotation at each of the 4 hierarchical levels (l0, l1, l2, and l3), each level has annotation (with cell types only from the target level) for all cells
#' @export pre.selected.markersM: a list, for all the intermediate parent types. For each type: expressed, not_expressed
getPreMarkers <- function(p2, ann.all.type) {
  ann.by.parentM <- ann.all.type$ann.by.parentM
  de.info.per.parentM <- ann.by.parentM %>%
    pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))

  pre.selected.markersM <- de.info.per.parentM %>% lapply(function(dfs)
    list(
      positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
      negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
    )
  ) #postive markers: Z>0.001; negative markers: Z<-0.001; here are pre-selected markers

  return(pre.selected.markersM)
}



#' Select markers from pre-markers
#' @param pre.selected.markersM a list, for all the intermediate parent types. For each type: expressed, not_expressed
#' @param cm.norm the tf-idf data. row: genes; col: cells
#' @param ann.all.type contains ann.by.parentM and ann.by.all.type, among which Mann.by.parentM is what we will use in this function: comparing to ann_by_all_typeM, ann_by_parentM lacks the annotations for the lowest-level clusters.
#' @export
selectMarkerAllTypes <- function(pre.selected.markersM, cm.norm, ann.all.type, max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10) {
  ann.by.parentM <- ann.all.type$ann.by.parentM
  marker.infoM <- names(pre.selected.markersM) %>% setNames(., .) %>% #This is a time consuming step
    lapply(function(n) selectMarkersPerType(cm.norm, ann.by.parentM[[n]], pre.selected.markersM[[n]],
            parent=n, max.iters=max.iters, max.uncertainty=max.uncertainty,
              log.step=log.step, verbose=verbose, n.cores=n.cores))
  #Automatic marker selection (for all the clusters, cell subtypes and types from the prior manual annotation information_                               ann.by.parent) from the DE   genes. All the used 22 genes are DE genes, in the last step you will also find that all these 22                          genes are selected as markers.

  marker.listM <- setNames(marker.infoM, NULL) %>% unlist(recursive=F)
  marker.listM %<>% .[sapply(., function(x) length(x$expressed) >= 0)] #markers for all the cell types (including root, intermediate parent types and bottom clusters)
  return(marker.listM)
}



#' calculating uncertainty by cell type per level # this function may be not used, try scoreCellUncertaintyPerLevel after re-annotation instead
#' @inheritParams getCellTypeScoreInfo ?
#' @param marker.listM marker list for all cell types (including intermediated parent types and bottom clusters) names(marker_listM)
#' @param cm.norm the tf-idf data. row: genes; col: cells
#' @param ann.all.type contains ann.by.parentM and ann.by.all.type, among which ann.by.parentM is what we will use in this function: comparing to ann_by_all_typeM, ann_by_parentM lacks the annotations for the lowest-level clusters.
#' @export mean.unc.by.types a vector, with names: parent cell types; value: uncertainty
meanUncByTypes <- function(marker.listM, cm.norm, ann.all.type, ann.by.levelM){
  ann.by.all.typeM <- ann.all.type$ann.by.all.typeM

  c.scores <- lapply(marker.listM, function(n) CellAnnotatoR:::getCellTypeScoreInfo(n, t(as.matrix(cm.norm)))) %>%
    lapply(`[[`, "scores")

  types.by.level <- lapply(ann.by.levelM, function(n) intersect(unique(n), c.scores%>%names()))
  scores.per.level <- lapply(types.by.level, function(n) c.scores[n]%>% as.data.frame(optional=T)) %<>% normalizeScores()

  confidence <- scores.per.level %>% lapply(function(x) mapply(function(m,n) x[m,n], lapply(ann.by.all.typeM[names(x)],names), colnames(x)))

  unc.by.cells.per.type <- lapply(confidence,function(m) lapply(m, function(n) 1-n))
  n.per.type <- lapply(confidence, function(m) sapply(m,length))
  mean.unc.by.types <- lapply(confidence, function(m) 1-sapply(m,mean))
  return (list(mean.unc.by.types=mean.unc.by.types, unc.by.cells.per.type=unc.by.cells.per.type,
               n.per.type=n.per.type, marker.list=marker.listM))
}



#' Prepare the "trait" file for plotUncHierarchy() from the positive, negative and coverage uncertainty data (output from scoreCellUncertaintyPerLevel())
#' @param unc.per.cell a list, each element represents the annotation of all cells at on hierarchical level
#' @param ann.by.level a list of lists, each main element is one hierarchical level; at each hierarchical level, there is a sublist, which contains 3 elements: positive, negative and coverage uncertainty.
#' @export ann.unc a dataframe, there are 4 columns:"type" #cell types   "positive" "negative" "coverage" #uncertainty
uncToTreeTrait <- function(unc.per.cell, ann.by.level) {
  ann.unc <- list()
  if (names(ann.by.level[1]) == "annotation") {
    for (i in names(unc.per.cell) %>% setNames(.,.)) {
      ann.unc[[i]] <- merge(as.data.frame(ann.by.level$annotation[[i]]), as.data.frame(unc.per.cell[[i]]),
                            by="row.names", all=TRUE)
      names(ann.unc[[i]])[2] <- "type"
      ann.unc[[i]] <- aggregate(ann.unc[[i]][,3:5], ann.unc[[i]][,2,drop=FALSE], mean)
    }
  } else {
    for (i in names(unc.per.cell) %>% setNames(.,.)) {
      ann.unc[[i]] <- merge(as.data.frame(ann.by.level[[i]]), as.data.frame(unc.per.cell[[i]]),
                            by="row.names", all=TRUE)
      names(ann.unc[[i]])[2] <- "type"
      ann.unc[[i]] <- aggregate(ann.unc[[i]][,3:5], ann.unc[[i]][,2,drop=FALSE], mean)
    }
  }
  
  rbinded <- rbind(ann.unc[[1]], ann.unc[[2]])
  L <- names(ann.unc) %>% length()
  if (L>2) {
    for (l in 3:L){
      rbinded <- rbind(rbinded, ann.unc[[l]])
    }
  }

  ann.unc <- rbinded
  ann.unc %<>% .[!duplicated(ann.unc$type),]
  return (ann.unc)
}



#' Plot uncertainty (positive, negative, coverage) on hierarchical trees
#' @param c.data output file from re-annotation (by assignCellsByScores())
#' @param trait.file contain the mean positive, negative, coverage uncertainties for each cell type (including root, node and tips)
#' @param unc one of c("positive", "negative", "coverage")
#' @export
plotTypeHierarchy <- function(c.data, trait.file, layout="slanted", xlims=NULL, font.size=3, ...) {
  if (!requireNamespace("ggtree", quietly=T))
    stop("You need to install package 'ggtree' to be able to plot hierarchical tree. ",
         "`Try devtools::install_github('YuLab-SMU/ggtree')`")

  c.df <- classificationTreeToDf(c.data$classification.tree)
  cg <- c.df %$% data.frame(parent=Parent, node=Node) %>% ape::as.phylo() #as.phylo(): we got edge, Nnode, node.label, tip.label (c.df only has parent and node)
  cg$edge.length <- rep(1, nrow(cg$edge))

  if (is.null(xlims)) {
    xlims <- c(0, max(c.df$PathLen) + 0.5)
  }

  tf <- trait.file
  names(tf)[1] <- "labels"
  tf$labels <- gsub(" ", "", tf[,1])

  tipid <- cg$tip.label %>% data.frame("labels"=., "node"=ggtree::nodeid(cg, .)) #get tip id
  nodeid <- cg$node.label %>% data.frame("labels"=., "node"=ggtree::nodeid(cg, .)) #get label id
  id <- rbind(tipid, nodeid) %>% as.data.frame()

  tfid <- merge(tf, id, by="labels")
  tfid$node <- as.numeric(tfid$node)
  cgtree <- dplyr::full_join(cg, tfid, by="node")

  a <- ggtree::ggtree(cgtree, aes(color=positive), layout = layout, ...) +
      scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"))+
      ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
      ggplot2::xlim(xlims) #the presence or absence of [node] makes no difference here
  b <- ggtree::ggtree(cgtree, aes(color=negative), layout = layout, ...) +
      scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"))+
      ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
      ggplot2::xlim(xlims)
  c <- ggtree::ggtree(cgtree, aes(color=coverage), layout = layout, ...) +
      scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"))+
      ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
      ggplot2::xlim(xlims)
    cowplot::plot_grid(a, b, c, nrow=1)
}
