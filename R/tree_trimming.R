markerListToMarkupLi <- function(marker.list, file="", top.parents, group.by.parent=T) {
  if (!group.by.parent) {
    markup.text <- names(marker.list) %>% sort() %>%
      lapply(function(n) markersToMarkupLi(marker.list[[n]], top.parents, n)) %>% paste(collapse="")
  } else {
    for (n in names(marker.list)) {
      if (length(marker.list[[n]]$parent) == 0) {
        marker.list[[n]]$parent <- "root"
      }
    }

    ml.per.parent <- marker.list %>% split(sapply(., `[[`, "parent")) %$%
      c(list(root="root"), .[names(.) != "root"])
    markup.text <- names(ml.per.parent) %>% lapply(function(pn)
      paste0("## ", pn, "\n\n", markerListToMarkupLi(ml.per.parent[[pn]], group.by.parent=F))) %>%
      paste0(collapse="")
  }

  if (!is.null(file) && nchar(file) > 0) {
    cat(markup.text, file=file)
    return(file)
  }

  return(markup.text)
}


markersToMarkupLi <- function(markers, root, name) {
  if (!is.null(markers$parent) && (markers$parent %in% root)) {
    markers$parent <- NULL
  }

  expr <- paste0("expressed: ", paste0(markers$expressed, collapse=", "), "\n")
  if (!is.null(markers$not_expressed) && length(markers$not_expressed) > 0) {
    not.expr <- paste0("not expressed: ", paste0(markers$not_expressed, collapse=", "), "\n")
  } else {
    not.expr <- ""
  }

  parent <- if (!is.null(markers$parent)) paste0("subtype of: ", markers$parent[1], "\n") else ""

  return(paste0("> ", name, "\n", expr, not.expr, parent, "\n"))
}



# Typically used for preparing the tree output from re-annotation bu removing singleton types (only has single subtype), so that it can be used for tree pruning
treePrunePrep <- function(ann.by.level, min.cluster.size=6){
  if (ann.by.level[[1]]%>%unique%>%length>1)
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))
                           %>%setNames(names(ann.by.level[[1]]))), ann.by.level)

  parent.singleton <-  lapply(1:(length(ann.by.level)-1), function(i)
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>%
    Reduce(c, .) %>% .[sapply(., function(x) length(unique(x)) == 1)] %>%  #remove %>% .[sapply(., length) > 5] after Reduce
    names %>% setdiff(ann.by.level[[length(ann.by.level)]]%>%unique)

  depth <- length(ann.by.level)
  message("Remove singletons ...")
  for (i in 2:depth){
    singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)
    while (length(singletons)!=0){
      for (j in singletons){
        j.cells <- ann.by.level[[i]] %>% .[.==j] %>% names
        for (k in i:(depth-1))
          ann.by.level[[k]][j.cells] <- ann.by.level[[k+1]][j.cells]}
      singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)}}

  message("Remove small clusters ...")
  #Remove small clusters (size of 5 cells or less):
  for (i in 2:depth){
    small.clusters <- ann.by.level[[i]]%>%split(.,.)%>%.[sapply(.,function(x) length(x)<=min.cluster.size)]%>%names
    if (length(small.clusters)!=0)
      for (s in small.clusters)
        ann.by.level <- rmSmallCluster(s, ann.by.level, i)}

  #Remove repeated layers (at lower level)
  message("Remove repeated layers ...")
  ann.by.level <- rmRepeatLayer(ann.by.level)

  return(ann.by.level)
}


rmSmallCluster <- function(small.clust, ann.by.level, layer.index){ #here layer.index>1
  depth <- length(ann.by.level)
  cells <- ann.by.level[[layer.index]] %>% .[.==small.clust] %>% names
  parent <- ann.by.level[[layer.index-1]][cells] %>% unique
  cells <- ann.by.level[[layer.index-1]] %>% .[.==parent] %>% names
  for (k in depth:layer.index)
    ann.by.level[[k]][cells] <- ann.by.level[[layer.index-1]][cells]

  return(ann.by.level)
}


# rmSingletons <- function(singleton, ann.by.level, layer.index){
#   depth <- length(ann.by.level)
#   cells <- ann.by.level[[layer.index]] %>% .[.==singleton] %>% names
#   # keep <- layer.index+1
#   # for (k in layer.index:(depth-2))
#   #   if (length(unique(ann.by.level[[k+2]][cells]))==1) keep=k+2
#   #
#   # mapply(function(i,j) ann.by.level[[i]][cells]<-ann.by.level[[j]][cells], layer.index:(depth-keep+layer.index), keep:depth)
#   # for (u in (depth-keep+layer.index+1):(depth-1))
#   #   ann.by.level[[u]][cells] <- ann.by.level[[depth]][cells]
#   for (k in layer.index:(depth-1))
#     ann.by.level[[k]][cells] <- ann.by.level[[k+1]][cells]
#   return(ann.by.level)
# }
#
# for (i in 2:depth)
#   #Remove singletons:
#   singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)
# while (length(singletons)>0)
#   for (j in singletons)
#     ann.by.level <- rmSingletons(j, ann.by.level, i)
# singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)


rmRepeatLayer <- function(ann.by.level){
  tree.depth=length(ann.by.level)
  for (i in tree.depth:2) {
    tip.layer.index <- i
    if (length(ann.by.level[[i]]%>%unique %>% setdiff(ann.by.level[[i-1]]%>%unique))==0) {
      tip.layer.index <- i-1
    } else {break}
  }
  ann.by.level %<>% .[1:tip.layer.index]

  return(ann.by.level)
}


getMarkersForLayers23 <- function(p2, ann.by.level, start.layer.index){
  ann.by.parent <- lapply(start.layer.index:(start.layer.index+1), function(i)
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>%
  Reduce(c, .) %>% .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]

  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  if (length(ann.by.parent)>0) {
  #pre-select markers
    de.info.per.parent <- ann.by.parent %>%
     pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T)) #CHANGE FROM F TO T 7.22.

    pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
      list(
        positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
        negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
      )
    )

    #Select markers
    marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
      selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
                         parent=n,max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))

    marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)

    return(list(marker.list=marker.list, ann.by.parent=ann.by.parent))
  } else{
    return(list(marker.list=NULL, ann.by.parent=NULL))
  }
}


getMarkersForNextLayer <- function(p2, ann.by.level, start.layer.index){
  ann.by.parent <- split(ann.by.level[[start.layer.index+1]], ann.by.level[[start.layer.index]])%>%
    .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]
  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()

  if (length(ann.by.parent)>0) {
    #pre-select markers
    de.info.per.parent <- ann.by.parent %>%
      pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))

    pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
      list(
        positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
        negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
      )
    )

    #Select markers
    marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
      selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
                           parent=n, max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))

    marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)

    message("\ngetMarkersForNextLayer is finished!")
    return(marker.list)
  } else{
    message("\ngetMarkersForNextLayer is finished!")
    return(NULL)
  }
}


liftTree <- function(parent, ann.by.level, layer.index){
  depth <- names(ann.by.level)%>%length()
  cells.p <- ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==parent] %>% names
  for (i in (layer.index+1):(depth-1)) {
    ann.by.level[[i]][cells.p] <- ann.by.level[[i+1]][cells.p]}
  message("\nNow the liftTree is finished!")
  return(ann.by.level)
}


#remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
removeParentNoPos <- function(p2, ann.by.level, ann.by.level.sub, layer.index){
  message("\nGet markers for parents under the grandparent ", ann.by.level.sub[[layer.index]]%>%unique)
  markers.subtypes <- getMarkersForNextLayer(p2, ann.by.level.sub, layer.index) #what about the other un-considered parents have no markers?

  types.no.pos <- markers.subtypes[sapply(markers.subtypes, function(x) length(x$expressed) == 0)] %>% names()
  message("\nThe cell types that have no positive markers are ", types.no.pos)
  parents.no.pos <- c()
  for (type in types.no.pos) {
    if (ann.by.level.sub[[layer.index+2]][ann.by.level.sub[[layer.index+1]]%>%.[.==type]%>%names]%>%unique%>%length>1) {
      parents.no.pos<-c(parents.no.pos, type)
    }
  }

  message("\nParent types that have no positive markers are ", parents.no.pos)
  #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
  while (length(parents.no.pos) > 0) {
    for (parent in parents.no.pos) {
      message("The no-positive-marker parent to be removed is ", parent)
      ann.by.level.sub <-
        liftTree(parent, ann.by.level.sub, layer.index)
      ann.by.level <- liftTree(parent, ann.by.level, layer.index)
    }
    message("\nAfter removed the no-positive-marker parent, now we need to get markers for the new parent layer")
    markers.subtypes <- getMarkersForNextLayer(p2, ann.by.level.sub, layer.index)
    types.no.pos <- markers.subtypes[sapply(markers.subtypes, function(x) length(x$expressed) == 0)] %>% names()
    message("\nThe new cell types without positive marker are ", types.no.pos)
    parents.no.pos <- c()
    for (type in types.no.pos) {
      if (ann.by.level.sub[[layer.index + 2]][ann.by.level.sub[[layer.index + 1]] %>% .[. == type] %>% names] %>% unique %>% length > 1) {
        parents.no.pos <- c(parents.no.pos, type)
      }
    }
    message("\nThe new parents without positive marker are ", parents.no.pos)
  }
  message("\nNow, the removeParentNoPos is finished!")
  return(list(ann.by.level=ann.by.level, ann.by.level.sub=ann.by.level.sub, markers.subtypes=markers.subtypes))
}


getScore <- function(ann.by.level, marker.list, cm.norm, layer.index){
  types.no.pos <- marker.list[sapply(marker.list, function(x) length(x$expressed) == 0)] %>% names() #Here the types could be singleton or parent types
  #nonparent.no.pos <- sapply(types.no.pos, function(n) if (ann.by.level[[layer.index+2]][ann.by.level[[layer.index+1]]%>%.[.==n]%>%names]%>%unique%>%length==1) n)

  top.parent <- ann.by.level[[layer.index]] %>% unique()
  marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>%
    markerListToMarkupLi(top.parent, file="~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list.txt",
                         group.by.parent=F)

  clf.data <- getClassificationData(cm.norm, "~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list.txt")

  score.info <- getMarkerScoreInfo(clf.data) #For the tip clusters withough positive markers, we assign a score of 0 to them
  if (length(types.no.pos)!=0){ #
    for (type in types.no.pos) {
      cells.all <- score.info[[1]]$scores.raw %>% names()
      score.info[[type]]$scores.raw <- rep(0, length(cells.all)) %>% setNames(., cells.all)
      score.info[[type]]$mult <- rep(0, length(cells.all)) %>% setNames(., cells.all)
      score.info[[type]]$max.positive <- rep(0, length(cells.all)) %>% setNames(., cells.all)
      score.info[[type]]$scores <- rep(0, length(cells.all)) %>% setNames(., cells.all)
    }
  }
  message("\ngetScore is finished!")
  return(score.info)
}


#In this version, the t-test is removed
#When St (before replacement) < 0.25, St (after replacement) < St(before)+0.05: the replacement shouldn't be accepted.
#When St (before) > 0.25, as long as St(after) >= S(before), the replacement is accepted
#Here the St is normalized, and the normalization is based on the cell types under the same grandparents (the peers oof the parent on focus)
treePruneSingleLayer <- function(p2, ann.by.level, ann.by.level.sub, layer, wei=1){
  k=0
  layer.index <- match(layer, names(ann.by.level))
  message("\nThe current layer.index is ", layer.index)
  if (ann.by.level.sub[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, ann.by.level.sub=ann.by.level.sub, k=0))}
  depth <- names(ann.by.level)%>%length()
  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()

  #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
  message("\nNow remove the parents that have no positive markers...")
  parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, ann.by.level.sub, layer.index)
  ann.by.level.sub <- parent.no.pos.removed$ann.by.level.sub
  ann.by.level <- parent.no.pos.removed$ann.by.level
  message("\nGet markers for the parents under the grandparent type ", ann.by.level.sub[[layer.index]]%>%unique)
  markers.parents <- parent.no.pos.removed$markers.subtypes

  #Only cell types under the same grandparent are considered for St and its normalization
  #markers.parents <- getMarkersForNextLayer(p2, ann.by.level.sub, layer.index)
  message("\nGet St scores for the parents under the grandparent type ", ann.by.level.sub[[layer.index]]%>%unique)
  score.info <- getScore(ann.by.level.sub, markers.parents, cm.norm, layer.index)
  message("\nNormalize the St scores for the parents under the grandparent type ", ann.by.level.sub[[layer.index]]%>%unique)
  score.info %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores() #Normalize St# Here, the normalization is only based on the cell types that share the same grandparent as the parent on focus, not all the cell types at each layer, therefore, the score will differ from the what we get from getMarkerScoresPerCellType(), where all the cell types from each layer are considered.
  message("\nGet each cell's normalized St scores for the cell types that the cells are annotated with")
  score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
    cells<-ann.by.level.sub[[layer.index+1]][ann.by.level.sub[[layer.index+1]]==n]%>%names;
    score.info[cells,n]}) %>% unlist

  parent.types <- ann.by.level.sub[[layer.index+1]] %>% unique #WHAT ABOUT parent.types is empty
  sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
  cat("\nThe parent types under the grandparent type", ann.by.level.sub[[layer.index]]%>%unique, "are", parent.types)
  sink()

  message("\nThe parent types under the grandparent type ", ann.by.level.sub[[layer.index]]%>%unique, " are ", parent.types)
  #There will be as many markers.sub.replaced as the subtypes.
  for (parent in parent.types){
    message("\nThe parent to be tested at this step is ", parent)
    if (!(parent %in% (ann.by.level[[ann.by.level.sub %>% length()]] %>% unique))) { #sub not in the last layer
      sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
      cat("\nTemporarily replace the parent", parent, "by its kids, so we can compare with that without this replacement...")
      sink()

      message("\nTemporarily replace the parent ", parent, " by its kids, so we can compare with that without this replacement...")
      ann.by.level.replaced <- ann.by.level.sub #temporary replace, only replace the parent layer, not the others
      parent.cells <- ann.by.level.sub[[layer.index+1]] %>% .[.==parent] %>% names
      ann.by.level.replaced[[layer.index+1]][parent.cells] <- ann.by.level.replaced[[layer.index+2]][parent.cells]
      message ("\nMarker secltion for all types under the grandparent ", ann.by.level.sub[[layer.index]]%>%unique, " while the parent ", parent, " being replaced by its kids ...")
      markers.parent.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
      if (sum(sapply(markers.parent.replaced, function(x) length(x$expressed) == 0))==0) { #*****''
        message ("\nGet normalized St scores for all the types under the grandparent ", ann.by.level.sub[[layer.index]]%>%unique, " while the parent ", parent, " being replaced by its kids ...")
        score.info.replaced <- getScore(ann.by.level.replaced, markers.parent.replaced, cm.norm, layer.index)
        score.info.replaced %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores()#Normalize St
        score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
          cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
          score.info.replaced[cells,n]}) %>% unlist
        message("\nTest and replace if necessary")
        message("\nIf the grandparent layer is the root...")

        sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
        cat("\nThe mean St score before replacement: ", mean(score.info))
        cat("\nThe mean St score after replacement: ", mean(score.info.replaced))
        sink()

        if ((layer.index==1) && (sum(score.info)<0.75*sum(score.info.replaced))) { #For layer 1
          sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
          cat("\n(layer.index==1) && (sum(score.info)<0.75*sum(score.info.replaced)), so we accept the replacement.")
          sink()

          ann.by.level <- liftTree(parent, ann.by.level, layer.index);
          ann.by.level.sub <- liftTree(parent, ann.by.level.sub, layer.index);
          k <- k+1; score.info <- score.info.replaced}
        if (layer.index>1){#For the other layers:
            message("\nIf the grandparent is not the root...")
            # when mean(score.info <= 0.25)
            if (mean(score.info)<=0.25 && mean(score.info.replaced)>(mean(score.info)+0.05)) {
              sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
              cat("\nFor layer.index>1:")
              cat("\nmean(score.info)<=0.25 && mean(score.info.replaced)>(mean(score.info)+0.05), therefore we accept the replacement.")
              sink()

              message("\nNow we remove parent ", parent, " and lift up the part of the tree below it")
              ann.by.level <- liftTree(parent, ann.by.level, layer.index);
              ann.by.level.sub <- liftTree(parent, ann.by.level.sub, layer.index);
              k <- k+1;
              score.info <- score.info.replaced}

            if (mean(score.info)>0.25 && mean(score.info)<=mean(wei*score.info.replaced))  { # when mean(score.info >0.25)
              sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
              cat("\nmean(score.info)>0.25 && mean(score.info)<=mean(wei*score.info.replaced), therefore we accept the replacement.")
              sink()

              message("\nNow we remove parent ", parent, " and lift up the part of the tree below it")
              ann.by.level <- liftTree(parent, ann.by.level, layer.index);
              ann.by.level.sub <- liftTree(parent, ann.by.level.sub, layer.index);
              k <- k+1;
              score.info <- score.info.replaced}} # For the rest of layers, wei=1
      }
    }
  }
  return(list(ann.by.level=ann.by.level, ann.by.level.sub=ann.by.level.sub, k=k))
}


treePrunesSingleLayer <- function(p2, ann.by.level, ann.by.level.sub, layer, wei=1){
  message("\nNow we are starting to trim the grandparent type ", ann.by.level.sub[[layer]]%>%unique, "s family")
  ann.by.level <- treePruneSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
  k <- ann.by.level$k
  message("\nThe k is ", k)
  ann.by.level.sub <- ann.by.level$ann.by.level.sub
  ann.by.level <- ann.by.level$ann.by.level
  while (k>0){
    sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
    cat("\nBecause k>0, we need to rerun the trimming for the grandparent type", ann.by.level.sub[[layer]]%>%unique)
    sink()

    message("\nBecause k>0, we need to rerun the trimming for the grandparent type ", ann.by.level.sub[[layer]]%>%unique)
    ann.by.level.2 <- treePruneSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
    ann.by.level <- ann.by.level.2$ann.by.level
    ann.by.level.sub <- ann.by.level.2$ann.by.level.sub
    k = ann.by.level.2$k
    message("\nThe k is ", k)
  }
  return(ann.by.level)
}


treePrunes <- function(p2, ann.by.level, wei=1, start.layer.index=1) {
  #ann.by.level <- treePrunePrep(ann.by.level) #remove parent types with single subtype
  if (ann.by.level[[1]]%>%unique%>%length!=1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))%>%setNames(ann.by.level[[1]]%>%names)),ann.by.level)
  }
  sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out") #for debugging
  cat("Let's start")
  sink()

  tree.depth <- length(ann.by.level)
  for (layer in names(ann.by.level)[start.layer.index:(tree.depth-2)]) { #the LAYERS compared actually the two below the starting layer, not the starting layer and the one below it.
    sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
    cat("\n\nThe present layer is", layer)
    sink()

    message("\nThe present layer is ", layer)
    layer.types <- ann.by.level[[layer]] %>% unique()
    for (type in layer.types) {
      sink("/home/rstudio/CellAnnotatoR-dev/notebooks/log.out", append=TRUE) #for debugging
      cat("\nThe present grandparent cell type on focus is", type, " at layer ", layer)
      sink()

      message("\nThe present grandparent cell type on focus is ", type, " at layer ", layer)
      type.cells <- ann.by.level[[layer]] %>% .[.==type] %>% names
      ann.by.level.sub <- list()
      for (i in 1:tree.depth){
        ann.by.level.sub[[names(ann.by.level[i])]] <- ann.by.level[[i]][type.cells]}
      ann.by.level <- treePrunesSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
      message("\nThe trimming for the grandparent ", type, " is now finished!")
    }
  }

  for (i in tree.depth:2) {
    tip.layer.index <- i
    if (length(ann.by.level[[i]]%>%unique %>% setdiff(ann.by.level[[i-1]]%>%unique))==0) {
      tip.layer.index <- i-1
    } else {break}
  }
  ann.by.level %<>% .[1:tip.layer.index]
  return(ann.by.level)
}



# annToTreeDf <- function(ann.by.level) {
#   if (ann.by.level[[1]]%>%unique%>%length>1){
#     ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))
#                            %>%setNames(names(ann.by.level[[1]]))), ann.by.level)
#   }
#
#   depth <- length(ann.by.level)
#   ann.split <- list()
#   for (i in 1:(depth-1)){
#     ann.split[[paste0("l",i)]] <- ann.by.level[[i+1]] %>% split(ann.by.level[[i]]) %>% lapply(unique)
#   }
#
#   parent <- c()
#   node <- c()
#   path.len <- c()
#   for (i in 1:length(ann.split)) {
#     for (j in names(ann.split[[i]])){
#       for (l in ann.split[[i]][[j]]){
#         if (j != l){
#           path.len <- c(path.len,i)
#           parent <- c(parent, j)
#           node <- c(node, l)
#         }
#       }
#     }
#   }
#
#   root <- ann.split[[1]] %>% names()
#   parent[parent==root] <- "root"
#
#   c.df <- tibble(
#     Parent = parent,
#     Node = node,
#     PathLen = as.numeric(path.len)
#   )
#
#   return(c.df)
# }

# # TODO: rename it to getClfDataInferringMarkers
# getClfData <- function(p2, ann.by.level, name, outfile.path){
#   if ((ann.by.level[[1]] %>% unique %>% length) != 1) {
#     ann.by.level <- list(root=rep("root", length(ann.by.level[[1]])) %>%
#                            setNames(names(ann.by.level[[1]]))) %>% c(ann.by.level)
#   }
#
#   ann.by.parent <-  lapply(1:(length(ann.by.level)-1), function(i)
#     split(ann.by.level[[i+1]], ann.by.level[[i]])) %>%
#     Reduce(c, .) %>% .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]
#
#   message("Selecting markers ...")
#   de.info.per.parent <- ann.by.parent %>% #ann: under one parent, cell annotation to its kids for cells only belonging to the given parent.
#     pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))
#
#   pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
#     list(
#       positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
#       negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
#     )
#   )
#
#   cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
#   marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
#     selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
#                          parent=n,max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))
#
#   marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)
#
#   message("Writing marker list to markup file ...")
#   #Check whether all clusters have markers:
#   if ((marker.list %>% .[sapply(., function(x) length(x$expressed) == 0)] %>% length) >0) {return(NULL)}
#
#   marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>% markerListToMarkup(file=paste0(outfile.path, "marker_list_", name, ".txt"),
#                        group.by.parent=F) #I may need to update this function in the future
#
#   clf.data <- getClassificationData(cm.norm, paste0(outfile.path, "marker_list_", name, ".txt"))
#   return(clf.data)
# }


#name: marker file name
# #outfile.path: where to save the marker file
# reAnnoPip <- function(p2, ann.by.level, name, outfile.path, graph=NULL, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
#   if (is.null(clf.data)){clf.data <- getClfData(p2, ann.by.level, name, outfile.path)}
#   message("Re-annotation ...")
#   if (is.null(clf.data)) {return(NULL)}
#   ann.by.level <-
#     assignCellsByScores(graph=NULL, clf.data, #The graph is removed here
#                         clusters = NULL, uncertainty.thresholds=uncertainty.thresholds) #we had clusters=p2$clusters$PCA$leiden, but lots of tip clusters dropped after re-annotation, therefore, I reset it to NULL.
#   #saveRDS(ann.by.level, paste0(outfile.path, "ann_by_level_", name, ".rds"))
#
#   # message("Checking the data ...")
#   # ann.by.level$annotation[[length(ann.by.level$annotation)]]%>%unique()
#   # ann.by.level$annotation %>% sapply(unique) %>% Reduce(c,.) %>% unique() %>% length()
#   #
#   # message("Preparing file for plotting hierarchy ...")
#   # c.df <- annToTreeDf(ann.by.level$annotation)
#   # saveRDS(c.df, file=paste0(outfile.path, "c_df_", name, ".rds")) #For plotting hierarchy
#
#   return(list(ann.by.level=ann.by.level, clf.data=clf.data))
# }



plotUncertaintyPerClustByLevel <- function(cell.uncertainty.per.level, annotation.per.level, layer){
  unc.per.clust <- scoreClusterUncertaintyPerLevel(cell.uncertainty.per.level, annotation.per.level[[layer]])
  plotUncertaintyPerClust(unc.per.clust[[layer]], annotation.per.level[[layer]], n.col=3, text.angle=60)}


#Slice tree layers and save them to files, so that we can use the files for self-projection analyses in python
sliceTreeLayers <- function(ann.by.level, outfile.path, name.base) {
  depth <- length(ann.by.level)
  for (l in 1:depth){
    ge.obs <- as.data.frame(cbind(ann.by.level[[l]] %>% names, ann.by.level[[l]] %>% Reduce(c,.)));
    colnames(ge.obs) <- c("cells", "clusters");
    write.csv(ge.obs, file = paste0(outfile.path, name.base, l, ".csv"), row.names = FALSE)}
}


#For self-projection
#Function to find the corresponding layers from the two compared trees based on the unique cell types (those shared the most cell types):
findLayerMatches <- function(layer, tree){
  matched = 0
  match.rate = 0
  for (k in 1:length(tree)) {
    rate <- sum(!is.na(match(unique(tree[[k]]), unique(layer))))/length(union(unique(tree[[k]]), unique(layer)))
    if (rate>match.rate) {
      matched=k
      match.rate=rate}
  }
  return(matched)
}


#Function for calculating the self-projection accuracy for each layer of the simpler one of two compared trees
getSelfProjAcc <- function(tree1, tree2){
  len1 <- if (length(unique(tree1[[1]]))==1) length(tree1)-1 else length(tree1)
  len2 <- if (length(unique(tree2[[1]]))==1) length(tree2)-1 else length(tree2)
  acc <- c()
  if (len1>=len2) {
    message("The second tree is simpler.")
    for (m in tree2){
      matched=findLayerMatches(m, tree1)
      acc <- c(acc, sum(m==tree1[[matched]])/length(m))} #== only apply to two vectors with equal lengths
    return(acc)}

  for (m in tree1){
    meassage("The first tree is simpler.")
    matched=findLayerMatches(m, tree2)
    acc <- c(acc, sum(m==matched)/length(m))}
  return(acc)
}


#Calculating accuracy for re-annotation results (comparing the annotation before and after the re-annotation)
#The annotation before the re-annotation is considered to be the groundtruth annotation
#Here, tree1 is always the original annotation tree (assumed to be truth), while tree2 is always the simpler one (less layers)
measureQ <- function(len, start, tree1, tree2, verbose=F){
  conm <- lapply (start:len%>%setNames(paste0("l", start:len)), function(m) {
    matched.index=findLayerMatches(tree2[[m]], tree1) #find the matched layer of tree2 in tree1
    matched <- tree1[[matched.index]]
    sapply(unique(matched), function(c){
      positive <- matched[matched==c] %>% names #The cells belong to cluster c in the orginal annotation-tree tree1 (supposed to be the truth)
      for (parent0 in split(tree1[[m+1]], tree1[[m]])) {
        if (c %in% unique(parent0)) mother=parent0} #positive and negative should be under the same parent.

      negative <- mother[!(mother==c)] %>% names  #The negative ones are those under the same parent but in different clusters in the original tree
      if ((length(negative)==0) && (verbose==T)) warning("There are no negative cells for cluster ", c, " because its parent is a singleton! ")
      predicted <- tree2[[m]][tree2[[m]]==c] %>% names
      tpr <- sum(!is.na(match(positive, predicted)))/length(positive) #acc = TP rate = TP/(TP+FN)
      fpr<-sum(!is.na(match(negative, predicted)))/length(negative) #Be careful, match is really tricky, when there is duplicates in the compared vectors
      #Calculate the false positive (The truth is Negative but is predicted to be Positive) rate=FP/(FP+TN)
      precision <- sum(!is.na(match(positive, predicted)))/(sum(!is.na(match(positive, predicted))) + sum(!is.na(match(negative, predicted))))
      cr <- sum(!is.na(match(negative, predicted)))/sum(!is.na(match(positive, predicted)))
      list(precision=precision, recall=tpr, false.positive.rate=fpr, confusion.rate=cr)
    })
  })
  return(conm)
}


# #Function for calculating the per-cluster self-projection accuracy for each layer of the simpler one of two compared trees
# getSelfProjQcPerCluster <- function(tree1, tree2){ #Please provide the original annotation-tree as tree1 and the predicted annotation-tree as tree2
#   if (length(unique(tree1[[1]]))>1) {tree1 <- c(list(root=rep("root", length(tree1[[1]]))%>%setNames(., names(tree1[[1]]))), tree1)}
#
#   len1 <- length(tree1)-1
#   len2 <- if (length(unique(tree2[[1]]))==1) length(tree2)-1 else length(tree2)
#
#   if (len1>=len2){
#     message("The second tree is simpler.")
#     start <- if (len2<length(tree2)) 2 else 1
#     return(measureQ(len2, start, tree1, tree2))
#   }
#   #When len1<len2:
#   message("The original tree is simpler.")
#   start <- if (len1<length(tree1)) 2 else 1
#   return(measureQ(len1, start, tree2, tree1))
# }
#
#




#Here, we use correlation coefficient of two cells/clusters for measuring their cells
#Since snRNA-seq does not follow normal distribution, therefore we will use non-parametric spearman correlation, rather than pearson correlation.
getHcHierTrunc <- function(p2, trunc.quantile=0.985, min.cluster.size=6){
  feature.matrix=p2$reductions$PCA
  c.dist <- (1 - cor(Matrix::t(feature.matrix), method="spearman")) %>% as.dist() #as.dist(), only oen side 	of the diagonal is maintained
  #for cor(), the default method="pearson", which is suitable for variables that are both normally distributed, which is however not the case for snRNA-seq, therefore, the non-parametric methods "spearman" should be more proper
  #Distance between cells*************
  clusts <- hclust(c.dist, method="ward.D2") #The default method = "complete": Minimize the maximum distance between two points in different clusters (farthest neighbour), other method options:  METHODS <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
  #I wanted to use "centroid", however, it is not leading to a monotone distance measure, or equivalently the resulting dendrograms can have so called inversions or reversals which are hard to interpret.
  #Therefore, I will choose "ward.D2", Ward's minimum variance ( Minimize the squared euclidian distances to the cluster means) method aims at finding compact, spherical clusters
  #ward.D2 implement Ward's (1963) clustering criterion (while "ward" does not)
  #Distance between clusters******************
  #trunc.quantile: where to the truncate the tree, e.g. 0.985 quantile
  cut <- quantile(clusts$height, trunc.quantile)
  h <- clusts$height[clusts$height>cut[[1]]] %>% rev()
  clusts.hc <- cutree(clusts, h=h) %>% as.matrix() #get a tree with all the layers I want to keep for now.

  # Rename each layer as l1, l2 ..., and rename the clusters at each layer as l1_1, l1_2, l1_3 ....
  for (i in 1:ncol(clusts.hc)) {
    clusts.hc[,i] %<>% paste0("l", i, "_", .)
  }

  # At each cutting line, this program assign all the cells to the available clusters at the cutting lines.
  # Here the output is like the annotation at different layers, though in the format of a matrix
  clusts.hc %<>% cbind(rownames(.)) %>% `colnames<-`(paste0("l", 1:ncol(.))) %>% 	tibble::as_tibble()

  message("Change the matrix/Table to list ...")
  hierarchy.hc <- splitClusteringDf(clusts.hc) %>% simplifyHierarchy()

  # Get ann.by.level
  ann.by.level <- list()
  layers.no <- ncol(clusts.hc)
  ann.by.level[["root"]] <- rep("root", length(clusts.hc[,1][[1]])) %>% setNames(., 	clusts.hc[,layers.no][[1]])
  for (layer in 2:(layers.no-1)) {
    name <- colnames(clusts.hc[,layer])
    ann.by.level[[name]] <- clusts.hc[,layer][[1]] %>% setNames(., 	clusts.hc[,layers.no][[1]])
  }
  message("Now the ann.by.level has been created ...")
  message("Simply ann.by.level ...")
  ann.by.level <- treePrunePrep(ann.by.level, min.cluster.size=min.cluster.size)

  return(ann.by.level)
}

#Read in the "self.proj.acc", "conf.rate.norm" results from the parentfortrimkid.txt file
parseConfFile <- function(file.path){
  lines <- readLines(file.path)
  acc.conf <- list()
  n=0
  for (i in 1:length(lines)) {
    if (strsplit(lines[i], " ", fixed=T) %>% sapply(`[[`,1)=="#"){
      if (n==1){
        pn <- length(block)/5
        lp <- list()
        for (j in 1:pn) {
          parent <- block[5*j-4]
          lp <- c(lp, list(setNames(block[(5*j-1):(5*j)]%>%as.numeric, block[(5*j-3):(5*j-2)]))%>%setNames(parent))}
        acc.conf = c(acc.conf, list(lp)%>%setNames(layer))
      }
      layer <- strsplit(lines[i], " ", fixed=T) %>% sapply(`[[`,2);
      n=1;
      block <- c()} else{
        block <- c(block, strsplit(lines[i], " ", fixed=T)[[1]])}
    if (i==length(lines)){
      pn <- length(block)/5
      lp <- list()
      for (j in 1:pn) {
        parent=block[5*j-4];
        lp <- c(lp, list(setNames(block[(5*j-1):(5*j)]%>%as.numeric, block[(5*j-3):(5*j-2)]))%>%setNames(parent))}
      acc.conf = c(acc.conf, list(lp)%>%setNames(layer))}
  }
  return(acc.conf)
}

#Cut kids and the tree below off if normalized confusion rate of a parent >0.3
treeTruncate <- function(parent, ann.by.level, layer.index){
  depth <- names(ann.by.level)%>%length()
  layer <- ann.by.level[[layer.index]]
  cells.p <- layer[layer==parent] %>% names
  for (i in (layer.index+1):depth) {
    ann.by.level[[i]][cells.p] <- layer[cells.p]}
  return(ann.by.level)
}

#Repeat treeCut through the entire tree
treeTruncates <- function(ann.by.level, acc.conf, conf.thresh){
  for (l in 1:(length(ann.by.level)-1)){
    for (k in unique(ann.by.level[[l]])){
      if ((k %in% names(acc.conf)) && acc.conf[[k]][2]>conf.thresh){
        ann.by.level <- treeTruncate(k, ann.by.level, l)}
    }
  }
  tree.depth <- length(ann.by.level)
  for (i in tree.depth:2) { #remove repeated layers
    tip.layer.index <- i
    if (length(ann.by.level[[i]]%>%unique %>% setdiff(ann.by.level[[i-1]]%>%unique))==0) {
      tip.layer.index <- i-1
    } else {break}
  }
  ann.by.level %<>% .[1:tip.layer.index]
  return(ann.by.level)
}


#withdraw each quality parameters from the output of measureQ() for each clusters (unique) of one entire tree
getSPQcEle <- function(self.proj.qc, row.index){
  len <- length(self.proj.qc)
  ele <- 1:(len-1) %>% sapply(function(l) {
    ltypes <- setdiff(self.proj.qc[[l]]%>%colnames, self.proj.qc[[len]]%>%colnames)
    self.proj.qc[[l]][row.index, ltypes]%>%Reduce(c,.)})
  self.proj.qc[[len]][row.index, ]%>%Reduce(c,.)
  return(c(ele%>%Reduce(c,.), self.proj.qc[[len]][row.index, ]%>%Reduce(c,.)))
}


# #Self Projection Quality Control
# #annotation.before <- p4$clusters$PCA$leiden #annotation before re-annotation (reAnnoPip)
# #annotation.after <- reann.p4$ann.by.level$annotation$l1 #annotation after re-annotation (reAnnoPip)
# #This is used for the self projection quaility control after re-annotating the leiden-clustering-based hierarchy
# selfProjQCperCluster <- function(annotation.before, annotation.after) {
#   #Get self projection accuracy:
#   acc <- sum(annotation.before == annotation.after)/length(annotation.before) #The annotation before and after need to have the same length.
#   #Get confusion rate (FP/TP): the cluster sample sized is balanced by adjusting weights inversely proportional to class frequencies in the input data
#   clusters <- unique(annotation.before)
#   conf.fpr <- sapply(clusters, function(x) {
#     wei.x <- length(annotation.before)/table(annotation.before)[x]
#     wei.other <- 1/(1-(1/wei.x))
#     positive <- annotation.before%>%.[.==x] %>% names
#     negative <- names(annotation.before) %>% setdiff(positive)
#     predicted <- annotation.after%>%.[.==x] %>% names
#     true.positive <- positive[!is.na(match(positive, predicted))]
#     false.positive <- negative[!is.na(match(negative, predicted))]
#     #Balanced confusion rate
#     conf <- (length(false.positive)*wei.other)/(length(true.positive)*wei.x)
#     #False positive rate
#     fpr <- length(false.positive)/length(negative)
#     c(conf=conf, fpr=fpr)
#   })
#   colnames(conf.fpr) <- clusters
#   rownames(conf.fpr) <- c("conf", "fpr")
#   return(list(accuracy=acc, confusion.rate.false.pos.rate=conf.fpr))
# }

