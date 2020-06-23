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
treePrunePrep <- function(ann.by.level){
  if (ann.by.level[[1]]%>%unique%>%length>1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))
                           %>%setNames(names(ann.by.level[[1]]))), ann.by.level)
  }
  
  parent.singleton <-  lapply(1:(length(ann.by.level)-1), function(i) 
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>%
    Reduce(c, .) %>% .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) == 1)] %>% 
    names %>% setdiff(ann.by.level[[length(ann.by.level)]]%>%unique)
  
  depth <- length(ann.by.level)
  for (i in 1:depth){
    singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)
    while (length(singletons)!=0){
      for (j in singletons){
        j.cells <- ann.by.level[[i]] %>% .[.==j] %>% names
        for (k in i:(depth-1))
          ann.by.level[[k]][j.cells] <- ann.by.level[[k+1]][j.cells]
      }
      singletons <- intersect(parent.singleton, ann.by.level[[i]] %>% unique)
    }
  }
  
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
  
  cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
  if (length(ann.by.parent)>0) {
  #pre-select markers
    de.info.per.parent <- ann.by.parent %>%
     pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=F))
  
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






treePruneSingleLayer2 <- function(p2, ann.by.level, layer, wei=1){
  n=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, n=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, layer.index)
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    score.info <- getScore(ann.by.level, parent.no.pos.removed$markers.subtypes, cm.norm, layer.index)
    # score.info <- sapply(names(score.info), function(n) {
    #   cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
    #   score.info[[n]]$scores[cells]%>%sum})
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[[n]]$scores[cells]})   
    
    subtypes <- ann.by.level[[layer.index+1]] %>% unique 
    ann.by.level.trimmed <- ann.by.level
    #There will be as many markers.sub.replaced as the subtypes.
    for (sub in subtypes){
      if (!(sub %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) {
        ann.by.level.replaced <- ann.by.level
        sub.cells <- ann.by.level[[layer.index+1]] %>% .[.==sub] %>% names
        
        sub.subtypes <- ann.by.level[[layer.index+2]][sub.cells] %>% unique
        ann.by.level.replaced[[layer.index+1]][sub.cells] <- ann.by.level.replaced[[layer.index+2]][sub.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", sub, " being replaced by its subtypes ...")
        markers.sub.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        #Only take the markers for the subtypes on focus!
        markers.sub.replaced %<>% .[sub.subtypes]
        # if (sum(sapply(markers.sub.replaced, function(x) length(x$expressed) > 0))>1) {
        #   score.info.replaced <- getScore(ann.by.level.replaced, markers.sub.replaced, cm.norm, layer.index)
        #   score.info.replaced <- sapply(sub.subtypes, function(n) {
        #     cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names; 
        #     score.info.replaced[[n]]$scores[cells]%>%sum}) %>% sum
        if (sum(sapply(markers.sub.replaced, function(x) length(x$expressed) > 0))>1) {
          score.info.replaced <- getScore(ann.by.level.replaced, markers.sub.replaced, cm.norm, layer.index)
          score.info.replaced <- lapply(sub.subtypes%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names; 
            score.info.replaced[[n]]$scores[cells]}) %>% unlist     
          #Test and replace if necessary
          t <- t.test(score.info[[sub]], score.info.replaced)
          #For the second layer, we want to keep the big clusters, so set wei to 0.8 here
          if ((layer.index==1) && (t$p.value<0.05)&&(t$estimate[[1]]<0.8*t$estimate[[2]])) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1} else {
          if ((layer.index>1) &&(t$p.value<0.05)&&(t$estimate[[1]]<wei*t$estimate[[2]])) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1}} # For the rest of layers, wei=1
          #if (score.info[[sub]] <= score.info.replaced) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1}        
        }
      }
    }
    return(list(ann.by.level=ann.by.level.trimmed, n=n))
  }
}





treePruneSingleLayer5 <- function(p2, ann.by.level, layer, wei=1){
  n=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, n=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, layer.index)
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    score.info <- getScore(ann.by.level, parent.no.pos.removed$markers.subtypes, cm.norm, layer.index)
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[[n]]$scores[cells]}) %>% unlist
    
    parent.types <- ann.by.level[[layer.index+1]] %>% unique 
    #There will be as many markers.sub.replaced as the subtypes.
    for (parent in parent.types){
      if (!(parent %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) { #sub not in the last layer
        ann.by.level.replaced <- ann.by.level
        parent.cells <- ann.by.level[[layer.index+1]] %>% .[.==parent] %>% names
        ann.by.level.replaced[[layer.index+1]][parent.cells] <- ann.by.level.replaced[[layer.index+2]][parent.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", parent, " being replaced by its subtypes ...")
        markers.parent.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        if (sum(sapply(markers.parent.replaced, function(x) length(x$expressed) == 0))==0) {
          score.info.replaced <- getScore(ann.by.level.replaced, markers.parent.replaced, cm.norm, layer.index)
          score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
            score.info.replaced[[n]]$scores[cells]}) %>% unlist
          #Test and replace if necessary
          t <- t.test(score.info, score.info.replaced)
          if ((layer.index==1) && (t$p.value<0.05)&&(t$estimate[[1]]<0.75*t$estimate[[2]])) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); n <- n+1; score.info <- score.info.replaced} else {
            #For the second layer, we want to keep the big clusters, so set wei to 0.8 here
            if ((layer.index>1) && (t$p.value<0.05)&&(t$estimate[[1]]<wei*t$estimate[[2]])) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); n <- n+1; score.info <- score.info.replaced}} # For the rest of layers, wei=1
        }
      }
    }
    return(list(ann.by.level=ann.by.level, n=n))
  }
}



treePruneSingleLayer4 <- function(p2, ann.by.level, layer, wei=1){
  n=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, n=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, layer.index)
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    score.info <- getScore(ann.by.level, parent.no.pos.removed$markers.subtypes, cm.norm, layer.index)
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[[n]]$scores[cells]}) %>% unlist
    
    subtypes <- ann.by.level[[layer.index+1]] %>% unique 
    ann.by.level.trimmed <- ann.by.level
    #There will be as many markers.sub.replaced as the subtypes.
    for (sub in subtypes){
      if (!(sub %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) {
        ann.by.level.replaced <- ann.by.level
        sub.cells <- ann.by.level[[layer.index+1]] %>% .[.==sub] %>% names
        sub.subtypes <- ann.by.level[[layer.index+2]][sub.cells] %>% unique
        ann.by.level.replaced[[layer.index+1]][sub.cells] <- ann.by.level.replaced[[layer.index+2]][sub.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", sub, " being replaced by its subtypes ...")
        markers.sub.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        if (sum(sapply(markers.sub.replaced, function(x) length(x$expressed) > 0))>1) {
          score.info.replaced <- getScore(ann.by.level.replaced, markers.sub.replaced, cm.norm, layer.index)
          score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
            score.info.replaced[[n]]$scores[cells]}) %>% unlist
          #Test and replace if necessary
          t <- t.test(score.info, score.info.replaced)
          if ((layer.index==1) && (t$p.value<0.05)&&(t$estimate[[1]]<0.75*t$estimate[[2]])) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1} else {
          #For the second layer, we want to keep the big clusters, so set wei to 0.8 here
            if ((layer.index>1) && (t$p.value<0.05)&&(t$estimate[[1]]<wei*t$estimate[[2]])) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1}} # For the rest of layers, wei=1
          #if (score.info[[sub]] <= score.info.replaced) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1}        
        }
      }
    }
    return(list(ann.by.level=ann.by.level.trimmed, n=n))
  }
}



treePruneSingleLayer3 <- function(p2, ann.by.level, layer, wei=1){
  n=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, n=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, layer.index)
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    score.info <- getScore(ann.by.level, parent.no.pos.removed$markers.subtypes, cm.norm, layer.index)
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[[n]]$scores[cells]}) %>% unlist
    
    subtypes <- ann.by.level[[layer.index+1]] %>% unique 
    ann.by.level.trimmed <- ann.by.level
    #There will be as many markers.sub.replaced as the subtypes.
    for (sub in subtypes){
      if (!(sub %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) {
        ann.by.level.replaced <- ann.by.level
        sub.cells <- ann.by.level[[layer.index+1]] %>% .[.==sub] %>% names
        sub.subtypes <- ann.by.level[[layer.index+2]][sub.cells] %>% unique
        ann.by.level.replaced[[layer.index+1]][sub.cells] <- ann.by.level.replaced[[layer.index+2]][sub.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", sub, " being replaced by its subtypes ...")
        markers.sub.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        if (sum(sapply(markers.sub.replaced, function(x) length(x$expressed) > 0))>1) {
          score.info.replaced <- getScore(ann.by.level.replaced, markers.sub.replaced, cm.norm, layer.index)
          score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
            score.info.replaced[[n]]$scores[cells]}) %>% unlist
          #Test and replace if necessary
          t <- t.test(score.info, score.info.replaced)
          #For the second layer, we want to keep the big clusters, so set wei to 0.8 here
          if ((t$p.value<0.05)&&(t$estimate[[1]]<wei*t$estimate[[2]])) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1} # For the rest of layers, wei=1
          #if (score.info[[sub]] <= score.info.replaced) {ann.by.level.trimmed <- liftTree(sub, ann.by.level.trimmed, layer.index); n <- n+1}        
        }
      }
    }
    return(list(ann.by.level=ann.by.level.trimmed, n=n))
  }
}




treePruneSingleLayer6 <- function(p2, ann.by.level, layer, wei=1){
  n=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, n=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, layer.index)
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    score.info <- getScore(ann.by.level, parent.no.pos.removed$markers.subtypes, cm.norm, layer.index)
    score.info %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores() #Normalize St# Here, the normalization is only based on the cell types that are considered at each layer, not all the cell types at each layer, therefore, the score will differ from the what we get from getMarkerScoresPerCellType(), where all the cell types from each layer are considered.
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[n][cells,]}) %>% unlist
    
    parent.types <- ann.by.level[[layer.index+1]] %>% unique 
    #There will be as many markers.sub.replaced as the subtypes.
    for (parent in parent.types){
      if (!(parent %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) { #sub not in the last layer
        ann.by.level.replaced <- ann.by.level
        parent.cells <- ann.by.level[[layer.index+1]] %>% .[.==parent] %>% names
        ann.by.level.replaced[[layer.index+1]][parent.cells] <- ann.by.level.replaced[[layer.index+2]][parent.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", parent, " being replaced by its subtypes ...")
        markers.parent.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        if (sum(sapply(markers.parent.replaced, function(x) length(x$expressed) == 0))==0) {
          score.info.replaced <- getScore(ann.by.level.replaced, markers.parent.replaced, cm.norm, layer.index)
          score.info.replaced %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores(min.val=1e-10)#Normalize St
          score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
            score.info.replaced[n][cells,]}) %>% unlist
          #Test and replace if necessary
          #t <- t.test(score.info, score.info.replaced)
          if ((layer.index==1) && (sum(score.info)<0.75*sum(score.info.replaced))) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); n <- n+1; score.info <- score.info.replaced} else {
            #For the other layers:
            if ((layer.index>1) && (sum(score.info)<=wei*sum(score.info.replaced))) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); n <- n+1; score.info <- score.info.replaced}} # For the rest of layers, wei=1
        }
      }
    }
    return(list(ann.by.level=ann.by.level, n=n))
  }
}



getMarkersForNextLayer <- function(p2, ann.by.level, start.layer.index){
  ann.by.parent <- split(ann.by.level[[start.layer.index+1]], ann.by.level[[start.layer.index]])%>% 
    .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]
  cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
  
  if (length(ann.by.parent)>0) {
    #pre-select markers
    de.info.per.parent <- ann.by.parent %>%
      pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=F))
    
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
    
    return(marker.list)
  } else{
    return(NULL)
  }
}



liftTree <- function(parent, ann.by.level, layer.index){
  depth <- names(ann.by.level)%>%length()
  cells.p <- ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==parent] %>% names 
  for (i in (layer.index+1):(depth-1)) {
    ann.by.level[[i]][cells.p] <- ann.by.level[[i+1]][cells.p]}
  return(ann.by.level)
}



#remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
removeParentNoPos <- function(p2, ann.by.level, ann.by.level.sub, layer.index){
  markers.subtypes <- getMarkersForNextLayer(p2, ann.by.level.sub, layer.index) #what about the other un-considered parents have no markers?
  types.no.pos <- markers.subtypes[sapply(markers.subtypes, function(x) length(x$expressed) == 0)] %>% names()
  parents.no.pos <- c()
  for (type in types.no.pos) {
    if (ann.by.level.sub[[layer.index+2]][ann.by.level.sub[[layer.index+1]]%>%.[.==type]%>%names]%>%unique%>%length>1) {
      parents.no.pos<-c(parents.no.pos, type)}}
  #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
  while (length(parents.no.pos)>0){ 
    for (parent in parents.no.pos) {ann.by.level.sub<-liftTree(parent, ann.by.level.sub, layer.index); ann.by.level<-liftTree(parent, ann.by.level, layer.index)}
    markers.subtypes <- getMarkersForNextLayer(p2, ann.by.level.sub, layer.index)
    types.no.pos <- markers.subtypes[sapply(markers.subtypes, function(x) length(x$expressed) == 0)] %>% names()
    parents.no.pos <- c()
    for (type in types.no.pos) {
      if (ann.by.level.sub[[layer.index+2]][ann.by.level.sub[[layer.index+1]]%>%.[.==type]%>%names]%>%unique%>%length>1) {
        parents.no.pos<-c(parents.no.pos, type)}}
  } 
  
  return(list(ann.by.level=ann.by.level, ann.by.level.sub=ann.by.level.sub, markers.subtypes=markers.subtypes))
}





getScore <- function(ann.by.level, marker.list, cm.norm, layer.index){
  types.no.pos <- marker.list[sapply(marker.list, function(x) length(x$expressed) == 0)] %>% names() #Here the types could be singleton or parent types
  #nonparent.no.pos <- sapply(types.no.pos, function(n) if (ann.by.level[[layer.index+2]][ann.by.level[[layer.index+1]]%>%.[.==n]%>%names]%>%unique%>%length==1) n)
  
  top.parent <- ann.by.level[[layer.index]] %>% unique()
  marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>% 
    markerListToMarkupLi(top.parent, file="~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list.txt", 
                         group.by.parent=F)     
  
  clf.data <- getClassificationData(cm.norm, "~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list.txt", 
                                    prenormalized=T)
  
  score.info <- getMarkerScoreInfo(clf.data) #For the tip clusters withough positive markers, we assign a score of 0 to them
  if (length(types.no.pos)!=0){
    for (type in types.no.pos)
      cells.all <- score.info[[1]]$scores.raw %>% names()
    score.info[[type]]$scores.raw <- rep(0, length(cells.all)) %>% setNames(., cells.all)
    score.info[[type]]$mult <- rep(0, length(cells.all)) %>% setNames(., cells.all)
    score.info[[type]]$max.positive <- rep(0, length(cells.all)) %>% setNames(., cells.all)
    score.info[[type]]$scores <- rep(0, length(cells.all)) %>% setNames(., cells.all)
  }
  
  return(score.info)
}



treePruneSingleLayer <- function(p2, ann.by.level, ann.by.level.sub, layer, wei=1){
  k=0
  layer.index <- match(layer, names(ann.by.level))
  if (ann.by.level.sub[[layer.index+1]]%>%unique%>%length==1) {return(list(ann.by.level=ann.by.level, k=0))} else {
    depth <- names(ann.by.level)%>%length()  
    cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
    
    #remove the parent type at layer.index+1 if no marker identified for it, and then lift the entire related lineage up.
    parent.no.pos.removed <- removeParentNoPos(p2, ann.by.level, ann.by.level.sub, layer.index)
    ann.by.level.sub <- parent.no.pos.removed$ann.by.level.sub
    ann.by.level <- parent.no.pos.removed$ann.by.level
    
    markers.parents <- getMarkersForNextLayer(p2, ann.by.level, layer.index)
    score.info <- getScore(ann.by.level, markers.parents, cm.norm, layer.index)
    score.info %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores() #Normalize St# Here, the normalization is only based on the cell types that are considered at each layer, not all the cell types at each layer, therefore, the score will differ from the what we get from getMarkerScoresPerCellType(), where all the cell types from each layer are considered.
    score.info <- lapply(names(score.info)%>%setNames(.,.), function(n) {
      cells<-ann.by.level[[layer.index+1]][ann.by.level[[layer.index+1]]==n]%>%names;
      score.info[cells,n]}) %>% unlist
    
    parent.types <- ann.by.level.sub[[layer.index+1]] %>% unique 
    #There will be as many markers.sub.replaced as the subtypes.
    for (parent in parent.types){
      if (!(parent %in% (ann.by.level[[ann.by.level %>% length()]] %>% unique))) { #sub not in the last layer
        ann.by.level.replaced <- ann.by.level #******
        parent.cells <- ann.by.level[[layer.index+1]] %>% .[.==parent] %>% names
        ann.by.level.replaced[[layer.index+1]][parent.cells] <- ann.by.level.replaced[[layer.index+2]][parent.cells]
        message ("\nMarker secltion for the layer layer.index+1 with ", parent, " being replaced by its subtypes ...")
        markers.parent.replaced <- getMarkersForNextLayer(p2, ann.by.level.replaced, layer.index)
        if (sum(sapply(markers.parent.replaced, function(x) length(x$expressed) == 0))==0) { #*****''
          score.info.replaced <- getScore(ann.by.level.replaced, markers.parent.replaced, cm.norm, layer.index)
          score.info.replaced %<>% lapply(function(n) n$scores) %>% as.data.frame() %>% normalizeScores()#Normalize St
          score.info.replaced <- lapply(names(score.info.replaced)%>%setNames(.,.), function(n) {
            cells<-ann.by.level.replaced[[layer.index+1]][ann.by.level.replaced[[layer.index+1]]==n]%>%names;
            score.info.replaced[cells,n]}) %>% unlist
          #Test and replace if necessary
          if ((layer.index==1) && (sum(score.info)<0.75*sum(score.info.replaced))) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); ann.by.level.sub <- liftTree(parent, ann.by.level.sub, layer.index); k <- k+1; score.info <- score.info.replaced} else {
            #For the other layers:
            t <- t.test(score.info, wei*score.info.replaced)
            if ((layer.index>1) && (t$p.value<0.05) && (t$estimate[[1]]<t$estimate[[2]])) {ann.by.level <- liftTree(parent, ann.by.level, layer.index); ann.by.level.sub <- liftTree(parent, ann.by.level.sub, layer.index); k <- k+1; score.info <- score.info.replaced}} # For the rest of layers, wei=1
        }
      }
    }
    return(list(ann.by.level=ann.by.level, ann.by.level.sub=ann.by.level.sub, k=k))
  }
}



treePrunesSingleLayer <- function(p2, ann.by.level, ann.by.level.sub, layer, wei=1){
  ann.by.level <- treePruneSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
  k <- ann.by.level$k
  ann.by.level.sub <- ann.by.level$ann.by.level.sub
  ann.by.level <- ann.by.level$ann.by.level
  while (k>0){
    ann.by.level.2 <- treePruneSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
    ann.by.level <- ann.by.level.2$ann.by.level
    ann.by.level.sub <- ann.by.level.2$ann.by.level.sub
    k = ann.by.level.2$k
  }
  return(ann.by.level)
}


treePrunes <- function(p2, ann.by.level, wei=1, start.layer.index=1) {
  #ann.by.level <- treePrunePrep(ann.by.level) #remove parent types with single subtype
  if (ann.by.level[[1]]%>%unique%>%length!=1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))%>%setNames(ann.by.level[[1]]%>%names)),ann.by.level)
  }
  
  tree.depth <- length(ann.by.level)
  for (layer in names(ann.by.level)[start.layer.index:(tree.depth-2)]) { #the LAYERS compared actually the two below the starting layer, not the starting layer and the one below it.
    layer.types <- ann.by.level[[layer]] %>% unique()
    for (type in layer.types) {
      type.cells <- ann.by.level[[layer]] %>% .[.==type] %>% names
      ann.by.level.sub <- list()
      for (i in 1:tree.depth){
        ann.by.level.sub[[names(ann.by.level[i])]] <- ann.by.level[[i]][type.cells]}
      ann.by.level <- treePrunesSingleLayer(p2, ann.by.level, ann.by.level.sub, layer, wei)
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



treePrunesSingleLayer2 <- function(p2, ann.by.level, layer, wei=1){
  ann.by.level <- treePruneSingleLayer(p2, ann.by.level, layer, wei)
  n <- ann.by.level$n
  ann.by.level.2 <- ann.by.level$ann.by.level
  while (n>0){
    ann.by.level.3 <- treePruneSingleLayer(p2, ann.by.level.2, layer, wei)
    ann.by.level.2 <- ann.by.level.3$ann.by.level
    n = ann.by.level.3$n
  }
  return(ann.by.level.2)
}

#ann.by.level.1 <- tree.prunes.single.layer(ann.by.level, layer)
#ann.by.level.1 %>% lapply(unique)
treePrunes2 <- function(p2, ann.by.level, wei=1, start.layer.index=1) {
  #ann.by.level <- treePrunePrep(ann.by.level) #remove parent types with single subtype
  if (ann.by.level[[1]]%>%unique%>%length!=1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))%>%setNames(ann.by.level[[1]]%>%names)),ann.by.level)
  }
  
  tree.depth <- length(ann.by.level)
  for (layer in names(ann.by.level)[start.layer.index:(tree.depth-2)]) { #the LAYERS compared actually the two below the starting layer, not the starting layer and the one below it.
    ann.by.level2 <- list()
    layer.types <- ann.by.level[[layer]] %>% unique()
    for (type in layer.types) {
      type.cells <- ann.by.level[[layer]] %>% .[.==type] %>% names
      ann.by.level.sub <- list()
      for (i in 1:tree.depth){
        ann.by.level.sub[[names(ann.by.level[i])]] <- ann.by.level[[i]][type.cells]}
      ann.by.level.sub <- treePrunesSingleLayer(p2, ann.by.level.sub, layer, wei)
      if (length(ann.by.level2) == 0){
        ann.by.level2<-ann.by.level.sub}else{
        for (i in 1:tree.depth){
          ann.by.level2[[names(ann.by.level2[i])]] <- c(ann.by.level2[[i]], ann.by.level.sub[[i]])}
      }
    }
    ann.by.level <- ann.by.level2
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

#ann.by.level <- ann.by.level.l113
#ann.by.level.out <- tree.prunes(ann.by.level)


annToTreeDf <- function(ann.by.level) {
  if (ann.by.level[[1]]%>%unique%>%length>1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))
                           %>%setNames(names(ann.by.level[[1]]))), ann.by.level)
  }
  
  depth <- length(ann.by.level)
  ann.split <- list()
  for (i in 1:(depth-1)){
    ann.split[[paste0("l",i)]] <- ann.by.level[[i+1]] %>% split(ann.by.level[[i]]) %>% lapply(unique)
  }
  
  parent <- c()
  node <- c()
  path.len <- c()
  for (i in 1:length(ann.split)) {
    for (j in names(ann.split[[i]])){
      for (l in ann.split[[i]][[j]]){
        if (j != l){
          path.len <- c(path.len,i)
          parent <- c(parent, j)
          node <- c(node, l)
        }
      }
    }
  }
  
  root <- ann.split[[1]] %>% names()
  parent[parent==root] <- "root"
  
  c.df <- tibble(
    Parent = parent,
    Node = node,
    PathLen = as.numeric(path.len)
  )
  
  return(c.df)
}


getClfData <- function(p2, ann.by.level, name){
  if (ann.by.level[[1]]%>%unique%>%length!=1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]]))%>%setNames(., names(ann.by.level[[1]]))), 
                      ann.by.level)
  }
  
  ann.by.parent <-  lapply(1:(length(ann.by.level)-1), function(i) 
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>%
    Reduce(c, .) %>% .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]
  
  message("Selecting markers ...")
  de.info.per.parent <- ann.by.parent %>% 
    pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))
  
  pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
    list(
      positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
      negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
    )
  )
  
  cm.norm <- p2$counts %>% Matrix::t() %>% normalizeTfIdfWithFeatures() 
  marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
    selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
                         parent=n,max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))
  
  marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)
  
  message("Writing marker list to markup file ...")
  marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>% 
    markerListToMarkup(file=paste0("~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list_", name, ".txt"), 
                       group.by.parent=F) #I may need to update this function in the future
  
  clf.data <- getClassificationData(cm.norm, paste0("~/CellAnnotatoR-dev/notebooks/12hclust_trimming/marker_list_", name, ".txt"), 
                                    prenormalized=T)  
  return(clf.data)
}



reAnnoPip <- function(p2, ann.by.level, name){

  clf.data <- getClfData <- function(p2, ann.by.level, name)
  message("Re-annotation ...")
  ann.by.level <-
    assignCellsByScores(p2$graphs$PCA, clf.data,
                        clusters = NULL) #we had clusters=p2$clusters$PCA$leiden, but lots of tip clusters dropped after re-annotation, therefore, I reset it to NULL.
  saveRDS(ann.by.level, paste0("~/CellAnnotatoR-dev/notebooks/12hclust_trimming/ann_by_level_", name, ".rds"))
  
  message("Checking the data ...")
  ann.by.level$annotation[[length(ann.by.level$annotation)]]%>%unique() #49
  ann.by.level$annotation %>% sapply(unique) %>% Reduce(c,.) %>% unique() %>% length() #96 clusters
  
  message("Preparing file for plotting hierarchy ...")
  c.df <- annToTreeDf(ann.by.level$annotation)
  saveRDS(c.df, file=paste0("~/CellAnnotatoR-dev/notebooks/12hclust_trimming/c.df.", name, ".rds")) #For plotting hierarchy
  
  return(list(ann.by.level=ann.by.level, clf.data=clf.data))
}
