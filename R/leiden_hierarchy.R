#TODO: sort the functions
# Building leiden clustering based hirarchy start here

#Self Projection Quality Control
#annotation.before <- p4$clusters$PCA$leiden #annotation before re-annotation (reAnnoPip)
#annotation.after <- reann.p4$ann.by.level$annotation$l1 #annotation after re-annotation (reAnnoPip)
#This is used for the self projection quaility control after re-annotating the leiden-clustering-based hierarchy
selfProjQCperCluster <- function(annotation.before, annotation.after) {
  #Get self projection accuracy:
  acc <- sum(annotation.before == annotation.after)/length(annotation.before)
  #The annotation before and after need to have the same length.
  #Get confusion rate (FP/TP):
  #the cluster sample size is balanced by adjusting weights inversely proportional to class frequencies in the input data
  clusters <- unique(annotation.before)

  conf.fpr <- sapply(clusters, function(x) {
    wei.x <- length(annotation.before)/table(annotation.before)[x]
    wei.other <- 1/(1-(1/wei.x))
    positive <- annotation.before%>%.[.==x] %>% names
    negative <- names(annotation.before) %>% setdiff(positive)
    predicted <- annotation.after%>%.[.==x] %>% names
    true.positive <- positive[!is.na(match(positive, predicted))]
    false.positive <- negative[!is.na(match(negative, predicted))]
    #Balanced confusion rate
    conf <- (length(false.positive)*wei.other)/(length(true.positive)*wei.x)
    #False positive rate
    fpr <- length(false.positive)/length(negative)
    c(conf=conf, fpr=fpr)
  })

  colnames(conf.fpr) <- clusters
  rownames(conf.fpr) <- c("conf", "fpr")

  return(list(accuracy=acc, confusion.rate.false.pos.rate=conf.fpr))
}


# ann.by.level, one layer, can't be singleton!!! otherwise error: subscript out of bound
getClfDataInferringMarkers <- function(p2, ann.by.level, name, outfile.path){
  #Check singleton
  if (length(ann.by.level)==1 && (ann.by.level[[1]] %>% unique %>% length)==1) return(NULL)
  if ((ann.by.level[[1]] %>% unique %>% length) != 1) {
    ann.by.level <- list(root=rep("root", length(ann.by.level[[1]])) %>%
                    setNames(names(ann.by.level[[1]]))) %>% c(ann.by.level)
  }
  # #Check small clusters (< 6 cells)
  # if ((lapply(1:(length(ann.by.level)-1), function(i)
  #   split(ann.by.level[[i+1]], ann.by.level[[i]])) %>% Reduce(c, .) %>%
  #   .[sapply(., length) < 6] %>% length)>0) {return(NULL)}

  ann.by.parent <-  lapply(1:(length(ann.by.level)-1), function(i)
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>% Reduce(c, .) %>%
    .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]

  if ((ann.by.parent %>% .[sapply(., function(c) length(unique(c))==1)] %>% length) > 0) {return(NULL)}

  marker.list <- selectMarkersByParent(p2, ann.by.parent, z.threshold=0, append.auc=T,
        max.iters=50, max.uncertainty=0.25, max.pos.markers=2, log.step=5, verbose=1, n.cores=10)

  if (is.null(marker.list)) {return(NULL)}

  message("Check whether all clusters have markers: ")
  if ((marker.list %>% .[sapply(., function(x) length(x$expressed) == 0)] %>% length) >0){
                             return(NULL)
    }
  marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>%
    markerListToMarkup(file=paste0(outfile.path, "marker_list_", name, ".txt"), group.by.parent=F)

  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  clf.data <- getClassificationData(cm.norm, paste0(outfile.path, "marker_list_", name, ".txt"))

  return(clf.data)
}




selectMarkersByParent <- function(p2, ann.by.parent, z.threshold=0, append.auc=T, max.iters=50,
                         max.uncertainty=0.25, max.pos.markers=2, log.step=5, verbose=1, n.cores=10){
  message("Selecting markers for clusters under the same parents ...")
  de.info.per.parent <- ann.by.parent %>%
    pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))

  # if (sum(de.info.per.parent %>% sapply(., function(c) c%>%length<2)) > 0) {
  #   return(NULL)
  # }

  pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
    list(positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
         negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
    )
  )
  #TODO: Error: (list) object cannot be coerced to type 'double' FIX THE BUG
  if (sum(pre.selected.markers %>% sapply(., function(c) c$positive%>%lapply(length))==0) > 0) {
    return(NULL)
  }
  message("Debug: Error in rowSums(.) : 'x' must be an array of at least two dimensions")
  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
    selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
    parent=n, max.iters=50, max.uncertainty=0.25, max.pos.markers=2, log.step=5, verbose=1, n.cores=10))

  marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)

  return(marker.list)
}



#Outfile.path: where to save the marker file
reAnnoPip <- function(p2, ann.by.level, name, outfile.path, graph=NULL, clf.data=NULL,
          clusters = NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
  #Check small clusters (< 6 cells)
  if (sum((ann.by.level[[1]] %>% table) < 6) > 0) {return(NULL)} #only one layer allowed
  #Check singletons
  if ((ann.by.level[[1]] %>% unique %>% length)==1) {return(NULL)}

  if (is.null(clf.data)) {clf.data <- getClfDataInferringMarkers(p2, ann.by.level, name, outfile.path)}

  if (is.null(clf.data) || (clf.data$marker.list %>% length)==1) {return(NULL)} #Check whether the new clf.data is still NULL(i.e. some clusters have no markers).

  message("Re-annotation ...")
  ann.by.level <-
    assignCellsByScores(graph=NULL, clf.data=clf.data, score.info=NULL,
                                    clusters=NULL, verbose=0, uncertainty.thresholds, max.depth=NULL)

  return(list(ann.by.level=ann.by.level, clf.data=clf.data))
}




###########################################
reAnnotationQC <- function(reann, p, sp.thresholds){
  if (!is.null(reann)){#cluster size before reann>5 cells, all clusters have markers
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){#no cluster dropping during reannotation (indicating bad markers)
      if (sum(sp$confusion.rate.false.pos.rate[1,] > sp.thresholds$confusion.rate)==0 &&
          sum((reann$ann.by.level$annotation$l1 %>% table)<6)==0){#confusion rate not above the threshold
        passQC <- TRUE
        return(passQC)
      }
    }
  }
  passQC <- FALSE
  return(passQC)
}


#direction=c("left", "right"),
getBigClusters <- function(p.left, p.right, p.right.old, res.left, res.right, res.right.old,
            reann.right.old, out.name, outfile.path, graph=NULL, clf.data=NULL, clusters=NULL,
            uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), type="PCA",
                           method=conos::leiden.community, n.iterations=50, name="leiden"){
  # if (direction=="left"){
  n = 0
  while (sum((p.right$clusters$PCA$leiden %>% table) < 6) > 0){
    if ((res.right-res.left) <= 0.05){
      break
    }
    n = n+1
    p.right.old <- p.right
    p.right$getKnnClusters(type="PCA", method=conos::leiden.community,#decrease res.right, get new p.right
                     resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
    res.right.old = res.right
    res.right = res.left+(res.right-res.left)/2
  }
  if (n>0){#update the reann.right before decreasing res
    if (length(p.right.old$clusters$PCA$leiden %>% table)==1){#Singleton can't be brought in, but may be produced inside.
      reann.right.old <- list(ann.by.level=list(annotation=list(l1=p.right.old$clusters$PCA$leiden)),
                         clf.data=list(marker.list=list()))
    } else{
        reann.right.old <- reAnnoPip(p.right.old, list(l=p.right.old$clusters$PCA$leiden),
                           out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
    }
  }
  return(list(p.left=p.left, p.right=p.right, p.right.old=p.right.old, res.left=res.left,
                      res.right=res.right, res.right.old=res.right.old, reann.right.old))
  # }
#
#   if (direction=="right"){
#     while (sum((p$clusters$PCA$leiden %>% table) < 6) > 1){
#       p$getKnnClusters(type="PCA", method=conos::leiden.community,
#                        resolution=res.right-(res.right-res.left)/2, n.iterations=50, name="leiden")
#       res.left.old = res.left
#       res.lef = res.right-(res.right-res.left)/2
#     }
#     return(list(p=p, res.left=res.left, res.left.old=res.left.old))
#   }
}



#Get the first layer of leiden clusters (>1 clusters) for the leiden clustering based hierarchy  #TESTED!
#p2: dataset in pagoda form.
getLayer1KnnClusters <- function(p2, res.start=0.01, res.step=0.01){
  p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start,
                    n.iterations=50, name="leiden")
  clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  n = 0
  while (clusters.n==1){
    n = n+1
    p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start+n*res.step,
                      n.iterations=50, name="leiden")
    clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  }
  return(list(p2=p2, res=res.start+n*res.step))
}
#e.g. l2 <- getLayer1Clusters(p2, res.step=0.01)
#l2$res
#l2$p2$clusters$PCA$leiden %>% table



#To solve the problem that the given res.max is too small.
#If KNN clustering gets a singleton,
#not sure whether it is a real singleton or just because the resolution is too low,
#we try to increase the res to get a non-singleton clustering,
#and see wheteher the new clustering fullfil requirements below or not later on.
#reann.left: the reannotation result at resolution, res.min
tailorResRange <- function(p.left, reann.left, res.min, res.max, res.switch, out.name, outfile.path,
                           graph=NULL, res.max.update=1, clf.data=NULL, clusters = NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  message("sp.thresholds ", sp.thresholds)
  if ((res.max-res.min) <= res.switch) {res.max = (res.max+1)}

  check.singleton <- getLayer1KnnClusters(p.left, res.start=res.max, res.step=1)#Increase res until not singleton
  p.right <- check.singleton$p2
  res.max <- check.singleton$res #Update max.res.middle

  if (sum((p.right$clusters$PCA$leiden %>% table) < 6) > 0){ #All clusters should have >= 6 cells, otherwise not pass
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left,
                reann.right=NULL, res.min=res.min, res.max=res.max))
  }

  reann.right <- reAnnoPip(p.right, list(l=p.right$clusters$PCA$leiden), out.name, outfile.path,
                           graph, clf.data, clusters, uncertainty.thresholds)

  if  (!is.null(reann.right)){
    sp <- selfProjQCperCluster(p.right$clusters$PCA$leiden, reann.right$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
      if (sum((reann.right$ann.by.level$annotation$l1 %>% table)<6)==0 &&
          sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
        res.min <- res.max
        p.left <- p.right
        reann.left <- reann.right
        res.max <- res.max+res.max.update
        return(tailorResRange(p.left, reann.left, res.min, res.max, res.switch,
                              out.name, outfile.path, graph, res.max.update, clf.data,
                              clusters, uncertainty.thresholds, sp.thresholds))
      }
    }
  } #At res.max, the clustering shouldn't pass quality control
  return(list(p.left=p.left, p.right=p.right, reann.left=reann.left,
              reann.right=reann.right, res.min=res.min, res.max=res.max))
}

################################

#####Main function 2
#res.max.update must be larger than res.switch, otherwise, you will be in trouble for the first step of findBestRes
getNextLayerKnnClusters <- function(p2, annotation, out.name, outfile.path, min.res.start=0,
                           graph=NULL, res.max.update=1, clf.data=NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           max.res.middle=1, res.switch=0.1, type="PCA",
                           method=conos::leiden.community, n.iterations=10, name="leiden",
                           sp.thresholds=list(accuracy=0.85, confusion.rate=3/17),
                           clustering.type=NULL, embeding.type=NULL){
  #TODO: message to be removed
  message("sp.thresholds", sp.thresholds)
  ann.by.parents <- list()
  marker.list <- list()
  #resol.last.layer <- list()
  ann.by.cluster <- split(annotation, annotation)
  for (cl in names(ann.by.cluster)){
    message("Now we will get next layer for cluster ", cl)
    p <- p2$misc$rawCounts[ann.by.cluster[cl][[1]]%>%names,] %>% Matrix::t() %>%
      vpscutils::GetPagoda(clustering.type, embeding.type, n.cores=1) #Withdraw data subset for cl.
    message("Get Knn clusters at res.min, which should pass QC")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=min.res.start,
                     n.iterations=10, name="leiden") #Knn clustering
    #Because min.res.start=0, so the Knn cluster must be a singleton that can't do reannotation.
    p.left <- p #has the clustering at res.min
    reann.left <- list(ann.by.level=list(annotation=list(l1=p.left$clusters$PCA$leiden)),
                  clf.data=list(marker.list=list())) #Build reann.left manually
    #Tailor the max.res (and maybe min.res as well) for cl:
    #p input here doesn't matter
    new.res.range <- tailorResRange(p.left, reann.left, min.res.start, max.res.middle, res.switch, out.name,
          outfile.path, graph, res.max.update, clf.data, clusters, uncertainty.thresholds, sp.thresholds)
    p.left <- new.res.range$p.left
    p.right <- new.res.range$p.right #Match reann, not reann.left
    reann.left <- new.res.range$reann.left #Must pass self projection QC, #match p.left
    reann.right <- new.res.range$reann.right #Can't pass self projection QC
    res.left <- new.res.range$res.min
    res.right <- new.res.range$res.max

    best.res <- findBestRes(p.left, p.right, reann.left, reann.right, res.left, res.right, out.name,
             outfile.path, res.switch, type="PCA",sp.thresholds, method, n.iterations, name, graph,
             clf.data, uncertainty.thresholds)
    reann.left <- best.res$reann.left
    #TODO: to be removed
    #res.right <- best.res$res.right
    message("Now we got next layer for ", cl)

    if ((reann.left$ann.by.level$annotation$l1 %>% table %>% length) > 1){#Do not want to export singletons
    #if ((reann.left$clf.data$marker.list %>% length) > 0){
      annotation <- reann.left$ann.by.level$annotation$l1 %>% sapply(function(n) paste0(cl,"_", n)) %>%
                                              setNames(reann.left$ann.by.level$annotation$l1 %>% names)
      ann.by.parents <- c(ann.by.parents, list(annotation) %>% setNames(cl))
      marker.list.cl <- reann.left$clf.data$marker.list %>% lapply(function(sub.clust) {sub.clust$"parent"=cl
                                                                                        sub.clust})
      marker.list <- c(marker.list, list(marker.list.cl) %>% setNames(cl))
      #marker.list <- c(marker.list, list(reann.left$clf.data$marker.list) %>% setNames(cl))
      #TODO: to be removed
      #resol.last.layer <- c(resol.last.layer, list(res.right) %>% setNames(cl))
    }
  }
  return(list(ann.by.parents=ann.by.parents, marker.list=marker.list))
  #TODO: to be removed
  #return(list(ann.by.parents=ann.by.parents, marker.list=marker.list, resol.last.layer=resol.last.layer))
}


#####Main function 3
#Look for the maximum resolution (e.g. the maximum number of predicted clusters) that have its predicted clusters all pass quality control (self projection accuracy via re-annotation: confusion rate)
#TODO: ADD reann.left
#res.switch: only affect the starting of findBestRes, or the switch between findBestRes and findSmallestRes, when finding the max res that fullfill all the classifyign requirements.
findBestRes <- function(p.left, p.right, reann.left, reann.right, res.left, res.right, out.name, outfile.path,
               res.switch=0.1, type="PCA", sp.thresholds=list(accuracy=0.85, confusion.rate=3/17),
               method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
               clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
  message("sp.thresholds, res.switch ", sp.thresholds, " ", "res.switch")
  if ((res.right-res.left)<=res.switch) {
    return(list(p.left=p.left, reann.left=reann.left, res.left=res.left, res.right=res.right))
  }
  message("Now, we are heading to the left for the next smallest resolution ...")
  to.left <- findNextSmallestRes(p.left, p.right, p.right.old=p.right,
             reann.left, reann.right, reann.right.old=reann.right, res.left, res.right,
             res.right.old=res.right, out.name, outfile.path, res.switch, type="PCA", method,
             n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds)
  p.left <-  to.left$p.left  #pass QC
  p.right <- to.left$p.right    #not pass QC
  reann.left <- to.left$reann.left #This must be a good one. Pass QC
  reann.right <- to.left$reann.right   #not pass QC
  res.left <- to.left$res.left
  res.right <- to.left$res.right
  if (length(reann.left$ann.by.level$annotation$l1 %>% table)==1) {
    return(list(p.left=p.left, reann.left=reann.left, res.left=res.left, res.right=res.right))
    }
  if ((res.right-res.left)>res.switch){
    message("Now, we are heading to the right for the next biggest resolution ...")
    to.right <- findNextBiggestRes(p.left, p.left.old=p.left, p.right, reann.left,
                reann.left.old=reann.left, reann.right, res.left, res.left.old=res.left,
                res.right, out.name, outfile.path, res.switch, type="PCA", method, n.iterations,
                name, graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds)
    p.left <- to.right$p.left
    p.right <- to.right$p.right
    res.left <- to.right$res.left #the highest res that work
    res.right <- to.right$res.right  #the lowest res that doesn't work
    reann.left <- to.right$reann.left
    reann.right <- to.right$reann.right #In case we want to turn left again
  }
  return(findBestRes(p.left, p.right, reann.left, reann.right, res.left, res.right, out.name, outfile.path, res.switch, type, sp.thresholds,
                     method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
                     clf.data=NULL, uncertainty.thresholds))
}


#####Main function 4
#To left
#There will be problem for this function, if the smallest res is very close (diff < 0.0001) to the original res.left, then the run may break.
#The error will look like: Error in ann.by.level[[i + 1]] : subscript out of bounds (reAnnPip: ann.by.level only has one cluster! singleton)
#Possible solution: increase relax the sp.threshold
findNextSmallestRes <- function(p.left, p.right, p.right.old, reann.left, reann.right, reann.right.old,
            res.left, res.right,res.right.old, out.name, outfile.path, res.switch=0.1, type="PCA",
            method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL,
            clusters = NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
            sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  #TODO: message to be removed.
  message("res.right, res.left, res.switch, sp.thresholds ", res.right, res.left, res.switch,
          " ", sp.thresholds)

  #TODO: message to be removed
  message("reAnnotation QC ...")
  #reann and p should match.
  pass.reAnnotation.QC <- reAnnotationQC(reann.right, p.right, sp.thresholds)

  if (pass.reAnnotation.QC){ #We look to the left for a res at which clustering pass QC #Pass
    return(list(p.left=p.right, p.right=p.right.old, reann.left=reann.right,
                reann.right=reann.right.old, res.left=res.right, res.right=res.right.old))
  }

  if ((res.right-res.left)<=res.switch && !pass.reAnnotation.QC){ #Stop #This can't be introduced!
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left, reann.right=reann.right,
                res.left=res.left, res.right=res.right))
  }

  message("Now we will decrease resolution to next middle ...")
  #Decrease res.right until no cluster has less than 6 cells..
  p.right.old = p.right
  p.right$getKnnClusters(type="PCA", method=conos::leiden.community,
                   resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
  res.right.old = res.right
  res.right = res.left+(res.right-res.left)/2

  reann.right.old <- reann.right
  #Make sure all cluster >5 cells
  big.clusters <- getBigClusters(p.left, p.right, p.right.old, res.left, res.right, res.right.old,
            reann.right.old, out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds,
            type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden")
  p.left <- big.clusters$p.left
  p.right <- big.clusters$p.right
  p.right.old <- big.clusters$p.right.old
  res.left <- big.clusters$res.left
  res.right <- big.clusters$res.right
  res.right.old <- big.clusters$res.right.old
  reann.right.old <- big.clusters$reann.right.old

  if (length(p.right$clusters$PCA$leiden %>% table) == 1){#Singleton can't do reannotation. #Considered as pass
    reann.right <- list(ann.by.level=list(annotation=list(l1=p.right$clusters$PCA$leiden)),#Build reann manually.
             clf.data=list(marker.list=list())) #Singleton's reann should pass all self projection QC.
    return(list(p.left=p.right, p.right=p.right.old, reann.left=reann.right, reann.right=reann.right.old,
                res.left=res.right, res.right=res.right.old))
  }

  # if (sum((p$clusters$PCA$leiden %>% table) < 6) > 1){ #All clusters should >= 6 cells, otherwise not pass
  #   return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
  #               res.left=res.left, res.left.old=res.left.old, res.right=res.right))
  # }
  reann.right <- reAnnoPip(p.right, list(l=p.right$clusters$PCA$leiden),
                     out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
  message("Now the new res.right and res.left become ", res.right, " ", res.left)

  return(findNextSmallestRes(p.left, p.right, p.right.old, reann.left, reann.right, reann.right.old,
         res.left, res.right, res.right.old, out.name, outfile.path, res.switch, type="PCA", method,
         n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds))
}



#####Main function 5
#The starting reann should fufill the 3 if conditions in the function!!!
#Look to the right for the maximum separable resolution (i.e. the 1st one that doesn't pass)
#i.e. find the first resolution, at which the clustering results does not pass self projection QC.
#reann.old: no matter the initial input or the generated one from below will all be fine.
findNextBiggestRes <- function(p.left, p.left.old=p.left, p.right, reann.left, reann.left.old=reann.left,
                           reann.right, res.left, res.left.old=res.left, res.right, out.name, outfile.path,
                           res.switch=0.1, type="PCA", method=conos::leiden.community, n.iterations=50,
                           name="leiden", graph=NULL, clf.data= NULL, clusters = NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  #TODO: message to be removed...
  message("1.res.switch, res.right, res.left, sp.thresholds ", res.switch, res.right, res.left, " ", sp.thresholds)

  pass.reAnnotation.QC <- reAnnotationQC(reann.left, p.left, sp.thresholds)

  if ((res.right-res.left)<=res.switch && pass.reAnnotation.QC){
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left, reann.right=reann.right,
                res.left=res.left, res.right=res.right))}

  if (pass.reAnnotation.QC){
    message("Now we will increase resolution to next middle ...")
    p.left.old <- p.left #Pass QC
    reann.left.old <- reann.left #Pass QC
    #Increase res.left until no cluster has less than 6 cells
    p.left$getKnnClusters(type="PCA", method=conos::leiden.community, #get new p.left clustering
                     resolution=res.right-(res.right-res.left)/2, n.iterations=50, name="leiden")
    res.left.old = res.left
    res.left = res.right-(res.right-res.left)/2
    # big.clusters <- getBigClusters(p, res.left, res.right, direction="right", #All clusters >= 6 cells
    #                         type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden")
    # p <- big.clusters$p
    # res.left <- big.clusters$res.left
    # res.left.old <- big.clusters$res.left.old
    # if (sum((p$clusters$PCA$leiden %>% table) < 6) > 0){ #All clusters should >= 6 cells, otherwise not pass
    #   return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
    #               res.left=res.left, res.left.old=res.left.old, res.right=res.right))
    # }

    if (length(p.left$clusters$PCA$leiden %>% table)==1){#Singleton can't be brought in, but may be produced inside.
      reann.left <- list(ann.by.level=list(annotation=list(l1=p.left$clusters$PCA$leiden)), #Counted as pass
                    clf.data=list(marker.list=list()))
      return(list(p.left=p.left, p.right=p.right, reann.left=reann.left, reann.right=reann.right,
                  res.left=res.left, res.right=res.right))
    }#reAnnPip will output NULL on singletons, one has to create reann manually.

    reann.left <- reAnnoPip(p.left, list(l=p.left$clusters$PCA$leiden),
                       out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
    message("Now the new res.left and res.right become ", res.left, " ", res.right)

    return(findNextBiggestRes(p.left, p.left.old, p.right, reann.left, reann.left.old, reann.right,
                              res.left, res.left.old, res.right, out.name, outfile.path, res.switch,
                              type="PCA", method, n.iterationss, name,
                              graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds))
  }

  return(list(p.left=p.left.old, p.right=p.left, reann.left=reann.left.old, reann.right=reann.left,
              res.left=res.left.old, res.right=res.left))
}


#####Main function 6
getNextLayersKnnClusters <- function(p2, annotation, out.name, outfile.path, layer.n,
                                     min.res.start=0, graph=NULL, res.max.update=1, clf.data=NULL,
                                     uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                                     max.res.middle=1, max.res.increase=1, res.switch=0.05, type="PCA",
                                     method=conos::leiden.community, n.iterations=50, name="leiden",
                                     sp.thresholds=list(accuracy=0.8, confusion.rate=0.25),
                                     clustering.type=NULL, embeding.type=NULL){
  next.layers <- list()
  for (l in 2:layer.n){
    if (length(annotation) > 0){
      message("Now we are going to get layer ", l, "...")
      next.layer <- getNextLayerKnnClusters(p2, annotation, out.name, outfile.path, min.res.start, graph,
                    res.max.update, clf.data, uncertainty.thresholds, max.res.middle, res.switch, type="PCA",
                    method, n.iterations, name, sp.thresholds, clustering.type, embeding.type)
      saveRDS(next.layer, file="~/CellAnnotatoR-clean/notebooks/25files/next.layer.rds")

      if (length(next.layer$ann.by.parents)==0) {annotation=c()} else{
        annotation <- next.layer$ann.by.parents %>% Reduce(c,.)
        annotation <- annotation[annotation %in% (table(annotation) %>% .[.>=12] %>% names)]
        #Collect all non-singleton clusters
        # annotation <- split(next.layer$ann.by.parents %>% Reduce(c,.),
        #                     annotation[names(next.layer$ann.by.parents %>% Reduce(c,.))]) %>%
        #                     .[sapply(.,function(k) (table(k)%>%length)>1)] %>% Reduce(c,.)
        # annotation <- annotation[annotation %in% (table(annotation) %>% .[.>=12] %>% names)]
      }
      max.res.middle <- max.res.middle + max.res.increase
      #max.res.middle <- (next.layer$resol.last.layer %>% sapply(function(n) n[1]) %>% max) + max.res.increase
      message("max.res.increase: ", max.res.increase)
      next.layers <- c(next.layers, list(next.layer) %>% setNames(paste0("l", l)))
    }
  }
  return(next.layers)
}

#####Main function 7
gatherAnnByLevelMarkerList <- function(reann.layer1, ann.layer1, ann.next.layers){
  #Gather ann.by.level
  ann.by.level <- list()
  marker.list <- list()
  ann.by.level <- c(ann.by.level, list(ann.layer1)%>%setNames("l1"))
  #Take out ann.by.parents from other layers
  ann.by.level.r <- ann.next.layers%>% lapply(function(n) n$ann.by.parents %>% Reduce(c,.))
  ann.by.level <- c(ann.by.level, ann.by.level.r)
  message("Make all layers have the same cell order ...")
  for (i in 1:(length(ann.by.level)-1)){#In order to use split, all layers should have the same cell order.
    if (length(ann.by.level[[i+1]]) != length(ann.by.level[[i]])){
      ann.by.level[[i+1]] <- c(ann.by.level[[i+1]],
            ann.by.level[[i]][names(ann.by.level[[i]]) %>% setdiff(names(ann.by.level[[i+1]]))])
    }
    ann.by.level[[i+1]] %<>% .[ann.by.level[[i]] %>% names]
  }
  #Check if the last two layers are the same, if yes, remove the last layer
  depth <- length(ann.by.level)
  if (sum(ann.by.level[[depth]]==ann.by.level[[depth-1]])==length(ann.by.level[[depth]])) {
    ann.by.level[depth] = NULL
  }
  #Gather marker.list,
  # reann.layer1$clf.data$marker.list %<>% lapply(function(t) {
  #   t$parent="l"
  #   t
  # })
  marker.list <- c(marker.list, list(l1=list(reann.layer1$clf.data$marker.list)%>%setNames("root"))) # For layer 1
  marker.list <- c(marker.list, ann.next.layers%>% lapply(function(n) n$marker.list) %>%
                     setNames(names(ann.next.layers))) #For the rest layers
  return(list(ann.by.level=ann.by.level, marker.list=marker.list))
}


#####Main function8:
#Get the leiden clustring based hierarchy:
#At the same time, we will also do re-annotation, so the outcome will be a ready hierarchy with re-annotation based on the selected markers
#layer.n: the number of layers that the users want to get. The more layers, the more clusters one will get, though one will not get endless number of clusters, becasue a. the feature genes will eventually become exhausted, b. the size of cluster can't be smaller than 6, or c. the cells within clusters become homogeneous eventually.
# If encountering error "Error in ann.by.level[[i + 1]] : subscript out of bounds", that means when running findBiggestRes(...), res.right approaches res.left too closely. Solution: set min.res.update=FALSE, and at the same time set res.switch=0
#res.max.update must be larger than res.switch, otherwise, you will be in trouble for the first step of findBestRes
##2. Self-projection quality control for re-annotation (one way to predict cell types based on markers)
#1. reann shouldn't be NA (as long as not all clusters have predicted markers, then reann will be NA)
#2. sp shouldn't have NaN or Inf
#3. All the clusters should hava a minimum size of 6 after re-annotation
#4. The self projection accuracy by re-annotation has to be >=0.6 (or 0.8) and the confusion rate has to be <=2/3 (0.25)
getLeidenHierarchy <- function(p2, out.name, outfile.path, layer.n=3, res.step.layer1=0.01, min.res=0,
                               max.res.layer2=1, max.res.increase=1, res.switch=0.05, clustering.type=NULL,
                               embeding.type=NULL, res.max.update=1, clf.data=NULL, clusters = NULL,
                               uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                               sp.thresholds=list(accuracy=0.8, confusion.rate=0.25), graph=NULL,
                               type="PCA", method=conos::leiden.community, n.iterations=50,
                               name="leiden"){
  if (res.max.update <= res.switch){
    stop("res.max.update must be larger than res.switch")}
  message("Get Knn clusters for layer 1 ...")
  #Get a big picutre in layer1, we will accept the first non-singletong KNN clustering result here.
  layer1 <- getLayer1KnnClusters(p2, res.start=0.01, res.step=res.step.layer1)
  p2 <- layer1$p2
  reann.layer1 <- reAnnoPip(p2, list(l1=p2$clusters$PCA$leiden), out.name, outfile.path, graph,
                            clf.data, clusters=NULL, uncertainty.thresholds) #Not singleton!

  message("Get Knn clusters for the rest layers ...")
  pass.reAnnotation.QC <- reAnnotationQC(reann.layer1, p2, sp.thresholds)
  if (pass.reAnnotation.QC){
    annotation.layer1 <- reann.layer1$ann.by.level$annotation$l1 %>% sapply(function(n) paste0("l_",n)) %>%
                         setNames(reann.layer1$ann.by.level$annotation$l1 %>% names())
    next.layers <- getNextLayersKnnClusters(p2, annotation=annotation.layer1, out.name, outfile.path,
                 layer.n, min.res.start=min.res, graph, res.max.update, clf.data, uncertainty.thresholds,
                 max.res.middle=max.res.layer2, max.res.increase, res.switch, type, method, n.iterations, name,
                 sp.thresholds, clustering.type, embeding.type)
  } else {stop("There is no subclusters in the cell population, or adjust sp.thresholds (confusion.rate = 1/accuracy -1), or try more gene features...")}

  return(gatherAnnByLevelMarkerList(reann.layer1, annotation.layer1, next.layers))
}


# #Check if the last 2 layers have the same clustering, if it is, drop the last layer
# if (length(ann.by.level %>%.[[length(.)]]%>%unique)==length(ann.by.level%>%.[[length(.)-1]]%>%unique)){
#   ann.by.level[length(ann.by.level)] <- NULL}






#For plotting hierarchy:
annToTreeDf <- function(ann.by.level) {
  if (ann.by.level[[1]]%>%unique%>%length>1){
    ann.by.level <- c(list(root=rep("root", length(ann.by.level[[1]])) %>%
                           setNames(names(ann.by.level[[1]]))), ann.by.level)
  }

  depth <- length(ann.by.level)
  ann.split <- list()
  for (i in 1:(depth-1)){
    ann.split[[paste0("l",i)]] <- ann.by.level[[i+1]] %>%
      split(ann.by.level[[i]]) %>% lapply(unique)
  }

  parent <- c()
  node <- c()
  path.len <- c()
  for (i in 1:length(ann.split)){
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

  c.df <- tibble( Parent = parent, Node = node, PathLen = as.numeric(path.len))

  return(c.df)
}


#' Add color argument
#' @param classificaiton should be either classificaiton.tree or a classification dataframe
#' @export
plotTypeHierarchyLi <- function(classification, layout="slanted",
                                xlims=NULL, font.size=3, col="black", ...) { #modified from Viktor's function
  if (!requireNamespace("ggtree", quietly=T))
    stop("You need to install package 'ggtree' to be able to plot hierarchical tree. ",
         "`Try devtools::install_github('YuLab-SMU/ggtree')`")
  if (!is.data.frame(classification)){
    c.df <- classificationTreeToDf(classification)
  } else{
    c.df <- classification
    cg <- c.df %$% data.frame(parent=Parent, node=Node) %>% ape::as.phylo()
    cg$edge.length <- rep(1, nrow(cg$edge))

    if (is.null(xlims)) {
      xlims <- c(0, max(c.df$PathLen) + 0.5)
    }

    ggtree::ggtree(cg, layout = layout, color=col,  ...) +
      ggtree::geom_rootpoint() +
      ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]),
                          size=font.size, color=col) + ggplot2::xlim(xlims)
  }
}


#' Prepare the "trait" file for plotUncHierarchy() from the positive, negative and coverage uncertainty data (output from scoreCellUncertaintyPerLevel())
#' @param unc.per.cell a list, each element represents the annotation of all cells at on hierarchical level
#' @param ann.by.level a list of lists, each main element is one hierarchical level; at each hierarchical level, there is a sublist, which contains 3 elements: positive, negative and coverage uncertainty.
#' @export ann.unc a dataframe, there are 4 columns:"type" #cell types   "positive" "negative" "coverage" #uncertainty
uncToTreeTrait <- function(unc.per.cell, ann.by.level) {
  ann.unc <- list()
  for (i in names(unc.per.cell) %>% setNames(.,.)) {
    ann.unc[[i]] <- merge(as.data.frame(ann.by.level$annotation[[i]]),
                          as.data.frame(unc.per.cell[[i]]), by="row.names", all=TRUE)
    names(ann.unc[[i]])[2] <- "type"
    ann.unc[[i]] <- aggregate(ann.unc[[i]][,3:5], ann.unc[[i]][,2,drop=FALSE], mean)
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


#Plot uncertainty on trees
#' Plot uncertainty (positive, negative, coverage) on hierarchical trees
#' @param c.data output file from re-annotation (by assignCellsByScores())
#' @param trait.file contain the mean positive, negative, coverage uncertainties for each cell type (including root, node and tips)
#' @param unc one of c("positive", "negative", "coverage")
#' @export
plotUcHierarchy <- function(c.data, trait.file, layout="slanted",
                            xlims=NULL, font.size=3, col.range=c(0,1), ...) {
  if (!requireNamespace("ggtree", quietly=T))
    stop("You need to install package 'ggtree' to be able to plot hierarchical tree. ",
         "`Try devtools::install_github('YuLab-SMU/ggtree')`")

  c.df <- classificationTreeToDf(c.data$classification.tree)
  cg <- c.df %$% data.frame(parent=Parent, node=Node) %>% ape::as.phylo()
  #as.phylo(): we got edge, Nnode, node.label, tip.label (c.df only has parent and node)
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

  pos.unc <- ggtree::ggtree(cgtree, aes(color=positive), layout = layout, ...) +
    scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"), limits=col.range)+
    ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
    ggplot2::xlim(xlims) #the presence or absence of [node] makes no difference here
  neg.unc <- ggtree::ggtree(cgtree, aes(color=negative), layout = layout, ...) +
    scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"), limits=col.range)+
    ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
    ggplot2::xlim(xlims)
  cov.unc <- ggtree::ggtree(cgtree, aes(color=coverage), layout = layout, ...) +
    scale_color_gradientn(colours=c("dark blue", "dark green", "dark orange", "red"), limits=col.range)+
    ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size) +
    ggplot2::xlim(xlims)

  cowplot::plot_grid(pos.unc, neg.unc, cov.unc, nrow=1)
}


#output from getLeidenHierarchy
getUncTree <- function(ann.by.level, marker.list){
  #Get clf.data
  m <- Reduce(c, marker.list) %>% Reduce(c,.)
  marker.list <- list()
  for (t in 1:length(m)) {
    if (m[[t]]$parent=="root"){
      marker.list <- c(marker.list, m[t] %>% setNames(paste0("l_", names(m[t]))))
    } else{
      marker.list <- c(marker.list, m[t] %>% setNames(paste0(m[[t]]$parent, "_", names(m[t]))))
    }
  }

  # marker.list <- 1:length(m) %>% lapply(function(n) m[n] %>%
  #                   setNames(paste0(m[[n]]$parent, "_", names(m[n])))) %>% Reduce(c,.)
  clf.data <- list()
  clf.data$marker.list <- marker.list
  clf.data$cm.norm <- p2$counts# %>% Matrix::t()
  clf.data$classification.tree <- createClassificationTree(marker.list)
  #Get score info
  #marker.score.per.cluster <- getMarkerScoresPerCellType(clf.data, score.info=NULL, aggr=T,
                                                         #underParent=T)
  score.info<- getMarkerScoreInfo(clf.data)
  score.info.per.cluster <- getMarkerScoresPerCellType(clf.data, score.info)
  score.info.by.layer <-
    ann.by.level %>% lapply(unique) %>% lapply(function(n)
      score.info.per.cluster[, n])
  #Get uncertainty
  unc.info <-
    names(ann.by.level)[1:length(ann.by.level)] %>% setNames(., .) %>%
    plapply(function(n) scorePerCellUncertainty(ann.by.level[[n]],
                                      score.info.by.layer[[n]], score.info), verbose=T)

  #Prepare the uncertainty data (as "trait")
  unc.trait <- uncToTreeTrait(unc.info, ann.by.level)
  #saveRDS(unc.trait, file="~/CellAnnotatoR-clean/notebooks/26files/unc.trait.rds")

  #Plot uncertainty
  plotUncHierarchy(clf.data, unc.trait)
}


