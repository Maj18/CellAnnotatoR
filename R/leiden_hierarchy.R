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



getClfDataInferringMarkers <- function(p2, ann.by.level, name, outfile.path){
  if ((ann.by.level[[1]] %>% unique %>% length) != 1) {#here, ann.by.level can't be singleton!!!
    ann.by.level <- list(root=rep("root", length(ann.by.level[[1]])) %>%
                    setNames(names(ann.by.level[[1]]))) %>% c(ann.by.level)
  }
  ann.by.parent <-  lapply(1:(length(ann.by.level)-1), function(i)
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>% Reduce(c, .) %>%
    .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]

  marker.list <- selectMarkersByParent(p2, ann.by.parent, z.threshold=0, append.auc=T,
                 max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10)

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
                                  max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10){
  message("Selecting markers for clusters under the same parents ...")
  de.info.per.parent <- ann.by.parent %>%
    pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))

  pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
    list(positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
         negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
    )
  )

  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
    selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
               parent=n, max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))

  marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)

  return(marker.list)
}



#Outfile.path: where to save the marker file
#TODO: CHECK THIS FUNCITON WHY THE OUTPUT IS NULL??
reAnnoPip <- function(p2, ann.by.level, name, outfile.path, graph=NULL, clf.data=NULL,
          clusters = NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
  if (is.null(clf.data)) {clf.data <- getClfDataInferringMarkers(p2, ann.by.level, name, outfile.path)}

  if (is.null(clf.data)) {return(NULL)} #Check whether the new clf.data is still NULL.

  message("Re-annotation ...")
  ann.by.level <-
    assignCellsByScores(graph=NULL, clf.data=clf.data, score.info=NULL,
                                    clusters=NULL, verbose=0, uncertainty.thresholds, max.depth=NULL)

  return(list(ann.by.level=ann.by.level, clf.data=clf.data))
}




###########################################
reAnnotationQC <- function(reann, p, sp.thresholds){
  if (!is.null(reann)){
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
      if (sum(sp$confusion.rate.false.pos.rate[1,] > sp.thresholds$confusion.rate)==0 &&
          sum((reann$ann.by.level$annotation$l1 %>% table)<6)==0){
        passQC <- TRUE
        return(passQC)
      }
    }
  }
  passQC <- FALSE
  return(passQC)
}

#direction=c("left", "right"),
getBigClusters <- function(p, res.left, res.right, res.right.old, type="PCA",
                           method=conos::leiden.community, n.iterations=50, name="leiden"){
  # if (direction=="left"){
  while (sum((p$clusters$PCA$leiden %>% table) < 6) > 0){
    p$getKnnClusters(type="PCA", method=conos::leiden.community,
                     resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
    res.right.old = res.right
    res.right = res.left+(res.right-res.left)/2
  }
  return(list(p=p, res.right=res.right, res.right.old=res.right.old))
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
tailorResRange <- function(p, reann.left, res.min, res.max, res.switch, out.name, outfile.path,
                           graph=NULL, res.max.update=1, clf.data=NULL, clusters = NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  message("sp.thresholds ", sp.thresholds)
  if ((res.max-res.min) <= res.switch) {res.max = (res.max+1)}

  check.singleton <- getLayer1KnnClusters(p, res.start=res.max, res.step=1)#Increase res until no singleton
  p <- check.singleton$p2
  res.max <- check.singleton$res #Update max.res.middle

  reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph, clf.data,
                     clusters, uncertainty.thresholds)

  if  (!is.null(reann)){
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
      if (sum((reann$ann.by.level$annotation$l1 %>% table)<6)==0 &&
          sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
        res.min <- res.max
        res.max <- res.max+res.max.update
        return(tailorResRange(p, reann.left=reann, res.min, res.max, res.switch,
                              out.name, outfile.path, graph, res.max.update, clf.data,
                              clusters, uncertainty.thresholds, sp.thresholds))
      }
    }
  } #At res.max, the clustering shouldn't pass quality control
  return(list(p=p, reann.left=reann.left, reann=reann, res.min=res.min, res.max=res.max))
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
    p.left <- p
    reann.left <- list(ann.by.level=list(annotation=list(l1=p$clusters$PCA$leiden)),
                  clf.data=list(marker.list=list())) #Build reann.left manually
    #Tailor the max.res (and maybe min.res as well) for cl:
    #p input here doesn't matter
    new.res.range <- tailorResRange(p, reann.left, min.res.start, max.res.middle, res.switch, out.name,
          outfile.path, graph, res.max.update, clf.data, clusters, uncertainty.thresholds, sp.thresholds)
    p <- new.res.range$p #Match reann, not reann.left
    reann <- new.res.range$reann #Can't pass self projection QC
    reann.left <- new.res.range$reann.left #Must pass self projection QC, #match p.left
    res.left <- new.res.range$res.min
    res.right <- new.res.range$res.max

    best.res <- findBestRes(p, p.left, reann, reann.left, out.name, res.left, res.right, outfile.path, res.switch,
             type="PCA",sp.thresholds, method, n.iterations, name, graph, clf.data, uncertainty.thresholds)
    reann.left <- best.res$reann.left
    #TODO: to be removed
    #res.right <- best.res$res.right
    message("Now we got next layer for ", cl)

    if ((reann.left$clf.data$marker.list %>% length) > 0){#Do not want to export singletons
      annotation <- reann.left$ann.by.level$annotation$l1 %>% sapply(function(n) paste0(cl,"_", n)) %>%
                                              setNames(reann.left$ann.by.level$annotation$l1 %>% names)
      ann.by.parents <- c(ann.by.parents, list(annotation) %>% setNames(cl))
      marker.list <- c(marker.list, list(reann.left$clf.data$marker.list) %>% setNames(cl))
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
findBestRes <- function(p, p.left, reann, reann.left, out.name, res.left, res.right, outfile.path,
               res.switch=0.1, type="PCA", sp.thresholds=list(accuracy=0.85, confusion.rate=3/17),
               method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
               clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
  message("sp.thresholds, res.switch ", sp.thresholds, " ", "res.switch")
  if ((res.right-res.left)<=res.switch) {
    return(list(p.left=p.left, reann.left=reann.left, res.left=res.left, res.right=res.right))
  }
  message("Now, we are heading to the left for the next smallest resolution ...")
  to.left <- findNextSmallestRes(p, p.left, reann, reann.left, out.name, outfile.path, res.left, res.right,
             res.right.old=res.right, res.switch, type="PCA", method, n.iterations, name, graph, clf.data,
             clusters, uncertainty.thresholds, sp.thresholds)
  p.left <-  to.left$p
  reann.left <- to.left$reann #This must be a good one. Pass QC
  res.left <- to.left$res.right
  res.right <- to.left$res.right.old
  if (length(reann$ann.by.level$annotation$l1 %>% table)==1) {
    return(list(p.left=p.left, reann.left=reann.left, res.left=res.left, res.right=res.right))
    }
  if ((res.right-res.left)>res.switch){
    message("Now, we are heading to the right for the next biggest resolution ...")
    to.right <- findNextBiggestRes(p, p.old=p, out.name, outfile.path, reann=reann.left, reann.old=reann.left,
                res.left, res.left.old=res.left, res.right, res.switch, type="PCA", method, n.iterations,
                name, graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds)
    p <- to.right$p
    p.left <- to.right$p.old
    res.left <- to.right$res.left.old #the highest res that work
    res.right <- to.right$res.left #the lowest res that doesn't work
    reann.left <- to.right$reann.old
    reann <- to.right$reann #In case we want to turn left again
  }
  return(findBestRes(p, p.left, reann, reann.left, out.name,
                     res.left, res.right, outfile.path, res.switch, type, sp.thresholds,
                     method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
                     clf.data=NULL, uncertainty.thresholds))
}


#####Main function 4
#To left
#There will be problem for this function, if the smallest res is very close (diff < 0.0001) to the original res.left, then the run may break.
#The error will look like: Error in ann.by.level[[i + 1]] : subscript out of bounds (reAnnPip: ann.by.level only has one cluster! singleton)
#Possible solution: increase relax the sp.threshold
findNextSmallestRes <- function(p, p.left, reann, reann.left, out.name, outfile.path, res.left, res.right,
            res.right.old, res.switch=0.1, type="PCA", method=conos::leiden.community,
            n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, clusters = NULL,
            uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
            sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  #TODO: message to be removed.
  message("res.right, res.left, res.switch, sp.thresholds ", res.right, res.left, res.switch,
          " ", sp.thresholds)
  if ((res.right-res.left)<=res.switch){ #Stop and res.left considered not pass #This can't be introduced!
    return(list(p=p.left, reann=reann.left, res.left=res.left, res.right.old=res.right, res.right=res.left))
  }

  #TODO: message to be removed
  message("reAnnotation QC ...")
  #reann and p should match.
  pass.reAnnotation.QC <- reAnnotationQC(reann, p, sp.thresholds)
  if (pass.reAnnotation.QC){ #We look to the left for a res at which clustering pass QC #Pass
    return(list(p=p, reann=reann, res.left=res.left, res.right.old=res.right.old, res.right=res.right))
  }

  message("Now we will decrease resolution to next middle ...")
  #Decrease res.right until no cluster has less than 6 cells..
  p$getKnnClusters(type="PCA", method=conos::leiden.community,
                   resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
  res.right.old = res.right
  res.right = res.left+(res.right-res.left)/2

  big.clusters <- getBigClusters(p, res.left, res.right, res.right.old, #direction="left",#Make sure all cluster >5 cells
                  type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden")
  p <- big.clusters$p
  res.right.old <- big.clusters$res.right.old
  res.right <- big.clusters$res.right

  if (length(p$clusters$PCA$leiden %>% table) == 1){#Singleton can't do reannotation.
    reann <- list(ann.by.level=list(annotation=list(l1=p$clusters$PCA$leiden)),#Build reann manually.
             clf.data=list(marker.list=list())) #Singleton's reann should pass all self projection QC.
    return(list(p=p, reann=reann, res.left=res.left,
                          res.right.old=res.right.old, res.right=res.right))
  }

  # if (sum((p$clusters$PCA$leiden %>% table) < 6) > 1){ #All clusters should >= 6 cells, otherwise not pass
  #   return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
  #               res.left=res.left, res.left.old=res.left.old, res.right=res.right))
  # }

  reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden),
                     out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
  message("Now the new res.right and res.left become ", res.right, " ", res.left)

  return(findNextSmallestRes(p, p.left, reann, reann.left, out.name,
         outfile.path, res.left, res.right, res.right.old, res.switch, type="PCA", method,
         n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds))
}



#####Main function 5
#The starting reann should fufill the 3 if conditions in the function!!!
#Look to the right for the maximum separable resolution (i.e. the 1st one that doesn't pass)
#i.e. find the first resolution, at which the clustering results does not pass self projection QC.
#reann.old: no matter the initial input or the generated one from below will all be fine.
findNextBiggestRes <- function(p, p.old=p, out.name, outfile.path, reann, reann.old=reann, res.left,
                           res.left.old=res.left, res.right, res.switch=0.1, type="PCA",
                           method=conos::leiden.community, n.iterations=50, name="leiden",
                           graph=NULL, clf.data= NULL, clusters = NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           sp.thresholds=list(accuracy=0.85, confusion.rate=3/17)){
  #TODO: message to be removed...
  message("1.res.switch, res.right, res.left, sp.thresholds ", res.switch, res.right, res.left, " ", sp.thresholds)
  #Stop and res.left considered not pass #shouldn't be introduced
  if ((res.right-res.left)<=res.switch) {return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
                          res.left=res.left, res.left.old=res.left.old, res.right=res.right))} #counted as Not pass

  pass.reAnnotation.QC <- reAnnotationQC(reann, p, sp.thresholds)
  if (pass.reAnnotation.QC){
    message("Now we will increase resolution to next middle ...")
    p.old <- p #Pass QC
    reann.old <- reann #Pass QC
    #Increase res.left until no cluster has less than 6 cells
    p$getKnnClusters(type="PCA", method=conos::leiden.community,
                     resolution=res.right-(res.right-res.left)/2, n.iterations=50, name="leiden")
    res.left.old = res.left
    res.left = res.right-(res.right-res.left)/2
    # big.clusters <- getBigClusters(p, res.left, res.right, direction="right", #All clusters >= 6 cells
    #                         type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden")
    # p <- big.clusters$p
    # res.left <- big.clusters$res.left
    # res.left.old <- big.clusters$res.left.old
    if (sum((p$clusters$PCA$leiden %>% table) < 6) > 0){ #All clusters should >= 6 cells, otherwise not pass
      return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
                  res.left=res.left, res.left.old=res.left.old, res.right=res.right))
    }

    if (length(p$clusters$PCA$leiden %>% table)==1){#Singleton can't be brought in, but may be produced inside.
      reann <- list(ann.by.level=list(annotation=list(l1=p$clusters$PCA$leiden)), #Counted as pass
                    clf.data=list(marker.list=list()))
      return(list(p=p, p.old=p, reann=reann, reann.old=reann,
                  res.left.old=res.left, res.left=res.right, res.right=res.right))
    }#reAnnPip doesn't work on singletons, one has to create reann manually.
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden),
                       out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
    message("Now the new res.left and res.right become ", res.left, " ", res.right)

    return(findNextBiggestRes(p, p.old, out.name, outfile.path, reann, reann.old, res.left, res.left.old,
                              res.right, res.switch, type="PCA", method, n.iterationss, name,
                              graph, clf.data, clusters, uncertainty.thresholds, sp.thresholds))
  }

  return(list(p=p, p.old=p.old, reann=reann, reann.old=reann.old,
              res.left=res.left, res.left.old=res.left.old, res.right=res.right))
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
  for (i in 1:(length(ann.by.level)-1)){ #In order to use split, all layers should have the same cell order.
    if (length(ann.by.level[[i+1]]) != length(ann.by.level[[i]])){
      ann.by.level[[i+1]] <- c(ann.by.level[[i+1]], ann.by.level[[i]][names(ann.by.level[[i]]) %>%
                                                             setdiff(names(ann.by.level[[i+1]]))])
    }
    ann.by.level[[i+1]] %<>% .[ann.by.level[[i]] %>% names]
  }
  #Gather marker.list
  marker.list <- c(marker.list, list(reann.layer1$clf.data$marker.list)%>%setNames("l1")) # For layer 1
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
  } else {stop("There is no subclusters in the cell population, or try more gene features...")}

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
    ann.split[[paste0("l",i)]] <- ann.by.level[[i+1]] %>% split(ann.by.level[[i]]) %>% lapply(unique)
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


#' @param classificaiton should be either classificaiton.tree or a classification dataframe
#' @export
plotTypeHierarchyLi <- function(classification, layout="slanted", xlims=NULL, font.size=3, col="black", ...) { #modified from Viktor's function
  if (!requireNamespace("ggtree", quietly=T)){
    stop("You need to install package 'ggtree' to be able to plot hierarchical tree. ",
         "`Try devtools::install_github('YuLab-SMU/ggtree')`")
  }
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
      ggtree::geom_label2(ggplot2::aes(label=c(cg$tip.label, cg$node.label)[node]), size=font.size, color=col) +
      ggplot2::xlim(xlims)
  }
}

