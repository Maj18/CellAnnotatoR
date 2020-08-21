
#' Self Projection Quality Control
#'
#' @description Used for the self projection quaility control after re-annotating
#'  the leiden-clustering-based hierarchy.
#'  Confusion rate: incorrectly predicted cell number that is normalized as being
#'  divided by the corrected predicted cells.
#'  The cluster sample size is balanced by adjusting weights inversely
#'  proportional to class frequencies in the input data
#' @param annotation.before,annotation.after The original and predicted
#' (based on selected markers) cell annotation, they need to have the same length.
#' @return A list including one vector with self projection accuracy for each parent,
#' and one matrix with columns being clusters, while rows being "conf", "fpr".
selfProjQCperCluster <- function(annotation.before, annotation.after) {
  #Get self projection accuracy:
  acc <- sum(annotation.before == annotation.after)/length(annotation.before)
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




#' Get classification data and marker list
#'
#' @description The markers are selected from differentially expressed genes and
#' the markers are selected under parents, which means the selected markers are supposed
#' to distinguish one cluster from all the rest that belong to the same parent cluster.
#' @param p2
#' @param ann.by.level
#' @param name
#' @param outfile.path
#' @param stringent
#' @return
getClfDataInferringMarkers <- function(p2, ann.by.level, name, outfile.path, stringent=T){
  #Check singleton
  if (length(ann.by.level)==1 && (ann.by.level[[1]] %>% unique %>% length)==1 && stringent)
    return(NULL)

  if ((ann.by.level[[1]] %>% unique %>% length) != 1) {
    ann.by.level <- list(root=rep("root", length(ann.by.level[[1]])) %>%
                    setNames(names(ann.by.level[[1]]))) %>% c(ann.by.level)
  }

  ann.by.parent <-  lapply(1:(length(ann.by.level)-1), function(i)
    split(ann.by.level[[i+1]], ann.by.level[[i]])) %>% Reduce(c, .) %>%
    .[sapply(., length) > 5] %>% .[sapply(., function(x) length(unique(x)) > 1)]

  if ((ann.by.parent %>% .[sapply(., function(c) length(unique(c))==1)] %>% length)>0
      && stringent) {return(NULL)}

  marker.list <- selectMarkersByParent(p2, ann.by.parent, z.threshold=0, append.auc=T,
        max.iters=50, max.uncertainty=0.25, max.pos.markers=2, log.step=5, verbose=1,
        n.cores=10, stringent=stringent)

  if (is.null(marker.list) && stringent) {return(NULL)}

  message("Check whether all clusters have markers: ")
  if ((marker.list %>% .[sapply(., function(x) length(x$expressed) == 0)] %>% length)>0 &&
      stringent) {return(NULL)
  }

  marker.list %>% .[sapply(., function(x) length(x$expressed) > 0)] %>%
    markerListToMarkup(file=paste0(outfile.path, "marker_list_", name, ".txt"), group.by.parent=F)

  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  clf.data <- getClassificationData(cm.norm, paste0(outfile.path, "marker_list_", name, ".txt"))

  return(clf.data)
}


selectMarkersByParent <- function(p2, ann.by.parent, z.threshold=0, append.auc=T, max.iters=50,
                         max.uncertainty=0.25, max.pos.markers=10, log.step=5,
                         verbose=1, n.cores=1, stringent=T){
  message("Selecting markers for clusters under the same parents ...")
  de.info.per.parent <- ann.by.parent %>%
    pblapply(function(ann) p2$getDifferentialGenes(groups=ann, z.threshold=0, append.auc=T))

  pre.selected.markers <- de.info.per.parent %>% lapply(function(dfs)
    list(positive=lapply(dfs, function(df) as.character(df$Gene[df$Z > 0.001])),
         negative=lapply(dfs, function(df) as.character(df$Gene[df$Z < -0.001]))
    )
  )
  #Check whether all clusters have pre selected positive markers
  if (sum(pre.selected.markers %>% sapply(., function(c) c$positive%>%sapply(length)==0)%>%Reduce(c,.))>0
      && stringent) {return(NULL)}

  cm.norm <- p2$counts %>% Matrix::t() #%>% normalizeTfIdfWithFeatures()
  marker.info <- names(de.info.per.parent) %>% setNames(., .) %>% lapply(function(n)
    selectMarkersPerType(cm.norm, ann.by.parent[[n]], pre.selected.markers[[n]],
    parent=n, max.iters=50, max.uncertainty=0.25, max.pos.markers=10, log.step=5, verbose=1, n.cores=1))

  marker.list <- setNames(marker.info, NULL) %>% unlist(recursive=F)

  return(marker.list)
}






#Outfile.path: where to save the marker file
reAnnoPip <- function(p2, ann.by.level, name, outfile.path, graph=NULL, clf.data=NULL,
          clusters = NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
          stringent=T){
  #Check small clusters (< 6 cells)
  if (length(ann.by.level)==1 && sum((ann.by.level[[1]] %>% table) < 6)>0
      && stringent) {return(NULL)} #only one layer allowed
  #Check singletons
  if (length(ann.by.level)==1 && (ann.by.level[[1]] %>% unique %>% length)==1 &&
       stringent) {return(NULL)}

  if (is.null(clf.data))
    {clf.data <- getClfDataInferringMarkers(p2, ann.by.level, name, outfile.path, stringent)}

  if (is.null(clf.data) || (clf.data$marker.list %>% length)==1){
    if (stringent) {return(NULL)}
  } #Check whether the new clf.data is still NULL(i.e. some clusters have no markers).

  message("Re-annotation ...")
  ann.by.level <-
    assignCellsByScores(graph=NULL, clf.data=clf.data, score.info=NULL,
                      clusters=NULL, verbose=0, uncertainty.thresholds, max.depth=NULL)

  return(list(ann.by.level=ann.by.level, clf.data=clf.data))
}







reAnnotationQC <- function(reann, p, confusion.rate.threshold, certainty.threshold=0.5){
  if (!is.null(reann)){#cluster size before reann>5 cells, all clusters have markers
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    #no cluster dropping during reannotation (indicating bad markers)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
      if (sum(sp$confusion.rate.false.pos.rate[1,] > confusion.rate.threshold)==0 &&
          #confusion rate not above the threshold
          sum((reann$ann.by.level$annotation$l1 %>% table)<6)==0){
        cells.by.clusters <- reann$ann.by.level$annotation$l1 %>% split(.,.) %>% lapply(names)
        if (length(cells.by.clusters)==1) {#Check singletons
          passQC <- TRUE
          return(passQC)
        }

        mean.st.per.cluster <- sapply(cells.by.clusters%>%names, function(c){
          (reann$ann.by.level$scores$l1)[cells.by.clusters[[c]], c] %>% mean})
        if (sum(mean.st.per.cluster <= certainty.threshold)==0){
          passQC <- TRUE
          return(passQC)
        }
      }
    }
  }
  passQC <- FALSE
  return(passQC)
}







getBigClusters <- function(sub.graph=NULL, p.left, p.right, p.right.old, res.left, res.right,
            res.right.old, reann.right.old, out.name, outfile.path, graph=NULL, clf.data=NULL,
            clusters=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
            type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden"){
  n = 0
  while (sum((p.right$clusters$PCA$leiden %>% table) < 6) > 0){
    if ((res.right-res.left) <= 0.05){
      break
    }
    n = n+1
    p.right.old <- p.right
    if (is.null(sub.graph)){
      p.right$getKnnClusters(type="PCA", method=conos::leiden.community,#decrease res.right, get new p.right
                             resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
    } else {
      sub.clusters <- conos::leiden.community(sub.graph, resolution=res.left+(res.right-res.left)/2, n.iterations=50)
      p.right$clusters$PCA$leiden <- sub.clusters$membership
    }

    res.right.old = res.right
    res.right = res.left+(res.right-res.left)/2
  }

  if (n>0){#update the reann.right before decreasing res
    #Singleton can't be brought in, but may be produced inside.
    if (length(p.right.old$clusters$PCA$leiden %>% table)==1){
      reann.right.old <- list(ann.by.level=list(annotation=list(l1=p.right.old$clusters$PCA$leiden)),
                         clf.data=list(marker.list=list()))
    } else{
        reann.right.old <- reAnnoPip(p.right.old, list(l=p.right.old$clusters$PCA$leiden),
                           out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
    }
  }
  return(list(p.left=p.left, p.right=p.right, p.right.old=p.right.old, res.left=res.left,
                      res.right=res.right, res.right.old=res.right.old, reann.right.old=reann.right.old))
}







#Get the first layer of leiden clusters (>1 clusters) for the leiden clustering based hierarchy
#p2: dataset is in pagoda form.
#sub.graph=NULL: equal to data.splitting="dataset"
getLayer1KnnClusters <- function(sub.graph=NULL, p2, res.start=0.01, res.step=0.01){
  if (!is.null(sub.graph)){
    #clustering based on subgraph.
    sub.clusters <- conos::leiden.community(sub.graph, resolution=res.start, n.iterations=50)
    p2$clusters$PCA$leiden <- sub.clusters$membership
  }else{
    p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start,
                      n.iterations=50, name="leiden")
  }

  clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  n = 0
  while (clusters.n==1){
    n = n+1
    if (!is.null(sub.graph)){
      #clustering based on subgraph.
      sub.clusters <- conos::leiden.community(sub.graph, resolution=res.start+n*res.step, n.iterations=50)
      p2$clusters$PCA$leiden <- sub.clusters$membership
    }else{
      p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start+n*res.step,
                        n.iterations=50, name="leiden")
    }
    clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  }
  return(list(p2=p2, res=res.start+n*res.step))
}







#To solve the problem that the given res.max is too small.
#If KNN clustering gets a singleton,
#not sure whether it is a real singleton or just because the resolution is too low,
#we try to increase the res to get a non-singleton clustering,
#and see wheteher the new clustering fullfil requirements below or not later on.
#sub.graph=NULL: equal to data.splitting="dataset"
tailorResRange <- function(sub.graph, p.left, reann.left, res.min, res.max, res.switch, out.name,
                           outfile.path, graph=NULL, res.max.update=1, clf.data=NULL, clusters = NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           confusion.rate.threshold=3/17, certainty.threshold=0.5){
  message("confusion.rate.threshold ", confusion.rate.threshold)
  if ((res.max-res.min) <= res.switch) {res.max = (res.max+1)}
  #Increase res until not singleton
  check.singleton <- getLayer1KnnClusters(sub.graph=sub.graph, p2=p.left, res.start=res.max, res.step=1)
  p.right <- check.singleton$p2
  res.max <- check.singleton$res #Update max.res.middle
  #All clusters should have >= 6 cells, otherwise not pass
  if (sum((p.right$clusters$PCA$leiden %>% table) < 6) > 0){
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left,
                reann.right=NULL, res.min=res.min, res.max=res.max))
  }

  reann.right <- reAnnoPip(p.right, list(l=p.right$clusters$PCA$leiden), out.name, outfile.path,
                           graph, clf.data, clusters, uncertainty.thresholds)

  pass.reAnnotation.QC <-
    reAnnotationQC(reann.right, p.right, confusion.rate.threshold, certainty.threshold)
  if (pass.reAnnotation.QC){
    res.min <- res.max
    p.left <- p.right
    reann.left <- reann.right
    res.max <- res.max+res.max.update
    return(tailorResRange(sub.graph, p.left, reann.left, res.min, res.max, res.switch,
                          out.name, outfile.path, graph, res.max.update, clf.data,
                          clusters, uncertainty.thresholds, confusion.rate.threshold,
                          certainty.threshold))
  }

  return(list(p.left=p.left, p.right=p.right, reann.left=reann.left,
              reann.right=reann.right, res.min=res.min, res.max=res.max))
}







#Main function 2
getNextLayerKnnClusters <- function(p2, annotation, out.name, outfile.path, min.res.start=0,
                           graph=NULL, res.max.update=1, clf.data=NULL,
                           uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                           max.res.middle=1, res.switch=0.1, type="PCA",
                           method=conos::leiden.community, n.iterations=10, name="leiden",
                           confusion.rate.threshold=3/17, clustering.type=NULL, embeding.type=NULL,
                           data.splitting=c("dataset", "graph"), reannotation=FALSE,
                           certainty.threshold=0.5){

  ann.by.parents <- list()
  marker.list <- list()
  ann.by.cluster <- split(annotation, annotation)
  for (cl in names(ann.by.cluster)){
    message("Now we will get next layer for cluster ", cl)
    p <- p2$misc$rawCounts[ann.by.cluster[cl][[1]]%>%names,] %>% Matrix::t() %>%
      vpscutils::GetPagoda(clustering.type, embeding.type, n.cores=1) #Withdraw data subset for cl.

    if (data.splitting=="graph"){
      sub.graph <- igraph::induced_subgraph(p2$graphs$PCA, ann.by.cluster[cl][[1]]%>%names)
      sub.clusters <- conos::leiden.community(sub.graph, resolution=min.res.start, n.iterations=50)
      p$clusters$PCA$leiden <- sub.clusters$membership
    }

    if (data.splitting=="dataset"){
      #leiden clustering
      p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=min.res.start,
                       n.iterations=10, name="leiden") #Knn clustering
      sub.graph <- NULL
    }
    #Because min.res.start=0, so the Knn cluster must be a singleton that can't do reannotation.
    p.left <- p #has the clustering at res.min
    reann.left <- list(ann.by.level=list(annotation=list(l1=p.left$clusters$PCA$leiden)),
                  clf.data=list(marker.list=list())) #Build reann.left manually
    #Tailor the max.res (and maybe min.res as well) for cl:
    new.res.range <- tailorResRange(sub.graph, p.left, reann.left, res.min=min.res.start,
                     res.max=max.res.middle, res.switch, out.name,outfile.path, graph,
                     res.max.update, clf.data, clusters, uncertainty.thresholds,
                     confusion.rate.threshold, certainty.threshold)
    p.left <- new.res.range$p.left
    p.right <- new.res.range$p.right
    reann.left <- new.res.range$reann.left #Must pass self projection QC, #match p.left
    reann.right <- new.res.range$reann.right #Can't pass self projection QC
    res.left <- new.res.range$res.min
    res.right <- new.res.range$res.max

    best.res <- findBestRes(sub.graph, p.left, p.right, reann.left, reann.right, res.left, res.right,
                            out.name, outfile.path, res.switch, type="PCA", confusion.rate.threshold,
                            method, n.iterations, name, graph, clf.data, uncertainty.thresholds,
                            certainty.threshold)
    reann.left <- best.res$reann.left
    p.left <- best.res$p.left
    message("Now we got next layer for ", cl)

    if (reannotation){
      if ((reann.left$ann.by.level$annotation$l1 %>% table %>% length) > 1){#Do not want to export singletons
        annotation <- reann.left$ann.by.level$annotation$l1 %>% sapply(function(n) paste0(cl,"_", n)) %>%
                                                setNames(reann.left$ann.by.level$annotation$l1 %>% names)
        ann.by.parents <- c(ann.by.parents, list(annotation) %>% setNames(cl))
        marker.list.cl <- reann.left$clf.data$marker.list %>% lapply(function(sub.clust) {
          sub.clust$"parent"=cl
          sub.clust})
        marker.list <- c(marker.list, list(marker.list.cl) %>% setNames(cl))
      }
    }

    if (!reannotation){
      if ((p.left$clusters$PCA$leiden %>% table %>% length) > 1) #&& sum((p.left$clusters$PCA$leiden %>% table)<6)==0)
        {annotation <- p.left$clusters$PCA$leiden %>% sapply(function(n) paste0(cl,"_", n)) %>%
                                                  setNames(p.left$clusters$PCA$leiden %>% names)
          ann.by.parents <- c(ann.by.parents, list(annotation) %>% setNames(cl))
          marker.list.cl <- reann.left$clf.data$marker.list %>% lapply(function(sub.clust) {
            sub.clust$"parent"=cl
            sub.clust})
          marker.list <- c(marker.list, list(marker.list.cl) %>% setNames(cl))
        }
    }

  }
  return(list(ann.by.parents=ann.by.parents, marker.list=marker.list))
}





#Main function 3
#Look for the maximum resolution (e.g. the maximum number of predicted clusters)
#that have its predicted clusters all pass quality control (self projection accuracy via re-annotation: confusion rate)
findBestRes <- function(sub.graph, p.left, p.right, reann.left, reann.right, res.left,
               res.right, out.name, outfile.path, res.switch=0.1, type="PCA",
               confusion.rate.threshold=3/17, method=conos::leiden.community,
               n.iterations=10, name="leiden", graph=NULL, clf.data=NULL,
               uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
               certainty.threshold=0.5){

  if ((res.right-res.left)<=res.switch) {
    return(list(p.left=p.left, reann.left=reann.left, res.left=res.left, res.right=res.right))
  }

  pass.reAnnotation.QC <-
    reAnnotationQC(reann.right, p.right, confusion.rate.threshold, certainty.threshold)
  if (pass.reAnnotation.QC) {
    stop("Something is wrong here, reann.right must not pass reAnnotation QC....")}
  message("Now, we are heading to the left for the next smallest resolution ...")
  to.left <- findNextSmallestRes(sub.graph, p.left, p.right, p.right.old=p.right,
             reann.left, reann.right, reann.right.old=reann.right, res.left, res.right,
             res.right.old=res.right, out.name, outfile.path, res.switch, type="PCA", method,
             n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds,
             confusion.rate.threshold, certainty.threshold)
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
    to.right <- findNextBiggestRes(sub.graph, p.left, p.left.old=p.left, p.right, reann.left,
                                   reann.left.old=reann.left, reann.right, res.left, res.left.old=res.left,
                                   res.right, out.name, outfile.path, res.switch, type="PCA", method,
                                   n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds,
                                   confusion.rate.threshold, certainty.threshold)
    p.left <- to.right$p.left
    p.right <- to.right$p.right
    res.left <- to.right$res.left #the highest res that work
    res.right <- to.right$res.right  #the lowest res that doesn't work
    reann.left <- to.right$reann.left
    reann.right <- to.right$reann.right #In case we want to turn left again
  }
  return(findBestRes(sub.graph, p.left, p.right, reann.left, reann.right, res.left, res.right,
                     out.name, outfile.path, res.switch, type, confusion.rate.threshold,
                     method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
                     clf.data=NULL, uncertainty.thresholds, certainty.threshold))
}






#Main function 4
findNextSmallestRes <- function(sub.graph, p.left, p.right, p.right.old, reann.left, reann.right,
            reann.right.old, res.left, res.right,res.right.old, out.name, outfile.path, res.switch=0.1,
            type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL,
            clf.data= NULL, clusters = NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
            confusion.rate.threshold=3/17, certainty.threshold=0.5){

  message("reAnnotation QC ...")
  pass.reAnnotation.QC <-
    reAnnotationQC(reann.right, p.right, confusion.rate.threshold, certainty.threshold)
  message("confusion.rate.threshold ", confusion.rate.threshold)
  if (pass.reAnnotation.QC){ #We look to the left for a res at which clustering pass QC #Pass
    return(list(p.left=p.right, p.right=p.right.old, reann.left=reann.right,
                reann.right=reann.right.old, res.left=res.right, res.right=res.right.old))
  }

  if ((res.right-res.left)<=res.switch && !pass.reAnnotation.QC){ #Stop #This can't be introduced!
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left, reann.right=reann.right,
                res.left=res.left, res.right=res.right))
  }

  message("Now we will decrease resolution to next middle ...")
  p.right.old = p.right
  if (is.null(sub.graph)){
    p.right$getKnnClusters(type="PCA", method=conos::leiden.community,
                  resolution=res.left+(res.right-res.left)/2, n.iterations=50, name="leiden")
  }
  if (!is.null(sub.graph)){
    sub.clusters <-
      conos::leiden.community(sub.graph, resolution=res.left+(res.right-res.left)/2, n.iterations=50)
    p.right$clusters$PCA$leiden <- sub.clusters$membership
  }

  res.right.old = res.right
  res.right = res.left+(res.right-res.left)/2

  reann.right.old <- reann.right
  big.clusters <- getBigClusters(sub.graph, p.left, p.right, p.right.old, res.left, res.right, res.right.old,
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

  reann.right <- reAnnoPip(p.right, list(l=p.right$clusters$PCA$leiden),
                     out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
  message("Now the new res.right and res.left become ", res.right, " ", res.left)

  return(findNextSmallestRes(sub.graph, p.left, p.right, p.right.old, reann.left, reann.right, reann.right.old,
                             res.left, res.right, res.right.old, out.name, outfile.path, res.switch, type="PCA",
                             method, n.iterations, name, graph, clf.data, clusters, uncertainty.thresholds,
                             confusion.rate.threshold, certainty.threshold))
}






#####Main function 5
#The starting reann should fufill the 3 if conditions in the function!!!
#Look to the right for the maximum separable resolution (i.e. the 1st one that doesn't pass)
#i.e. find the first resolution, at which the clustering results does not pass self projection QC.
#reann.old: no matter the initial input or the generated one from below will all be fine.
findNextBiggestRes <- function(sub.graph, p.left, p.left.old=p.left, p.right, reann.left,
                               reann.left.old=reann.left, reann.right, res.left, res.left.old=res.left, res.right,
                               out.name, outfile.path, res.switch=0.1, type="PCA", method=conos::leiden.community,
                               n.iterations=50, name="leiden", graph=NULL, clf.data= NULL, clusters = NULL,
                               uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                               confusion.rate.threshold=3/17, certainty.threshold=0.5){

  pass.reAnnotation.QC <-
    reAnnotationQC(reann.left, p.left, confusion.rate.threshold, certainty.threshold)

  if ((res.right-res.left)<=res.switch && pass.reAnnotation.QC){
    return(list(p.left=p.left, p.right=p.right, reann.left=reann.left, reann.right=reann.right,
                res.left=res.left, res.right=res.right))}

  if (pass.reAnnotation.QC){
    message("Now we will increase resolution to next middle ...")
    p.left.old <- p.left #Pass QC
    reann.left.old <- reann.left #Pass QC
    #Increase res.left until no cluster has less than 6 cells
    if (is.null(sub.graph)) {
      p.left$getKnnClusters(type="PCA", method=conos::leiden.community, #get new p.left clustering
                      resolution=res.right-(res.right-res.left)/2, n.iterations=50, name="leiden")
    }
    if (!is.null(sub.graph)){
      sub.clusters <-
        conos::leiden.community(sub.graph, resolution=res.right-(res.right-res.left)/2, n.iterations=50)
      p.left$clusters$PCA$leiden <- sub.clusters$membership
    }

    res.left.old = res.left
    res.left = res.right-(res.right-res.left)/2

    if (length(p.left$clusters$PCA$leiden %>% table)==1){#Singleton can't be brought in, but may be produced inside.
      reann.left <- list(ann.by.level=list(annotation=list(l1=p.left$clusters$PCA$leiden)), #Counted as pass
                    clf.data=list(marker.list=list()))
    } else {#reAnnPip will output NULL on singletons, one has to create reann manually.
        reann.left <- reAnnoPip(p.left, list(l=p.left$clusters$PCA$leiden),
                       out.name, outfile.path, graph, clf.data, clusters, uncertainty.thresholds)
    }
    message("Now the new res.left and res.right become ", res.left, " ", res.right)

    return(findNextBiggestRes(sub.graph, p.left, p.left.old, p.right, reann.left,
                              reann.left.old, reann.right, res.left, res.left.old, res.right, out.name,
                              outfile.path, res.switch, type="PCA", method, n.iterationss, name,
                              graph, clf.data, clusters, uncertainty.thresholds,
                              confusion.rate.threshold, certainty.threshold))
  }

  return(list(p.left=p.left.old, p.right=p.left, reann.left=reann.left.old, reann.right=reann.left,
              res.left=res.left.old, res.right=res.left))
}


#Main function 6
getNextLayersKnnClusters <- function(p2, annotation, out.name, outfile.path, layer.n, layer.no.factor=c(0, 2/3),
                                     min.res.start=0, graph=NULL, res.max.update=1, clf.data=NULL,
                                     uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                                     max.res.middle=1, max.res.increase=1, res.switch=0.05, type="PCA",
                                     method=conos::leiden.community, n.iterations=50, name="leiden",
                                     confusion.rate.threshold=1/9, confusion.rate.maximum=3/7,
                                     clustering.type=NULL, embeding.type=NULL, data.splitting="dataset",
                                     reannotation=FALSE, certainty.threshold=0.5){
  next.layers <- list()
  for (l in 2:layer.n){
    if (length(annotation) > 0){
      message("Now we are going to get layer ", l, "...")
      next.layer <- getNextLayerKnnClusters(p2, annotation, out.name, outfile.path, min.res.start, graph,
                    res.max.update, clf.data, uncertainty.thresholds, max.res.middle, res.switch, type="PCA",
                    method, n.iterations, name, confusion.rate.threshold, clustering.type, embeding.type,
                    data.splitting, reannotation, certainty.threshold)

      saveRDS(next.layer, file=paste0(outfile.path, "next.layer.rds"))

      if (length(next.layer$ann.by.parents)==0) {annotation=c()} else{
        annotation <- next.layer$ann.by.parents %>% Reduce(c,.)
        annotation <- annotation[annotation %in% (table(annotation) %>% .[.>=12] %>% names)]
      }

      confusion.rate.threshold <-
        pmin((confusion.rate.threshold+layer.no.factor), confusion.rate.maximum)
      max.res.middle <- max.res.middle + max.res.increase

      next.layers <- c(next.layers, list(next.layer) %>% setNames(paste0("l", l)))
    }
  }
  return(next.layers)
}






#Main function 7
gatherAnnByLevelMarkerList <- function(reann.layer1, ann.layer1, ann.next.layers){
  #Gather ann.by.level
  ann.by.level <- list()
  marker.list <- list()
  ann.by.level <- c(ann.by.level, list(ann.layer1)%>%setNames("l1"))
  #Take out ann.by.parents from other layers
  ann.by.level.r <- ann.next.layers%>% lapply(function(n) n$ann.by.parents %>% Reduce(c,.))
  ann.by.level <- c(ann.by.level, ann.by.level.r)
  message("Make all layers have the same cell order ...")
  for (i in 1:(length(ann.by.level)-1)){
    #In order to use split, all layers should have the same cell order.
    if (length(ann.by.level[[i+1]]) != length(ann.by.level[[i]])){
      ann.by.level[[i+1]] <- c(ann.by.level[[i+1]],
        ann.by.level[[i]][names(ann.by.level[[i]]) %>% setdiff(names(ann.by.level[[i+1]]))])
    }
    ann.by.level[[i+1]] %<>% .[ann.by.level[[i]] %>% names]
  }
  #Check if the last two layers are the same, if yes, remove the last layer
  depth <- length(ann.by.level)
  if (sum(ann.by.level[[depth]]==ann.by.level[[depth-1]])==length(ann.by.level[[depth]])){
    ann.by.level[depth] = NULL
  }
  #Gather marker.list
  reann.layer1$clf.data$marker.list %<>% lapply(function(t) {
    t$parent="root"
    t
  })
  # For layer 1
  marker.list <- c(marker.list, list(l1=list(reann.layer1$clf.data$marker.list)%>%setNames("root")))
  marker.list <- c(marker.list, ann.next.layers%>% lapply(function(n) n$marker.list) %>%
                     setNames(names(ann.next.layers))) #For the rest layers
  return(list(ann.by.level=ann.by.level, marker.list=marker.list))
}





#####Main function8:
#Get the leiden clustring based hierarchy:
#layer.n: the number of layers that the users want to get.
#The more layers, the more clusters one will get, though one will not get endless number of clusters,
#becasue a. the feature genes will eventually become exhausted, b. the size of cluster can't be smaller than 6,
#or c. the cells within clusters become homogeneous eventually.
#1. reann shouldn't be NA (as long as not all clusters have predicted markers, then reann will be NA)
#2. sp shouldn't have NaN or Inf
#3. All the clusters should hava a minimum size of 6 after re-annotation
#4. The self projection accuracy by re-annotation has to be >=0.6 (or 0.8) and the confusion rate has to be <=2/3 (0.25)
getLeidenHierarchy <- function(p2, out.name, outfile.path, layer.n=10, layer.no.factor=0.2,
                               res.step.layer1=0.01, min.res=0, max.res.layer2=1, max.res.increase=1,
                               res.switch=0.05, clustering.type=NULL,embeding.type=NULL,
                               res.max.update=1, clf.data=NULL, clusters = NULL,
                               uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75),
                               confusion.rate.threshold=1/9, confusion.rate.maximum=3/7, graph=NULL,
                               type="PCA", method=conos::leiden.community, n.iterations=50,
                               name="leiden", data.splitting="dataset", reannotation=FALSE
                               , certainty.threshold=0.5){
  if (res.max.update <= res.switch){
    stop("res.max.update must be larger than res.switch")}
  message("Get Knn clusters for layer 1 ...")
  #Get a big picutre in layer1, we will accept the first non-singletong KNN clustering result here.
  layer1 <- getLayer1KnnClusters(sub.graph=NULL, p2, res.start=0.01, res.step=res.step.layer1)
  p2 <- layer1$p2
  reann.layer1 <- reAnnoPip(p2, list(l1=p2$clusters$PCA$leiden), out.name, outfile.path, graph,
                            clf.data, clusters=NULL, uncertainty.thresholds) #Not singleton!

  message("Get Knn clusters for the rest layers ...")
  pass.reAnnotation.QC <-
    reAnnotationQC(reann.layer1, p2, confusion.rate.threshold, certainty.threshold)
  if (pass.reAnnotation.QC){
    if (reannotation){
      annotation.layer1 <- reann.layer1$ann.by.level$annotation$l1 %>% sapply(function(n) paste0("l_",n)) %>%
        setNames(reann.layer1$ann.by.level$annotation$l1 %>% names())
    } else {
      annotation.layer1 <- p2$clusters$PCA$leiden %>% sapply(function(n) paste0("l_",n)) %>%
                setNames(reann.layer1$ann.by.level$annotation$l1 %>% names())
    }

    next.layers <- getNextLayersKnnClusters(p2, annotation=annotation.layer1, out.name, outfile.path,
                   layer.n, layer.no.factor, min.res.start=min.res, graph, res.max.update, clf.data,
                   uncertainty.thresholds, max.res.middle=max.res.layer2, max.res.increase, res.switch,
                   type, method, n.iterations, name, confusion.rate.threshold, confusion.rate.maximum,
                   clustering.type, embeding.type, data.splitting, reannotation, certainty.threshold)
  } else {
    stop("There is no subclusters in the cell population,
         or adjust confusion.rate.threshold (confusion.rate = 1/accuracy -1), or try more gene features...")}

  return(gatherAnnByLevelMarkerList(reann.layer1, annotation.layer1, next.layers))
}





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



#' Plot cell type hierarchy
#'
#' @description modified from Viktor's function
#' @param classificaiton should be either classificaiton.tree or a classification dataframe
#' @inheritParams ggtree::ggtree
#' @inheritParams ggtree::geom_rootpoint
#' @inheritParams ggtree::geom_label2
#' @inheritParams ggplot2::xlim
#' @inheritParams ggplot2::aes
#' @inheritParams ape::as.phylo
#' @export
plotTypeHierarchyLi <- function(classification, layout="slanted",
                                xlims=NULL, font.size=3, col="black", ...) {
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




#' Uncertainty to tree traits
#'
#' @description Prepare the "trait" file for plotUncHierarchy() from the positive, negative and coverage uncertainty data (output from scoreCellUncertaintyPerLevel())
#' @param unc.per.cell a list, each element represents the annotation of all cells at on hierarchical level
#' @param ann.by.level a list of lists, each main element is one hierarchical level; at each hierarchical level, there is a sublist, which contains 3 elements: positive, negative and coverage uncertainty.
#' @return ann.unc a dataframe, there are 4 columns:"type" #cell types   "positive" "negative" "coverage" #uncertainty
#' @export
uncToTreeTrait <- function(unc.per.cell, ann.by.level) {
  ann.unc <- list()
  for (i in names(unc.per.cell) %>% setNames(.,.)) {
    ann.unc[[i]] <- merge(as.data.frame(ann.by.level[[i]]),
    #ann.unc[[i]] <- merge(as.data.frame(ann.by.level$annotation[[i]]),
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




#' Plot uncertainty on trees
#'
#' @description  Plot uncertainty (positive, negative, coverage) on hierarchical trees
#' @param c.data output file from re-annotation (by assignCellsByScores())
#' @param trait.file contain the mean positive, negative, coverage uncertainties for each cell type (including root, node and tips)
#' @param unc one of c("positive", "negative", "coverage")
#' @inheritParams dplyr::full_join
#' @inheritParams ggtree::geom_label2
#' @inheritParams ggplot2::xlim
#' @inheritParams cowplot::plot_grid
#' @inheritParams ggtree::nodeid
#' @inheritParams ape::as.phylo
#' @inheritParams dplyr::full_join
#' @inheritParams ggtree::ggtree
#' @export
plotUncHierarchy <- function(c.data, trait.file, layout="slanted",
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

  clf.data <- list()
  clf.data$marker.list <- marker.list
  clf.data$cm.norm <- p2$counts# %>% Matrix::t()
  clf.data$classification.tree <- createClassificationTree(marker.list)

  score.info<- getMarkerScoreInfo(clf.data)
  score.info.per.cluster <- getMarkerScoresPerCellType(clf.data, score.info)
  score.info.by.layer <-
    ann.by.level %>% lapply(unique) %>% lapply(function(n){
      n <- intersect(n, score.info.per.cluster%>%colnames)
      score.info.per.cluster[, n]})
  #Get uncertainty
  unc.info <-
    names(ann.by.level)[1:length(ann.by.level)] %>% setNames(., .) %>%
    plapply(function(n) scorePerCellUncertainty(ann.by.level[[n]],
                                      score.info.by.layer[[n]], score.info), verbose=T)

  #Prepare the uncertainty data (as "trait")
  unc.trait <- uncToTreeTrait(unc.info, ann.by.level)
  #saveRDS(unc.trait, file="~/CellAnnotatoR-clean/notebooks/26files/unc.trait.rds")

  #Plot uncertainty
  #saveRDS(clf.data, file="~/CellAnnotatoR-clean/notebooks/28files/clf.data.rds")
  #saveRDS(unc.trait, file="~/CellAnnotatoR-clean/notebooks/28files/unc.trait.rds")
  #clf.data <- readRDS("/Users/yuanli/Documents/Degree_project/Hierarchical_annotation/bin/CellAnnotatoR-clean/notebooks/28files/clf.data.rds")
  #unc.trait <- readRDS("/Users/yuanli/Documents/Degree_project/Hierarchical_annotation/bin/CellAnnotatoR-clean/notebooks/28files/unc.trait.rds")
  plotUncHierarchy(clf.data, unc.trait)
}



compareAdjustedRandIndex <- function(ann.marker.level.constant.acc, ann.by.level.manual){
  adj.rand.index <- lapply(ann.marker.level.constant.acc, function(acc) {
    abl <- acc$ann.by.level
    iter <- min(length(abl), length(ann.by.level.manual))
    adj.rand.index <- sapply(1:iter, function(l){
      mclust::adjustedRandIndex(abl[[l]], ann.by.level.manual[[l]])})
    if (length(ann.by.level.manual) > length(abl)){
      adj.rand.index2 <- sapply((iter+1):length(ann.by.level.manual), function(ll){
        mclust::adjustedRandIndex(abl[[iter]], ann.by.level.manual[[ll]])})
      adj.rand.index <- c(adj.rand.index, adj.rand.index2)
    } else if (length(ann.by.level.manual) < length(abl)) {
      adj.rand.index2 <- sapply((iter+1):length(abl), function(lll){
        mclust::adjustedRandIndex(abl[[lll]], ann.by.level.manual[[iter]])})
      adj.rand.index <- c(adj.rand.index, adj.rand.index2)
    }
    adj.rand.index
  })
  return(adj.rand.index)
}



#How many layers and how clusters per layer, per tree
getTreeFeature <- function(tree){
  layer.no <- length(tree)
  cluster.no.per.layer <- sapply(tree, function(l) unique(l)%>%length)
  return(list(layer.no=layer.no, cluster.no.per.layer=cluster.no.per.layer))
}


#How many layers and how clusters per layer
getTreeFeatures <- function(ann.marker.level.constant.acc){
  layer.no <- c()
  cluster.no.per.layer  <- c()
  tf <- lapply(ann.marker.level.constant.acc, function(acc){
    abl <- acc$ann.by.level
    tf <- getTreeFeature(abl)
  })
  layer.no <- sapply(tf, function(acc) acc$layer.no)
  cluster.no.per.layer <- sapply(tf, function(acc) acc$cluster.no.per.layer)
  return(list(layer.no=layer.no, cluster.no.per.layer=cluster.no.per.layer))
}

#For Hierarchy result presentation
getClusterNo <- function(df) {
  lapply(df, function(cr) {
    sapply(cr$ann.by.level, function(l){
      split(l,l) %>% length
    })
  })
}

#getClusterNo(ann.marker.level.constant.acc4)

#gather the confusion rate, cluster number per layer and the layer names information
gatherConfRClustNo <- function(df, conf.rates){
  confusion.rate <- c()
  cluster.no <- c()
  layers <- c()
  clust.no.list <- getClusterNo(df)
  for (r in 1:length(clust.no.list)) {
    confusion.rate <- c(confusion.rate, rep(conf.rates[r], length(clust.no.list[[r]])))
    cluster.no <- c(cluster.no, clust.no.list[[r]])
    layers <- c(layers, clust.no.list[[r]]%>%names)
  }
  return(list(confusion.rate=confusion.rate, cluster.no=cluster.no, layers=layers))
}



#Correlation measurement: proportionality
getProportionalitySingleLayer <- function(layer1, layer2, prop.matrix){
  propr.single.layer <-lapply(layer1, function(cluster.1) {
    sapply(layer2, function(cluster.2){
      prop.matrix[cluster.1, cluster.2]%>%mean
    })
  })
  return(propr.single.layer)
}



getProportionality <- function(ann.by.level1, ann.by.level2,
                               raw.counts, metric = c("rho", "phi", "phs", "cor", "vlr"),
                               ivar = "clr", symmetrize = FALSE,  p = 100){
  prop <- propr::propr(raw.counts%>%t(), metric = metric, ivar = "clr", symmetrize = FALSE, p = 100)
  depth.min <- min(length(ann.by.level1), length(ann.by.level2))
  ann.by.level1.cluster <- lapply(ann.by.level1, function(l) split(l,l) %>% lapply(names))
  ann.by.level2.cluster <- lapply(ann.by.level2, function(l) split(l,l) %>% lapply(names))

  props <- mapply(function(layer1, layer2) {
    lapply(layer1, function(cluster.1) {
      sapply(layer2, function(cluster.2){
        prop@matrix[cluster.1, cluster.2] %>% mean
      })
    })
  }, ann.by.level1.cluster[1:depth.min], ann.by.level2.cluster[1:depth.min])

  return(props %>% lapply(function(l) data.frame(l)))
}
