
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







# Building leiden clustering based hirarchy start here

#Get the first layer of leiden clusters (>1 clusters) for the leiden clustering based hierarchy  #TESTED!
getLayer1KnnClusters <- function(p2, res.start=0.01, res.step=0.01){ 
  p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start, n.iterations=50, name="leiden")
  clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  n = 0
  while (clusters.n==1){
    n = n+1
    p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.start+n*res.step, n.iterations=50, name="leiden")
    clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  }
  return(list(p2=p2, res=res.start+n*res.step))
}
#e.g. l2 <- getLayer1Clusters(p2, res.step=0.01)
#l2$res
#l2$p2$clusters$PCA$leiden %>% table

# res.max.update muste be larger than 
updateResRange <- function(p, res.min, res.max, out.name, outfile.path, graph, res.max.update=1, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){ #To solve the problem that the givem res.max is too small.
  #If KNN clustering gets a singleton, not sure whether it is a real singleton or just because the resolution is too low, we try to increase the res to get a non-singleton clustering, and see wheteher the new clustering fullfil requirements below or not later on.
  check.singleton <- getLayer1KnnClusters(p, res.start=res.max, res.step=1)
  p <- check.singleton$p2
  res.max <- check.singleton$res #Update max.res.middle 
  
  reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph, clf.data=NULL, uncertainty.thresholds)
  
  if (!is.null(reann)){
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
      if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)==0 && sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
        res.min <- res.max
        res.max <- res.max+res.max.update
        return(updateResRange(p, res.min, res.max))}}}
  return(list(p=p, reann=reann, res.min=res.min, res.max=res.max))
}


#Get the last layer of leiden clusters (>1 clusters) for the leiden clustering based hierarchy
#res.max.update must be larger than res.switch, otherwise, you will be in trouble for the first step of findBestRes.
getNextLayerKnnClusters <- function(p2, annotation, out.name, outfile.path, min.res.start, graph=NULL, res.max.update=1, clf.data=NULL, 
                                    uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), max.res.middle=1, res.switch=0.05, 
                                    type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", 
                                    sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){ 
  ann.by.cluster <- split(annotation, annotation)
  ann.by.parents <- list()
  marker.list <- list()
  resol.last.layer <- list()
  for (cl in names(ann.by.cluster)){
    message("Now we will get next layer for ", cl)
    #sub <- igraph::induced_subgraph(p2$graphs$PCA, ann.by.cluster[cl][[1]]%>%names) #Split the graph
    #igraph::V(sub)  # Show information about the graph
    #clusters <- conos::leiden.community(sub, resolution=max.res.middle, n.iterations=10) #Build clusters based on the split graph
    p <- p2$misc$rawCounts[ann.by.cluster[cl][[1]]%>%names,] %>% Matrix::t() %>% 
      vpscutils::GetPagoda(clustering.type=clustering.type, embeding.type=embeding.type, n.cores=1)
    new.res.range <- updateResRange(p, min.res.start, max.res.middle, out.name, outfile.path, graph, res.max.update, clf.data, uncertainty.thresholds, sp.thresholds) #Tailor the min.res and max.res for cluster cl.
    p <- new.res.range$p
    reann <- new.res.range$reann
    res.left <- new.res.range$res.min
    res.right <- new.res.range$res.max

    if (is.null(reann)){ #Some clusters do not have markers
      best.res <- findBestRes(p, reann, out.name, res.left, res.right, res.switch=0.05, outfile.path, type, sp.thresholds=sp.thresholds, method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75))
    } else{
      sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
      if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){ #Some TP=0
        best.res <- findBestRes(p, reann, out.name, res.left, res.right, res.switch=0.05, outfile.path, type, sp.thresholds=sp.thresholds, method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75))
      } else if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0 || sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)>0){
        best.res <- findBestRes(p, reann, out.name, res.left, res.right, res.switch=0.05, outfile.path, type, sp.thresholds=sp.thresholds, method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75))}}
    p <- best.res$p
    reann <- best.res$reann
    res.left <- best.res$res.left
    res.right <- best.res$res.right

    message("Now we got next layer for ", cl)
    #What if reann=NULL??????
    annotation <- reann$ann.by.level$annotation$l1 %>% sapply(function(n) paste0(cl,"_", n)) %>% setNames(reann$ann.by.level$annotation$l1%>%names())
    ann.by.parents <- c(ann.by.parents, list(annotation)%>%setNames(cl))
    marker.list <- c(marker.list, list(reann$clf.data$marker.list)%>%setNames(cl))
    resol.last.layer <- c(resol.last.layer, list(c(res.left, res.right))%>%setNames(cl))
  }
  return(list(ann.by.parents=ann.by.parents, marker.list=marker.list, resol.last.layer=resol.last.layer))
}


#Look for the maximum resolution (e.g. the maximum number of predicted clusters) that have its predicted clusters all pass quality control (self projection accuracy via re-annotation: confusion rate)
findBestRes <- function(p, reann, out.name, res.left, res.right, res.switch=0.05, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", type="PCA", sp.thresholds=list(accuracy=0.8, confusion.rate=0.25), method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75)){
  #res.switch: only affect the starting of findBestRes, or the switch between findBestRes and findSmallestRes, when finding the max res that fullfill all the classifyign requirements.
  if ((res.right-res.left)<=res.switch) { 
    return(list(p=p, reann=reann, res.left=res.left, res.right=res.right))} 
  message("Now, we are heading to the left for the minimum working resolution ...")
  to.left <- findSmallestRes(p, reann, out.name, outfile.path, res.left, res.right, res.right.old=res.right, res.switch=0.05)
  p <-  to.left$p
  reann <- to.left$reann
  res.left <- to.left$res.right 
  res.right <- to.left$res.right.old
  if (length(reann$ann.by.level$annotation$l1 %>% table)==1) {return(list(p=p, reann=reann, res.left=res.left, res.right=res.right))}
  message("Now we check res.switch, res.right, res.left: ", res.switch, res.right, res.left)
  if ((res.right-res.left)>res.switch){
    message("Now, we are heading to the right for the maximum working resolution ...")
    to.right <- findBiggestRes(p, out.name, outfile.path, reann, reann.old=reann, res.left, res.left.old=res.left, res.right, res.switch=0.05, type, sp.thresholds=sp.thresholds)
    p <- to.right$p
    res.left <- to.right$res.left.old #the highest res that work
    res.right <- to.right$res.left #the lowest res that doesn't work
    reann <- to.right$reann.old}
  return(findBestRes(p, reann, out.name, res.left, res.right, res.switch=0.05))
}  

#There will be problem for this function, if the smallest res is very close (diff < 0.0001) to the original res.left, then the run may break.
#The error will look like: Error in ann.by.level[[i + 1]] : subscript out of bounds (reAnnPip: ann.by.level only has one cluster! singleton)
#Possible solution: increase relax the sp.threshold
findSmallestRes <- function(p, reann, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", res.left, res.right, res.right.old=res.right, res.switch=0.05, type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){ #At least one of the if has to be fullfiled by reann!!!!.
  #To avoid the boundary the maximum working resolution
  if ((res.right-res.left)<=res.switch) {#In this case, we are not sure whether the res.right works or not, so for safe side, we are going to consider it as not working.
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      reann$ann.by.level$annotation$l1 %<>% sapply(function(n) n=1)}
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph=NULL, clf.data=NULL, uncertainty.thresholds)
    return(list(p=p, reann=reann, res.left=res.left, res.right=res.left, res.right.old=res.right))} 
  if (is.null(reann)) {
    message("Now we will decrease resolution to half ...")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      #Here, reann is empty, we need to build a reann mually:
      reann <- list(ann.by.level=list(annotation=list(l1=p$clusters$PCA$leiden)), clf.data=list(marker.list=list())) #res.left not necessarily be singleton, so we'd better keep this this way
      return(list(p=p, reann=reann, res.left=res.left, res.right=res.left+(res.right-res.left)/2, res.right.old=res.right))}
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph, clf.data=NULL, uncertainty.thresholds)
    res.right.old=res.right
    res.right=res.left+(res.right-res.left)/2
    message("Now the new res.right becomes ", res.right)
    return(findSmallestRes(p, reann, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", res.left, res.right, res.right.old, res.switch=0.05))}
  sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
  if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){#To make sure no cluster dropping during re-annotation, which is either becasue bad marker selection (annPip requires all clusters have to have makers!!!)
    message("Now we will decrease resolution to half ...")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      reann$ann.by.level$annotation$l1 %<>% sapply(function(n) n=1)
      return(list(p=p, reann=reann, res.left=res.left, res.right=res.left+(res.right-res.left)/2, res.right.old=res.right))}
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph, clf.data=NULL, uncertainty.thresholds)
    res.right.old=res.right
    res.right=res.left+(res.right-res.left)/2
    message("Now the new res.right becomes ", res.right)
    return(findSmallestRes(p, reann, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", res.left, res.right, res.right.old, res.switch=0.05))}
  if (sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)>0){
    message("Now we will decrease resolution to half ...")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      reann$ann.by.level$annotation$l1 %<>% sapply(function(n) n=1)
      return(list(p=p, reann=reann, res.left=res.left, res.right=res.left+(res.right-res.left)/2, res.right.old=res.right))}
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph, clf.data=NULL, uncertainty.thresholds)
    res.right.old=res.right
    res.right=res.left+(res.right-res.left)/2
    message("Now the new res.right becomes ", res.right)
    return(findSmallestRes(p, reann, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", res.left, res.right, res.right.old, res.switch=0.05)) }
  if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0){
    message("Now we will decrease resolution to half ...")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      reann$ann.by.level$annotation$l1 %<>% sapply(function(n) n=1)
      return(list(p=p, reann=reann, res.left=res.left, res.right=res.left+(res.right-res.left)/2, res.right.old=res.right))}
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph, clf.data=NULL, uncertainty.thresholds)
    res.right.old=res.right
    res.right=res.left+(res.right-res.left)/2
    message("Now the new res.right becomes ", res.right)
    return(findSmallestRes(p, reann, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", res.left, res.right, res.right.old, res.switch=0.05))}
  return(list(p=p, reann=reann, res.left=res.left, res.right=res.right, res.right.old=res.right.old))
}      
 

#The starting reann should fufill the 3 if conditions in the function!!!
#Look to the right for the maximum separable resolution
findBiggestRes <- function(p, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", reann, reann.old=reann, res.left, res.left.old=res.left, res.right, res.switch=0.05, type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){
  message("Now we check res.switch, res.right, res.left: ", res.switch, res.right, res.left)
  if ((res.right-res.left)<=res.switch) {#In this case, we are not sure whether the res.right works or not, so for safe side, we are going to consider it as not working.
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left, n.iterations=10, name="leiden")
    if (length(p$clusters$PCA$leiden %>% table)==1){ #If the parent cluster is a singleton, stop!
      reann$ann.by.level$annotation$l1[] <- 1
    }
    
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", graph=NULL, clf.data=NULL, uncertainty.thresholds)
    return(list(p=p, reann=reann, reann.old=reann, res.left.old=res.left, res.left=res.right,  res.right=res.right))}
  if (!is.null(reann)){
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){#To make sure no cluster dropping during re-annotation, which is either becasue bad marker selection (annPip requires all clusters have to have makers!!!)
      if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)==0 && sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
        message("Now we will increase resolution to half ...")
        p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.right-(res.right-res.left)/2, n.iterations=10, name="leiden")
        #The findBestRes() function makes sure that no singleton should go in to this function, we don't have to test it.
        reann.old <- reann
        reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph, clf.data=NULL, uncertainty.thresholds)
        res.left.old=res.left
        res.left=res.right-(res.right-res.left)/2
        message("Now the new res.left becomes ", res.left)
        return(findBiggestRes(p, out.name, outfile.path="~/CellAnnotatoR-dev/notebooks/24files/", reann, reann.old, res.left, res.left.old, res.right, res.switch=0.05, type=type, n.iterations=n.iterations, name=name, graph=graph, clf.data= clf.data, uncertainty.thresholds=uncertainty.thresholds))
      } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
    } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
  } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
}     
  

#Get the leiden clustring based hierarchy:
#At the same time, we will also do re-annotation, so the outcome will be a ready hierarchy with re-annotation based on the selected markers
#layer.n: the number of layers that the users want to get. The more layers, the more clusters one will get, though one will not get endless number of clusters, becasue a. the feature genes will eventually become exhausted, b. the size of cluster can't be smaller than 6, or c. the cells within clusters become homogeneous eventually.
# If encountering error "Error in ann.by.level[[i + 1]] : subscript out of bounds", that means when running findBiggestRes(...), res.right approaches res.left too closely. Solution: set min.res.update=FALSE, and at the same time set res.switch=0
#res.max.update must be larger than res.switch, otherwise, you will be in trouble for the first step of findBestRes.
getLeidenHierarchy <- function(p2, out.name, outfile.path, layer.n=3, res.step.layer1=0.01, min.res=0, max.res.layer2=1, max.res.increase=1.5, res.switch=0.05, clustering.type=NULL, embeding.type=NULL, graph=NULL, res.max.update=1, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){
  if (res.max.update <= res.switch) {Error("res.max.update must be larger than res.switch")}
  message("Get Knn clusters for layer 1 ...")
  #Predict cluster for layer 1, the point is to get a big picture, so therefore we will accept the first non-singletong KNN clustering result here.
  message("Now we are going to get layer 1 ...")
  l1 <- getLayer1KnnClusters(p2, res.start=0.01, res.step=res.step.layer1)
  #Quality control of layer1 clustering by SELF PROJECTION USING re-annotation based on the selected markers
  ##1. Re-annotation
  message("Reannotation for layer 1 ...")
  reann <- reAnnoPip(l1$p2, list(l1=l1$p2$clusters$PCA$leiden), out.name, outfile.path, graph=NULL, clf.data=NULL, uncertainty.thresholds=uncertainty.thresholds)
  ##2. Self-projection quality control for re-annotation (one way to predict cell types based on markers)
  if (is.null(reann)) {stop("There is no subclusters in the cell population, or try more gene features...")} #1. reann shouldn't be NA (as long as not all clusters have predicted markers, then reann will be NA)
  sp <- selfProjQCperCluster(l1$p2$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
  if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){ #2. sp shouldn't have NaN or Inf
    stop("There is no subclusters in the cell population, or try more gene features...")}
  #3. All the clusters should hava a minimum size of 6 after re-annotation
  #4. The self projection accuracy by re-annotation has to be >=0.6 (or 0.8) and the confusion rate has to be <=2/3 (0.25)
  if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0 || sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)>0){
    stop("There is no subclusters in the cell population, or try more gene features...")}
  #The 4 requirements for accepting clustering results also apply to other layers.
  # To predict clusters for the rest of layers:
  annotation1 <- reann$ann.by.level$annotation$l1 %>% sapply(function(n) paste0("l",n)) %>% setNames(reann$ann.by.level$annotation$l1%>%names())
  annotation <- annotation1 #In then end, remove annotation1
  next.layers <- list()
  max.res.middle = max.res.layer2
  for (l in 2:layer.n){
    if (length(annotation)>0){ #This (together with the annotation assignment below) makes sure that the program will not continue clustering if the last two layers have the same clustering.
      message("Now we are going to get layer ", l, "...")
      next.layer <- getNextLayerKnnClusters(p2, annotation, out.name, outfile.path, min.res.start=min.res, graph=NULL, res.max.update=1, clf.data=NULL, 
                    uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), max.res.middle, res.switch, type="PCA", 
                    method=conos::leiden.community, n.iterations=10, name="leiden", sp.thresholds)
      #Collect all non-singleton clusters and continue with next round of clustering
      annotation <- split(next.layer$ann.by.parents %>% Reduce(c,.), annotation[names(next.layer$ann.by.parents %>% Reduce(c,.))]) %>%
                      .[sapply(.,function(k) (table(k)%>%length)>1)] %>% Reduce(c,.)
      annotation <- annotation[annotation %in% (table(annotation)%>%.[.>=12]%>%names)]
      max.res.middle <- (next.layer$resol.last.layer%>%sapply(function(n) n[2]) %>% max) + max.res.increase 
      next.layers <- c(next.layers, list(next.layer)%>%setNames(paste0("l", l)))
    }
  }
  ann.by.level <- list()
  marker.list <- list()
  ann.by.level <- c(ann.by.level, list(annotation1)%>%setNames("l1"))
  #Take out ann.by.parents from other layers
  ann.by.level.r <- next.layers%>% lapply(function(n) n$ann.by.parents %>% Reduce(c,.))
  ann.by.level <- c(ann.by.level, ann.by.level.r)
  for (i in 1:(length(ann.by.level)-1)){ #In order to use split, all the layers should have the same cell order
    ann.by.level[[i+1]] <- c(ann.by.level[[i+1]], ann.by.level[[i]][names(ann.by.level[[i]]) %>% setdiff(names(ann.by.level[[i+1]]))])}
  
  #Check if the last 2 layers have the same clustering, if it is, drop the last layer
  if (length(ann.by.level %>%.[[length(.)]]%>%unique)==length(ann.by.level%>%.[[length(.)-1]]%>%unique)){
    ann.by.level[length(ann.by.level)] <- NULL}
  
  #Make marker.list
  marker.list <- c(marker.list, list(reann1$clf.data$marker.list)%>%setNames("l1")) # For the 1st layer
  marker.list <- c(marker.list, next.layers%>% lapply(function(n) n$marker.list) %>% setNames(names(next.layers))) #For the rest layers

  return(list(ann.by.level=ann.by.level, marker.list=marker.list))
}








#Self Projection Quality Control
#annotation.before <- p4$clusters$PCA$leiden #annotation before re-annotation (reAnnoPip)
#annotation.after <- reann.p4$ann.by.level$annotation$l1 #annotation after re-annotation (reAnnoPip)
#This is used for the self projection quaility control after re-annotating the leiden-clustering-based hierarchy
selfProjQCperCluster <- function(annotation.before, annotation.after) {
  #Get self projection accuracy:
  acc <- sum(annotation.before == annotation.after)/length(annotation.before) #The annotation before and after need to have the same length.
  #Get confusion rate (FP/TP): the cluster sample sized is balanced by adjusting weights inversely proportional to class frequencies in the input data
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


