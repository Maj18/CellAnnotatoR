
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
getLayer1KnnClusters <- function(p2, res.step=0.01){ 
  p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=0.01, n.iterations=50, name="leiden")
  clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  n = 0
  while (clusters.n==1){
    n = n+1
    p2$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=0.01+n*res.step, n.iterations=50, name="leiden")
    clusters.n <- p2$clusters$PCA$leiden %>% unique %>% length
  }
  return(list(p2=p2, res=0.01+n*res.step))
}
#e.g. l2 <- getLayer1Clusters(p2, res.step=0.01)
#l2$res
#l2$p2$clusters$PCA$leiden %>% table

#Get the last layer of leiden clusters (>1 clusters) for the leiden clustering based hierarchy
getNextLayerKnnClusters <- function(p2, annotation, out.name, outfile.path, min.res, graph=NULL, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), max.res.middle=1, res.switch=0.5, type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){ 
  ann.by.cluster <- split(annotation, annotation)
  ann.by.parents <- list()
  marker.list <- list()
  resol.last.layer <- list()
  for (cl in names(ann.by.cluster)){
    message("Now we will get next layer for ", cl)
    p <- vpscutils::GetPagoda(p2$misc$rawCounts[ann.by.cluster[cl][[1]]%>%names,]%>%Matrix::t(), clustering.type=clustering.type, embeding.type=embeding.type)
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=max.res.middle, n.iterations=10, name="leiden")
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph=graph, clf.data= clf.data, uncertainty.thresholds=uncertainty.thresholds)
   
    if (!is.null(reann)){
      sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
      if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
        if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)==0 && sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
          stop("The maximum resolution/max.res.middle is too low, suggest to increase max.res.middle if it is for layer 2, otherwise try to increase max.res.increase!")}}}
    res.left <- min.res
    res.right <- max.res.middle
    if (is.null(reann)){
      best.res <- findBestRes(p, reann, res.left, res.right, res.switch, type, sp.thresholds=sp.thresholds)
    } else{
      sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
      if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){
        best.res <- findBestRes(p, reann, res.left, res.right, res.switch, type, sp.thresholds=sp.thresholds)
      } else if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0 || sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)>0){
        best.res <- findBestRes(p, reann, res.left, res.right, res.switch, type, sp.thresholds=sp.thresholds)}}
    p <- best.res$p
    reann <- best.res$reann
    res.left <- best.res$res.left
    res.right <- best.res$res.right

    message("Now we got next layer for ", cl)
    annotation <- reann$ann.by.level$annotation$l1 %>% sapply(function(n) paste0(cl,"_",n)) %>% setNames(reann$ann.by.level$annotation$l1%>%names())
    ann.by.parents <- c(ann.by.parents, list(annotation)%>%setNames(cl))
    marker.list <- c(marker.list, list(reann$clf.data$marker.list)%>%setNames(cl))
    resol.last.layer <- c(resol.last.layer, list(c(res.left, res.right))%>%setNames(cl))
  }
  return(list(ann.by.parents=ann.by.parents, maker.list=marker.list, resol.last.layer=resol.last.layer))
}

#Look for the maximum resolution (e.g. the maximum number of predicted clusters) that have its predicted clusters all pass quality control (self projection accuracy via re-annotation: confusion rate)
findBestRes <- function(p, reann, res.left, res.right, res.switch, type="PCA", sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){
  #res.switch: only affect the starting of findBestRes, or the switch between findBestRes and findSmallestRes, when finding the max res that fullfill all the classifyign requirements.
  if ((res.right-res.left)<=res.switch) {return(list(p=p, reann=reann, res.left=res.left, res.right=res.right))}
  if ((res.right-res.left)>res.switch){
    message("Now, we are heading to the left for the minimum working resolution ...")
    to.left <- findSmallestRes(p, reann, out.name, outfile.path, res.left, res.right, res.right.old=res.right, sp.thresholds=sp.thresholds)
    p <-  to.left$p
    reann <- to.left$reann
    res.left <- to.left$res.right 
    res.right <- to.left$res.right.old
    if ((res.right-res.left)>res.switch){
      message("Now, we are heading to the right for the maximum working resolution ...")
      to.right <- findBiggestRes(p, out.name, outfile.path, reann, reann.old=reann, res.left, res.left.old=res.left, res.right, type, sp.thresholds=sp.thresholds)
      p <- to.right$p
      res.left <- to.right$res.left.old #the highest res that work
      res.right <- to.right$res.left #the lowest res that doesn't work
      reann <- to.right$reann.old}
    return(findBestRes(p, reann, res.left, res.right, res.switch))}
}  

#There will be problem for this function, if the smallest res is very close (diff < 0.0001) to the original res.left, then the run may break.
#The error will look like: Error in ann.by.level[[i + 1]] : subscript out of bounds 
#Possible reasons for this is that he sp.threshold is too stringent (e.g. when acc=0.9 and conf.rate=1/9)
#Need to find a way to solve this problem????????????????????????????????????????
#Possible solution: increase relax the sp.threshold
findSmallestRes <- function(p, reann, out.name, outfile.path, res.left, res.right, res.right.old=res.right, type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){ #At least one of the if has to be fullfiled by reann!!!!.
  if (is.null(reann)) {
    message("Now we will decrease resolution to half ...")
    p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
    reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph=graph, clf.data= clf.data, uncertainty.thresholds=uncertainty.thresholds)
    res.right.old=res.right
    res.right=res.left+(res.right-res.left)/2
    message("Now the new res.right becomes ", res.right)
    return(findSmallestRes(p, reann, out.name, outfile.path, res.left, res.right, res.right.old)) 
  } else{
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){#To make sure no cluster dropping during re-annotation, which is either becasue bad marker selection (annPip requires all clusters have to have makers!!!)
      message("Now we will decrease resolution to half ...")
      p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
      reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph, clf.data, uncertainty.thresholds)
      res.right.old=res.right
      res.right=res.left+(res.right-res.left)/2
      message("Now the new res.right becomes ", res.right)
      return(findSmallestRes(p, reann, out.name, outfile.path, res.left, res.right, res.right.old))}
    else if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0 || sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)>0){
      message("Now we will decrease resolution to half ...")
        p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.left+(res.right-res.left)/2, n.iterations=10, name="leiden")
        reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph=graph, clf.data, uncertainty.thresholds)
        res.right.old=res.right
        res.right=res.left+(res.right-res.left)/2
        message("Now the new res.right becomes ", res.right)
        return(findSmallestRes(p, reann, out.name, outfile.path, res.left, res.right, res.right.old)) }
  } 
  return(list(p=p, reann=reann, res.left=res.left, res.right=res.right, res.right.old=res.right.old))
}      
 

#The starting reann should fufill the 3 if conditions in the function!!!
findBiggestRes <- function(p, out.name, outfile.path, reann, reann.old=reann, res.left, res.left.old=res.left, res.right, type="PCA", method=conos::leiden.community, n.iterations=50, name="leiden", graph=NULL, clf.data= NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){
  if (!is.null(reann)){
    sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
    if (!is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){#To make sure no cluster dropping during re-annotation, which is either becasue bad marker selection (annPip requires all clusters have to have makers!!!)
      if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)==0 && sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
        message("Now we will increase resolution to half ...")
        p$getKnnClusters(type="PCA", method=conos::leiden.community, resolution=res.right-(res.right-res.left)/2, n.iterations=10, name="leiden")
        reann.old <- reann
        reann <- reAnnoPip(p, list(l=p$clusters$PCA$leiden), out.name, outfile.path, graph=graph, clf.data= clf.data, uncertainty.thresholds=uncertainty.thresholds)
        res.left.old=res.left
        res.left=res.right-(res.right-res.left)/2
        message("Now the new res.left becomes ", res.left)
        return(findBiggestRes(p, out.name, outfile.path, reann, reann.old, res.left, res.left.old, res.right, type=type, n.iterations=n.iterations, name=name, graph=graph, clf.data= clf.data, uncertainty.thresholds=uncertainty.thresholds))
      } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
    } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
  } else {return(list(p=p, reann=reann, reann.old=reann.old, res.left=res.left, res.left.old=res.left.old, res.right=res.right))}
}     
  
    

#Get the leiden clustring based hierarchy:
#At the same time, we will also do re-annotation, so the outcome will be a ready hierarchy with re-annotation based on the selected markers
#layer.n: the number of layers that the users want to get. The more layers, the more clusters one will get, though one will not get endless number of clusters, becasue a. the feature genes will eventually become exhausted, b. the size of cluster can't be smaller than 6, or c. the cells within clusters become homogeneous eventually.
# If encountering error "Error in ann.by.level[[i + 1]] : subscript out of bounds", that means when running findBiggestRes(...), res.right approaches res.left too closely. Solution: set min.res.update=FALSE, and at the same time set res.switch=0
getLeidenHierarchy <- function(p2, out.name, outfile.path, layer.n=3, res.step.layer1=0.01, min.res=0, min.res.update=FALSE, max.res.layer2=1, max.res.increase=4, res.switch=0.1, clustering.type=NULL, embeding.type=NULL, graph=NULL, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), sp.thresholds=list(accuracy=0.8, confusion.rate=0.25)){
  message("Get Knn clusters for layer 1 ...")
  #Predict cluster for layer 1, the point is to get a big picture, so therefore we will accept the first non-singletong KNN clustering result here.
  message("Now we are going to get layer 1 ...")
  l1 <- getLayer1KnnClusters(p2, res.step=res.step.layer1)
  #Quality control of layer1 clustering by SELF PROJECTION USING re-annotation based on the selected markers
  ##1. Re-annotation
  message("Reannotation for layer 1 ...")
  reann1 <- reAnnoPip(l1$p2, list(l1=l1$p2$clusters$PCA$leiden), out.name, outfile.path, graph=NULL, clf.data=NULL, uncertainty.thresholds=uncertainty.thresholds)
  ##2. Self-projection quality control for re-annotation (one way to predict cell types based on markers)
  if (is.null(reann)) {stop("There is no subclusters in the cell population, or try more gene features...")} #1. reann shouldn't be NA (as long as not all clusters have predicted markers, then reann will be NA)
  sp <- selfProjQCperCluster(p$clusters$PCA$leiden, reann$ann.by.level$annotation$l1)
  if (is.nan(sum(sp$confusion.rate.false.pos.rate[1,]))){ #2. sp shouldn't have NaN or Inf
    stop("There is no subclusters in the cell population, or try more gene features...")}
  #3. All the clusters should hava a minimum size of 6 after re-annotation
  #4. The self projection accuracy by re-annotation has to be >=0.6 (or 0.8) and the confusion rate has to be <=2/3 (0.25)
  if (sum(reann$ann.by.level$annotation$l1 %>% table()<6)>0 || sum(sp$confusion.rate.false.pos.rate[1,]>sp.thresholds$confusion.rate)==0){
    stop("There is no subclusters in the cell population, or try more gene features...")}
  #The 4 requirements for accepting clustering results also apply to other layers.
  # To predict clusters for the rest of layers:
  annotation1 <- reann1$ann.by.level$annotation$l1 %>% sapply(function(n) paste0("l_",n)) %>% setNames(reann1$ann.by.level$annotation$l1%>%names())
  annotation <- annotation1 #In then end, remove annotation1
  next.layers <- list()
  min.res <- if (min.res.update) l1$res else 0
  max.res.middle = max.res.layer2
  for (l in 2:layer.n){
    message("Now we are going to get layer ", l, "...")
    next.layer <- getNextLayerKnnClusters(p2, annotation, out.name, outfile.path, min.res, graph=NULL, clf.data=NULL, uncertainty.thresholds=c(coverage=0.5, negative=0.5, positive=0.75), max.res.middle, res.switch, type="PCA", method=conos::leiden.community, n.iterations=10, name="leiden", sp.thresholds)
    annotation <- next.layer$ann.by.parents %>% Reduce(c,.)
    min.res <- if (min.res.update) next.layer$resol.last.layer%>%sapply(function(n) n[1]) %>% min/2 else 0
    max.res.middle <- (next.layer$resol.last.layer%>%sapply(function(n) n[2]) %>% max) + max.res.increase 
    res.switch <- 0.5*res.switch
    #sp.thresholds <- list(accuracy=pmax(0.6, sp.thresholds$accuracy-0.1), confusion.rate=pmin(2/3, sp.thresholds$confusion.rate+1/9))
    next.layers <- c(next.layers, list(next.layer)%>%setNames(paste0("l_", l)))
  }

  #Make ann.by.level for entire leiden hierarchy:
  ann.by.level <- list()
  marker.list <- list()
  ann.by.level <- c(ann.by.level, list(annotation1)%>%setNames("l1"))
  #For other layers: ann.by.parents=ann.by.parents, maker.list=marker.list, resol.last.layer=resol.last.layer
  ann.by.parents.layers <- next.layers%>% lapply(function(n) lapply(names(n$ann.by.parents), function(p) sapply(n$ann.by.parents[[p]], function(m) 
            paste0(p, "_", m))%>%setNames(n$ann.by.parents[[p]]%>%names)))%>% 
            setNames(names(next.layers)) #a list of layers, each layer has a list of parents: Cells, named; each parent, named.
  
  for (l in names(ann.by.parents.layers)){
    annotation.l <- ann.by.parents.layers[[l]] %>% Reduce(c,.)  #in the end, annotation.l should also be removed
    ann.by.level <- c(ann.by.level, list(annotation.l)%>%setNames(l))
  }
  
  #Make marker.list
  marker.list <- c(marker.list, list(reann1$clf.data$marker.list)%>%setNames("l1")) # For the 1st layer
  marker.list <- c(marker.list, next.layers%>% lapply(function(n) n$marker.list) %>% setNames(names(next.layers))) #For the rest layers
  return(lsit(ann.by.level=ann.by.level, marker.list=marker.list))
}

# Building leiden clustering based hirarchy finish here  







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


