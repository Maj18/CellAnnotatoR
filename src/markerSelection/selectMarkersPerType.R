#' @export
# I have modified the selectMarkersPerType(): replace cm.norm() with t(cm.norm()) in mean.unc.per.type <- (1 - getMeanConfidencePerType(marker.list, t(cm.norm), annotation)).
# I have also replaced unique(unlist(markers.per.type)) with unique(c(unlist(markers.per.type$positive),unlist(markers.per.type$negative)))
selectMarkersPerType <- function(cm.norm, annotation, markers.per.type, marker.list=NULL, max.iters=ncol(cm.norm), parent=NULL,
                                 max.uncertainty=0.25, verbose=0, min.pos.markers=1, max.pos.markers=10, log.step=1, n.cores=1, refinement.period=10, return.all=F) {
  if (verbose > 0) message("Running marker selection for parent type '", parent, "'")
  
  for (n in names(marker.list)) {
    marker.list[[n]]$locked <- T
  }
  
  marker.list %<>% prepareInitialMarkerList(names(markers.per.type$positive), parent=parent)
##names(ann_by_parent3[["l1_1"]])
##marker.list %<>% prepareInitialMarkerList(names(markers.per.type$positive), parent=parent)
# markers.per.type <- preSelectMarkerCandidates(de_per_parent3[["l1_15"]], blacklist=NULL)
#The if condition has been added to take into considertion of those cell type that has no marker genes identified at all
  if (is.null(unique(c(unlist(markers.per.type$positive),unlist(markers.per.type$negative))))){
    cm.norm <- matrix(,ncol=1, nrow=length(names(annotation)))
    rownames(cm.norm) <- names(annotation)} else{
      cm.norm %<>% .[names(annotation), unique(c(unlist(markers.per.type$positive),unlist(markers.per.type$negative)))] %>% as.matrix()
    }
#old:cm.norm %<>% .[unique(c(unlist(markers.per.type$positive),unlist(markers.per.type$negative))), names(annotation)]
##cm.norm <- matrix(,nrow=1, ncol=length(names(ann_by_parent3[["l1_1"]])))
##colnames(cm.norm) <- names(ann_by_parent3[["l1_1"]])
##cm_norm_annotation %<>% .[names(ann_by_parent3[["l1_15"]]), unique(c(as.character(unlist(markers.per.type$positive)),as.character(unlist(markers.per.type$negative))))] %>% as.matrix()
  mean.unc.per.type <- (1 - getMeanConfidencePerType(marker.list, cm.norm, annotation))
#Here, the original script is "mean.unc.per.type <- (1 - getMeanConfidencePerType(marker.list, cm.norm, annotation))" 
#This step is get the mean uncertainty (1-s(t)) score for each cell type**********
##mean.unc.per.type <- (1 - getMeanConfidencePerType(marker.list, t(cm.norm), ann_by_parent3[["l1_1"]]))
##max.iters = ncol(cm_norm_annotation)  # to get the number of DE genes
  did.refinement <- F
  for (i in 1:max.iters) {
    n.markers.per.cell <- sapply(marker.list, function(x) length(x$expressed))[names(mean.unc.per.type)]
    type.mask <- (n.markers.per.cell < max.pos.markers) & (sapply(markers.per.type$positive, length)[names(mean.unc.per.type)] > 0)
    if (type.mask== 0)
      break
    cell.type <- mean.unc.per.type[type.mask] %>% which.max() %>% names()
    m.update <- getNextMarkers(cell.type, cm.norm, annotation, marker.list=marker.list, markers.per.type=markers.per.type,
                               verbose=(verbose > 1), n.cores=n.cores)
    #m.update <- getNextMarkers(cell.type, t(cm_norm_annotation), ann_by_parent3[["l1_15"]], marker.list=marker.list, markers.per.type=markers.per.type,verbose=(verbose > 1), n.cores=n.cores)
    marker.list.new <- marker.list
    marker.list.new[[cell.type]]$expressed %<>% union(m.update$expressed)
    marker.list.new[[cell.type]]$not_expressed %<>% union(m.update$not_expressed)
    
    markers.per.type %<>% updateMarkersPerType(marker.list=setNames(list(m.update), cell.type))
    
    mean.unc.per.type.new <- (1 - getMeanConfidencePerType(marker.list.new, cm.norm, annotation))
    if (mean(mean.unc.per.type.new) < mean(mean.unc.per.type)) {
      mean.unc.per.type <- mean.unc.per.type.new
      marker.list <- marker.list.new
      did.refinement <- F
    }
    
    if (verbose && (log.step > 0) && (i %% log.step == 0)) {
      message("Iteration ", i, ". Max uncertainty: ", round(max(mean.unc.per.type), 3), ", mean uncertainty: ", round(mean(mean.unc.per.type), 3),
              ". Target type: ", cell.type, ", gene: ", m.update$expressed)
    }
    
    if ((max(mean.unc.per.type) < max.uncertainty) && (sapply(marker.list, function(x) length(x$expressed)) >= min.pos.markers))
      break
    
    if ((refinement.period > 0) && (i %% refinement.period == 0) && !did.refinement) {
      if (verbose) message("Refine markers...")
      marker.list %<>% filterMarkerListByScore(cm.norm, annotation, verbose=(verbose > 1), n.cores=n.cores)
      did.refinement <- T
    }
  }
  
  if ((refinement.period != 0) && !did.refinement) {
    if (verbose) message("Refine markers...")
    marker.list %<>% filterMarkerListByScore(cm.norm, annotation, verbose=(verbose > 1), n.cores=n.cores)
  }
  
  for (n in names(marker.list)) {
    marker.list[[n]]$locked <- NULL
  }
  
  if (return.all)
    return(list(marker.list=marker.list, markers.per.type=markers.per.type))
  
  return(marker.list)
}