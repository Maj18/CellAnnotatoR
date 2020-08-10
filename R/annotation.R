#' @useDynLib CellAnnotatoR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %<>% %$% %>%
NULL

## Correct unloading of the library
.onUnload <- function (libpath) {
  library.dynam.unload("CellAnnotatoR", libpath)
}

expandAnnotationToClusters <- function(scores, clusters) {
  clusters <- droplevels(clusters[rownames(scores)])
  ann.update <- split(1:length(clusters), clusters) %>%
    sapply(function(cids) colSums(scores[cids, ,drop=F]) %>% which.max() %>% names()) %>%
    .[clusters] %>% setNames(names(clusters))

  return(ann.update)
}

annotationFromScores <- function(scores, clusters=NULL) {
  if (!is.null(clusters))
    return(expandAnnotationToClusters(scores, clusters))

  return(colnames(scores)[apply(scores, 1, which.max)] %>% setNames(rownames(scores))) #apply(scores, 1, which.max): in the data.frame scores, for each row, which is the biggest one.
} #For the two kids "l22_7"  "l24_18" under parent l13_7, all cells have been assigned tol24_18, that's why l22_7 was dropped!

#' Diffuse Score per Type
#' @inheritDotParams diffuseGraph fading fading.const verbose tol score.fixing.threshold
diffuseScorePerType <- function(scores.per.type, graph, parents, cbs.per.type, verbose, n.cores=1, ...) {
  plapply(parents, function(p)
    diffuseGraph(igraph::induced_subgraph(graph, cbs.per.type[[p]]),
                 scores=scores.per.type[[p]], verbose=verbose, ...),
    verbose=(verbose > 0), n.cores=n.cores)
}

#' Diffuse Graph
#' @description Run diffusion on graph
#' @param graph graph to diffuse on
#' @param scores table of scores
#' @param fading fading level for graph diffusion
#' @param fading.const constant in exponent for graph diffusion
#' @param score.fixing.threshold threshold for a label to be considered certain. Such labels can't be changed during diffusion.
#' @param verbose print progress bar
#' @param tol tolerance for diffusion stopping
diffuseGraph <- function(graph, scores, fading=10, fading.const=0.5, score.fixing.threshold=0.8,
                         verbose=FALSE, max.iters=1000, tol=1e-3) {
  cbs <- igraph::V(graph)$name
  if (length(cbs) == 0)
    return(NULL)

  scores <- as.matrix(scores[cbs, ])
  scores[rowSums(scores) < 1e-5, ] <- 1
  scores %<>% `/`(rowSums(.))

  if (any(is.na(scores)))
    stop("NAs in scores")

  edges <- igraph::as_edgelist(graph)

  is.fixed <- (apply(scores, 1, max) > score.fixing.threshold)

  if (nrow(edges) == 0)
    return(scores)

  edge.weights <- igraph::edge.attributes(graph)$weight
  res <- conos:::smoothMatrixOnGraph(edges, edge.weights, scores, is.label.fixed=is.fixed, max_n_iters=max.iters,
                                     diffusion_fading=fading, diffusion_fading_const=fading.const, verbose=verbose,
                                     tol=tol, normalize=T)
  return(res)
}

#' Assign Cells By Scores
#' If we are only to reannoate one layer of a hierarchy, or the St_norm got over entire layer, we can still use this function.
#' If we want to reannotate the entire hierarchy, the St_norm of which acquired under the same parents, then this funciton is NOT PROPER ANYMORE.
#' @description Assign cell types for each cell based on type scores. Optionally uses `clusters` to expand annotation.
#'
#' @param graph cell graph from Seurat, Pagoda2 or some other tool. Can be either in igraph or adjacency matrix format.
#'    Use `graph=NULL` to skip graph diffusion step and get raw score annotation (useful when debug marker genes).
#' @param score.info cell type scores from `getMarkerScoreInfo` function. Re-estimated if NULL
#' @param clusters vector with cluster labels named by cell ids. Used to expand annotation on these clusters.
#' @param verbose verbosity level (from 0 to 2)
#' @param max.depth maximal depth for which annotation is done. Useful during manual marker selection
#' @inheritParams getMarkerScoreInfo
#' @inheritDotParams diffuseScorePerType
#' @return list with parameters:\itemize{
#'   \item{annotation: annotation per level}
#'   \item{scores: assignment scores per level}
#'   \item{annotation.filt: the same as annotation, but cells, which don't pass QC are assigned to NA class}
#' }
#'
#' @examples
#'   clf_data <- getClassificationData(cm, marker_path)
#'   ann_by_level <- assignCellsByScores(graph, clf_data, clusters=clusters)
#'   Li: set default coverage=0.93, rather than 0.5, this can make sure min(St)=0.1$$$$$$$$$$$$$$$$$$
#' @export
assignCellsByScores <- function(graph, clf.data, score.info=NULL, clusters=NULL, verbose=0, uncertainty.thresholds=c(coverage=0.05, negative=0.5, positive=0.75), max.depth=NULL, ...) { #What this function does, is to compare the normalized St scores among the kid clusters under the same parent, and for each cell the kid with the highest score will be assigned to that cells, that is how re-annotation will be down, and filtering is completely based on the uncertainty threshold given, if larger, then ti will be turned own (NA). The recommended threshold for coverage is 0.93 that corresponds to a St of 0.1 (non-normalized), which we assume to be minimum to be considered as the cell expressing the marker genes.
  if (!is.null(clusters)) {
    clusters <- as.factor(clusters)
  }

  if (is.null(score.info)) {
    score.info <- getMarkerScoreInfo(clf.data)
  }

  scores <- getMarkerScoresPerCellType(clf.data, score.info=score.info) #here, St are normalzed under parents #len=90
  if (!is.null(graph)) {
    if ((is(graph, "Matrix") || is(graph, "matrix")) && ncol(graph) == nrow(graph)) {
      graph <- igraph::graph_from_adjacency_matrix(graph, weighted=T)
    } else if (!is(graph, "igraph")) {
      stop("Unknown graph format. Only adjacency matrix or igraph are supported")
    }

    if (length(setdiff(igraph::V(graph)$name, rownames(scores))) > 0)
      stop("Not all cells from the graph are presented in clf.data")

    if (length(setdiff(rownames(scores), igraph::V(graph)$name)) > 0) {
      warning("Not all cells from the clf.data are presented in the graph. Omitting ",
              nrow(scores) - length(igraph::V(graph)$name), " cells")
      scores %<>% .[igraph::V(graph)$name,]
    }
  }

  subtypes.per.depth.level <- classificationTreeToDf(clf.data$classification.tree) %$% #len=90
    split(., PathLen) %>% lapply(function(df) split(df$Node, df$Parent)) %>% .[order(names(.))] #Parent-kid relationship

  max.depth <- if (is.null(max.depth)) length(subtypes.per.depth.level) else min(max.depth, length(subtypes.per.depth.level))
  c.ann <- rep("root", nrow(scores)) %>% setNames(rownames(scores))
  c.ann.filt <- c.ann
  ann.by.level <- list()
  ann.filt.by.level <- list()
  scores.by.level <- list()
  scores.posterior <- scores

  possible.ann.levels <- c()
  for (pl in 1:max.depth) {
    if (verbose > 0) message("Level ", pl, "...")

    c.subtypes.per.parent <- subtypes.per.depth.level[[pl]] #parent-kids
    possible.ann.levels %<>% c(unlist(c.subtypes.per.parent)) %>% unique() #kid types

    c.parents <- names(c.subtypes.per.parent) %>% setNames(., .)
    cbs.per.type <- lapply(c.parents, function(p) names(c.ann)[c.ann == p]) %>% .[sapply(., length) > 0]
    #c.ann comes from the last cycle, so the annotation has been updated cycle by cycle, however the scores remains the same for the original tree, don't we need to update the score as well?
       #all cells for each type at this level

    if (length(cbs.per.type) == 0) # some subtypes have deeper level of annotation, but non of them is found in the dataset
      break

    c.parents %<>% .[names(cbs.per.type)]

    scores.per.type <- lapply(c.parents, function(p) scores[cbs.per.type[[p]], c.subtypes.per.parent[[p]], drop=F])
#scores is a dataframe, and there is no "root", therefore we use the kids to locate it. Here we get scores per type, the cells are also right.
    if (!is.null(graph)) {
      scores.per.type %<>% diffuseScorePerType(graph, c.parents, cbs.per.type, verbose=(verbose > 1), ...)

      for (cs in scores.per.type) {
        scores.posterior[rownames(cs), colnames(cs)] <- cs
      }
    }
    #Here, the name of scores.per.type is still "root", but the scores shown are for it's kids
    res.ann <- lapply(scores.per.type, annotationFromScores, clusters) %>% Reduce(c, .) #this is an important step, here we use St to reannotate the cells: which kid cluster each cell belongs to $$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #HERE, THE CELL TYPES DROPPED?????
    #For the two kids "l22_7"  "l24_18" under parent l13_7, all cells have been assigned tol24_18, that's why l22_7 was dropped!
    res.ann.filt <- filterAnnotationByUncertainty(res.ann, scores.posterior[,possible.ann.levels], score.info=score.info,
                                                  cur.types=unique(res.ann), clusters=clusters, thresholds=uncertainty.thresholds)
    #Here, we have quality check.
    #scores.posterior[,possible.ann.levels]: The kids' St scores (normalized)
    #score.info=score.info: the St (not normalized), Sp, Sn scores
    #cur.types: the kid types
    #res.ann: all cells (could be less, depends on what parent type is)
    #scores: each type has all cells, haven't withdrawed cells for types.

    c.ann[names(res.ann)] <- res.ann
    c.ann.filt[names(res.ann.filt)] %<>% is.na() %>% ifelse(NA, res.ann.filt)

    level.name <- paste0("l", pl)
    ann.by.level[[level.name]] <- c.ann
    ann.filt.by.level[[level.name]] <- c.ann.filt
    scores.by.level[[level.name]] <- scores.posterior[,possible.ann.levels] #Those are the normalized scores
    if (verbose > 0) message("Done")
  }

  return(list(annotation=ann.by.level, scores=scores.by.level, annotation.filt=ann.filt.by.level))
}


## Score assignment

#' Normalize Total Count with Features
#' @description Normalize `cm` matrix using total counts and then column-wise min-max scaling
#'
#' @param max.quantile quantile to be used for max estimation during scaling
#' @return Normalized matrix of the same shape as `cm`
#' @inheritParams getClassificationData
#'
#' @export
normalizeTCWithFeatures.ORIG <- function(cm, max.quantile=0.95, max.smooth=1e-10, transpose=T) {
  cm %<>% as("dgCMatrix") %>% Matrix::drop0()

  # Total count normalization (i.e. TF-step) (divided by all genes' expression sum within a cell)
  cm@x <- cm@x / rep(Matrix::colSums(cm), diff(cm@p))
  cm <- Matrix::t(cm) #Now rows become cells, col becomes genes

  # Factors for min-max gene normalization
  max.vals <- split(cm@x, rep(1:(length(cm@p)-1), diff(cm@p)))[paste0(1:ncol(cm))] # a list of vectors, each vector is for a gene, vector lengths differ between genes.
  max.vals[is.null(max.vals)] <- c()
  max.vals %<>% sapply(quantile, max.quantile) %>% `+`(max.smooth) # Robust alternative to maximum #The 95% quantile of all non-0 values/expressions for each gene. #maximum could be an outlier.

  cm@x <- cm@x / rep(max.vals, diff(cm@p)) # fast way to do columnwise normalization for sparse matrices
  cm@x %<>% pmin(1.0)

  if (!transpose)
    return(cm)

  return(Matrix::t(cm))
}

#' Normalize TF-IDF with Features
#' @description Normalize `cm` matrix using total counts, column-wise min-max scaling and then IDF normalization
#' @inheritDotParams normalizeTCWithFeatures max.quantile max.smooth
#'
#' @export
normalizeTfIdfWithFeatures.ORIG <- function(cm, ...) {
  tf.idf <- normalizeTCWithFeatures.ORIG(cm, transpose=F, ...) #HERE .ORIG IS ADDED
  # IDF-factors: log(1 + fraction of expressing cells)
  idf.weights <- log(1 + nrow(tf.idf) / (Matrix::colSums(tf.idf > 0) + 1)) #row: cells, col: genes
  tf.idf@x <- tf.idf@x * rep(idf.weights, diff(tf.idf@p))
  return(Matrix::t(tf.idf))
}

#' Get Marker Score Info
#' @description estimate info, neccessary for scoring of cells by cell types
#'
#' @param clf.data classification data from `getClassificationData`
#' @inheritDotParams getCellTypeScoreInfo aggr
#' @return List of score info for each cell type. See `CellAnnotatoR:::getCellTypeScoreInfo` for more info
#' @export
getMarkerScoreInfo <- function(clf.data, ...) {
  lapply(clf.data$marker.list, getCellTypeScoreInfo, clf.data$cm, ...)
}

#' Merge Scores
#' @param score.name type of score to merge
#' @param score.infos scoring info per cell per cell type
mergeScores <- function(score.name, score.infos, aggr.func=c) {
  names(score.infos[[1]][[score.name]]) %>% setNames(., .) %>% lapply(function(n)
    lapply(score.infos, `[[`, c(score.name, n)) %>% Reduce(aggr.func, .))
}

#' Merge Score Infos
#' @description merge score information from multiple datasets
#'
#' @inheritParams mergeScores
#' @inheritParams plapply
#' @inherit getMarkerScoreInfo return
#' @examples
#'   clf_datas <- lapply(cms, getClassificationData, marker_path)
#'   score_infos <- lapply(clf_datas, getMarkerScoreInfo)
#'   all_score_info <- mergeScoreInfos(score_infos, verbose=T)
#'
#' @export
mergeScoreInfos <- function(score.infos, verbose=F) {
  names(score.infos[[1]]) %>% setNames(., .) %>% plapply(mergeScores, score.infos, verbose=verbose)
}

mergeAnnotationInfos <- function(ann.infos) {
  return(list(
    annotation=mergeScores("annotation", ann.infos),
    scores=mergeScores("scores", ann.infos, aggr.func=rbind),
    annotation.filt=mergeScores("annotation.filt", ann.infos)
  ))
}

#' Get Cell Type Score Info
#' @description estimate info, neccessary for scoring of cells by cell types for the specified cell type
#'
#' @param markers element of marker list. List with fields `expressed` and `not_expressed`
#' @param tf.idf TF-IDF normalized matrix. Can be obtained with `normalizeTfIdfWithFeatures`
#' @param aggr if scores must be aggregated for the whole cell type or returned for each marker separately
#' @return list with score info:\itemize{
#'   \item{scores.raw: scores from positive markers}
#'   \item{mult: score multiplier, estimated with negative markers}
#'   \item{max.positive: maximal expression of positive markers. Used for estimation of negative scores}
#'   \item{scores: final scores. Equal to `scores * score.mult`}
#' }
getCellTypeScoreInfo <- function(markers, tf.idf, aggr=T) {
  expressed.genes <- intersect(markers$expressed, colnames(tf.idf))
  if (length(expressed.genes) == 0) {
    if (aggr) {
      scores <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf)) #row: cells # give each cell a score of 0
    } else {
      scores <- matrix(0, nrow=nrow(tf.idf), ncol=0, dimnames=list(rownames(tf.idf), character())) #Here it is an empty matrix, only has row (cell) names.
    }

    max.positive.expr <- setNames(rep(0, nrow(tf.idf)), rownames(tf.idf)) #now the same as the scores above, a list of cells, each cell has a score of 0
  } else { #marker gene is not empty
    c.submat <- tf.idf[, expressed.genes, drop=F]
    c.submat.t <- Matrix::t(c.submat) #here, row becomes genes and col becomes cells
    scores <- if (aggr) Matrix::colSums(c.submat.t) else c.submat
    max.positive.expr <- apply(c.submat.t, 2, max) #col max, for a given cell, which positive marker has max tf.idf?
  }

  not.expressed.genes <- intersect(markers$not_expressed, colnames(tf.idf))
  if (length(not.expressed.genes) == 0) {
    score.mult <- setNames(rep(1, nrow(tf.idf)), rownames(tf.idf)) #No need to consider aggr?
  } else {
    max.negative.expr <- apply(tf.idf[, not.expressed.genes, drop=F], 1, max) #For a given cell, which negative marker has the max tf.idf?
    # max.negative.expr <- sparseRowMax(tf.idf[, not.expressed.genes, drop=F])
    score.mult <- pmax(max.positive.expr - max.negative.expr, 0) / max.positive.expr
    score.mult[is.na(score.mult)] <- 0
  }

  res <- list(scores.raw=scores, mult=score.mult, max.positive=max.positive.expr, scores=(scores * score.mult))
  return(res)
}

normalizeScores <- function(scores, min.val=1e-10) {
  scores[rowSums(scores) < 1e-10,] <- 1
  scores  %<>% `/`(rowSums(.))
  return(scores)
}

#' Li: Return initial scores of each cell type for each cell
#' Li: In the old version, the St nromalization is on all the cell types at each layer, not on cell types under the same parent.
#' Li: Now, the St nromalization is on all the cell types under the same parent.
#' Li: @param score.info pre-calculated score info from `getMarkerScoreInfo`
#' Li: @param aggr should individual gene scores be aggregated per cell type? If `FALSE`,
#' Li: returns list of data.frames, showing scores of each gene for each cell.
#' Useful for debugging list of markers.
#' @inheritParams getMarkerScoreInfo
#' @return data.frame with rows corresponding to cells and columns corresponding to cell types.
#'   Values are cell type scores, normalized per level of hierarchy
getMarkerScoresPerCellType <- function(clf.data, score.info=NULL, aggr=T, underParent=T) {
  if (is.null(score.info)) {
    score.info <- getMarkerScoreInfo(clf.data, aggr=aggr)
  }

  scores <- lapply(score.info, `[[`, "scores")
  if (!aggr)
    return(lapply(scores, as.matrix) %>% lapply(as.data.frame, optional=T))

  clf.nodes <- classificationTreeToDf(clf.data$classification.tree)
  scores %<>% as.data.frame(optional=T)

  if (underParent) {
    parent <- split(clf.nodes$Parent, clf.nodes$PathLen)
    node <- split(clf.nodes$Node, clf.nodes$PathLen)
    family <- lapply(1:length(parent), function(n) split(node[[n]], parent[[n]])) %>% unlist(recursive=F)
    #TODO: Li added drop=F, may need check! Removed!
    message("Normalized St under parents ...")
    scores <- lapply(family, function(n) scores[, n]%<>% normalizeScores()) # Li: normlization: types under the same parent
    coln<- scores %>% sapply(colnames) %>% Reduce(c,.)
    scores %<>% as.data.frame
    colnames(scores) <- coln
  } else{
    for (nodes in split(clf.nodes$Node, clf.nodes$PathLen)) {
      #TODO: Li added drop=F, may need check! Removed
      message("Normalized St over the entire layer ...")
      scores[, nodes]  %<>% normalizeScores()}
  }

  return(scores)
}



#' Li: The algorithm used here looks like HMM full probabilistic model (P(S,π|HMM,θ); S: an observed sequence, π: a state path, HMM with parameters θ, The probability P(S,π|HMM,θ) that an HMM with parameters θ generates a state path π and an observed sequence S is the product of all the emission probabilities and transition probabilities that were used).
#' 1) each branch/cell type chain/lineage = sequence (made of symbol, here symbol=cell type), 2) state=layer, 3) the emission probability at each state = the normalized St score for each cell type, however at each state, the set of emission probabilities will vary from lineage to lineage, because we normalize scores under parent types, rather than for the entire layer, which means that the kid type probabilities under the same parent will sum up to one, and the cell types from all other parents of the same layer will equal to 0.  4) the transition probability from one layer to the next is 1, self transition probability is 0.
#' HMM (https://www.nature.com/articles/nbt1004-1315)
#' Li: This is a hurrican script
#' Extract Cell Type Probabilities
#' @description Re-normalizes scores to estimate full probability for each cell to belong to a specific cell type
#' @return data.frame of annotation probabilities with cell by rows and cell types by columns
#' @export
extractCellTypeProbs <- function(annotation.scores, clf.tree, ann.level=NULL) {
  extractLevelProbs <- function(scores, tree, cur.node) {
    sub.tree <- if (is.null(cur.node)) tree else tree[[cur.node]]

    if (is.null(names(sub.tree))) {
      return(scores[sub.tree] / pmax(rowSums(scores[sub.tree]), 1e-5))
    }

    if ((length(sub.tree) == 1)) {
      return(extractLevelProbs(scores, sub.tree, names(sub.tree)[[1]]))
    }

    if (is.null(sub.tree))
      stop("Error: is.null(tree[[cur.node]]), ", cur.node)

    return(names(sub.tree) %>% setNames(., .) %>% lapply(function(n) as_tibble(extractLevelProbs(scores, sub.tree, n)) * scores[[n]]) %>% Reduce(cbind, .))
  }

  probs <- extractLevelProbs(annotation.scores, classificationTreeToHierarhy(clf.tree, max.depth=ann.level), NULL) %>%
    as.data.frame() %>% set_rownames(rownames(annotation.scores))
  probs[probs < 1e-5] <- 0.0
  probs <- probs / pmax(rowSums(probs), 1e-5)
  probs[rowSums(probs) < 0.1,] <- 1 / ncol(probs)
  return(probs)
}
