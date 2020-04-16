getMeanConfidencePerType <- function(marker.list, cm.norm, annotation) {
  c.scores <- lapply(marker.list, CellAnnotatoR:::getCellTypeScoreInfo, cm.norm) %>%
    lapply(`[[`, "scores") %>% as.data.frame(optional=T) %>% normalizeScores()
  confidence <- CellAnnotatoR:::getAnnotationConfidence(annotation, c.scores)
  return(confidence %>% split(annotation[names(.)]) %>% sapply(mean))
}
#getMeanConfidencePerType(marker.list, t(cm.norm), ann_by_parent3[["l1_1"]])
#print(confidence %>% split(ann_by_parent3[["l1_1"]][names(.)]) %>% sapply(mean))
#confidence <- CellAnnotatoR:::getAnnotationConfidence(ann_by_parent3[["l1_1"]], c.scores)
#Here, getCellTypeScoreInfo is changed to CellAnnotatoR:::getCellTypeScoreInfo
#here, getAnnotationConfidence is changed to CellAnnotatoR:::getAnnotationConfidence
# c.scores <- lapply(marker.list, CellAnnotatoR:::getCellTypeScoreInfo, t(cm.norm)) %>%
#   lapply(`[[`, "scores") %>% as.data.frame(optional=T) %>% normalizeScores()
# getAnnotationConfidence <- function(annotation, scores) {
#   mapply(function(i,j) scores[i, j], 1:nrow(scores), match(annotation, colnames(scores))) %>%
#     setNames(rownames(scores)) %>% unlist()
# }
## match(annotation, colnames(scores)): pick up the annotation inf for particular cell types
# here, mapply(function(i,j) scores[i, j], 1:nrow(scores), match(annotation, colnames(scores))
# does the work as mapply(function(i,j) scores[i, j], i cell list, j cell type list)

