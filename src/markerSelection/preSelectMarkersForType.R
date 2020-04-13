#Here I have modified the function "preSelectMarkersForType", to take the "NULL" (no DE genes) conditions into consideration:
preSelectMarkersForType <- function(de.df, whitelist=NULL, blacklist=NULL, min.pos.markers=5, max.pos.markers=100,
                                    min.pos.specificity=0.2, min.pos.expression.frac=0.1,
                                    min.pos.markers.soft=as.integer(round(mean(c(min.pos.markers, max.pos.markers)))),
                                    min.pos.specificity.soft=0.75, min.pos.expression.frac.soft=0.25,
                                    pos.expression.frac.weight=0.2, max.neg.expression.frac=0.1) {
  if (!is.null(whitelist)) {
    de.df %<>% dplyr::filter(Gene %in% whitelist)
  }
  if (!is.null(blacklist)) {
    de.df %<>% dplyr::filter(!(Gene %in% blacklist))
  }
  if (length(de.df) == 0) { #To handle those cells without DE genes
    neg.markers <- NA
    pos.markers <- NA} else {
      de.pos <- de.df[de.df$Z > 0, ]
      if (sum(de.pos$Specificity > min.pos.specificity) < min.pos.markers) {
        pos.markers <- de.pos$Gene[order(de.pos$Specificity, decreasing=T) %>% .[1:min(min.pos.markers, length(.))]]
      } else {
        de.pos %<>% .[.$Specificity > min.pos.specificity,]
        if (sum(de.pos$ExpressionFraction > min.pos.expression.frac) < min.pos.markers) {
          pos.markers <- de.pos$Gene[order(de.pos$ExpressionFraction, decreasing=T) %>% .[1:min(min.pos.markers, length(.))]]
        } else {
          de.pos %<>% .[.$ExpressionFraction > min.pos.expression.frac,] %>%
            .[order(.$Specificity + .$ExpressionFraction * pos.expression.frac.weight, decreasing=T),]
          soft.mask <- (de.pos$ExpressionFraction > min.pos.expression.frac.soft) & (de.pos$Specificity > min.pos.specificity.soft)
          if (sum(soft.mask) > min.pos.markers.soft) {
            de.pos %<>% .[soft.mask, ]
          } else {
            de.pos %<>% .[1:min(nrow(.), min.pos.markers.soft), ]
          }
          pos.markers <- de.pos %>% .$Gene %>% .[1:min(length(.), max.pos.markers)]
        }
      }
      # Negative: ExpressionFraction < 0.1 && Z < 0 && top by specificity (or > 0.95)
      neg.markers <- de.df %>% .[(.$Z < 0) & (.$ExpressionFraction < max.neg.expression.frac), ] %>% .$Gene
    }
  return(list(positive=pos.markers[!is.na(pos.markers)], negative=neg.markers[!is.na(neg.markers)]))
}    
