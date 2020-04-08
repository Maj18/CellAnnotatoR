# Analysis of marmoset data

This file shows how to analyze the visual cortex (V1) data from marmoset  using CellAnnotatoR [CellAnnotatoR tutorial](https://github.com/khodosevichlab/CellAnnotatoR/blob/master/vignettes/mca_marker_selection.md) 

```{r}
library(ggplot2)
library(magrittr)
library(Matrix)
library(pbapply)
library(dplyr)
library(pagoda2)
library(CellAnnotatoR)

theme_set(theme_bw())
``

The visual cortex (V1) data from marmoset has been re-processed

## Import the visual cortex (V1) data from marmoset

```{r}
p3 <- readRDS("~/data/marmoset/marm022_v1_p2.rds")
names(p3)
```

## Run clustering and visualize.

```{r}
p3$getKnnClusters(type="PCA", method=conos::leiden.community, n.iterations=10, 
                  resolution=6, name="leiden")
conos::embeddingPlot(p3$embeddings$PCA$UMAP_graph, groups=p3$clusters$PCA$leiden)
```

## Define hierarchy

The automated algorithm is rather simple and will be improved in future. And currently, it can not be used on the aligned data.

```{r} 
hierarchy3 <- deriveHierarchy(p3$reductions$PCA, p3$clusters$PCA$leiden, max.depth=3)
clf_tree3 <- hierarchyToClassificationTree(hierarchy3)
plotTypeHierarchy(clf_tree3)
```

Annotation by level:

```{r}
ann_by_level3 <- mergeAnnotationByLevels(p3$clusters$PCA$leiden, clf_tree3)
plotAnnotationByLevels(p3$embeddings$PCA, ann_by_level3, font.size=c(2, 3), n.col=3, shuffle.colors=T)
```

## Marker selection

First, we need to find DE genes for each sub-brunch of the hierarchy:

```{r}
ann_by_parent3 <- getAnnotationPerParent(clf_tree3, p3$clusters$PCA$leiden)
as.character(ann_by_parent3)
```

Get DE genes for each "parent" cell type

```{r}
de_info_per_parent3 <- ann_by_parent3 %>% 
  pblapply(function(ann) p3$getDifferentialGenes(groups=ann, z.threshold=0))
de_info_per_parent3 <- pbmapply(function(de, ann) prepareDeInfo(de, ann, cm.raw=p3$misc$rawCounts, n.cores=10), 
                               de_info_per_parent3, ann_by_parent3)
```

Next, do TF-IDF normalization of the matrix:

```{r}
cm_norm <- p3$counts %>% normalizeTfIdfWithFeatures() %>% .[names(p3$clusters$PCA$leiden), ]
```
To be continued