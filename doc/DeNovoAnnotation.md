# De novo annotation of cell types:

You can do completely de-novo annotation. The goal here is to learn how to build a hierarchy and to pick markers. So you can do it by yourself without any automated code. The process is described in the "de-novo annotation" part of the documentation: https://github.com/khodosevichlab/CellAnnotatoR#de-novo-annotation .

It seems the empty markers for some cell type causes lots of problem????

```{r}
library(ggplot2)
library(magrittr)
library(Matrix)
library(pbapply)
library(dplyr)
library(pagoda2)
library(CellAnnotatoR)
library(parallel)

theme_set(theme_bw())
```

## 1. Build a hierarchy

KNN Clustering (unsupervised) to get cell clusters, the number of clusters may vary from run to run! I The clusters here are flat. the hierarchical levels are give by the classification tree.

```{r}
p3$getKnnClusters(type="PCA", method=conos::leiden.community, n.iterations=10,resolution=6, name="leiden")
```

## 2. Find marker candidates either with differential expression or using prior knowledge

```{r}
de_per_parent3 <- ann_by_parent3 %>% pblapply(function(ann) p3$getDifferentialGenes(groups=ann, z.threshold=0)) 
MarkerCandidates <- names(de_per_parent3) %>% setNames(., .) %>% lapply(function(n) lapply(de_per_parent3[[n]], preSelectMarkersForType))
```

```
Or, to filter marker candidates and add them to the marker.list, we should run the following command instead of the one above:

marker_info<- names(de_info_per_parent) %>% setNames(., .) %>% lapply(function(n)
  selectMarkersPerType(cm_norm, ann_by_parent[[n]], 
                        preSelectMarkerCandidates(de_info_per_parent[[n]], blacklist=marker_blacklist), 
                        marker.list=marker_list_prior,
                        parent=n,max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))

Here, the marker.list=marker_list_prior is where we bring in the prior marker list knowledge. However, for the moment, the prior marker list has predefined cell type labels, while for our KNN clustering above, we got random cell type labels, how could we connect them two together, will the command above automatically do it for us?

marker_list<- setNames(marker_info, NULL) %>% unlist(recursive=F)
```

```{t}
source("~/src/markerSelection/preSelectMarkersForType.R") #This version took cells types withought DE identified into consideration
```

de_per_parent3, store DE gene between subtypes within each parent cell type + DE genes between parent cell types under root.

saveRDS(MarkerCandidates, file = "~/data/marmoset/denovoannotation/MarkerCandidates.rds")

MarkerCandidates <- readRDS(file="~/data/marmoset/denovoannotation/MarkerCandidates.rds")

## 3. Plot the markers on your data and ensure that they suit your case

Below I will use cell type l1_11 (see "~/data/marmoset/hierarchy_auto/marmosetv1_hierarchy.pdf") as an example

Seurat has its own functions for plotting gene expression, but for general case CellAnnotatoR provides the functionn plotGeneExpression(genes, embedding, cm, ...). It returns list of plots for individual genes. Note: matrix cm must be transposed, i.e. have genes as columns and cells as rows.

```{r}
p3 <- readRDS(file="~/data/marmoset/p3.rds")
markerlistl1_11[1:12] %>% plotGeneExpression(p3$embeddings$PCA$UMAP, p3$counts)
conos::embeddingPlot(p3$embeddings$PCA$UMAP_graph, groups=p3$clusters$PCA$leiden)
```

Put gene expression map and cell type map together for comparison:

```{r}
pdf("~/result/marmoset/Denovo/markerGexpressionMap/l1_11$26Gene1-12.pdf")
cowplot::plot_grid(markerlistl1_11[1:12] %>% plotGeneExpression(p3$embeddings$PCA$UMAP, p3$counts), conos::embeddingPlot(p3$embeddings$PCA$UMAP_graph, groups=p3$clusters$PCA$leiden))
dev.off()
```

If you want to use panel of violinplots instead, you can use plotExpressionViolinMap(genes, cm, annotation). It suits better for large panels of markers:

```{r}
markerlistl1_11 %>% plotExpressionViolinMap(p3$counts, p3$clusters$PCA$leiden)
c("LOC100385341") %>% plotExpressionViolinMap(p3$counts, p3$clusters$PCA$leiden) #cell types 49, 51
```

## 4. Add the markers to the file and re-run the classification

Adding markers is trivial and described in the Garnett specification

Export the marker data to a plain txt file

```{r}
sink("~/data/marmoset/denovoannotation/marker.file.txt")
for (i in names(MarkerCandidates)){
  for (j in names(MarkerCandidates[[i]])) {
    for (k in names(MarkerCandidates[[i]][[j]])){
      if (k == "positive"){
        kk = "expressed"
      } else if (k == "negative") {
        kk = "not expressed"
      }
      cat(paste0(kk,": "))
      if (!identical(MarkerCandidates[[i]][[j]][[k]], character(0))){
        cat(paste0(strsplit(MarkerCandidates[[i]][[j]][[k]], split="\t"), collapse=", "), "\n")
      } else {
        cat("NULL", "\n")
      }
    }
  }
}
sink()
```

### First, getClassificationData performs TF-IDF normalization inside, which doesn't depend on the marker list. So we can do it only once:

```{r}
cm_norm <- p3$misc$rawCounts %>% normalizeTfIdfWithFeatures() %>% t()
```

saveRDS(cm_norm, file="~/data/marmoset/denovoannotation/cm_norm.rds")

cm_norm <- readRDS("~/data/marmoset/denovoannotation/cm_norm.rds")

```{r}
markers <- parseMarkerFile("~/data/marmoset/denovoannotation/marker.file.txt")
```

clf_data <- getClassificationData(cm_norm, markers, prenormalized=TRUE) #the following 4 commands do the same thing as this command line

```{r}
gi <- CellAnnotatoR:::unifyGeneIds(cm_norm, data.gene.id.type="SYMBOL", marker.gene.id.type="SYMBOL", db=NULL, verbose=F)
classification.tree <- createClassificationTree(markers) #This tree is the same as what we got from "hierarchy3 <- deriveHierarchy(p3$reductions$PCA, p3$clusters$PCA$leiden, max.depth=3)", "clf_tree3 <- hierarchyToClassificationTree(hierarchy3)".
res <- list(cm=gi$cm, classification.tree=classification.tree, gene.table=gi$gene.table, marker.list=markers)
clf_data <- res
```

saveRDS(clf_data,file="~/data/marmoset/denovoannotation/clf_data.rds")

clf_data <- readRDS("~/data/marmoset/denovoannotation/clf_data.rds")

saveRDS(classification.tree,file="~/data/marmoset/denovoannotation/classification.tree.rds")

saveRDS(clf_tree3,file="~/data/marmoset/denovoannotation/clf_tree3.rds")

### Second, the most time-consuming step for classification is label propagation on graph. It improves classification quality in many cases, but for getting approximate results we can avoid this. To do so, it's enough to pass NULL instead of the graph object:

```{r}
clusters <- p3$clusters$PCA$leiden
ann_by_level <- assignCellsByScores(NULL, clf_data, clusters=clusters)
```

names(ann_by_level)
[1] "annotation"      "scores"          "annotation.filt"

saveRDS(ann_by_level, file = "~/data/marmoset/ann_by_level.rds")
ann_by_level <- readRDS("~/data/marmoset/ann_by_level.rds")

```
(So, during marker selection it's recommended to re-run the following codes every time you update markers:
clf_data <- getClassificationData(cm_norm, marker_path, prenormalized=T)
ann_by_level <- assignCellsByScores(NULL, clf_data, clusters=clusters)
ann_by_level <- assignCellsByScores(NULL, clf_data, clusters=p3$clusters$PCA$leiden)
)
```

## 5. Check the result

Plot annotation

```{r}
png("~/result/marmoset/Denovo/markerGexpressionMap/OldNewClusters.png")
cowplot::plot_grid(conos::embeddingPlot(p3$embeddings$PCA$UMAP_graph, groups=p3$clusters$PCA$leiden), conos::embeddingPlot(p3$embeddings$PCA$UMAP_graph, groups=as.factor(ann_by_level$annotation$l3)))
dev.off()
```

Or one can plot by levels:

```{r}
png("~/result/marmoset/Denovo/markerGexpressionMap/origvsnew_bylevel.png")
ann_by_level3 <- mergeAnnotationByLevels(p3$clusters$PCA$leiden, clf_tree3)
cowplot::plot_grid(plotAnnotationByLevels(p3$embeddings$PCA, ann_by_level3, n.col=3, font.size=c(2, 3)), plotAnnotationByLevels(p3$embeddings$PCA, ann_by_level$annotation, n.col=3, font.size=c(2, 3)))
dev.off()
```

Uncertainty plots per cell

```{r}
score_info <- CellAnnotatoR:::getMarkerScoreInfo(clf_data)
unc_per_cell <- CellAnnotatoR:::scoreCellUncertaintyPerLevel(ann_by_level, score_info)
lapply(unc_per_cell, function(unc) CellAnnotatoR:::plotUncertaintyPerCell(p3$embeddings$PCA, unc)) %>% 
  cowplot::plot_grid(plotlist=., ncol=1, labels=names(unc_per_cell), label_y=0.9, label_x=0.01)
```

Level 1 assignment scores shown on cell embedding:

```{r}
plotAssignmentScores(p3$embeddings$PCA, ann_by_level$scores$l1, clf_data$classification.tree, parent="root", build.panel=T, size=0.1, n.col=5)
```

Positive markers violin plot

```{r}
marker_list[sort(unique(ann_by_level$l1))] %>% lapply(`[[`, "expressed") %>% unlist() %>% unique() %>% plotExpressionViolinMap(p3$counts, ann_by_level$l1, gene.order=T)
```

Save to file:

```{r}
markerListToMarkup(marker_list, file="mca_markers_auto.txt") # With this funciton, I don't need the self-written script above
```
