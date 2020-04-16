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
```

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

### Automated 

The automated algorithm is rather simple and will be improved in future. And currently, it can not be used on the aligned data.

```{r} 
hierarchy3 <- deriveHierarchy(p3$reductions$PCA, p3$clusters$PCA$leiden, max.depth=3)
clf_tree3 <- hierarchyToClassificationTree(hierarchy3)
plotTypeHierarchy(clf_tree3)
```
### Manual

1. For macaque_v1

```{r}
hierarchy3mv1 <- list('Neurons'=list('Inhibitory'=list('Id2'=c('Id2_Pde7b', 'Id2_Trpc3', 'Id2_Adamtsl1'), 'Pvalb'=list('Pvalb_Slit2'=c('Pvalb_Slit2_Sulf1', 'Pvalb_Slit2_Epha6'), 'Pvalb_Fstl5 '), 'Sst'=c(' Sst_Calb1', ' Sst_Slit2 ', ' Sst_Myo16'), 'Vip'=list('Vip_Nrg1' , 'Vip_Grm1'=c('Vim_Grm1_Sema3c', 'Vim_Grm1_Kiaa1217'))), 'Excitatory'=list('L2_4_Cux2'=list('L2_Cux2_Adcyap1' , 'L2_Cux2_Fbn2', 'L3_4_Rorb'=c('L3_Rorb_Ca10', 'L4_Rorb_Kcnh8')), 'L5_6_B3galt2'=list('L6_B3galt2_Pcp4'=c('L6_B3galt2_Pcp4_Htr2c', 'L6_B3galt2_Pcp4_Pcsk5', 'L6_B3galt2_Pcp4_Sema3e', 'L6_B3galt2_Pcp4_Tafa1', 'L6_B3galt2_Pcp4_Gria4'), 'L5_B3galt2_Atp2b4'=c('L5_B3galt2_Atp2b4_Dach1', 'L5_B3galt2_Atp2b4_Ptger3', 'L5_B3galt2_Atp2b4_Themis', 'L5_B3galt2_Atp2b4_Znf804b')))))
clf_tree3mv1 <- hierarchyToClassificationTree(hierarchy3mv1)
plotTypeHierarchy(clf_tree3mv1)
``

2. For macaque_v4

```{r}
hierarchy3mv4<- list('Neurons'=list('Inhibitory'=list('Id2'=list('Id2_Lamp5'=c('Id2_Lamp5_Nos1' , 'Id2_Lamp5_Pbx3' , 'Id2_Lamp5_Nmbr'), 'Id2_Nckap5' , 'Id2_Pax6'),'Pvalb'=list('Pvalb_Slit2'=c('Pvalb_Slit2_Sulf1' , 'Pvalb_Slit2_Col5a2' , 'Pvalb_Slit2_L3mbtl4' , 'Pvalb_Slit2_Glra3'), 'Pvalb_Fstl5'), 'Sst'=c("Sst_Calb1", "Sst_Slit2", "Sst_Pcsk5", "Sst_Reln", "Sst_Alcam"), 'Vip'=list('Vip_Nrg1'=c('Vip_Nrg1_Ttn' , 'Vip_Nrg1_Sema5a'), 'Vip_Grm1' , 'Vip_Arhgap18' , 'Vip_Sema3e' )), 'Excitatory'=list('L2_4_Cux2'=list('L2_Cux2_Adcyap1' , 'L2_Cux2_Fbn2' , 'L3_4_Rorb'=c('L3_Rorb_Cdh9' , 'L3_4_Rorb_Syt2' , 'L4_Rorb_Shisa9')), 'L5_6_B3galt2'=list('L6_B3galt2_Nr4a2', 'L5_6_B3galt2_Tle4'=c('L5_6_B3galt2_Tle4_Htr2c', 'L5_6_B3galt2_Tle4_Pcsk5', 'L5_6_B3galt2_Tle4_Adamts19'), 'L5_6_B3galt2_Lrrk1'=c('L5_6_B3galt2_Lrrk1_Lama4' , 'L5_6_B3galt2_Lrrk1_Spon1' , 'L5_6_B3galt2_Lrrk1_Prkg1' , 'L5_6_B3galt2_Lrrk1_Grin3a' , 'L5_6_B3galt2_Lrrk1_Tll1'), 'L5_6_B3galt2_Themis' , 'L5_B3galt2_Pde3a' ))), 'Astrocytes', 'Microglia', 'Oligodendrocytes', 'Oligodendrocyte_Precursors', 'Vascular') 
clf_tree3m <- hierarchyToClassificationTree(hierarchy3mv4)
plotTypeHierarchy(clf_tree3mv4)
```

3. For human_v2

```{r}
hierarchy3hv2 <- list('neuron'=list('Inhibitory'=list('Id2'=list('Id2_Lamp5'=c('Id2_Lamp5_Nos1', 'Id2_Lamp5_Crh', 'Id2_Lamp5_Nmbr'), 'Id2_Nckap5', 'Id2_Pax6'), 'Pvalb'=c('Pvalb_Nos1', 'Pvalb_Sulf1', 'Pvalb_Lgr5', 'Pvalb_Crh'), 'Sst'=c('Sst_Th', 'Sst_Tac3', 'Sst_Tac1', 'Sst_Stk32a', 'Sst_Nos1', 'Sst_Calb1', 'Sst_Npy', 'Sst_Isoc1'), 'Vip'=c('Vip_Abi3bp', 'Vip_Tyr', 'Vip_Crh', 'Vip_Sstr1', 'Vip_Sema3c', 'Vip_Nrg1', 'Vip_Cbln1')), 'Excitatory'=list('L2_3_Cux2'=list('L2_Cux2_Lamp5'=c('L3_Cux2_Prss12', 'L3_Cux2_Prss12'), 'L2_3_Cux2_Frem3', 'L3_Cux2_Prss12'), 'L4_Rorb'=c('L4_Rorb_Mme', 'L4_Rorb_Met', 'L4_Rorb_Arhgap15'), 'L5_6_Themis'=c('L5_6_Themis_Sema3a', 'L5_6_Themis_Ntng2'), 'L5_6_Fezf2'=c(' L5_6_Fezf2_Lrrk1', ' L5_6_Fezf2_Tle4', ' L5_6_Fezf2_Adra1a'))))
clf_tree3hv2 <- hierarchyToClassificationTree(hierarchy3hv2)
plotTypeHierarchy(clf_tree3h2)
```

## Marker selection

Annotation by level:

```{r}
ann_by_level3 <- mergeAnnotationByLevels(p3$clusters$PCA$leiden, clf_tree3)
plotAnnotationByLevels(p3$embeddings$PCA, ann_by_level3, font.size=c(2, 3), n.col=3, shuffle.colors=T)
```

First, we need to find DE genes for each sub-brunch of the hierarchy:

```{r}
ann_by_parent3 <- getAnnotationPerParent(clf_tree3, p3$clusters$PCA$leiden)
```

Get DE genes for each "parent" cell type

```{r}
de_per_parent3 <- ann_by_parent3 %>% 
  pblapply(function(ann) p3$getDifferentialGenes(groups=ann, z.threshold=0))
names(de_per_parent3)

de_info_per_parent3 <- pbmapply(function(de, ann) prepareDeInfo(de, ann, cm.raw=p3$misc$rawCounts, n.cores=10), 
                               de_per_parent3, ann_by_parent3)
```

Next, do TF-IDF normalization of the matrix:

```{r}
cm_norm <- p3$misc$rawCounts %>% normalizeTfIdfWithFeatures() 
cm_norm_annotation <- cm_norm[,names(p3$clusters$PCA$leiden)] %>% t(.)
```

Finally we run makrker selection for each of the sub-branches:

```{r}
marker_info <- names(de_per_parent3) %>% setNames(., .) %>% lapply(function(n)
  selectMarkersPerType(cm_norm_annotation, ann_by_parent3[[n]], 
                       preSelectMarkerCandidates(de_per_parent3[[n]], blacklist=NULL), 
                       parent=n, max.iters=50, max.uncertainty=0.25, log.step=5, verbose=1, n.cores=10))
           
```

To be continued