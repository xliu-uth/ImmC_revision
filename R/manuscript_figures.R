library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggalluvial)

data_folder <- 'C:\\Users/huang/Dropbox/AXL/Projects/ImmClassifier/Manuscript/data_for_figure_plots/'

#########################################################################
#                               Figure 2                                #
#########################################################################


######################   Figure 2A   ######################

bulk_immc <- fread(paste(data_folder, 'immc/bulk.output.txt', sep = "/"))
bulk_immc[, Known:=gsub("_[0-9]+$", "", Cell)]
bulk_immc[grepl("ERY1|GRAN1|PRE_BCELL2|GMP|CMP|MEP|MEGA1|HSC", Cell), Known:="CD34+"]


from_celltype <- c('DENDA1', 'DENDA2', 
                   paste0('TCELLA', c(1:4, 6:8)),
                   paste0('ERY', c(2:5)),
                   paste0('GRAN',2:3),
                   paste0('MONO',1:2),
                   paste0('BCELLA',1:4),
                   paste0('NKA',1:3))

to_celltype <- c('pDC', 'mDC', rep('CD8', 4), rep('CD4', 3),
                 rep('ERY', 4), 
                 rep('GRAN',2),
                 rep('MONO',2),
                 rep('BCELL', 4),
                 rep('NK', 3))

for (i in 1:length(from_celltype)){
  bulk_immc[Known==from_celltype[i], Known:=to_celltype[i]]  
}


bulk_immc[, annot:=ImmClassifier_prediction]
bulk_immc[grepl("CD4", ImmClassifier_prediction), annot:="CD4"]
bulk_immc[grepl("CD8", ImmClassifier_prediction), annot:="CD8"]
bulk_immc[grepl("CD34", ImmClassifier_prediction), annot:="CD34+"]
bulk_immc[grepl("Mega|Platelet", ImmClassifier_prediction), annot:="Mega/Platelet"]

ord1 <- c('L:B', 'CD4', 'CD8', 'L:NK', 'M:mDC', 'M:pDC', 
          'M:Mono', 'M:Neu','M:Ery', 'CD34+', 'Mega/Platelet', 
          'L:unconvT:MAIT','L:PC','M:Eos')
bulk_immc[!annot %in% ord1, annot:="Other"]


plot_heatmap(table(bulk_immc[, c(c('Known', 'annot'))]), title = "", palette= "RdPu", 
             minshow = 1, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1 = ord1,
             ord2 = c('BCELL', 'CD4', 'CD8', 'NK', 'mDC', 'pDC', 
                      'MONO', 'GRAN', 'ERY', 'CD34+','MEGA2', 'NKA4', 
                      'PRE_BCELL3',  'EOS2', 'BASO1'),
             measures = c('Recall', 'Precision', 'f1') )  




######################   Figure 2B   ######################



skcm_meta <- data.table(readRDS(paste(data_folder,"metadata/SKCM-NH-meta.rds", sep = "/")))
setkey(skcm_meta, Cell)
skcm_cluster_names <- data.table(readRDS(paste(data_folder,"metadata/SKCM-NH-clusternames.rds", sep = "/")))
setkey(skcm_cluster_names, clusterID )
skcm_fine_meta <- fread(paste(data_folder,"metadata/SKCM-NH-fineT-meta.txt", sep = "/"))
setkey(skcm_fine_meta, `Cell Name`)
skcm_fine_cluster_names <- c('Exhaustion/cell-cycle',
                             'Exhaustion/heat shock protein',
                             'Exhaustion',
                             'Memory/effector',
                             'Early activated cells',
                             'Memory/effector')
names(skcm_fine_cluster_names) <- c('CD8_1', 'CD8_2', 'CD8_3',
                                    'CD8_4', 'CD8_5', 'CD8_6')

skcm_immc <- fread(paste(data_folder, 'immc/skcm.output.txt', sep = "/"))
skcm_immc[, Known:=skcm_cluster_names[skcm_meta[Cell, 'Cluster'], clusterName]]
skcm_immc[, annot:=ImmClassifier_prediction]
skcm_immc <- skcm_immc[!is.na(Known),]

skcm_immc[grepl("ympho|Memory|Exh", Known), Known:="T cells"]
skcm_immc[grepl("Treg", annot), annot:="Treg"]
skcm_immc[grepl("DC", annot), annot:="DC"]
skcm_immc[grepl("L:T", annot), annot:="T cells"]
skcm_immc[grepl("Mono|Mac", annot), annot:="Mono|Mac"]

ord1 <- c('L:B', 'L:PC',  'T cells','Treg',  'Mono|Mac', 'DC')

skcm_immc[!annot %in% ord1, annot:='Other']

plot_heatmap(table(skcm_immc[, c(c('Known', 'annot'))]), title = "", 
             minshow = 19, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=ord1,
             ord2 = c('B', 'Plasma Cells', 'T cells',
                      'Treg', 'Monocytes/Macrophage', 'DC'),
             measures = c('Recall', 'Precision', 'f1') )  




skcm_immc <- fread(paste(data_folder, 'immc/skcm.output.txt', sep = "/"))
skcm_immc[, Known:=skcm_fine_cluster_names[skcm_fine_meta[Cell, Cluster]]]
skcm_immc[, annot:=ImmClassifier_prediction]
skcm_immc <- skcm_immc[!is.na(Known),]
#skcm_immc[grepl("xha", Known), Known:="exhausted T cells"]
#skcm_immc[grepl("L:T:CD4", annot), annot:="CD4 T cells"]
skcm_immc[grepl("L:T:CD8:CM$|L:T:CD8:EM$|L:T:CD8:TRM$", annot), annot:="CD8 memory T"]
skcm_immc[grepl("L:T:CD4:CM$|L:T:CD4:EM$|L:T:CD4:TRM$", annot), annot:="CD4 memory T"]

ord1 = c('L:T:CD8:Ex','L:T:CD8:EMRA','CD8 memory T','CD4 memory T')
skcm_immc[!annot %in% ord1, annot:='Other']


plot_heatmap(table(skcm_immc[, c(c('Known', 'annot'))]), title = "", palette='Blues',
             minshow = 15, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=ord1,#skcm_immc[, .N, by = annot][, annot],
             ord2 =  c('Exhaustion/cell-cycle',
                       'Exhaustion/heat shock protein','Exhaustion', 'Early activated cells','Memory/effector'),
             measures = c('Recall', 'Precision', 'f1') )  






######################   Figure 2C   ######################

hnscc_meta <- readRDS(paste(data_folder,'metadata/HNSCC.meta.rds', sep = "/"))
hnscc_immc <- fread(paste(data_folder,'immc/hnscc.output.txt', sep = "/"))

hnscc_immc[, Known:=hnscc_meta[Cell,'non_cancer_cell_type']]
hnscc_immc[, annot:=ImmClassifier_prediction]

hnscc_immc <- hnscc_immc[grepl("B|Den|Ma|T", Known), ]
hnscc_immc[grepl("L:T|L:unconv", annot), annot:="T cells"]
hnscc_immc[grepl("DC", annot), annot:="DC"]
ord1 <- c('L:PC', 'M:Mac', 'M:Mono', 'T cells', 'DC', 'M:Mast')
hnscc_immc[!annot %in% ord1, annot:='Other']


plot_heatmap(table(hnscc_immc[, c(c('Known', 'annot'))]), title = "", 
             minshow = 1, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=ord1,
             ord2 =  hnscc_immc[, .N, by = Known][, Known],
             measures = c('Recall', 'Precision', 'f1') )               



# plot legendcol

colfunc <- colorRampPalette(brewer.pal(9, 'RdPu'))
colfunc2 <- colorRampPalette(brewer.pal(9, 'Blues'))
barplot(rep(1, 10),col = c('white', colfunc(9)[c(3,5,7,9)],
                           'white', colfunc2(9)[c(3,5,7,9)]))



#########################################################################
#                               Figure 3                                #
#########################################################################



# Group the celltypes by levels

levels_celltype <- list(
  'Level 1'=c('CD34', 'L', 'M'),
  'Level 2'=c('L:B', 'L:T', 'L:NK',
              'DC', 'M:Mono',
              'M:Mac', 'M:Mast', 'M:Neu'),
  'Level 3'=c('CD4+T', 'CD8+T', 'mDC', 'pDC'),
  'Level 4'=c('CD4+Naive', 'CD4+CM', 'CD4+EM', 'CD4+Ex', 'CD4+Treg', 'CD4+Tfh',
              'CD8+Naive', 'CD8+CM', 'CD8+EM', 'CD8+Ex', 'MAIT')
)

############## brca3p dataset  ##############

brca3p_meta <- data.table(readRDS(paste(data_folder, 'metadata/BRCA-DP-meta.rds', sep = "/")))
brca3p_clsnames <- data.table(readRDS(paste(data_folder, 'metadata/BRCA-DP-clusternames.rds', sep = "/")))
setkey(brca3p_clsnames, ClusterID)
setkey(brca3p_meta, cellid)


brca3p_immc <- fread(paste(data_folder,'immc/brca3p.output.txt', sep = "/"))


setkey(brca3p_immc, Cell)
brca3p_immc[, Known:=brca3p_clsnames[brca3p_meta[Cell, cluster], Final_annotation]]


load(paste(data_folder, 'singleR/singler_brca3p.RData', sep = "/"))
brca3p.singler <- singler

brca3p.singler.df <- data.frame(singler=brca3p.singler$singler[[1]]$SingleR.single$labels, stringsAsFactors = F)
rownames(brca3p.singler.df ) <- brca3p.singler$singler[[1]]$SingleR.single$cell.names
brca3p_singler_dt <- data.table(cell=singler$singler[[1]]$SingleR.single$cell.names,
                                singler=singler$singler[[1]]$SingleR.single$labels)
brca3p_garnett_extend  <- readRDS(paste(data_folder,'garnett/brca.3p.garnett.extend.rds', sep = "/"))
brca3p_garnett_extend_dt <- data.table(cell=brca3p_garnett_extend$cellid, 
                                       garnett=brca3p_garnett_extend$cluster_ext_type)

rm(brca3p.singler)
rm(singler)

brca3p_umap <- fread(paste(data_folder,'umap/out-brca3p_umap/models//umap.pca_0500.cluster_25.k_25.res_1.0.umap_25/coordinates.txt', sep = "/"))

# Combine all prediction results and the original annotation into one data table
brca3p_dt <- data.table(
    merge(
      merge(
        merge(brca3p_umap[, .(Cell, UMAP1, UMAP2, Cluster)], 
              brca3p_immc, by.x='Cell', by.y='Cell', all.x = T),
        brca3p_singler_dt, by.x='Cell', by.y='cell', all.x = T),
      brca3p_garnett_extend_dt, by.x='Cell', by.y='cell', all.x = T)
    
)[, original_annotation:=brca3p_clsnames[brca3p_meta[Cell, cluster], Final_annotation]]

brca3p_dt[, Cluster:=as.character(Cluster)]
colnames(brca3p_dt)[grepl("ImmC", colnames(brca3p_dt))] <- "ImmC"
colnames(brca3p_dt)[grepl("singler", colnames(brca3p_dt))] <- "SingleR"
setkey(brca3p_dt, Cell)



# --------------------------- brca3p level 1 ------------------------------------ #
brca3p_dt[, original_annotation_L1:="not_specified"]
brca3p_dt[grepl("T:|NK:|B:|NKT", original_annotation), original_annotation_L1:='L']
brca3p_dt[grepl("DC|MONO|MA|NEU", original_annotation), original_annotation_L1:='M']

brca3p_dt[, ImmC_L1:="not_specified"]
brca3p_dt[grepl("L:", ImmC), ImmC_L1:='L']
brca3p_dt[grepl("M:", ImmC), ImmC_L1:='M']
brca3p_dt[grepl("CD34", ImmC), ImmC_L1:='CD34']


plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L1, annot=ImmC_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = ImmC_L1][, ImmC_L1],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 



brca3p_dt[, SR_L1:="not_specified"]
brca3p_dt[grepl("T_cell|NK|B_cell:|Pre-B|^B_cell$", SingleR), SR_L1:='L']
brca3p_dt[grepl("Mono|Neu|Mac|DC|Erythroblast", SingleR), SR_L1:='M']
brca3p_dt[grepl("GMP|CD34\\+|CMP", SingleR), SR_L1:='CD34']
brca3p_dt[grepl("Epi|iPS|Embryonic_stem_cells|Fib|Endo|Tissue_stem_cells|Chon|Smooth", SingleR), SR_L1:='Non-immune']

brca3p_dt[, .N, by = c('SingleR', 'SR_L1')]



plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L1, annot=SR_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = SR_L1][, SR_L1],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 


brca3p_dt[, GN_L1:=garnett]
brca3p_dt[grepl("T cell|NK|B", garnett), GN_L1:='L']
brca3p_dt[grepl("Mono|Dend", garnett), GN_L1:='M']
brca3p_dt[grepl("CD34", garnett), GN_L1:='CD34']


plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L1, annot=GN_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = GN_L1][, GN_L1],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 



# --------------------------- brca3p level 2 ------------------------------------ #

brca3p_dt[, original_annotation_L2:=gsub(":[^:]*$", "",original_annotation)]
brca3p_dt[grepl("DC", original_annotation), original_annotation_L2:="DC"]
brca3p_dt[, ImmC_L2:='not_specified']
brca3p_dt[grepl("L:T|unconvT", ImmC), ImmC_L2:='L:T']
brca3p_dt[grepl("DC", ImmC), ImmC_L2:='DC']
brca3p_dt[ImmC %in% c('L:B', 'L:NK','M:Mono', 'M:Mac', 'M:Mast', 'M:Neu'), ImmC_L2:=ImmC]

plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L2, annot=ImmC_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = ImmC_L2][, ImmC_L2],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 




brca3p_dt[, SR_L2:="not_specified"]
brca3p_dt[grepl("B_cell:|Pre-B|^B_cell$", SingleR), SR_L2:='B']
brca3p_dt[grepl("T_cell", SingleR), SR_L2:='T']
brca3p_dt[grepl("NK", SingleR), SR_L2:='NK']
brca3p_dt[grepl("Mono", SingleR), SR_L2:='Mono']
brca3p_dt[grepl("Mac", SingleR), SR_L2:='Mac']
brca3p_dt[grepl("Neu", SingleR), SR_L2:='Neu']
brca3p_dt[grepl("DC", SingleR), SR_L2:='mDC']
brca3p_dt[grepl("Ery", SingleR), SR_L2:='Ery']

brca3p_dt[grepl("GMP|CD34\\+|CMP|HSC", SingleR), SR_L2:='CD34']
brca3p_dt[grepl("Epi|iPS|Embryonic_stem_cells|Fib|Endo|Tissue_stem_cells|Chon|Smooth", SingleR), SR_L2:='Non-immune']

plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L2, annot=SR_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = SR_L2][, SR_L2],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 





brca3p_dt[, GN_L2:=garnett]
brca3p_dt[grepl("T cell", garnett), GN_L2:='T']



plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L2, annot=GN_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = GN_L2][, GN_L2],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 





# --------------------------- brca3p level 3 ------------------------------------ #

brca3p_dt[, original_annotation_L3:='not_specified']
brca3p_dt[grepl("CD4", original_annotation), original_annotation_L3:='CD4+T']
brca3p_dt[grepl("CD8", original_annotation), original_annotation_L3:='CD8+T']
brca3p_dt[grepl("mDC", original_annotation), original_annotation_L3:='mDC']
brca3p_dt[grepl("pDC", original_annotation), original_annotation_L3:='pDC']


brca3p_dt[, ImmC_L3:='not_specified']
brca3p_dt[grepl("L:T:CD4", ImmC), ImmC_L3:='CD4+T']
brca3p_dt[grepl("L:T:CD8", ImmC), ImmC_L3:='CD8+T']
brca3p_dt[grepl("mDC", ImmC), ImmC_L3:='mDC']
brca3p_dt[grepl("pDC", ImmC), ImmC_L3:='pDC']
plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L3, annot=ImmC_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = ImmC_L3][, ImmC_L3],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 



brca3p_dt[, SR_L3:="not_specified"]

brca3p_dt[grepl("CD4", SingleR), SR_L3:='CD4+T']
brca3p_dt[grepl("CD8", SingleR), SR_L3:='CD8+T']
brca3p_dt[grepl("DC", SingleR), SR_L3:='mDC']

plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L3, annot=SR_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = SR_L3][, SR_L3],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 




brca3p_dt[, GN_L3:=garnett]
brca3p_dt[grepl("CD4", garnett), GN_L3:='CD4+T']
brca3p_dt[grepl("CD8", garnett), GN_L3:='CD8+T']

plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L3, annot=GN_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = GN_L3][, GN_L3],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 



# --------------------------- brca3p level 4 ------------------------------------ #

brca3p_dt[, original_annotation_L4:=original_annotation]
brca3p_dt[!grepl("^T:", original_annotation), original_annotation_L4:='other']


brca3p_dt[, ImmC_L4:=ImmC]
brca3p_dt[!grepl("L:T:", ImmC_L4), ImmC_L4:='other']
brca3p_dt[, .N, by = ImmC_L4]


plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L4, annot=ImmC_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = ImmC_L4][, ImmC_L4],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




brca3p_dt[, SR_L4:=SingleR]
brca3p_dt[!grepl("T_cell", SingleR), SR_L4:='other']

plot_heatmap(table(brca3p_dt[, list(Known=original_annotation_L4, annot=SR_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca3p_dt[, .N, by = SR_L4][, SR_L4],
             ord2 =  brca3p_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




############## pbmc68k dataset  ##############


pbmc68k_meta <- data.table(readRDS(paste(data_folder, 'metadata/pbmc68k.meta.rds', sep = "/")))
setkey(pbmc68k_meta, barcodes)
pbmc68k_immc <- fread(paste(data_folder, 'immc/pbmc68k.output.txt', sep = "/"))
setkey(pbmc68k_immc, Cell)

load(paste(data_folder,'singleR/singler_pbmc68k.RData', sep = "/"))
pbmc68k_singler <- singler
pbmc68k_singler_dt <- data.table(Cell=pbmc68k_singler$singler[[1]]$SingleR.single$cell.names,
                                 singler=pbmc68k_singler$singler[[1]]$SingleR.single$labels)
setkey(pbmc68k_singler_dt, Cell)

pbmc68k_garnett_extend  <- readRDS(paste(data_folder,'garnett/pbmc68k.garnett.extend.rds', sep = "/"))
pbmc68k_garnett_extend_dt <- data.table(Cell = rownames(pbmc68k_garnett_extend), 
                                        garnett=pbmc68k_garnett_extend$cluster_ext_type)
setkey(pbmc68k_garnett_extend_dt, Cell)


pbmc68k_umap <- fread(paste(data_folder,'umap/out-pbmc_umap/models/umap.pca_1000.cluster_10.k_25.res_0.8.umap_25/coordinates.txt', sep = "/"))
setkey(pbmc68k_umap, Cell)

pbmc68k_dt <- data.table(
    merge(
      merge(
        merge(pbmc68k_umap[, .(Cell, UMAP1, UMAP2, Cluster)], 
              pbmc68k_immc, by.x='Cell', by.y='Cell', all.x = T),
        pbmc68k_singler_dt, by.x='Cell', by.y='Cell', all.x = T),
      pbmc68k_garnett_extend_dt, by.x='Cell', by.y='Cell', all.x = T))[, original_annotation:=pbmc68k_meta[Cell, 'celltype']]

pbmc68k_dt[, Cluster:=as.character(Cluster)]
colnames(pbmc68k_dt)[grepl("ImmC", colnames(pbmc68k_dt))] <- "ImmC"
colnames(pbmc68k_dt)[grepl("singler", colnames(pbmc68k_dt))] <- "SingleR"
setkey(pbmc68k_dt, Cell)

rm(singler)
rm(pbmc68k_singler)

# --------------------------- pbmc68k level 1 ------------------------------------ #
pbmc68k_dt[, original_annotation_L1:="not_specified"]

pbmc68k_dt[grepl("T|B|CD4|CD8|NK", original_annotation), original_annotation_L1:='L']
pbmc68k_dt[grepl("Den|Mono", original_annotation), original_annotation_L1:='M']
pbmc68k_dt[grepl("CD34", original_annotation), original_annotation_L1:='CD34']
pbmc68k_dt[, .N, by = original_annotation_L1]

pbmc68k_dt[, ImmC_L1:="not_specified"]
pbmc68k_dt[grepl("L:", ImmC), ImmC_L1:='L']
pbmc68k_dt[grepl("M:", ImmC), ImmC_L1:='M']
pbmc68k_dt[grepl("CD34", ImmC), ImmC_L1:='CD34']


plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L1, annot=ImmC_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = ImmC_L1][, ImmC_L1],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 




pbmc68k_dt[, SR_L1:=SingleR]
pbmc68k_dt[is.na(SR_L1), SR_L1:='Unknown']
pbmc68k_dt[grepl("T_cell|NK|B_cell:|Pre-B|^B_cell$",SingleR), SR_L1:='L']
pbmc68k_dt[grepl("Mono|Neu|Mac|DC|Erythroblast|Plate",SingleR), SR_L1:='M']
pbmc68k_dt[grepl("GMP|CD34\\+|CMP",SingleR), SR_L1:='CD34']
pbmc68k_dt[grepl("Epi|iPS|Embryonic_stem_cells|Fib|Endo|Tissue_stem_cells|Chon|Smooth|MEP",SingleR), SR_L1:='Non-immune']

plot_heatmap(table(pbmc68k_dt[!is.na(SR_L1), list(Known=original_annotation_L1, annot=SR_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = SR_L1][, SR_L1],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 


pbmc68k_dt[, GN_L1:=garnett]
pbmc68k_dt[grepl("T cell|NK|B", garnett), GN_L1:='L']
pbmc68k_dt[grepl("Mono|Dend", garnett), GN_L1:='M']
pbmc68k_dt[grepl("CD34", garnett), GN_L1:='CD34']

plot_heatmap(table(pbmc68k_dt[!is.na(SR_L1), list(Known=original_annotation_L1, annot=GN_L1)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = GN_L1][, GN_L1],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L1][, original_annotation_L1],
             measures = c('Recall', 'Precision', 'f1') ) 




# --------------------------- pbmc68k level 2 ------------------------------------ #


pbmc68k_dt[, original_annotation_L2:=original_annotation]

pbmc68k_dt[grepl("T|CD4", original_annotation), original_annotation_L2:='L:T']
pbmc68k_dt[grepl("CD34", original_annotation), original_annotation_L2:='CD34']
pbmc68k_dt[, .N, by = original_annotation_L2]



pbmc68k_dt[, ImmC_L2:=ImmC]
pbmc68k_dt[grepl("L:T", ImmC), ImmC_L2:='L:T']
pbmc68k_dt[grepl("DC", ImmC), ImmC_L2:='DC']
pbmc68k_dt[!ImmC_L2 %in% c('L:B', 'L:T', 'L:NK', 'DC', 'M:Mono','M:Mac', 'M:Mast', 'M:Neu'), ImmC_L2:='other']

plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L2, annot=ImmC_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = ImmC_L2][, ImmC_L2],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 




pbmc68k_dt[, SR_L2:="not_specified"]
pbmc68k_dt[grepl("B_cell:|Pre-B|^B_cell$", SingleR), SR_L2:='B']
pbmc68k_dt[grepl("T_cell", SingleR), SR_L2:='T']
pbmc68k_dt[grepl("NK", SingleR), SR_L2:='NK']
pbmc68k_dt[grepl("Mono", SingleR), SR_L2:='Mono']
pbmc68k_dt[grepl("Mac", SingleR), SR_L2:='Mac']
pbmc68k_dt[grepl("Neu", SingleR), SR_L2:='Neu']
pbmc68k_dt[grepl("DC", SingleR), SR_L2:='mDC']
pbmc68k_dt[grepl("Ery|Platelet", SingleR), SR_L2:='Ery/Mega']

pbmc68k_dt[grepl("GMP|CD34\\+|CMP", SingleR), SR_L2:='CD34']
pbmc68k_dt[grepl("Epi|iPS|Embryonic_stem_cells|Fib|Endo|Tissue_stem_cells|Chon|Smooth", SingleR), SR_L2:='Non-immune']

plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L2, annot=SR_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = SR_L2][, SR_L2],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 





pbmc68k_dt[, GN_L2:=garnett]
pbmc68k_dt[grepl("T cell", garnett), GN_L2:='T']



plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L2, annot=GN_L2)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = GN_L2][, GN_L2],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L2][, original_annotation_L2],
             measures = c('Recall', 'Precision', 'f1') ) 




# --------------------------- pbmc68k level 3 ------------------------------------ #

pbmc68k_dt[, original_annotation_L3:='not_specified']
pbmc68k_dt[grepl("CD4", original_annotation), original_annotation_L3:='CD4+T']
pbmc68k_dt[grepl("CD8", original_annotation), original_annotation_L3:='CD8+T']


pbmc68k_dt[, ImmC_L3:='not_specified']
pbmc68k_dt[grepl("L:T:CD4", ImmC), ImmC_L3:='CD4+T']
pbmc68k_dt[grepl("L:T:CD8", ImmC), ImmC_L3:='CD8+T']

plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L3, annot=ImmC_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = ImmC_L3][, ImmC_L3],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 



pbmc68k_dt[, SR_L3:="not_specified"]

pbmc68k_dt[grepl("CD4", SingleR), SR_L3:='CD4+T']
pbmc68k_dt[grepl("CD8", SingleR), SR_L3:='CD8+T']
pbmc68k_dt[grepl("DC", SingleR), SR_L3:='mDC']

plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L3, annot=SR_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = SR_L3][, SR_L3],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 




pbmc68k_dt[, GN_L3:=garnett]


plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L3, annot=GN_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = GN_L3][, GN_L3],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 






# --------------------------- pbmc68k level 4 ------------------------------------ #

pbmc68k_dt[, original_annotation_L4:=original_annotation]
pbmc68k_dt[!grepl("T|CD4", original_annotation), original_annotation_L4:='other']


pbmc68k_dt[, ImmC_L4:=ImmC]
pbmc68k_dt[!grepl("L:T:", ImmC_L4), ImmC_L4:='other']
pbmc68k_dt[, .N, by = ImmC_L4]


plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L4, annot=ImmC_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = ImmC_L4][, ImmC_L4],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




pbmc68k_dt[, SR_L4:=SingleR]
pbmc68k_dt[!grepl("T_cell", SingleR), SR_L4:='other']

plot_heatmap(table(pbmc68k_dt[, list(Known=original_annotation_L4, annot=SR_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=pbmc68k_dt[, .N, by = SR_L4][, SR_L4],
             ord2 =  pbmc68k_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




##################### HCC #############################

hcc_meta <- data.table(readRDS(paste(data_folder, 'metadata/HCC-ZZ-meta.rds', sep = "/")))
hcc_clusternames <- data.table(readRDS(paste(data_folder,'metadata/HCC-ZZ-clusternames.rds', sep = "/")))
setkey(hcc_meta, UniqueCell_ID)
setkey(hcc_clusternames, clusterID)
hcc_immc <- fread(paste(data_folder,'immc/hcc.output.txt', sep = "/"))


load(paste(data_folder,'singleR/singler_hcc.RData', sep = "/"))
hcc_singler_dt <- data.table(Cell=singler$singler[[1]]$SingleR.single$cell.names,
                             singler=singler$singler[[1]]$SingleR.single$labels)
setkey(hcc_singler_dt, Cell)


hcc_garnett_extend  <- readRDS(paste(data_folder,'garnett/hcc.garnett.extend.rds',sep = "/"))

hcc_garnett_extend_dt <- data.table(Cell = rownames(hcc_garnett_extend), 
                                    garnett=hcc_garnett_extend$cluster_ext_type)
setkey(hcc_garnett_extend_dt, Cell)


hcc_immc_ctrl <- fread(paste(data_folder, 'immc_ctrl/hcc.output.txt', sep = "/"))
colnames(hcc_immc_ctrl)[2] <- 'ctrl'
hcc_umap <- fread(paste(data_folder,'umap/out-hcc_all_umap/models/umap.pca_0500.cluster_25.k_25.res_1.0.umap_25/coordinates.txt', sep = "/"))
setkey(hcc_umap, Cell)

rm(singler)
rm(hcc_singler)


hcc_dt <- data.table(
    merge(
      merge(
        merge(hcc_umap[, .(Cell, UMAP1, UMAP2, Cluster)], 
              hcc_immc, by.x='Cell', by.y='Cell', all.x = T),
        hcc_singler_dt, by.x='Cell', by.y='Cell', all.x = T),
      hcc_garnett_extend_dt, by.x='Cell', by.y='Cell', all.x = T))[, original_annotation:=hcc_clusternames[hcc_meta[Cell, majorCluster], clusterName]]

hcc_dt[, Cluster:=as.character(Cluster)]
colnames(hcc_dt)[grepl("ImmC", colnames(hcc_dt))] <- "ImmC"
colnames(hcc_dt)[grepl("singler", colnames(hcc_dt))] <- "SingleR"
setkey(hcc_dt, Cell)


# --------------------------- hcc level 3 ------------------------------------ #

hcc_dt[, .N, by = original_annotation]
hcc_dt[, original_annotation_L3:=original_annotation]
hcc_dt[grepl("CD4|helper|reg", original_annotation), original_annotation_L3:='CD4+T']
hcc_dt[grepl("CD8", original_annotation), original_annotation_L3:='CD8+T']
hcc_dt[, .N, by = original_annotation_L3]

hcc_dt[, ImmC_L3:='not_specified']
hcc_dt[grepl("L:T:CD4", ImmC), ImmC_L3:='CD4+T']
hcc_dt[grepl("L:T:CD8", ImmC), ImmC_L3:='CD8+T']

plot_heatmap(table(hcc_dt[, list(Known=original_annotation_L3, annot=ImmC_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=hcc_dt[, .N, by = ImmC_L3][, ImmC_L3],
             ord2 =  hcc_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 



hcc_dt[, SR_L3:="not_specified"]

hcc_dt[grepl("CD4", SingleR), SR_L3:='CD4+T']
hcc_dt[grepl("CD8", SingleR), SR_L3:='CD8+T']
hcc_dt[grepl("DC", SingleR), SR_L3:='mDC']

plot_heatmap(table(hcc_dt[, list(Known=original_annotation_L3, annot=SR_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=hcc_dt[, .N, by = SR_L3][, SR_L3],
             ord2 =  hcc_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 




hcc_dt[, GN_L3:=garnett]


plot_heatmap(table(hcc_dt[, list(Known=original_annotation_L3, annot=GN_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=hcc_dt[, .N, by = GN_L3][, GN_L3],
             ord2 =  hcc_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 






# --------------------------- hcc level 4 ------------------------------------ #

hcc_dt[, .N, by = original_annotation]
hcc_dt[, original_annotation_L4:=original_annotation]
hcc_dt[grepl("Treg", original_annotation), original_annotation_L4:='Treg']

hcc_dt[, ImmC_L4:=ImmC]
hcc_dt[!grepl("L:T|MAIT", ImmC_L4), ImmC_L4:='other']
hcc_dt[, .N, by = ImmC_L4]


plot_heatmap(table(hcc_dt[, list(Known=original_annotation_L4, annot=ImmC_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=hcc_dt[, .N, by = ImmC_L4][, ImmC_L4],
             ord2 =  hcc_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




hcc_dt[, SR_L4:=SingleR]
hcc_dt[!grepl("T_cell", SingleR), SR_L4:='other']

plot_heatmap(table(hcc_dt[, list(Known=original_annotation_L4, annot=SR_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=hcc_dt[, .N, by = SR_L4][, SR_L4],
             ord2 =  hcc_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 



##################### brca5p #############################


brca5p_clusternames <- list(as.character(c(7,16,22,33)),
                            as.character(c(5,6,9,26,30)),
                            as.character(c(3,4)),
                            as.character(c(2,14,20,29,34)),
                            as.character(c(12,17,32,21)),
                            as.character(c(10,19)),
                            as.character(c(1,8,11,13,15,18,24,25,28,31,27)),
                            as.character(c(23)))

names(brca5p_clusternames) <- c( "T:CD4+Naive", "T:CD4+CM", "T:CD4+EM", "Treg",
                                 "T:CD8+Naive", "T:CD8+CM", "T:CD8+EM",
                                 "NKT")



brca5p_clusternames_dt <- data.table(clusterID=unlist(brca5p_clusternames), clusterName=gsub("[0-9]+$", "", names(unlist(brca5p_clusternames))))

setkey(brca5p_clusternames_dt, clusterID)


brca5p_immc <- fread(paste(data_folder, 'immc/brca5p.output.txt', sep = "/"))


load(paste(data_folder,"singleR/singler_brca5p.RData", sep = "/"))


brca5p_singler_dt <- data.table(Cell=singler$singler[[1]]$SingleR.single$cell.names,
                                SingleR=singler$singler[[1]]$SingleR.single$labels)
setkey(brca5p_singler_dt, Cell)

brca5p_garnett_extend  <- readRDS(paste(data_folder,'garnett/brca.5p.garnett.extend.rds', sep = "/"))
brca5p_garnett_extend_dt <- data.table(Cell = rownames(brca5p_garnett_extend), 
                                       garnett=brca5p_garnett_extend$cluster_ext_type)

setkey(brca5p_garnett_extend_dt, Cell)

brca5p_immc_ctrl <- fread(paste(data_folder, 'immc_ctrl/brca5p.output.txt', sep = "/"))
colnames(brca5p_immc_ctrl)[2] <- 'ctrl'

brca5p_umap <- fread(paste(data_folder,'umap/out-brca5p_umap/models//umap.pca_0500.cluster_25.k_25.res_1.0.umap_25/coordinates.txt', sep = "/"))
setkey(brca5p_umap, Cell)


brca5p_dt <- data.table(
merge(
      merge(
        merge(brca5p_umap[, .(Cell, UMAP1, UMAP2, Cluster)], 
              brca5p_immc, by.x='Cell', by.y='Cell', all.x = T),
        brca5p_singler_dt, by.x='Cell', by.y='Cell', all.x = T),
      brca5p_garnett_extend_dt, by.x='Cell', by.y='Cell', all.x = T))[, original_annotation:=brca5p_clusternames_dt[brca3p_meta[Cell, 'cluster'], 'clusterName']]

brca5p_dt[, Cluster:=as.character(Cluster)]
colnames(brca5p_dt)[grepl("ImmC", colnames(brca5p_dt))] <- "ImmC"
colnames(brca5p_dt)[grepl("SingleR", colnames(brca5p_dt))] <- "SingleR"
setkey(brca5p_dt, Cell)


rm(singler)


# --------------------------- brca5p level 3 ------------------------------------ #

brca5p_dt[, .N, by = original_annotation]
brca5p_dt[, original_annotation_L3:=original_annotation]
brca5p_dt[grepl("CD4|helper|reg", original_annotation), original_annotation_L3:='CD4+T']
brca5p_dt[grepl("CD8", original_annotation), original_annotation_L3:='CD8+T']
brca5p_dt[, .N, by = original_annotation_L3]

brca5p_dt[, ImmC_L3:='not_specified']
brca5p_dt[grepl("L:T:CD4", ImmC), ImmC_L3:='CD4+T']
brca5p_dt[grepl("L:T:CD8", ImmC), ImmC_L3:='CD8+T']

plot_heatmap(table(brca5p_dt[, list(Known=original_annotation_L3, annot=ImmC_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca5p_dt[, .N, by = ImmC_L3][, ImmC_L3],
             ord2 =  brca5p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 



brca5p_dt[, SR_L3:="not_specified"]

brca5p_dt[grepl("CD4", SingleR), SR_L3:='CD4+T']
brca5p_dt[grepl("CD8", SingleR), SR_L3:='CD8+T']


plot_heatmap(table(brca5p_dt[, list(Known=original_annotation_L3, annot=SR_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca5p_dt[, .N, by = SR_L3][, SR_L3],
             ord2 =  brca5p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 




brca5p_dt[, GN_L3:=garnett]


plot_heatmap(table(brca5p_dt[, list(Known=original_annotation_L3, annot=GN_L3)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca5p_dt[, .N, by = GN_L3][, GN_L3],
             ord2 =  brca5p_dt[, .N, by = original_annotation_L3][, original_annotation_L3],
             measures = c('Recall', 'Precision', 'f1') ) 






# --------------------------- brca5p level 4 ------------------------------------ #

brca5p_dt[, .N, by = original_annotation]
brca5p_dt[, original_annotation_L4:=original_annotation]

brca5p_dt[, ImmC_L4:=ImmC]
brca5p_dt[!grepl("L:T", ImmC_L4), ImmC_L4:='other']
brca5p_dt[, .N, by = ImmC_L4]


plot_heatmap(table(brca5p_dt[, list(Known=original_annotation_L4, annot=ImmC_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca5p_dt[, .N, by = ImmC_L4][, ImmC_L4],
             ord2 =  brca5p_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 




brca5p_dt[, SR_L4:=SingleR]
brca5p_dt[!grepl("T_cell", SingleR), SR_L4:='other']

plot_heatmap(table(brca5p_dt[, list(Known=original_annotation_L4, annot=SR_L4)]), title = "", 
             minshow = 0, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=brca5p_dt[, .N, by = SR_L4][, SR_L4],
             ord2 =  brca5p_dt[, .N, by = original_annotation_L4][, original_annotation_L4],
             measures = c('Recall', 'Precision', 'f1') ) 





# NA: not detected by predictor
# -1: not detected in original annotation
brca3p_immc_rec <- list(
  'Level 1' =c(-1, 98,97),
  'Level 2' =c(89, 97, 49, 56, 23, 93, 70, 75),
  'Level 3' =c(89, 39, 37, 79),
  'Level 4' =c(48, 10, 9, -1, 71, -1, 3, 0, 14, -1, -1)
)

brca3p_immc_prec <- list(
  'Level 1' =c(-1, 99,89),
  'Level 2' =c(98, 88, 88, 49, 66, 60, 61, 35),
  'Level 3' =c(69, 47, 59, 17),
  'Level 4' =c(87, 3, 0, -1, 7, -1, 20, 0, 16, -1, -1)
)


brca3p_immc_f1 <- list(
  'Level 1' =c(-1, 99, 93),
  'Level 2' =c(94, 92, 63, 52, 35, 73, 65, 47),
  'Level 3' =c(78, 43, 45, 28),
  'Level 4' =c(62, 5, 1, -1, 12, -1, 6, 0, 15, -1, -1)
)


brca3p_sr_rec <- list(
  'Level 1' =c(-1, 96, 89),
  'Level 2' =c(86, 92, 53, 3, 70, 76, NA, 74),
  'Level 3' =c(95, 24, 3, NA),
  'Level 4' =c(41, 50, 37, -1, 0, -1, 4, 7, 0, -1, -1))

brca3p_sr_prec <- list(
  'Level 1' =c(-1, 99, 98),
  'Level 2' =c(95, 87, 65, 22, 50, 70, NA, 80),
  'Level 3' =c(68, 45, 18, NA),
  'Level 4' =c(98, 7, 1, -1, 33, -1, 63, 1, 2, -1, -1)
)


brca3p_sr_f1 <- list(
  'Level 1' =c(-1, 97, 93),
  'Level 2' =c(90, 89, 59, 6, 58, 73, NA, 77),
  'Level 3' =c(79, 31, 5, NA),
  'Level 4' =c(58, 13, 1, -1, 1, -1, 7, 2, 0, -1, -1)
)




brca3p_gn_rec <- list(
  'Level 1' =c(-1, 83, 55),
  'Level 2' =c(94, 86, 50, 19, 42, NA, NA, NA),
  'Level 3' =c(31, 7, NA, NA),
  'Level 4' = c(NA, NA, NA, -1, NA, -1, NA, NA, NA, -1, -1)
)


brca3p_gn_prec <- list(
  'Level 1' =c(-1, 97, 91),
  'Level 2' =c(97, 93, 74, 41, 29, NA, NA, NA),
  'Level 3' =c(70, 65, NA, NA),
  'Level 4' = c(NA, NA, NA, -1, NA, -1, NA, NA, NA, -1, -1)
)

brca3p_gn_f1 <- list(
  'Level 1' =c(-1, 90, 68),
  'Level 2' =c(95, 89, 60, 26, 34, NA, NA, NA),
  'Level 3' =c(43, 13, NA, NA),
  'Level 4' = c(NA, NA, NA, -1, NA, -1, NA, NA, NA, -1, -1))





pbmc68k_immc_rec <- list(
  'Level 1' = c(14, 100, 97),
  'Level 2' = c(62, 99, 57, 40, 33, -1, -1, -1),
  'Level 3' = c(98, 37, -1, -1),
  'Level 4' = c(37, -1, -1, -1, 47, 0, 19, -1, -1, -1, -1))

pbmc68k_immc_prec <- list(
  'Level 1' = c(50, 100, 93),
  'Level 2' = c(100, 89, 97, 81, 55, -1, -1, -1),
  'Level 3' = c(30, 78, -1, -1),
  'Level 4' = c(7, -1, -1, -1, 31, 0, 68, -1, -1, -1, -1))

pbmc68k_immc_f1 <- list(
  'Level 1' = c(22, 100, 95),
  'Level 2' = c(76, 94, 72, 54, 42,  -1, -1, -1),
  'Level 3' = c(46, 50, -1, -1),
  'Level 4' = c(11,  -1, -1, -1, 37, 0, 29,  -1, -1, -1, -1))



pbmc68k_sr_rec <- list(
  'Level 1' = c(20, 100, 88),
  'Level 2' = c(64, 99, 49, NA, 93, -1, -1, -1),
  'Level 3' = c(99, 14, -1, -1),
  'Level 4' = c(64, -1, -1, -1, NA, NA, 0, -1, -1, -1, -1))

pbmc68k_sr_prec <- list(
  'Level 1' = c(77, 99, 96),
  'Level 2' = c(88, 88, 95, NA, 60, -1, -1, -1),
  'Level 3' = c(24, 54, -1, -1),
  'Level 4' = c(7, -1, -1, -1, NA, NA, 0, -1, -1, -1, -1))

pbmc68k_sr_f1 <- list(
  'Level 1' = c(32, 99, 92),
  'Level 2' = c(74, 93, 64, NA, 73, -1, -1, -1),
  'Level 3' = c(39, 22, -1, -1),
  'Level 4' = c(12, -1, -1, -1, NA, NA, 0, -1, -1, -1, -1 ))


pbmc68k_gn_rec <- list(
  'Level 1' = c(7, 100, 39),
  'Level 2' = c(65, 99, 51, 21, 38, -1, -1, -1),
  'Level 3' = c(63, 10, -1, -1),
  'Level 4' = c(NA, -1, -1, -1, NA, NA, NA, -1, -1, -1, -1 ))

pbmc68k_gn_prec <- list(
  'Level 1' = c(37, 98, 94),
  'Level 2' = c(97, 88, 80, 81, 73, -1, -1, -1),
  'Level 3' = c(55, 47, -1, -1),
  'Level 4' = c(NA, -1, -1, -1, NA, NA, NA, -1, -1, -1, -1 ))

pbmc68k_gn_f1 <- list(
  'Level 1' = c(12, 99, 55),
  'Level 2' = c(78, 93, 62, 33, 50, -1, -1, -1),
  'Level 3' = c(59, 16, -1, -1),
  'Level 4' = c(NA, -1, -1, -1, NA, NA, NA, -1, -1, -1, -1 ))



hcc_immc_rec <- list(
  'Level 3'=c(85, 90, -1, -1),
  'Level 4'=c(24, -1, -1, 37, 91, 1, 50, -1, 5, 4, 41)
)

hcc_immc_prec <- list(
  'Level 3'=c(76, 34, -1, -1),
  'Level 4'=c(43, -1, -1, 58, 52, 17, 58, -1, 3, 18, 88)
)
# cd8+em are mostly predicted to cd8+emra
hcc_immc_f1 <- list(
  'Level 3'=c(80, 50, -1, -1),
  'Level 4'=c(31, -1, -1, 45, 67, 1, 53, -1, 4, 7, 56)
)


hcc_sr_rec <- list(
  'Level 3'=c(94, 54, -1, -1),
  'Level 4'=c(20, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

hcc_sr_prec <- list(
  'Level 3'=c(61, 43, -1, -1),
  'Level 4'=c(49, -1, -1, NA, NA, NA, NA,-1,NA, NA, NA )
)

hcc_sr_f1 <- list(
  'Level 3'=c(74, 48, -1, -1),
  'Level 4'=c(29, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)



hcc_gn_rec <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

hcc_gn_prec <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

hcc_gn_f1 <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)


brca5p_immc_rec <- list(
  'Level 3'=c(85, 71, -1, -1),
  'Level 4'=c(10, 28, 11, -1, 76, -1, 1, 0, 20, -1, -1)
)

brca5p_immc_prec <- list(
  'Level 3'=c(86, 93, -1, -1),
  'Level 4'=c(17, 52, 33, -1, 45, -1, 6, 14, 58, -1, -1)
)

brca5p_immc_f1 <- list(
  'Level 3'=c(86, 81, -1, -1),
  'Level 4'=c(13, 36, 16, -1, 57, -1, 2, 0, 30, -1, -1)
)


brca5p_sr_rec <- list(
  'Level 3'=c(99, 29, -1, -1),
  'Level 4'=c(6, 75, 28, -1, 0, -1, 2, 0, 0, -1, -1)
)

brca5p_sr_prec <- list(
  'Level 3'=c(58, 96, -1, -1),
  'Level 4'=c(14, 27, 8, -1, 0, -1, 28, 1, 64, -1, -1)
)

brca5p_sr_f1 <- list(
  'Level 3'=c(73, 44, -1, -1),
  'Level 4'=c(9, 39, 13, -1, 0, -1, 4, 0, 1, -1, -1)
)



brca5p_gn_rec <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

brca5p_gn_prec <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

brca5p_gn_f1 <- list(
  'Level 3'=c(NA, NA, -1, -1),
  'Level 4'=c(NA, -1, -1, NA, NA, NA, NA, -1, NA, NA, NA)
)

# 
# i <- 4
# round(brca3p_immc_rec[[i]] * brca3p_immc_prec[[i]]/(brca3p_immc_rec[[i]] + brca3p_immc_prec[[i]]) *2,0)
# brca3p_immc_f1[[i]]
# 
# round(brca3p_sr_rec[[i]] * brca3p_sr_prec[[i]]/(brca3p_sr_rec[[i]] + brca3p_sr_prec[[i]]) *2, 0)
# brca3p_sr_f1[[i]]
# 
# round(brca3p_gn_rec[[i]] * brca3p_gn_prec[[i]]/(brca3p_gn_rec[[i]] + brca3p_gn_prec[[i]]) *2,0)
# brca3p_gn_f1[[i]]

#----------------- all performance ------------------#


perf <- NULL
for (i in 1:4){
  
  tmp <- data.table(celltype=levels_celltype[[paste('Level', i)]],
             level = i,
             immc_f1=brca3p_immc_f1[[paste('Level', i)]],
             sr_f1=brca3p_sr_f1[[paste('Level', i)]],
             gn_f1=brca3p_gn_f1[[paste('Level', i)]],
             immc_rec=brca3p_immc_rec[[paste('Level', i)]],
             sr_rec=brca3p_sr_rec[[paste('Level', i)]],
             gn_rec=brca3p_gn_rec[[paste('Level', i)]],
             immc_prec=brca3p_immc_prec[[paste('Level', i)]],
             sr_prec=brca3p_sr_prec[[paste('Level', i)]],
             gn_prec=brca3p_gn_prec[[paste('Level', i)]],
             dataset='brca3p'
             )
  
  perf <- rbind(perf, tmp)
}

for (i in 1:4){
  
  tmp <- data.table(celltype=levels_celltype[[paste('Level', i)]],
                    level = i,
                    immc_f1=pbmc68k_immc_f1[[paste('Level', i)]],
                    sr_f1=pbmc68k_sr_f1[[paste('Level', i)]],
                    gn_f1=pbmc68k_gn_f1[[paste('Level', i)]],
                    immc_rec=pbmc68k_immc_rec[[paste('Level', i)]],
                    sr_rec=pbmc68k_sr_rec[[paste('Level', i)]],
                    gn_rec=pbmc68k_gn_rec[[paste('Level', i)]],
                    immc_prec=pbmc68k_immc_prec[[paste('Level', i)]],
                    sr_prec=pbmc68k_sr_prec[[paste('Level', i)]],
                    gn_prec=pbmc68k_gn_prec[[paste('Level', i)]],
                    dataset='pbmc68k'
  )
  
  perf <- rbind(perf, tmp)
}

for (i in 3:4){
  
  tmp <- data.table(celltype=levels_celltype[[paste('Level', i)]],
                    level = i,
                    immc_f1=brca5p_immc_f1[[paste('Level', i)]],
                    sr_f1=brca5p_sr_f1[[paste('Level', i)]],
                    gn_f1=brca5p_gn_f1[[paste('Level', i)]],
                    immc_rec=brca5p_immc_rec[[paste('Level', i)]],
                    sr_rec=brca5p_sr_rec[[paste('Level', i)]],
                    gn_rec=brca5p_gn_rec[[paste('Level', i)]],
                    immc_prec=brca5p_immc_prec[[paste('Level', i)]],
                    sr_prec=brca5p_sr_prec[[paste('Level', i)]],
                    gn_prec=brca5p_gn_prec[[paste('Level', i)]],
                    dataset='brca5p'
  )
  
  perf <- rbind(perf, tmp)
}

for (i in 3:4){
  
  tmp <- data.table(celltype=levels_celltype[[paste('Level', i)]],
                    level = i,
                    immc_f1=hcc_immc_f1[[paste('Level', i)]],
                    sr_f1=hcc_sr_f1[[paste('Level', i)]],
                    gn_f1=hcc_gn_f1[[paste('Level', i)]],
                    immc_rec=hcc_immc_rec[[paste('Level', i)]],
                    sr_rec=hcc_sr_rec[[paste('Level', i)]],
                    gn_rec=hcc_gn_rec[[paste('Level', i)]],
                    immc_prec=hcc_immc_prec[[paste('Level', i)]],
                    sr_prec=hcc_sr_prec[[paste('Level', i)]],
                    gn_prec=hcc_gn_prec[[paste('Level', i)]],
                    dataset='hcc'
  )
  
  perf <- rbind(perf, tmp)
}

perf2 <- melt(perf, id.vars = c('celltype', 'level', 'dataset'),
              variable.factor = F, value.factor = F)
perf2[is.na(value), value:=0]

perf2[,method:=gsub("_.*", "", variable)]
perf2[, measure:=gsub(".*_", "", variable)]
ggplot(perf2[value!=-1,], 
       aes(as.character(level), value, fill=factor(method, levels = c('immc', 'sr', 'gn')))) + 
  geom_boxplot(outlier.shape = NA) +  
  facet_grid(factor(dataset, levels = c('pbmc68k', 'brca3p','brca5p', 'hcc'))~factor(measure, levels = c('rec', 'prec', 'f1')))+
  theme_classic() + 
  theme(axis.text.x = element_text(angle =0, size =12, hjust = 1,color= "black",), 
        axis.text.y = element_text(size = 12, hjust = 1,color = "black"),
        axis.title.x = element_text(size = 12, hjust = 1), 
        axis.title.y = element_text(size = 12, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.pos = "right")




# Fig3b

dcast(perf2[value!=-1 & dataset == 'brca3p' & measure == 'rec', mean(value), 
                 by = list(method, dataset, level, measure)], dataset+level~method)[, list(immc-gn, immc-sr)]


perf2[value!=-1 & dataset == 'brca3p' & variable == 'immc' & !celltype %in% c('L','M', 'DC'), mean(value)]
perf2[value!=-1 & dataset == 'brca3p' & variable == 'sr' & !celltype %in% c('L','M', 'DC'), mean(value)]
perf2[value!=-1 & dataset == 'brca3p' & variable == 'gn' & !celltype %in% c('L','M', 'DC'), mean(value)]




#########################################################################
#                               Figure 4                                #
#########################################################################



set.seed(100)
brca3p_subset_dt <- brca3p_dt[brca3p_dt[,sample(Cell, 200, replace = T) ,by=original_annotation][, V1],]

brca3p_celltypes <- c('MAST', 'B', 'mDC', 'NK', 'MONOCYTE', 
                      'CD4+T' ,'pDC', 'CD8+T' ,'NEUTROPHIL' ,'NKT', 
                      'MACROPHAGE' , 'Treg' , 'T', 'Other', 'Unassigned')


# rename Original annotation cell types
brca3p_subset_dt[grepl("^NKT", original_annotation), original_annotation:="NKT"]
brca3p_subset_dt[grepl("T:Reg", original_annotation), original_annotation:="Treg"]
brca3p_subset_dt[grepl("T:CD4", original_annotation), original_annotation:="CD4+T"]
brca3p_subset_dt[grepl("T:CD8", original_annotation), original_annotation:="CD8+T"]
brca3p_subset_dt[, original_annotation:=gsub(":[^:]+.*$", "", original_annotation)]
brca3p_subset_dt[, original_annotation:=gsub(":$", "", original_annotation)]
brca3p_subset_dt[, .N, by = original_annotation]
# rename ImmClassifier output cell types
brca3p_subset_dt[grepl("B", ImmC), ImmC:="B"]
brca3p_subset_dt[grepl("NK", ImmC), ImmC:="NK"]
brca3p_subset_dt[grepl("mDC", ImmC), ImmC:="mDC"]
brca3p_subset_dt[grepl("pDC", ImmC), ImmC:="pDC"]
brca3p_subset_dt[grepl("reg", ImmC), ImmC:="Treg"]
brca3p_subset_dt[grepl("CD4", ImmC), ImmC:="CD4+T"]
brca3p_subset_dt[grepl("CD8", ImmC), ImmC:="CD8+T"]
brca3p_subset_dt[grepl("Mono", ImmC), ImmC:="MONOCYTE"]
brca3p_subset_dt[grepl("Mac", ImmC), ImmC:="MACROPHAGE"]
brca3p_subset_dt[grepl("Mast", ImmC), ImmC:="MAST"]
brca3p_subset_dt[grepl("Neu", ImmC), ImmC:="NEUTROPHIL"]

brca3p_subset_dt[!ImmC %in% brca3p_celltypes, ImmC:="Other"]
brca3p_subset_dt[, .N, by = ImmC]


# rename SingleR output cell types
brca3p_subset_dt[is.na(SingleR), SingleR:="Unassigned"]
brca3p_subset_dt[grepl("B_cell:|^B_cell$|Pre-B", SingleR), SingleR:="B"]
brca3p_subset_dt[grepl("NK_cell", SingleR), SingleR:="NK"]
brca3p_subset_dt[grepl("Monocyte", SingleR), SingleR:="MONOCYTE"]
brca3p_subset_dt[grepl("DC:m", SingleR), SingleR:="mDC"]
brca3p_subset_dt[grepl("Macrophage", SingleR), SingleR:="MACROPHAGE"]
brca3p_subset_dt[grepl("Treg", SingleR), SingleR:="Treg"]
brca3p_subset_dt[grepl("T_cell:CD4", SingleR), SingleR:="CD4+T"]
brca3p_subset_dt[grepl("T_cell:CD8", SingleR), SingleR:="CD8+T"]
brca3p_subset_dt[grepl("Neutrophil", SingleR), SingleR:="NEUTROPHIL"]
brca3p_subset_dt[!SingleR %in% brca3p_celltypes, SingleR:="Other"]
brca3p_subset_dt[, .N, by = SingleR]
# rename Garnett output cell types
brca3p_subset_dt[grepl("B cells", garnett), garnett:="B"]
brca3p_subset_dt[grepl("CD4 T cells", garnett), garnett:="CD4+T"]
brca3p_subset_dt[grepl("CD8 T cells", garnett), garnett:="CD8+T"]
brca3p_subset_dt[grepl("Den", garnett), garnett:="mDC"]
brca3p_subset_dt[grepl("Mono", garnett), garnett:="MONOCYTE"]
brca3p_subset_dt[grepl("NK", garnett), garnett:="NK"]
brca3p_subset_dt[grepl("T cells", garnett), garnett:="T"]
brca3p_subset_dt[grepl("Unknown", garnett), garnett:="Unassigned"]
brca3p_subset_dt[!garnett %in% brca3p_celltypes, garnett:="Other"]


brca3p_subset_dt[, random:=brca3p_subset_dt$original_annotation[sample(1:nrow(brca3p_subset_dt), nrow(brca3p_subset_dt))]]
brca3p_subset_dt <- brca3p_subset_dt[, list(Cell, UMAP1, UMAP2, Cluster, original_annotation, random, ImmC, SingleR, garnett)]

order_methods <- c('ImmC', 'SingleR', 'garnett','original_annotation', 'random')

melt(brca3p_subset_dt,
            id.vars = c('Cell', 'UMAP1', 'UMAP2', 'Cluster')) %>%
  mutate(value = factor(value, levels=brca3p_celltypes)) %>%
  mutate(variable = factor(variable, levels=order_methods)) %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = value), size = .5, pch = 18) + 
  theme_minimal() + 
  theme(legend.pos = "bottom", 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = .1),
        axis.text.x = element_text(angle = 0, size = 8, hjust = 1,color= "black"), 
        axis.text.y = element_text(size = 8, hjust = 1,color = "black"),
        axis.title.x = element_text(angle = 0, size = 8, hjust = .5,color= "black"), 
        axis.title.y = element_text(angle = 90, size = 8, hjust = .5,color= "black"), 
        strip.text.x = element_text(size = 10, colour = "black", angle = 0)) +
  coord_equal () + 
  scale_color_manual(values = c(colorRampPalette(brewer.pal(12, "Paired"))(13), 'lightgrey', 'black'),guide = guide_legend(nrow=3,override.aes = list(size=4))) +
  facet_wrap(~variable, ncol = 3)



dist_plot(brca3p_subset_dt, palette='Paired', n1=12, n2=13)

brca3p_celltypes_alluvium <- c('mDC', 'pDC', 'MONOCYTE',  'MACROPHAGE' ,'NEUTROPHIL' ,'MAST',  
                               'B', 'CD4+T' ,'Treg',  'CD8+T' ,'NK','NKT', 
                              'Other', 'Unassigned')

brca3p_subset_dt[, .N, by = c('ImmC', 'SingleR', 'original_annotation')] %>%
  mutate(ImmC = factor(ImmC, levels = brca3p_celltypes_alluvium),
         SingleR = factor(SingleR, levels = brca3p_celltypes_alluvium),
         original_annotation = factor(original_annotation, levels = brca3p_celltypes_alluvium)
  ) %>%
  ggplot(aes( axis1 = ImmC, axis2 = original_annotation, axis3 = SingleR,y = N)) +
  scale_x_discrete(limits = c('ImmC', 'original_annotation', 'SingleR'), expand = c(.1, 0)) +
  scale_fill_manual(values =sample(colorRampPalette(brewer.pal(11, "Spectral"))(14)))+
  geom_alluvium(aes(fill = original_annotation)) +
  geom_stratum(width = .5, color = 'darkgrey') + 
  geom_text(stat = "stratum", infer.label = T)+
  #ggfittext::geom_fit_text(stat = "stratum", infer.label = T) +
  theme_minimal() + #coord_flip()+
  xlab("Method") +
  theme(legend.pos = "none") 



#########################################################################
#                               Figure 5                                #
#########################################################################

adt <- fread('datasets/pont/GSM4288255_Antibodies_count.csv')
count <- fread('datasets/pont/GSM4288255_Raw_UMI_genes.tsv')
umap <- fread('datasets/pont/GSM4288255_UMAP_coordinate.csv')
# 5559 T cells
singler_annot <- fread('datasets/pont/comparison/pont_singler_annot.txt')
colnames(singler_annot)[1] <- 'SingleR'
immc_annot <- fread(paste(data_folder, 'immc/pont.output.txt', sep = "/"))
colnames(immc_annot)[2] <- 'ImmC'



meta <- merge(merge(immc_annot,
                    merge(singler_annot, merge(umap,adt, by.x = 'id', by.y = 'id'), 
                          by.x = 'detail_cellid', by.y = 'id'),
                    by.y = 'detail_cellid', by.x = 'Cell'),
              count[, c('V1', 'CD4', 'CD8A', 'CD8B')],
              by.x = 'Cell', by.y = 'V1')

setkey(meta, Cell)
# set up manual gate

cd4_gate <- 3.25
cd8_gate <- 2.82

cd62_gate_cd4t <- 2.8 
cd45ra_gate_cd4t <- 2.6 


cd62_gate_cd8t <- 3.2 
cd45ra_gate_cd8t <- 3.75



cd19_gate <- 3.2
cd16_gate <- 2.65
delta <- 0#.5
cd3_gate <- 3


# set up truth label
meta[, truth:=NA]
meta[, truth:=ifelse(CD3_Ab>cd3_gate & !(UMAP_1>2.5 & UMAP_2<2.5) & (CD4_Ab>cd4_gate & CD8_Ab>cd8_gate), 'T:CD4+CD8+', truth)]
meta[, truth:=ifelse(CD3_Ab>cd3_gate & !(UMAP_1>2.5 & UMAP_2<2.5) & (CD4_Ab>cd4_gate & CD8_Ab<=cd8_gate), 'T:CD4+CD8-', truth)]
meta[, truth:=ifelse(CD3_Ab>cd3_gate & !(UMAP_1>2.5 & UMAP_2<2.5) & (CD4_Ab<cd4_gate & CD8_Ab>cd8_gate), 'T:CD4-CD8+', truth)]
meta[, truth:=ifelse(CD3_Ab>cd3_gate & !(UMAP_1>2.5 & UMAP_2<2.5) & (CD4_Ab<=cd4_gate & CD8_Ab<=cd8_gate), 'T:CD4-CD8-', truth)]

meta[, truth:=ifelse((!(CD3_Ab>cd3_gate& !(UMAP_1>2.5 & UMAP_2<2.5) & UMAP_1<4)) & CD19_Ab<cd19_gate & CD16_Ab>=cd16_gate, 'NK', truth)]
meta[, truth:=ifelse((!(CD3_Ab>cd3_gate& !(UMAP_1>2.5 & UMAP_2<2.5) & UMAP_1<4)) & CD19_Ab<cd19_gate & CD16_Ab<cd16_gate, 'Myeloid_cells', truth)]
meta[, truth:=ifelse((!(CD3_Ab>cd3_gate& !(UMAP_1>2.5 & UMAP_2<2.5) & UMAP_1<4)) & CD19_Ab>=cd19_gate & CD16_Ab<cd16_gate, 'B', truth)]

meta <- meta[!is.na(truth), ]
meta[, .N, by = truth]


# Figure 6a overview of predicted cell types

meta2 <- copy(meta)
meta2[, ImmC:=ifelse(grepl("L:T:CD4", ImmC), 'T:CD4+CD8-', ImmC)]
meta2[, ImmC:=ifelse(grepl("L:T:CD8", ImmC), 'T:CD4-CD8+', ImmC)]
meta2[, ImmC:=ifelse(grepl("M:|DC", ImmC), 'Myeloid_cells', ImmC)]
meta2[, ImmC:=ifelse(grepl("NK", ImmC), 'NK', ImmC)]
meta2[, ImmC:=ifelse(grepl("B|PC", ImmC), 'B', ImmC)]
meta2[, .N, by = ImmC]

meta2[, SingleR:=ifelse(grepl("B_cell",SingleR), "B", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK",SingleR), "NK", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD4\\+",SingleR), "T:CD4+CD8-", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD8\\+",SingleR), "T:CD4-CD8+", SingleR)]
meta2[, SingleR:=ifelse(grepl("T_cell",SingleR), "T", SingleR)]
meta2[, SingleR:=ifelse(grepl("Monocyte",SingleR), "Myeloid_cells", SingleR)]
meta2[, SingleR:=ifelse(grepl("HSC",SingleR), "CD34", SingleR)]
meta2[, .N, by = SingleR]


# Fig 5a
seed<- 29
set.seed(seed)
meta3 <- melt(meta2[,  c('Cell', 'UMAP_1', 'UMAP_2', 'truth', 'ImmC', 'SingleR')], 
              id.vars = c('Cell', 'UMAP_1', 'UMAP_2'), 
              variable.factor = F, value.factor = F)
ggplot(meta3, 
       aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = factor(value, levels = sample(meta3[, .N, by = value][,value]))), 
             size = .1, alpha = 1) +
  theme_classic() +
  coord_equal() +
  facet_wrap(~factor(variable, levels = c('truth', 'ImmC', 'SingleR'))) +
  theme(axis.text.x=element_text(size = 12, color = 'black'),
        axis.title.x=element_text(size = 12, color = 'black'),
        axis.text.y=element_text(size = 12, color = 'black'),
        axis.title.y=element_text(size = 12, color = 'black'),
        strip.text = element_text(size = 24),
        legend.pos = 'bottom') + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_brewer(palette="Paired")






# Figure 6b T cells CD4 vs. CD8 distribution

meta4 <- melt(meta2[grepl('T', truth),  c('Cell', 'CD4_Ab', 'CD8_Ab', 'truth', 'ImmC', 'SingleR')], 
              id.vars = c('Cell', 'CD4_Ab', 'CD8_Ab', 'truth'), 
              variable.factor = F, value.factor = F)
set.seed(seed)

ggplot(meta4[variable == 'SingleR',], 
       aes(CD4_Ab, CD8_Ab)) + 
  geom_point(aes(color = factor(value, levels = sample(meta4[, .N, by = value][,value]))), 
             size = 1, alpha = .8) +
  theme_classic() +
  coord_equal() + 
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  facet_grid(factor(value, levels = meta4[, .N, by = value][order(-N), value])~factor(truth,levels = meta4[, .N, by = truth][order(-N), truth]))+
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.pos = 'none') + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_brewer(palette="Dark2")

meta5 <- meta2[grepl('T', truth),]
meta5[, CD4_RNA:=CD4>=1]
meta5[, CD8_RNA:=CD8A+CD8B>=1]



ggplot(meta5, aes(CD4_Ab, CD8_Ab)) + 
  geom_point(color = ifelse(meta5$CD4_RNA & meta5$CD8_RNA, 
                            'orange', 
                            ifelse(meta5$CD4_RNA & ! meta5$CD8_RNA, 
                                   'red', 
                                   ifelse(!meta5$CD4_RNA & meta5$CD8_RNA, 
                                          'blue', '#616161'))),
             size = 1.5, alpha =1) +
  theme_classic() +
  #facet_grid(paste(CD4_RNA, CD8_RNA)~truth)+
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        strip.text.x = element_text(size = 16),
        legend.pos = 'right') + coord_equal()



meta5[, list(CD4_Ab>cd4_gate&CD8_Ab<cd8_gate, CD4_RNA)][, .N, by = c('V1', 'CD4_RNA')]

ggplot(meta5, aes(CD4_Ab, CD8_Ab)) + 
  geom_point(color = ifelse(grepl("CD4\\+", meta5$ImmC) ,  'red', 
                            ifelse(grepl("CD8\\+", meta5$ImmC), 'blue',  'black')),
             size = 1.5, alpha =1) +
  theme_classic() +
  #facet_grid(paste(CD4_RNA, CD8_RNA)~truth)+
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        strip.text.x = element_text(size = 16),
        legend.pos = 'right') + coord_equal()
meta5[, list(CD4_Ab>cd4_gate&CD8_Ab<cd8_gate, grepl("CD4\\+", ImmC))][, .N, by = c('V1', 'V2')]
meta5[, list(CD8_Ab>cd8_gate&CD4_Ab<cd4_gate, grepl("CD8\\+", ImmC))][, .N, by = c('V1', 'V2')]



ggplot(meta5, aes(CD4_Ab, CD8_Ab)) + 
  geom_point(color = ifelse(grepl("CD4\\+", meta5$SingleR) ,  'red', 
                            ifelse(grepl("CD8\\+", meta5$SingleR), 'blue',  'black')),
             size = 1.5, alpha =1) +
  theme_classic() +
  #facet_grid(paste(CD4_RNA, CD8_RNA)~truth)+
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        strip.text.x = element_text(size = 16),
        legend.pos = 'right') + coord_equal()
meta5[, list(CD4_Ab>cd4_gate&CD8_Ab<cd8_gate, grepl("CD4\\+", SingleR))][, .N, by = c('V1', 'V2')]
meta5[, list(CD8_Ab>cd8_gate&CD4_Ab<cd4_gate, grepl("CD8\\+", SingleR))][, .N, by = c('V1', 'V2')]


# Figure 6d Recall & Prec heatmap
source('utility.R')


#meta2 <- copy(meta[truth == 'T:CD4+CD8-' |  truth == 'T:CD4-CD8+', ])
meta2 <- copy(meta)
meta2[, truth:=ifelse(grepl("CD4\\+CD8-", truth), "CD4+T", truth)]
meta2[, truth:=ifelse(grepl("CD4-CD8\\+", truth), "CD8+T", truth)]
meta2[, .N, by = truth]
# add t phenotype 


meta2[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t) & grepl("CD4\\+T", truth), truth:=paste0(truth, "naive")]
meta2[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t) & grepl("CD4\\+T", truth), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t) & grepl("CD4\\+T", truth), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t)& grepl("CD4\\+T", truth), truth:=paste0(truth, "EM")]



meta2[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t) & grepl("CD8\\+T", truth), truth:=paste0(truth, "naive")]
meta2[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t) & grepl("CD8\\+T", truth), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t) & grepl("CD8\\+T", truth), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t)& grepl("CD8\\+T", truth), truth:=paste0(truth, "EM")]





meta2[, ImmC:=gsub("L:T:", "", ImmC)]
meta2[, ImmC:=gsub(":CM", "+TCM", ImmC)]
meta2[, ImmC:=gsub(":EM$", "+TEM", ImmC)]
meta2[, ImmC:=gsub(":EMRA", "+TEMRA", ImmC)]
meta2[, ImmC:=gsub(":Na.*$", "+Tnaive", ImmC)]

meta2[, SingleR:=gsub("_[Cc]entral_memory", "TCM", SingleR)]
meta2[, SingleR:=gsub("T_cell:", "", SingleR)]

meta2[, SingleR:=gsub("_naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_Naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_effector_memory$", "TEM", SingleR)]
meta2[, SingleR:=gsub("_effector_memory_RA", "TEMRA", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK", SingleR), "NK", SingleR)]



p1 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, ImmC)][, c('Known', 'annot')]),
                   ord1=meta2[, .N, by  = truth][, truth],
                   ord2=meta2[, .N, by  = truth][, truth],# measures = 'f1',
                   minshow = 0, size2 =6, minpct = 0,
                   measures = c('Recall', 'Precision'))


p2 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, SingleR)][, c('Known', 'annot')]),
                   ord1=meta2[, .N, by  = truth][, truth],
                   ord2=meta2[, .N, by  = truth][, truth], #measures = 'f1',
                   minshow = 0, size2 =6, minpct = 0,
                   measures = c('Recall', 'Precision'))

grid.arrange(p1, p2, ncol = 2)


meta3 <- copy(meta)


meta3[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t) & grepl("T", truth), truth:="naive"]
meta3[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t) & grepl("T", truth), truth:"EMRA"]
meta3[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t)& grepl("T", truth), truth:="CM"]
meta3[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t)& grepl("T", truth), truth:="EM"]



meta3[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t) & grepl("T", truth), truth:=paste0(truth, "naive")]
meta3[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t) & grepl("T", truth), truth:=paste0(truth, "EMRA")]
meta3[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t)& grepl("T", truth), truth:=paste0(truth, "CM")]
meta3[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t)& grepl("T", truth), truth:=paste0(truth, "EM")]


meta3[, .N, by = truth]



# merge CD4 and CD8
p1 <- plot_heatmap(table(meta3[, c('Known', 'annot'):=list(truth, ImmC)][, c('Known', 'annot')]),
                   ord1=meta3[, .N, by  = truth][, truth],
                   ord2=meta3[, .N, by  = truth][, truth],# measures = 'f1',
                   minshow = 0, size2 =6, minpct = 0,
                   measures = c('Recall', 'Precision'))


p2 <- plot_heatmap(table(meta3[, c('Known', 'annot'):=list(truth, SingleR)][, c('Known', 'annot')]),
                   ord1=meta3[, .N, by  = truth][, truth],
                   ord2=meta3[, .N, by  = truth][, truth], #measures = 'f1',
                   minshow = 0, size2 =6, minpct = 0,
                   measures = c('Recall', 'Precision'))

grid.arrange(p1, p2, ncol = 2)


t_celltypes <- c('CD4+Tnaive', 'CD4+TCM','CD4+TEM', 'CD4+TEMRA',
                 'CD8+Tnaive', 'CD8+TCM','CD8+TEM', 'CD8+TEMRA'
                 )

t_celltypes2 <- c('naive', 'CM','EM', 'EMRA')

immc_recall <- c(85,63,10,0,79,10,12,85)
immc_prec <- c(48,62,48,0,57,12,41,15)

immc2_recall <- c(86,60,11,82)
immc2_prec <- c(51,60,44,15)


singler_recall <- c(90,79,31,0,0,0,0,1)
singler_prec <- c(65,63,43,0,0,0,0,7)

singler2_recall <- c(85,75,20,1)
singler2_prec <- c(75,66,56,7)



ggplot(melt( data.table(celltype = t_celltypes,
                        immc_rec=immc_recall, 
                        singler_rec=singler_recall,
                        immc_prec=immc_prec, 
                        singler_prec=singler_prec), id.vars = c('celltype')),
       aes(celltype, variable)) + geom_tile(aes(fill = value), color = 'black') +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_classic() + 
  theme(axis.text.x= element_text(size = 16, color = 'black', angle=45, hjust =1 ),
        axis.text.y= element_text(size = 24, color = 'black') ) + coord_equal() 




ggplot(melt( data.table(celltype = t_celltypes2,
                        immc_rec=immc2_recall, 
                        singler_recall=singler2_recall,
                        immc_prec=immc2_prec, 
                        singler_prec=singler2_prec), id.vars = c('celltype')),
       aes(celltype, variable)) + geom_tile(aes(fill = value), color = 'black') +
  scale_fill_gradient2(low = "red", mid = "white", high = "springgreen4") +
  theme_classic() + 
  theme(axis.text.x= element_text(size = 16, color = 'black', angle=45, hjust =1 ),
        axis.text.y= element_text(size = 24, color = 'black') ) + coord_equal() 




# Figure 6c CD62L vs. CD45RA ImmC breakdown

ggplot(meta2[grepl("CD4\\+T|CD8\\+T", truth),], aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(size = .1) +
  theme_classic() + 
  geom_density_2d(lwd=.1) +
  facet_wrap(~grepl("CD4", truth)) + 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .1) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .1) + xlim(0, 6) + ylim(0,7) +
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .1 ) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .1)



meta2[,is_facet:=ifelse(grepl("T", ImmC), ImmC, 'other')]

ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(color = ifelse(grepl("CD4\\+T", meta2$truth), 'red', ifelse(grepl("CD8\\+T", meta2$truth), 'blue', ifelse(grepl("CD4\\+CD8\\+", meta2$truth), 'gold', ifelse(grepl("CD4-CD8-", meta2$truth), 'cyan', 'darkgrey')))), size = .1) +
  theme_classic() + facet_wrap(~factor(ImmC, levels = meta2[,.N, by = is_facet][order(-N), is_facet]), ncol = 4)+ 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .3) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .3) +
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .3) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .3)


meta2[ImmC=='CD4+TCM', list(ImmC, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+TEM', list(ImmC, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+Tnaive', list(ImmC, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+TEMRA', list(ImmC, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]





meta2[ImmC=='CD8+TCM', list(ImmC, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+TEM', list(ImmC, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+Tnaive', list(ImmC, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+TEMRA', list(ImmC, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]


meta2[,is_facet:=ifelse(grepl("T|CD8", SingleR), SingleR, 'other')]
ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  #geom_point(aes(color = gsub("\\+.*$","", truth)), size = .1) +
  #geom_point(color = ifelse(grepl("CD4", meta2$truth), 'red', 'blue'), size = .1) +
  geom_point(color = ifelse(grepl("CD4\\+T", meta2$truth), 'red', ifelse(grepl("CD8\\+T", meta2$truth), 'blue', ifelse(grepl("CD4\\+CD8\\+", meta2$truth), 'gold', ifelse(grepl("CD4-CD8-", meta2$truth), 'cyan', 'darkgrey')))), size = .1) +
  theme_classic() + facet_wrap(~factor(SingleR, levels = meta2[,.N, by = is_facet][order(-N), is_facet]), ncol = 4)+ 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .3) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .3) + 
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .3 ) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .3)


meta2[SingleR=='CD4+TCM', list(SingleR, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD4+TEM', list(SingleR, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD4+Tnaive', list(SingleR, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]





meta2[SingleR=='CD8+_Central_memory', list(SingleR, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+TEM', list(SingleR, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+Tnaive', list(SingleR, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+TEMRA', list(SingleR, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
