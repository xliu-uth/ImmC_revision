
data_folder <- 'C:\\Users/huang/Dropbox/AXL/Projects/ImmClassifier/Manuscript/data_for_figure_plots/'



celltype_match <- function(dt, from_celltype, to_celltype){

}




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
             minshow = 19, size =10, size2= 3, minpct=20, legend.pos = "none",
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
             minshow = 19, size =10, size2= 3, minpct=20, legend.pos = "none",
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

                                     

#########################################################################
#                               Figure 4                                #
#########################################################################

data_folder2 <- 'C:\\Users/huang/Dropbox/AXL/Projects/ImmClassifier/Manuscript/data_for_figure_plots/'
brca3p_meta <- data.table(readRDS(paste(data_folder, 'metadata/BRCA-DP-meta.rds', sep = "/")))
brca3p_clsnames <- data.table(readRDS(paste(data_folder, 'metadata/BRCA-DP-clusternames.rds', sep = "/")))
setkey(brca3p_clsnames, ClusterID)
setkey(brca3p_meta, cellid)

# load prediction results from ImmClassifier, SingleR and Garnett, baseline
brca3p_immc <- fread(paste(data_folder,'brca3p.output.txt', sep = "/"))
brca3p_dnn_stat <- fread('tensorflow/output/brca3p.deeplearning.ontotree.stats.txt', header = F)[, 1:39]
cnames <- c("Cell",ref_nodes[!grepl("X", ref_nodes)])
colnames(brca3p_dnn_stat) <- cnames
brca3p_dnn_stat[, Known:=brca3p_clsnames[brca3p_meta[Cell, cluster], Final_annotation]]


ggplot(melt(brca3p_dnn_stat, id.vars = c('Cell', 'Known')), aes(variable, value)) + 
  geom_boxplot() + facet_wrap(~Known) +
  theme(axis.text.x=element_text(angle = 90, size = 8))

brca3p_immc[, Known:=brca3p_clsnames[brca3p_meta[Cell, cluster], Final_annotation]]
brca3p_immc[, annot:=ImmClassifier_prediction]
table(brca3p_immc[, c('Known', 'annot')])  %>%
  plot_heatmap(title = "", palette= "RdPu", 
               minshow =5, size =16, size2= 5, minpct=0, measures = c('Recall', 'Precision'),
               ord1 = brca3p_immc[, .N, by = annot][, annot],
               ord2 =brca3p_immc[, .N, by = Known][, Known] )            

#########################################################################
#                               Figure 5                                #
#########################################################################
