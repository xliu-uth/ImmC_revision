


data_folder <- 'C:\\Users/huang/Dropbox/AXL/Projects/ImmClassifier/Manuscript/data_for_figure_plots/'




#########################################################################
#                               Figure 1                                #
#########################################################################


#########################################################################
#                               Figure 2                                #
#########################################################################


# Figure 2A

bulk_immc <- fread(paste(data_folder, 'immc/bulk.output.txt', sep = "/"))

rename_pairs <- list(from=c(),
                     to= c())



bulk_immc[, Known:=gsub("[0-9]*_[0-9]+$", "", Cell)]
bulk_immc[grepl("TCELLA[1-4]", Cell), Known:='CD8']
bulk_immc[grepl("TCELLA[6-8]", Cell), Known:='CD4']
bulk_immc[grepl("ERY1|GRAN1|PRE_BCELL2|GMP|CMP|MEP|MEGA1|HSC", Cell), Known:='CD34']
bulk_immc[grepl("DENDA1", Cell), Known:='pDC']
bulk_immc[grepl("DENDA2", Cell), Known:='mDC']
bulk_immc[grepl("NKA4", Cell), Known:='NKT']
bulk_immc[grepl("BCELL", Cell), Known:='BCELL']

bulk_immc[, annot:=ImmClassifier_prediction]
bulk_immc[grepl('CD4', ImmClassifier_prediction), annot:='CD4']
bulk_immc[grepl('CD8', ImmClassifier_prediction), annot:='CD8']
bulk_immc[grepl('mDC', ImmClassifier_prediction), annot:='mDC']
bulk_immc[grepl('pDC', ImmClassifier_prediction), annot:='pDC']
bulk_immc[grepl('L:B', ImmClassifier_prediction), annot:='B']
bulk_immc[grepl('Mega', ImmClassifier_prediction), annot:='Mega']
bulk_immc[grepl('Mono', ImmClassifier_prediction), annot:='Mono']
bulk_immc[grepl('Neu', ImmClassifier_prediction), annot:='Neu']
bulk_immc[grepl('NK', ImmClassifier_prediction), annot:='NK']
bulk_immc[grepl('Ery', ImmClassifier_prediction), annot:='Ery']
bulk_immc[grepl('Eos', ImmClassifier_prediction), annot:='Eos']
bulk_immc[grepl('CD34', ImmClassifier_prediction), annot:='CD34']
bulk_immc[, .N, by = c('Known', 'annot')]




plot_heatmap(table(bulk_immc[, c(c('Known', 'annot'))]), title = "", palette= "RdPu", 
             minshow = 19, size =10, size2= 3, minpct=0, legend.pos = "none",
             ord1 = bulk_immc[, .N, by = annot][, annot],
             ord2 =bulk_immc[, .N, by = Known][, Known], 
             measures = c('Recall', 'Precision') )  
             
             

#########################################################################
#                               Figure 3                                #
#########################################################################
skcm_immc <- fread(paste(data_folder, 'immc/skcm.output.txt', sep = "/"))
skcm_meta <- readRDS(paste(data_folder,"metadata/SKCM-NH-meta.rds", sep = "/"))
skcm_cluster_names <- readRDS(paste(data_folder,"metadata/SKCM-NH-clusternames.rds", sep = "/"))


skcm_immc[, Known:=skcm_cluster_names[skcm_meta[Cell, 'Cluster'], 'clusterName']]
skcm_immc[, annot:=ImmClassifier_prediction]
table(skcm_immc[, c('Known', 'annot')])  %>%
    plot_heatmap(title = "", palette= "RdPu", 
                 minshow =5, size =16, size2= 5, minpct=0, measures = c('Recall', 'Precision'),
                 ord1 = skcm_immc[, .N, by = annot][, annot],
                 ord2 =skcm_immc[, .N, by = Known][, Known], )                                                      
                                                

hnscc_meta <- readRDS(paste(data_folder,'metadata/HNSCC.meta.rds', sep = "/"))
hnscc_immc <- fread(paste(data_folder,'immc/hnscc.output.txt', sep = "/"))

hnscc_immc[, Known:=hnscc_meta[Cell,'non_cancer_cell_type']]
hnscc_immc[, annot:=ImmClassifier_prediction]
table(hnscc_immc[, c('Known', 'annot')])  %>%
  plot_heatmap(title = "", palette= "RdPu", 
               minshow =5, size =16, size2= 5, minpct=0, measures = c('Recall', 'Precision'),
               ord1 = hnscc_immc[, .N, by = annot][, annot],
               ord2 =hnscc_immc[, .N, by = Known][, Known], )                                                      

#########################################################################
#                               Figure 4                                #
#########################################################################


#########################################################################
#                               Figure 5                                #
#########################################################################
