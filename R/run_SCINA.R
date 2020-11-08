# local
library(SCINA)

## Figure R1: SCINA predicts bulk dataset
test_path <- 'SCINA/input_matrix/bulk.logrma.txt'

exp <-  fread(test_path)
exp_data <- data.frame(exp, stringsAsFactors = F)
rownames(exp_data) <- exp[,Gene]
sig1 <- preprocess.signatures('SCINA/Fig1b_kidney_signatures.csv')
sig1_celltypes <- unique(gsub("[0-9]+$", "", 
                              names(unlist(sig1)[unlist(sig1) %in% exp_data$Gene])))

bulk_scina <- data.table(Cell = colnames(exp_data)[-1],
                         scina=SCINA(exp_data[, -1], sig1[sig1_celltypes], 
                                     max_iter = 100, convergence_n = 10, 
                                     convergence_rate = 0.999, 
                                     sensitivity_cutoff = 0.9, 
                                     rm_overlap=TRUE, allow_unknown=TRUE, 
                                     log_file='SCINA.log')$cell_labels)

fwrite(bulk_scina, 'SCINA/output/bulk_scina.txt')

bulk_scina[, Known:=gsub("_[0-9]+$", "", Cell)]
bulk_scina[grepl("ERY1|GRAN1|PRE_BCELL2|GMP|CMP|MEP|MEGA1|HSC", Cell), Known:="CD34+"]


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
  bulk_scina[Known==from_celltype[i], Known:=to_celltype[i]]  
}


bulk_scina[, annot:=scina]
bulk_scina[grepl("NK", annot), annot:='NK']

ord1 <- c('B.cells', 'Tm.cells','Th1.cells', 'Th2.cells', 'Tfh.cells',
          'NK', 'Dendritic.cells',
          'Macrophages', 'Neutrophils', 'Mast.cells')
bulk_scina[!annot %in% ord1, annot:="Other"]


plot_heatmap(table(bulk_scina[, c(c('Known', 'annot'))]), title = "", palette= "RdPu", 
             minshow = 1, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1 = ord1,
             ord2 = c('BCELL', 'CD4', 'CD8', 'NK', 'mDC', 'pDC', 
                      'MONO', 'GRAN', 'ERY', 'CD34+','MEGA2', 'NKA4', 
                      'PRE_BCELL3',  'EOS2', 'BASO1'),
             measures = c('Recall', 'Precision') )  



## Figure R2: SCINA predicts skcm dataset 
# preprocess on gigantion
library(SCINA)
library(data.table)

x <- readRDS('skcm.nh.ltpm.rds')
y <- data.table(Cell=rownames(x), x[, -1])
z <- dcast(melt(y, id.var = 'Cell'), variable~Cell)

fwrite(z, '../processed_input_matrix/skcm_ltpm_scina_input.txt')


test_path <- 'SCINA/processed_input_matrix/skcm_ltpm_scina_input.txt'
exp <-  fread(test_path)
exp_data <- data.frame(exp, stringsAsFactors = F)
rownames(exp_data) <- exp[,variable]

sig1 <- preprocess.signatures('SCINA/signature_matrix/Fig1b_kidney_signatures.csv')
sig1_celltypes <- unique(gsub("[0-9]+$", "", 
                              names(unlist(sig1)[unlist(sig1) %in% exp_data$variable])))


skcm_scina <- data.table(Cell = colnames(exp_data)[-1],
                         scina=SCINA(exp_data[, -1], sig1[sig1_celltypes], 
                                     max_iter = 100, convergence_n = 10, 
                                     convergence_rate = 0.999, 
                                     sensitivity_cutoff = 0.9, 
                                     rm_overlap=TRUE, allow_unknown=TRUE, 
                                     log_file='SCINA.log')$cell_labels)

fwrite(skcm_scina, 'SCINA/scina_output/skcm_scina.txt')
### local

data_folder <- '../../workspace/ImmClassifier_manuscript/data'
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



skcm_scina <- fread('SCINA/output/skcm_scina.txt')
skcm_scina[, Known:=skcm_cluster_names[skcm_meta[Cell, 'Cluster'], clusterName]]
skcm_scina[, annot:=scina]
skcm_scina <- skcm_scina[!is.na(Known),]


skcm_scina[grepl("ympho|Memory|Exh", Known), Known:="T cells"]
skcm_scina[grepl("B.cells", annot), annot:="B"]
skcm_scina[grepl("Tfh|Th|Tm|CD8", annot), annot:="T cells"]


ord1 <- c('B', 'T cells', 'Treg.cells', 'Macrophages', 'Dendritic.cells', 
          'CD56dim.NK.cells', 'Neutrophils', 'Mast.cells', 'unknown')

skcm_scina[!annot %in% ord1, annot:='Other']

plot_heatmap(table(skcm_scina[, c(c('Known', 'annot'))]), title = "", 
             minshow = 19, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1=ord1,
             ord2 = c('B', 'Plasma Cells', 'T cells',
                      'Treg', 'Monocytes/Macrophage', 'DC'),
             measures = c('Recall', 'Precision') )  


## Figure R3: SCINA predicts hnscc dataset
x <- readRDS('hnscc.ltpm.rds')
y <- data.table(Cell=rownames(x), x[, -1])
z <- dcast(melt(y, id.var = 'Cell'), variable~Cell)

fwrite(z, '../processed_input_matrix/hnscc.ltpm_scina_input.txt')


test_path <- 'SCINA/processed_input_matrix/hnscc.ltpm_scina_input.txt'
exp <-  fread(test_path)
exp_data <- data.frame(exp, stringsAsFactors = F)
rownames(exp_data) <- exp[,variable]


hnscc_scina <- data.table(Cell = colnames(exp_data)[-1],
                          scina=SCINA(exp_data[, -1], sig1[sig1_celltypes], 
                                      max_iter = 100, convergence_n = 10, 
                                      convergence_rate = 0.999, 
                                      sensitivity_cutoff = 0.9, 
                                      rm_overlap=TRUE, allow_unknown=TRUE, 
                                      log_file='SCINA.log')$cell_labels)

fwrite(skcm_scina, 'SCINA/scina_output/hnscc_scina.txt')

### local

hnscc_meta <- readRDS(paste(data_folder,'metadata/HNSCC.meta.rds', sep = "/"))
hnscc_scina <- fread('SCINA/output/hnscc_scina.txt')

hnscc_scina[, Known:=hnscc_meta[Cell,'non_cancer_cell_type']]
hnscc_scina[, annot:=scina]

hnscc_scina <- hnscc_scina[grepl("B|Den|Ma|T", Known), ]
hnscc_scina[grepl("B.cells", annot), annot := 'B cell']
hnscc_scina[grepl("^T|CD8", annot), annot := 'T cell']

ord1 <- c('T cell', 'Macrophages', 'B cell', 'Mast.cells','Dendritic.cells',
          'Neutrophils',  'CD56dim.NK.cells', 'unknown')

plot_heatmap(table(hnscc_scina[, c(c('Known', 'annot'))]), title = "", 
             minshow = 19, size =10, size2= 3, minpct=20, legend.pos = "none",
             ord1= ord1, 
             ord2 = hnscc_scina[, .N, by = Known][, Known],
             measures = c('Recall', 'Precision') )  




## Figure R4: SCINA predicts brca3p dataset

# unicron prepare for input matrix for SCINA
x <- readRDS('../../../immune_pred/BRCA-DP/processed_data/brca.3p.count.rds')

y <- data.table(x[, -c(1:3)])
z <- dcast(melt(y, id.var = 'cellid'), variable~cellid)
fwrite(z, '../SCINA_input/brca3p_count_scina_input.txt')

test_path <- 'SCINA/input_matrix/brca3p_count_scina_input.txt'
exp <-  fread(test_path)
exp_data <- data.frame(exp, stringsAsFactors = F)
rownames(exp_data) <- exp[,variable]
exp_raw=log(exp_data[, -1]+1)
exp_norm[]=normalize.quantiles(exp_raw)


brca3p_scina <- data.table(Cell = colnames(exp_data)[-1],
                           scina=SCINA(exp_data[, -1], sig1[sig1_celltypes], 
                                       max_iter = 100, convergence_n = 10, 
                                       convergence_rate = 0.999, 
                                       sensitivity_cutoff = 0.9, 
                                       rm_overlap=TRUE, allow_unknown=TRUE, 
                                       log_file='SCINA.log')$cell_labels)




