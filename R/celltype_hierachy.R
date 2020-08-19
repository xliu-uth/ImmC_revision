library(data.table)

dnn_train_input <- data.table(readRDS('tensorflow/input/deeplearn.train.rds'))
dnn_train_input[grepl("Ref6.CD8..*[^M]$", clusterID), clusterID:="Ref6.CD8..gd"]

cell_hierarchy <- fread('annotation_data/dataset.txt')
setkey(cell_hierarchy, ClusterID)

dnn_train_input[, Hierarchy:=cell_hierarchy[clusterID, Hierarchy_version08112020]]
 

cluster_size <- dnn_train_input[, .N, by = Hierarchy]
cluster_size[, size:=ifelse(N<1000, 1000, ifelse(N<2000, 1200, ifelse(N<5000, 1400, 1600)))]
setkey(cluster_size, Hierarchy)


dnn_train_input_balanced <- dnn_train_input[Hierarchy!='X',.SD[sample(.N, cluster_size[Hierarchy, size], replace =T)],by = Hierarchy]

ref_nodes <-  unique(cover_set(cell_hierarchy$Hierarchy_version08112020))
ref_nodes <- ref_nodes[ref_nodes!='X']
ref_nodes
dnn_train_input_balanced[, Evopath:=sapply(dnn_train_input_balanced$Hierarchy, function(x) convert_to_bits(x, ref_nodes[ref_nodes!='X']))]
fwrite(data.table(dnn_train_input_balanced[, -1], Hierarchy=dnn_train_input_balanced$Hierarchy), 'tensorflow/input/deeplearn.train.balance.input08112020.txt')

