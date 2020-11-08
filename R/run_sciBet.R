# gigantion
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))

model <- readr::read_csv("sciBet/major_human_cell_types.csv")
proc_model <- dcast(melt(model, id.var = 'X1'), variable~X1)
rownames(proc_model) <- proc_model$variable
prd <- LoadModel_R(proc_model[, -1])



# bulk
exp <- fread('sciBet/scibet_input/bulk.logrma.txt')
exp_data <- dcast(melt(exp, id.var = 'Gene'), variable~Gene)

bulk_scibet <- data.table(Cell=exp_data[, variable], 
           scibet=prd(2^(exp_data[, -1])))

fwrite(bulk_scibet, 'sciBet/scibet_output/bulk_scibet.txt')

#skcm
exp <- readRDS('sciBet/scibet_input/skcm.nh.ltpm.rds')
skcm_scibet <- data.table(Cell=rownames(exp), 
                          scibet=prd(2^(exp)))
fwrite(skcm_scibet, 'sciBet/scibet_output/skcm_scibet.txt')

#hnscc

exp <- readRDS('sciBet/scibet_input/hnscc.ltpm.rds')
hnscc_scibet <- data.table(Cell=rownames(exp), 
                          scibet=prd(2^(exp)))
fwrite(hnscc_scibet, 'sciBet/scibet_output/hnscc_scibet.txt')


#brca3p

exp <- readRDS('sciBet/scibet_input/brca.3p.lcpm.rds')
brca3p_scibet <- data.table(Cell=rownames(exp), 
                           scibet=prd(2^(exp)))
fwrite(brca3p_scibet, 'sciBet/scibet_output/brca3p_scibet.txt')

#brca5p

exp <- readRDS('sciBet/scibet_input/brca.5p.lcpm.rds')
brca5p_scibet <- data.table(Cell=rownames(exp), 
                            scibet=prd(2^(exp)))
fwrite(brca5p_scibet, 'sciBet/scibet_output/brca5p_scibet.txt')

#pbmc68k

exp <- readRDS('sciBet/scibet_input/pbmc68k.lcpm.rds')
pbmc68k_scibet <- data.table(Cell=rownames(exp), 
                            scibet=prd(2^(exp)))
fwrite(pbmc68k_scibet, 'sciBet/scibet_output/pbmc68k_scibet.txt')


#hcc

exp <- readRDS('sciBet/scibet_input/hcc.lcpm.rds')
hcc_scibet <- data.table(Cell=rownames(exp), 
                             scibet=prd(2^(exp)))
fwrite(hcc_scibet, 'sciBet/scibet_output/hcc_scibet.txt')


# pont
exp <- fread('sciBet/scibet_input/immc_in.txt')
exp_data <- dcast(melt(exp, id.var = 'Gene'), variable~Gene)

pont_scibet <- data.table(Cell=exp_data[, variable], 
                          scibet=prd(2^(exp_data[, -1])))

fwrite(pont_scibet, 'sciBet/scibet_output/pont_scibet.txt')




