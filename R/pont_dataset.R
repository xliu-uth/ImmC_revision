library(data.table)
library(ggplot2)

library(RColorBrewer)
library(ggalluvial)

adt <- fread('datasets/pont/GSM4288255_Antibodies_count.csv')
count <- fread('datasets/pont/GSM4288255_Raw_UMI_genes.tsv')
umap <- fread('datasets/pont/GSM4288255_UMAP_coordinate.csv')
hto <- fread('datasets/pont/GSM4288255_HTO_identification.csv')


####################################
# run single-R
####################################

.libPaths("/data/xuanliu/R-3.6.0/library")
library(data.table)
library(SingleR)


prefix <- 'pont'

count <- fread(paste0(prefix, '_singler_in.txt'))
df <- data.frame(count, stringsAsFactors = F)
colnames(df)[1] <- 'Gene.Symbol'
df <- df[!duplicated(df$Gene.Symbol), ]
rownames(df) <- df$Gene.Symbol


singler = CreateSinglerSeuratObject(df[,-1], project.name=prefix,
                                    min.genes = 500, technology="10X", species = "Human", 
                                    normalize.gene.length = F, min.cells = 10, npca = 10,
                                    regress.out = "nUMI", reduce.seurat.object = T)

# The object can then be saved and uploaded to the SingleR web-app for further analysis and visualization or using functions available in the SingleR package (see vignette).
save(singler,file=paste0(prefix,'.RData'))





####################################
# Compare
####################################

# 5559 T cells
singler_annot <- fread('datasets/pont/comparison/pont_singler_annot.txt')
colnames(singler_annot)[1] <- 'SingleR'
immc_annot <- fread('datasets/pont/comparison/pont.output.txt')
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
meta2[, ImmC:=ifelse(grepl("^M|DC", ImmC), 'Myeloid_cells', ImmC)]
meta2[, ImmC:=ifelse(grepl("NK", ImmC), 'NK', ImmC)]
meta2[, ImmC:=ifelse(grepl("B", ImmC), 'B', ImmC)]


meta2[, SingleR:=ifelse(grepl("B_cell",SingleR), "B", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK",SingleR), "NK", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD4\\+",SingleR), "T:CD4+CD8-", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD8\\+",SingleR), "T:CD4-CD8+", SingleR)]
meta2[, SingleR:=ifelse(grepl("T_cell",SingleR), "T", SingleR)]
meta2[, SingleR:=ifelse(grepl("Monocyte",SingleR), "Myeloid_cells", SingleR)]
meta2[, SingleR:=ifelse(grepl("HSC",SingleR), "CD34", SingleR)]


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

#meta2 <- copy(meta[truth == 'CD4+T' |  truth == 'CD8+T', ])
meta2 <- copy(meta[truth == 'T:CD4+CD8-' |  truth == 'T:CD4-CD8+', ])
meta2[, truth:=ifelse(grepl("CD4\\+", truth), "CD4+T", truth)]
meta2[, truth:=ifelse(grepl("CD8\\+", truth), "CD8+T", truth)]
meta2[, .N, by = truth]
# add t phenotype 


meta2[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t) & grepl("CD4", truth), truth:=paste0(truth, "naive")]
meta2[CD45RA_Ab>(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t) & grepl("CD4", truth), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab>(cd62_gate_cd4t) & grepl("CD4", truth), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd4t) &  CD62L_Ab<=(cd62_gate_cd4t)& grepl("CD4", truth), truth:=paste0(truth, "EM")]



meta2[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t) & grepl("CD8", truth), truth:=paste0(truth, "naive")]
meta2[CD45RA_Ab>(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t) & grepl("CD8", truth), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab>(cd62_gate_cd8t) & grepl("CD8", truth), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate_cd8t) &  CD62L_Ab<=(cd62_gate_cd8t)& grepl("CD8", truth), truth:=paste0(truth, "EM")]





meta2[, ImmC:=gsub("L:T:", "", ImmC)]
meta2[, ImmC:=gsub(":CM", "+TCM", ImmC)]
meta2[, ImmC:=gsub(":EM$", "+TEM", ImmC)]
meta2[, ImmC:=gsub(":EMRA", "+TEMRA", ImmC)]
meta2[, ImmC:=gsub(":Na.*$", "+Tnaive", ImmC)]


meta2[, SingleR:=gsub("T_cell:", "", SingleR)]
meta2[, SingleR:=gsub("_central_memory", "TCM", SingleR)]
meta2[, SingleR:=gsub("_naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_Naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_effector_memory$", "TEM", SingleR)]
meta2[, SingleR:=gsub("_effector_memory_RA", "TEMRA", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK", SingleR), "NK", SingleR)]



p1 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, ImmC)][, c('Known', 'annot')]),
             ord1=meta2[, .N, by  = truth][, truth],
             ord2=meta2[, .N, by  = truth][, truth],# measures = 'f1',
             minshow = 0, size2 =6, minpct = 0)
             

p2 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, SingleR)][, c('Known', 'annot')]),
             ord1=meta2[, .N, by  = truth][, truth],
             ord2=meta2[, .N, by  = truth][, truth], #measures = 'f1',
             minshow = 0, size2 =6, minpct = 0)

grid.arrange(p1, p2, ncol = 2)


immc_f1 <- c(57, 79, 28, 59, 0, 67, 14 , 0)
singler_f1 <- c(73, 0, 39, 77, 1, 0, 0, 0)

immc_recall <- c(52, 87, 18, 83, 0, 84, 14, 0)
immc_prec <- c(63, 73, 57, 46, 0, 56, 15, 0)


singler_recall <- c(79, 0, 31, 90, 1, 0, 0, 0)
singler_prec <- c(67, 0, 50, 66, 8, 0, 0, 0)


ggplot(melt( data.table(celltype = meta2[, .N, by  = truth][, truth],
                        immc_rec=immc_recall, 
                        singler_rec=singler_recall,
                        immc_prec=immc_prec, 
                        singler_prec=singler_prec), id.vars = c('celltype')),
       aes(celltype, variable)) + geom_tile(aes(fill = value), color = 'black') +
  scale_fill_gradient2(low = "red", mid = "white", high = "darkblue") +
  theme_classic() + 
  theme(axis.text.x= element_text(size = 16, color = 'black', angle=45, hjust =1 ),
        axis.text.y= element_text(size = 24, color = 'black') ) + coord_equal()
  




# Figure 6c CD62L vs. CD45RA ImmC breakdown

ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(size = .1) +
  theme_classic() + 
  geom_density_2d(lwd=.1) +
  facet_wrap(~grepl("CD4", truth)) + 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .1) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .1) + xlim(0, 6) + ylim(0,7) +
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .1 ) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .1)





ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(aes(color = ImmC), size = .1) +
  theme_classic() + facet_wrap(~factor(ImmC, levels = meta2[,.N, by = ImmC][order(-N), ImmC]), ncol = 4)+ 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .1) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .1) + xlim(0, 6) + ylim(0,7) +
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .1 ) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .1)


meta2[ImmC=='CD4+TCM', list(ImmC, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+TEM', list(ImmC, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+Tnaive', list(ImmC, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD4+TEMRA', list(ImmC, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]





meta2[ImmC=='CD8+TCM', list(ImmC, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+TEM', list(ImmC, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+Tnaive', list(ImmC, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[ImmC=='CD8+TEMRA', list(ImmC, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]



ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(aes(color = SingleR), size = .1) +
  theme_classic() + facet_wrap(~factor(SingleR, levels = meta2[,.N, by = SingleR][order(-N), SingleR]), ncol = 4)+ 
  geom_vline(xintercept = cd62_gate_cd4t,  linetype="dashed", color = "red", size = .1) + 
  geom_hline(yintercept = cd45ra_gate_cd4t,  linetype="dashed", color = "red", size = .1) + xlim(0, 6) + ylim(0,7) +
  geom_vline(xintercept = cd62_gate_cd8t,  linetype="dashed", color = "blue", size = .1 ) + 
  geom_hline(yintercept = cd45ra_gate_cd8t,  linetype="dashed", color = "blue" , size = .1)


meta2[SingleR=='CD4+TCM', list(SingleR, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD4+TEM', list(SingleR, CD62L_Ab<=cd62_gate_cd4t & CD45RA_Ab<=cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD4+Tnaive', list(SingleR, CD62L_Ab>cd62_gate_cd4t & CD45RA_Ab>cd45ra_gate_cd4t)][, .N, by = V2][, list(V2,N/sum(N))]





meta2[SingleR=='CD8+_Central_memory', list(SingleR, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+TEM', list(SingleR, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab<=cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+Tnaive', list(SingleR, CD62L_Ab>cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
meta2[SingleR=='CD8+TEMRA', list(SingleR, CD62L_Ab<=cd62_gate_cd8t & CD45RA_Ab>cd45ra_gate_cd8t)][, .N, by = V2][, list(V2,N/sum(N))]
