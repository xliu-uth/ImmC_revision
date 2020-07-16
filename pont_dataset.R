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
cd62_gate <- 2.8 
cd45ra_gate <- 2.6
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

# Figure 1 overview of predicted cell types

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
meta2[, SingleR:=ifelse(grepl("HSC",SingleR), "CD34+", SingleR)]


seed<- 29
set.seed(seed)
meta3 <- melt(meta2[,  c('Cell', 'UMAP_1', 'UMAP_2', 'truth', 'ImmC', 'SingleR')], 
              id.vars = c('Cell', 'UMAP_1', 'UMAP_2'), 
              variable.factor = F, value.factor = F)
ggplot(meta3, 
       aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color = factor(value, levels = sample(meta3[, .N, by = value][,value]))), 
             size = 1, alpha = .8) +
  theme_classic() +
  coord_equal() +
  facet_wrap(~factor(variable, levels = c('truth', 'ImmC', 'SingleR'))) +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        legend.pos = 'bottom') + 
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_brewer(palette="Paired")



# Figure 2
# CD4 CD8


ggplot(meta2[grepl('T', truth),], 
       aes(CD4_Ab, CD8_Ab)) + 
  geom_point(color = 'black', alpha = .6, size = 1) +
  theme_classic() +
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        legend.pos = 'none') 
  

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
                            'purple', 
                            ifelse(meta5$CD4_RNA & ! meta5$CD8_RNA, 
                                   'red', 
                                   ifelse(!meta5$CD4_RNA & meta5$CD8_RNA, 
                                          'blue', '#616161'))),
             size = 1, alpha = .4) +
  theme_classic() +
  facet_grid(paste(CD4_RNA, CD8_RNA)~truth)+
  geom_vline(xintercept = cd4_gate,  linetype="dashed", color = "red" ) + 
  geom_hline(yintercept = cd8_gate, linetype="dashed", color = "red") +
  theme(axis.text.x=element_text(size = 20, color = 'black'),
        axis.title.x=element_text(size = 24, color = 'black'),
        axis.text.y=element_text(size = 20, color = 'black'),
        axis.title.y=element_text(size = 24, color = 'black'),
        strip.text.x = element_text(size = 16),
        legend.pos = 'right') + coord_equal()


###############################################

# Figure 3 Alluvial of broad celltypes



meta2 <- copy(meta)
meta2[, ImmC:=ifelse(grepl("^M:|DC",ImmC), 'Myeloid_cells',ImmC)]
meta2[, ImmC:=ifelse(grepl("CD34",ImmC), "CD34",  ImmC)]
meta2[, ImmC:=ifelse(grepl("CD4",ImmC), "CD4+", ifelse(grepl("CD8",ImmC), "CD8+", ImmC))]
meta2[, ImmC:=ifelse(grepl("B",ImmC), "B", ImmC)]
meta2[, ImmC:=ifelse(grepl("NK",ImmC), "NK", ImmC)]

meta2[, SingleR:=ifelse(grepl("B_cell:",SingleR), "B", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK",SingleR), "NK", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD4+",SingleR), "CD4+", ifelse(grepl("CD8+",SingleR), "CD8+", SingleR))]
meta2[, SingleR:=ifelse(grepl("Monocyte",SingleR), "Myeloid_cells", SingleR)]





ggplot(meta2, aes(CD3_Ab, CD14_Ab)) + geom_point(size = .1) + facet_wrap(~ImmC)


ggplot(meta2[grepl('CD4', truth), .N, by = c('ImmC',
                                             'truth', 
                                             'SingleR')][, cus_col:= paste0(ImmC, truth, SingleR)], 
       aes(axis1 = ImmC, axis2 = truth, axis3 = SingleR, y = N)) +
  scale_x_discrete(limits = c('ImmC', 'truth', 'SingleR'), 
                   expand = c(.1, 0)) +
  xlab("Method") +
  scale_fill_manual(values =colorRampPalette(brewer.pal(8, "Dark2"))(9))+
  geom_alluvium(aes(fill = truth))+
  geom_stratum(width = .3) + 
  geom_text(stat = "stratum", infer.label = T, size = 4)+
  theme_minimal() +
  theme(legend.pos = "none") 









ggplot(meta2, aes(CD62L_Ab, CD45RA_Ab)) + 
  geom_point(aes(color = ImmC), size = .1) +
  theme_classic() + facet_wrap(~ImmC)+ 
  geom_vline(xintercept = cd62_gate) + 
  geom_hline(yintercept = cd45ra_gate)



ggplot(meta2, aes(CD4_Ab, CD8_Ab)) + 
  geom_point(aes(color = ImmC), size = .1) +
  theme_classic() + facet_wrap(~factor(ImmC, levels = meta2[, .N, by = ImmC][order(-N),ImmC]))+ 
  geom_vline(xintercept = cd4_gate) + 
  geom_hline(yintercept = cd8_gate)



meta2[, SingleR:=ifelse(grepl("B",SingleR), "B", SingleR)]
meta2[, SingleR:=ifelse(grepl("NK",SingleR), "NK", SingleR)]
meta2[, SingleR:=ifelse(grepl("CD4+",SingleR), "CD4+", ifelse(grepl("CD8+",SingleR), "CD8+", SingleR))]


# cd4 cd8 distinction

ggplot(meta2[, .N, by = c('ImmC',
                          'truth', 
                          'SingleR')][, cus_col:= paste0(ImmC, truth, SingleR)], 
       aes(axis1 = ImmC, axis2 = truth, axis3 = SingleR, y = N)) +
  scale_x_discrete(limits = c('ImmC', 'truth', 'SingleR'), 
                   expand = c(.1, 0)) +
  xlab("Method") +
  scale_fill_manual(values =colorRampPalette(brewer.pal(8, "Dark2"))(8))+
  geom_alluvium(aes(fill = truth))+
  geom_stratum(width = .3) + 
  geom_text(stat = "stratum", infer.label = T, size = 6)+
  theme_minimal() +
theme(legend.pos = "none") 


ggplot(meta2, aes(CD4_Ab, CD8_Ab)) + geom_point(aes(color = ImmC), size = .1) +
  theme_classic() + facet_wrap(~ImmC) + 
  geom_vline(xintercept = cd4_gate) + 
  geom_hline(yintercept = cd8_gate)






# CD4+T 
meta2 <- meta[truth == 'CD4+T',]
meta2[, ImmC:=ifelse(grepl("CD8+",ImmC), "CD8+",  ImmC)]
meta2[, ImmC:=ifelse(grepl("B",ImmC), "B", ImmC)]
meta2[, ImmC:=ifelse(grepl("NK",ImmC), "NK", ImmC)]
meta2[, ImmC:=gsub("L:T:", "", ImmC)]
meta2[, SingleR:=gsub("T_cell:", "", SingleR)]


meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "N")]
meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EM")]


ggplot(meta2[grepl('CD4', truth), .N, by = c('ImmC',
                          'truth', 
                          'SingleR')][, cus_col:= paste0(ImmC, truth, SingleR)], 
       aes(axis1 = ImmC, axis2 = truth, axis3 = SingleR, y = N)) +
  scale_x_discrete(limits = c('ImmC', 'truth', 'SingleR'), 
                   expand = c(.1, 0)) +
  xlab("Method") +
  scale_fill_manual(values =colorRampPalette(brewer.pal(8, "Dark2"))(9))+
  geom_alluvium(aes(fill = truth))+
  geom_stratum(width = .3) + 
  geom_text(stat = "stratum", infer.label = T, size = 4)+
  theme_minimal() +
theme(legend.pos = "none") 




# CD8+T 
meta2 <- meta[truth == 'CD8+T' , ]
meta2[, ImmC:=ifelse(grepl("CD4+",ImmC), "CD4+",  ImmC)]
meta2[, ImmC:=ifelse(grepl("B",ImmC), "B", ImmC)]
meta2[, ImmC:=ifelse(grepl("NK",ImmC), "NK", ImmC)]
meta2[, SingleR:=ifelse(grepl("CD4+",SingleR), "CD4+",  SingleR)]
meta2[, ImmC:=gsub("L:T:", "", ImmC)]
meta2[, SingleR:=gsub("T_cell:", "", SingleR)]





meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "N")]
meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EM")]


ggplot(meta2[grepl('CD8', truth), .N, by = c('ImmC',
                          'truth', 
                          'SingleR')][, cus_col:= paste0(ImmC, truth, SingleR)], 
       aes(axis1 = ImmC, axis2 = truth, axis3 = SingleR, y = N)) +
  scale_x_discrete(limits = c('ImmC', 'truth', 'SingleR'), 
                   expand = c(.1, 0)) +
  xlab("Method") +
  scale_fill_manual(values =colorRampPalette(brewer.pal(8, "Dark2"))(8))+
  geom_alluvium(aes(fill = truth))+
  geom_stratum(width = .3) + 
  geom_text(stat = "stratum", infer.label = T, size = 4)+
  theme_minimal() +
theme(legend.pos = "none") 



source('utility.R')

meta2 <- copy(meta[truth == 'CD4+T' |  truth == 'CD8+T', ])

meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "naive")]
meta2[CD45RA_Ab>(cd45ra_gate+delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EMRA")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab>(cd62_gate+delta), truth:=paste0(truth, "CM")]
meta2[CD45RA_Ab<=(cd45ra_gate-delta) &  CD62L_Ab<=(cd62_gate-delta), truth:=paste0(truth, "EM")]


meta2[, ImmC:=gsub("L:T:", "", ImmC)]
meta2[, ImmC:=gsub(":CM", "+TCM", ImmC)]
meta2[, ImmC:=gsub(":EM", "+TEM", ImmC)]
meta2[, ImmC:=gsub(":Na.*$", "+Tnaive", ImmC)]


meta2[, SingleR:=gsub("T_cell:", "", SingleR)]
meta2[, SingleR:=gsub("_central_memory", "TCM", SingleR)]
meta2[, SingleR:=gsub("_naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_Naive", "Tnaive", SingleR)]
meta2[, SingleR:=gsub("_effector_memory$", "TEM", SingleR)]
meta2[, SingleR:=gsub("_effector_memory_RA", "TEMRA", SingleR)]



p1 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, ImmC)][, c('Known', 'annot')]),
             ord1=meta2[, .N, by  = truth][, truth],
             ord2=meta2[, .N, by  = truth][, truth],# measures = 'f1',
             minshow = 0, size2 =6, minpct = 0)
             

p2 <- plot_heatmap(table(meta2[, c('Known', 'annot'):=list(truth, SingleR)][, c('Known', 'annot')]),
             ord1=meta2[, .N, by  = truth][, truth],
             ord2=meta2[, .N, by  = truth][, truth], #measures = 'f1',
             minshow = 0, size2 =6, minpct = 0)

grid.arrange(p1, p2, ncol = 2)


immc_f1 <- c(58, 0, 50, 26, 42, 32, 6, 0)
singler_f1 <- c(72, 2, 0, 38, 76, 0, 0, 0)

immc_recall <- c(65, 0, 37, 17, 30, 31, 3 , 0)
immc_prec <- c(53, 0, 75, 52, 70, 33, 22, 0)


singler_recall <- c(79, 0, 0, 31, 90, 0, 0, 0)
singler_prec <- c(66, 46, 0, 41, 66, 0, 0, 0)

ggplot(data.table(celltype = meta2[, .N, by  = truth][, truth],
           immc=immc_recall, singler=singler_recall)[, delta:=immc-singler],
       aes(celltype,delta)) + geom_bar(stat="identity")


ggplot(melt(data.table(celltype = meta2[, .N, by  = truth][, truth],
                  immc=immc_recall, singler=singler_recall), id.vars = c('celltype')),
       aes(variable,value)) + 
  geom_boxplot(alpha = .8) + 
  geom_point(size = 4) + theme_classic()

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
  

