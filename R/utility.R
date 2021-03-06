library(dplyr)
library(ggplot2)
plot_heatmap <- function(data, title="",  ncol = 2, ord1 = NULL, ord2=NULL,
                         size = 10,size2=2, palette = 'RdPu',
                         minpct = 0, minshow=25, measures = c('f1', 'Precision', 'Recall'),
                         legend.pos='none') {
  
  TP <- data
  FN <- apply(data, 1, sum) - data
  FP <- t(apply(data, 2, sum) - t(data))
  TN <- sum(data) - TP -FN - FP
  
  accuracy <- (TP+TN)/(FN+FP+TP+TN)
  recall <- TP/(TP+FN)
  specificity <- TN/(TN+FP)
  precision <- TP/(TP+FP) # PPV
  
  f1 <- 2*(precision*recall)/(precision+recall)
  
  
  data <- rbind(data.frame(recall %>% melt, measure = 'Recall'),
                data.frame(precision %>% melt, measure = 'Precision'),
                data.frame(f1  %>% melt, measure = 'f1')
  ) %>% mutate(value = round(value,2)*100)
  
  
  if(!is.null(ord1)){
    data <- data %>%
      mutate(annot = as.character(annot)) %>%
      mutate(annot = ifelse(annot %in% ord1, annot, 'Other')) %>%
      mutate(annot = factor(annot, levels = c(ord1, 'Other')))
    
  }
  
  if(!is.null(ord2)){
    data <- data %>%
      mutate(Known = as.character(Known)) %>%
      mutate(Known = ifelse(Known %in% ord2, Known, 'Other')) %>%
      mutate(Known =  factor(Known, levels = c(ord2, 'Other')))
  }
  
  
  
  #colfunc <- colorRampPalette(c("white", "purple"))
  colfunc <- colorRampPalette(brewer.pal(9, palette))
  #zcut <- c(-1, 25, 50, 75, 100)
  zcut <- c(-1, 19, 40, 60, 80, 100)
  #zcut <- c(-1, 14, 30, 45, 60)
  
  data2 <- data %>% filter(measure %in% measures)
  data2[is.na(data2$value), 'value'] <- 0
  
  data2 %>%
    ggplot(aes(Known, annot,  fill = cut(value, zcut))) +
    geom_tile(color = "black") +
    #  scale_fill_brewer(palette = palette,  direction = 1) +
    theme_minimal()+
    theme(
      axis.text.x = element_text(angle = 45, size = size, hjust = 1,color= "black"),
      axis.text.y = element_text(size = size, hjust = 1,color = "black"),
      axis.title.x = element_text(size = size+4, hjust = 1),
      axis.title.y = element_text(size = size+4, hjust = 1),
      legend.text = element_text(size = size, angle = 0, hjust = 1),
      legend.title = element_text(size = size, hjust = 1),
      plot.title = element_text(size = 14, hjust = 1),
      strip.text.x = element_text(size = 14, colour = "black", angle = 0),
      legend.position=legend.pos) +
    coord_equal()+
    labs(title = title, x ="known", y = "predicted") +
    
    facet_wrap(~measure, ncol =  ncol)   +
    scale_fill_manual(values = c('white', colfunc(9)[c(3,5,7,9)])) +
    scale_colour_manual(values=c(rep("black", 3), rep("white",2))) +
    # scale_fill_manual(values = colfunc(20)) +
    geom_text(aes(label=ifelse(value>minshow, value, ""), color = cut(value, zcut)), size = size2) +
    guides(fill=guide_legend(title="percentage of row (recall), col(precision)"))
  
}




cover_set <- function(leaves){
  cset <- c()
  
  for(leaf in leaves){
    
    
    for (node in strsplit(leaf, ";")[[1]]){
      
      elements <- strsplit(node, ":")[[1]]
      i <- length(elements)
      
      current <- elements[1]
      
      cset <- c(cset, current)
      while (i > 1){
        current <- elements[1]
        for(j in 2:i){
          
          current <- paste(current, elements[j], sep = ":")
          
          cset <- c(cset, current)
        }
        i <- i-1
      }
      
      
      
      
    }
    
  }    
  return (cset)
}




# create child-parent-link

convert_to_bits <- function(nodes, ref.nodes){
  
  target.bits <- rep(0, length(ref.nodes))
  names(target.bits) <- ref.nodes
  
  
  for (node in strsplit(nodes, ";")[[1]]){
    
    elements <- strsplit(node, ":")[[1]]
    
    i <- 1
    current <- ""
    while(i <= length(elements)){
      
      if(i == 1){
        current <- paste0(current, elements[i])
      }else{
        current <- paste(current, elements[i], sep = ":")
      }
      target.bits[current] <- 1
      i <- i+1
    }
    
  }
  
  #print (paste0("return ",paste(target.bits, collapse = ",")))
  return(paste(target.bits, collapse = ","))
}




# match clusters between original and unsupervised umap

getNearestCentroid <- function(query, original.centroids){
  return(original.centroids[order(apply(original.centroids[, -1], 1, function(x) (x[1]-query[1])^2 + (x[2]-query[2])^2))[1], original_annotation])
}




dist_plot <- function(dt, palette, n1, n2){
  original_centroids <- dt[, {C1=mean(UMAP1);C2=mean(UMAP2); list(C1, C2)},by = original_annotation]
  setkey(original_centroids, original_annotation)
  
  immc_centroids <- dt[, {C1=mean(UMAP1);C2=mean(UMAP2); list(C1, C2)},by = ImmC]
  setkey(immc_centroids, ImmC)
  
  singler_centroids <- dt[, {C1=mean(UMAP1);C2=mean(UMAP2); list(C1, C2)},by = SingleR]
  setkey(singler_centroids, SingleR)
  
  garnett_centroids <- dt[, {C1=mean(UMAP1);C2=mean(UMAP2); list(C1, C2)},by = garnett]
  setkey(garnett_centroids, garnett)
  
  random_centroids <- dt[, {C1=mean(UMAP1);C2=mean(UMAP2); list(C1, C2)},by = random]
  setkey(random_centroids, random)
  

  
  
  immc_dist <-  data.table(celltype=immc_centroids$ImmC, dist = (immc_centroids[, -1] - original_centroids[immc_centroids$ImmC, -1])[, sqrt(C1^2+C2^2)], method = "ImmC")
  
  singler_dist <-  data.table(celltype=singler_centroids$SingleR,  dist=(singler_centroids[, -1] - original_centroids[singler_centroids$SingleR, -1])[, sqrt(C1^2+C2^2)], method = "SingleR")
  
  garnett_dist <-  data.table(celltype=garnett_centroids$garnett, dist = (garnett_centroids[, -1] - original_centroids[garnett_centroids$garnett, -1])[, sqrt(C1^2+C2^2)], method = "garnett")
  
  random_dist <-  data.table(celltype=random_centroids$random, dist = (random_centroids[, -1] - original_centroids[random_centroids$random, -1])[, sqrt(C1^2+C2^2)], method = "random")
 
  comb_dt <- rbind(immc_dist, singler_dist, garnett_dist, random_dist)

  ftmp <- comb_dt[!is.na(dist) & !celltype %in% c('Other', 'Unassigned'),]
  ftmp[,  method := factor(method, c('ImmC','SingleR', 'garnett','random'))]
  
  plt <- ggplot(ftmp, aes(x=method, y=dist)) +
    geom_boxplot(alpha = 0.2, color = "black", lwd = .3, outlier.shape = NA, width = .5) +
    geom_jitter(aes(fill = celltype), width = 0.2, size = 6, pch = 21, alpha = .8)+
    theme_bw() +
    theme(legend.pos = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, size = 32,hjust = 1, color= "black"),
          axis.text.y = element_text(size = 32, color = "black"),
          axis.title.x = element_text(angle = 0, size =0, color= "black"),
          axis.title.y = element_text(angle = 90, size = 0,color= "black")
    ) +
    scale_fill_manual(values = c(colorRampPalette(brewer.pal(n1, palette))(n2), 'lightgrey', 'black'),
                      guide = guide_legend(nrow=4,override.aes = list(size=4))) +
    stat_compare_means(data=ftmp,
                       mapping=aes(x=method, y=dist),
                       label = "p.signif",
                       method= "wilcox.test",
                       ref.group = 'random')
  
  return (plt)
}


plot_error <- function(ctrl.err, exp.err){
  data <- melt(ctrl.err - exp.err)
  
  max_abs <- max(abs(data$value)) * 1.1
  
  errorplt <- ggplot(data, aes(Var2, value)) +
    stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 1) +
    geom_jitter(aes(color = Var1), size = 2, alpha = .6, height = .1, width = .25)+
    theme_classic() + ylim (-max_abs, max_abs)
  
  return(errorplt)
}


