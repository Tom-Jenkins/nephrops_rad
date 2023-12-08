# =========================== #
#
# Nephrops Variant Analysis 2023
#
# Utility R functions
# 
# =========================== #

plot_pca <- function(pca_scores, percent, genlight_obj, rad_samples_meta, axes = c(1,2),
                     axis_lab = "PC", title = "Principal components analysis",
                     by = "site", cols = palette.colors(nPop(genlight_obj), "Set 3"),
                     labsize = 3.5) {
  
  # Create a data.frame containing individual coordinates
  ind_coords <- as.data.frame(pca_scores)
  
  # Rename columns of data.frame
  colnames(ind_coords) <- c("Axis1","Axis2","Axis3")
  
  # Add a column containing individuals
  ind_coords$Ind <- indNames(genlight_obj)
  
  # Add a column with the population IDs
  ind_coords$Site <- genlight_obj$pop
  
  # Subset metadata data.frame by individuals in genlight object
  rad_samples_meta_sub <- dplyr::filter(rad_samples_meta, Ind_ID %in% indNames(genlight_obj))
  
  # Add a column with male and female
  sex_strata <- c()
  for (i in ind_coords$Ind) {
    ind_sex <- filter(rad_samples_meta_sub, Ind_ID == i)$Sex
    sex_strata[i] <- ind_sex
  }
  ind_coords$sex <- sex_strata
  ind_coords$sex <- factor(ind_coords$sex, levels = c("M","F",""))
  
  # Calculate centroid (mean average) position for each population
  centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ Site,
                        data = ind_coords,
                        FUN = mean)
  
  # Add centroid coordinates to ind_coords dataframe
  ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
  
  # Custom x and y labels
  xlab <- paste0(axis_lab, " ", axes[1], " (", format(round(percent[axes[1]], 1), nsmall=1)," %)", sep="")
  ylab <- paste0(axis_lab, " ", axes[2], " (", format(round(percent[axes[2]], 1), nsmall=1)," %)", sep="")
  
  # Custom ggplot2 theme
  ggtheme <- theme(legend.title = element_blank(),
                   axis.text.y = element_text(colour="black", size=10),
                   axis.text.x = element_text(colour="black", size=10),
                   axis.title = element_text(colour="black", size=14),
                   legend.position = "right",
                   legend.text = element_text(size=15),
                   legend.key = element_rect(fill = NA),
                   legend.key.size = unit(0.7, "cm"),
                   legend.box.spacing = unit(0, "cm"),
                   panel.border = element_rect(colour="black", fill=NA, linewidth=1),
                   panel.background = element_blank(),
                   # title centered
                   plot.title = element_text(hjust=0.5, size=15) 
  )
  
  # Scatter plot coloured by site
  if (by == "site") {
    plt <- ggplot(data = ind_coords,
                  aes(x = !!parse_expr(paste0("Axis",axes[1])),
                      y = !!parse_expr(paste0("Axis",axes[2]))))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      # spider segments
      geom_segment(aes(xend = !!parse_expr(paste0("Axis",axes[1],".cen")),
                       yend = !!parse_expr(paste0("Axis",axes[2],".cen")),
                       colour = Site),
                    show.legend = FALSE, linewidth = 0.3)+
      # points
      geom_point(aes(fill = Site), shape = 21, size = 3,  show.legend = FALSE)+
      # centroids
      geom_label(data = centroid, aes(label = Site, fill = Site), size = labsize, alpha = 0.9,
                 show.legend = FALSE, label.padding = unit(0.1, "cm"))+
      # colouring
      scale_fill_manual(values = cols)+
      scale_colour_manual(values = cols)+
      # custom labels
      labs(x = xlab, y = ylab)+
      # title
      ggtitle(title)+
      # custom theme
      ggtheme
  }
  
  # Scatter plot coloured by sex
  if (by == "sex") {
    plt <- ggplot(data = ind_coords,
                  aes(x = !!parse_expr(paste0("Axis",axes[1])),
                      y = !!parse_expr(paste0("Axis",axes[2]))))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_point(aes(fill = sex), shape = 21, size = 3, show.legend = TRUE)+
      scale_fill_manual("Sex",
                        values = c("royalblue","#dd1c77","grey"),
                        labels = c("Male","Female","n/a"))+
      labs(x = xlab, y = ylab)+
      ggtitle("Nephrops PCA: coloured by sex")+
      ggtheme
  }
  
  return(plt)
}

plot_cross_entropy <- function(snmf_obj) {
  min_ce <- summary(snmf_obj)$crossEntropy["min", ]
  min_ce_df <- data.frame(K = 1:length(min_ce), cross_entropy = min_ce)
  ggplot(data = min_ce_df)+
    geom_point(aes(x=K, y=cross_entropy), shape = 19, size = 4, colour = "red")+
    scale_x_continuous(breaks = 1:length(min_ce))+
    ylab("Cross entropy")+
    xlab("Number of ancestral populations (K)")+
    ggtitle("SNMF: Cross entropy per K")
}

plot_admixture = function(snmf_obj, genlight_obj, K) {
  
  # Find and extract the lowest cross-entropy of all runs for K
  lowest_ce <- which.min(cross.entropy(snmf_obj, K = K))
  
  # Extract Q-matrix for the best run
  qmatrix <- as.data.frame(Q(snmf_obj, K = K, run = lowest_ce))
  
  # Label column names of qmatrix
  colnames(qmatrix) <- sapply(1:K, function(x) paste("Cluster", x))
  
  # Add individual IDs
  qmatrix$Ind = indNames(genlight_obj)
  
  # Add site IDs
  qmatrix$Site = as.character(genlight_obj$pop)
  
  # Convert data.frame to long format
  qmatrix <- pivot_longer(
    data = qmatrix,
    cols = starts_with("Cluster"),
    names_to = "Cluster",
    values_to = "Proportion"
  )
  
  # Plot admixture barplot
  barchart = ggplot(data=qmatrix, aes(x=Ind, y=Proportion, fill=Cluster))+
    geom_bar(stat = "identity", width=1, colour="black", linewidth=0.25)+
    coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    # scale_fill_manual(values = as.character(colours))+
    ylab("Admixture proportion")+
    # xlab("Individual")+
    theme(axis.text.y = element_text(face = "bold", colour = "black", size = 10),
          axis.text.x = element_text(angle = 360, colour = "black", size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.margin = margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")
    )
  return(barchart)
}

