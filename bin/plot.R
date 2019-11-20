## ---------------------------
## Script name: plot.R
## Purpose of script: plot
## Author: LXF
## Date Created: 2019-01-01
## ---------------------------
## Notes: some notes
## ---------------------------

setwd("sample_path")

library(ComplexHeatmap)
library(circlize)
library(plyr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggthemes)
library(ggrepel)

##### scatter with gene annotation #####
infile <- "path/file"
outfile <- "path/outfile"
pair_mat <- read.table(inflie, header = T, row.names = 1, sep = '\t', check.names = F)
pair_mat_plot <- subset(pair_mat, pair_mat$geneA != pair_mat$geneB)
pair_mat_sub <- pair_mat_plot[pair_mat_plot$pval > 3, ]
pair_mat_sub$pair_name <- ifelse(pair_mat_sub$pval > 3, paste(pair_mat_sub$geneA, "+", pair_mat_sub$geneB), "")
ggplot(pair_mat_plot, aes(x = rate, y = pval)) + 
  geom_point(color = "darkviolet", alpha = 0.5, size = 3) + 
  labs(y = expression(paste("-log"[10], " P value")), x = "Co-occurrence rate") +
  theme_base() +
  theme(axis.title.x = element_text(size = 30), axis.text.x = element_text(size=25),
        axis.title.y = element_text(size = 30), axis.text.y = element_text(size=25)) +
  geom_label_repel(aes(label = pair_name), size = 5, data = pair_mat_sub) 

##### heatmap with up-right triangle for p-value, down-left triangle for co-occurrence #####
reshape_tri <- function(x, v_names){
  x_wide <- reshape(x, v.names = v_names, timevar = "geneA", idvar = "geneB", direction = "wide")
  row.names(x_wide) <- x_wide[, 1]
  x_wide <- x_wide[, -1]
  x_wide[upper.tri(x_wide)] <- x_wide[lower.tri(x_wide)]
  colnames(x_wide) <- gsub(pattern = paste0(v_names, "."), "", colnames(x_wide))
  x_wide
}

pair_rate_wide <- reshape_tri(pair_mat, v_names = "rate")
pair_p_wide <- reshape_tri(pair_mat, v_names = "pval")

# color
col_rate <- circlize::colorRamp2(c(0, 1), c("white", "darkviolet"))
col_p0 <- circlize::colorRamp2(c(0, 3), c("white", "darkgreen"))

# legend
lgd01 <- Legend(col_fun = col_rate, title = "Co-occurrence rate", at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
               direction = "horizontal", legend_height = unit(2, "cm"), legend_width = unit(6, "cm"))
lgd02 <- Legend(col_fun = col_p0, title = expression(paste("-log"[10], " P value")), at = c(0, 1, 2, 3, 4), 
               direction = "horizontal", legend_height = unit(2, "cm"), legend_width = unit(6, "cm"))
lgd_list0 <- list(lgd01, lgd02)

pdf(outfile, width = 8, height = 8)
h0 <- Heatmap(pair_rate_wide, rect_gp = gpar(type = "none"), 
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(col = "grey", fill = NA))
               if(i >= j) {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = col_rate(pair_rate_wide[i, j])))
               } else {
                 grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = col_p0(pair_p_wide[i, j])))
               }
             }, cluster_rows = FALSE, cluster_columns = FALSE,
             show_heatmap_legend = FALSE
)
draw(h0, annotation_legend_list = lgd_list0, annotation_legend_side = "top")
dev.off()

##### heatmap with boxplot and line plot as top annotation #####
nUse <- 13
nUnuse <- 2
nNum <- 2

# data arrange
dat2 <- read.csv(infile, header = T, row.names = NULL, check.names = F)
meta_variant <- read.table(inflie, heade = T, sep = "\t", row.names = 1, check.names = F)
dat_new <- dat2[order(dat2$gene_name), ] 
dat_new$sample_num <- ncol(dat_new) - rowSums(is.na(dat_new[, 1:(ncol(dat_new)-nUnuse)])) - nUnuse
dat_new$RowMax <- apply(dat_new[, 1:nUse], 1, max, na.rm = T)
dat_new$RowMedian <- apply(dat_new[, 1:nUse], 1, median, na.rm = T)
dat_new <- subset(dat_new, dat_new$sample_num >= nNum | dat_new$RowMax >= 0.5)
write.table(dat_new, outfile, sep = ",", col.names = T, row.names = F)

dat_new$variant_num <- table(dat_new$gene_name)[dat_new$gene_name]
dat_sorted_gene <- dat_new[order(dat_new$variant_num, decreasing = T),]
gene_order <- dat_sorted_gene$gene_name
level_order <- match(unique(gene_order), levels(gene_order))
gene_order <- factor(gene_order, levels(gene_order)[level_order])

dat_mat <- data.matrix(dat_sorted_gene[, 1:(ncol(dat_sorted_gene)-6)])
dat_sorted_sample <- dat_mat[, match(rownames(meta), colnames(dat_mat))]
sample_median <- apply(dat_sorted_sample, 1, median, na.rm = T)
sample_median <- data.frame(dat_sorted_gene$gene_name, sample_median)

# Make the plot
ha_top <- HeatmapAnnotation(
  # line plot
  sample_num = anno_lines(dim(dat_sorted_sample)[2] - colSums(is.na(t(dat_sorted_sample))), 
                          gp = gpar(col = "darkviolet"), pt_gp = gpar(col = "violet"), pch = 16, 
                          add_points = TRUE, height = unit(2, "cm")),
  # boxplot for rate in different sample group
  sample_prop = anno_boxplot(t(dat_sorted_sample), gp = gpar(fill = "steelblue"), height = unit(4, "cm")))

col_fun_sample <- colorRamp2(c(0,5), c("white", "darkviolet"))

# right annotation with sample type
ha_right <- rowAnnotation(
  sample = meta$sample_id,
  type = meta$sample_type,
  col = list(type = c("sample1" = "lightslategray", "sample2" = "lightsteelblue", "sample3" = "steelblue")),
  annotation_legend_param = list(type = list(title = "Type", nrow = 1))
  )

col_fun_heatmap <- colorRamp2(c(0, 1), c("white", "darkgreen"))
dat_plot <- t(sqrt(dat_sorted_sample))
dat_plot[is.na(dat_plot)] <- 0

pdf(outflie, width = 14, height = 8)
h1 <- Heatmap(dat_plot, cluster_rows = FALSE, cluster_columns = FALSE,
        column_split = gene_order, cluster_column_slices = FALSE,
        column_order = order(sample_median$sample_median, decreasing = TRUE),
        row_split = meta$col_sep, cluster_row_slices = FALSE,
        col = col_fun_heatmap,
        show_column_names = F, column_title_rot = 90,
        border = T, heatmap_legend_param = list(direction = "horizontal"),
        top_annotation = ha_top,
        right_annotation = ha_right)
draw(h1,  merge_legend = TRUE, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off() 


##### plot like circos #####
# data arrange
dat_com <- read.csv(inflie, header = T, row.names = NULL, check.names = F)
dat_com <- arrange(dat_com, gene_order, observation, individual)
dat_com[,1] <- dat_com[,1] * 100

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
nObsType <- nlevels(as.factor(dat_com$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(dat_com$group)*nObsType, ncol(dat_com)) )
colnames(to_add) <- colnames(dat_com)
to_add$group <- rep(gene_order, each=empty_bar*nObsType )
to_add$gene_order <- rep(c(1:35), each=5)
dat_com <- rbind(dat_com, to_add)
dat_com <- dat_com %>% arrange(gene_order, individual)
dat_com$id <- rep( seq(1, nrow(dat_com)/nObsType) , each=nObsType)
dat_com <- dat_com[, c(5,2,4,1,6)]
space <- tail(dat_com, 5)
space$id <- 106
dat_com <- rbind(dat_com, space)

# prepare a data frame for base lines
base_data <- dat_com %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(dat_com) +      
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.7) +
  scale_fill_viridis(discrete=TRUE) +
  
  # Add a val=100/75/50/25 lines.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(dat_use$id),5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , color="grey", size = 6, angle = 0, fontface="bold", hjust=1) +
  ylim(-150, 100) +
  theme_minimal() + theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) + coord_polar() +
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.8 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)


