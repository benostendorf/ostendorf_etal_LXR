theme_dotplot <- 
  theme(
    axis.title = element_text(size = 6), 
    axis.text.x = element_text(size = 5), 
    axis.text.y = element_text(size = 5), 
    legend.text = element_text(size = 5), 
    legend.title = element_text(size = 5)
  )

## Plotting options
cols_RGX <- c("untreated" = "#A9A9A8", "treated" = "#0095FF")
cols_custom <- RColorBrewer::brewer.pal(3, "Set1")[c(2, 1)]

## Set ComplexHeatmap options
ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5), 
       "heatmap_column_names_gp" = gpar(fontsize = 5), 
       "heatmap_row_title_gp" = gpar(fontsize = 6), 
       "heatmap_column_title_gp" = gpar(fontsize = 7), 
       "legend_title_gp" = gpar(fontsize = 6), 
       "legend_labels_gp" = gpar(fontsize = 5), 
       "legend_grid_height" = unit(2, "mm"), 
       "legend_grid_width" = unit(2, "mm"))

get_max_coverage <- function(data, region){
  require(rtracklayer)
  require(GenomicRanges)
  ymax_ls <- list()
  for (i in seq_len(length(data))) {
    wig.data <- rtracklayer::import(data[i], format = "bigWig", selection=region)
    ymax_ls[i] <-  max(0, max(mcols(wig.data)[,1]))
  }
  return(ceiling(max(unlist(ymax_ls))))
}

custom_linewidth <- 0.2352539

theme_custom2 <-
  ggplot2::theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = custom_linewidth),
    axis.ticks = element_line(size = custom_linewidth),
    axis.ticks.length = unit(0.075, "cm"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5, colour = "black"),
    plot.title = element_text(size = 7, hjust = 0.5),
    panel.background = element_blank(),
    legend.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.key.size = unit(0.25, "line")
  )
