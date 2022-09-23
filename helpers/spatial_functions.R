


PlotEmbedding <- function(
  seurat_obj,
  group.by,
  reduction = 'umap',
  split.by = NULL,
  plot_under = FALSE,
  label = TRUE,
  point_size = 1,
  legend_point_size = 5,
  text_size = 3,
  raster = TRUE,
  order_points = "shuffle",
  x = NULL,
  y = NULL,
  raster_dpi = 400,
  raster_scale = 1,
  selected = NULL,
  plot_theme = NULL,
  color_df = NULL,
  plot_ratio=TRUE
){

  plot_df <- seurat_obj@meta.data
  if(!is.null(x) & !is.null(y)){
    plot_df$x <- plot_df[[x]]
    plot_df$y <- plot_df[[y]]
  } else if(reduction %in% names(seurat_obj@reductions)){
    plot_df$x <- seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$y <- seurat_obj@reductions[[reduction]]@cell.embeddings[,2]
  }


  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  if(label){
    for(cur_cluster in unique(plot_df[[group.by]])){
      cur_meta <- plot_df[plot_df[[group.by]] == cur_cluster,]
      df <- data.frame(
        cluster = cur_cluster,
        x = mean(cur_meta$x),
        y = mean(cur_meta$y)
      )
      centroid_df <- rbind(centroid_df, df)
    }
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x='x', y='y', color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name

  }

  # only show selected groups
  if(!is.null(selected)){
    plot_df[[group.by]][!(plot_df[[group.by]] %in% selected)] <- NA

    if(is.factor(plot_df[[group.by]])){
      plot_df[[group.by]] <- droplevels(plot_df[[group.by]])
    }

    if(label){
      centroid_df <- subset(centroid_df, cluster %in% selected)
    }
  }

  print(levels(plot_df[[group.by]]))

  # shuffle points:
  if(order_points == "shuffle"){
    plot_df <- plot_df[sample(nrow(plot_df)),]
  }

  # plot a single embedding
  if(is.null(split.by)){
    p <- .PlotSingleEmbedding(plot_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, plot_theme=plot_theme, plot_ratio=plot_ratio)
  } else{

    split_groups <- unique(plot_df[[split.by]])
    plot_list <- lapply(split_groups, function(cur_split){
      print(cur_split)
      cur_df <- plot_df[plot_df[[split.by]] == cur_split,]
      split_df <- plot_df[plot_df[[split.by]] != cur_split,]
      .PlotSingleEmbedding(cur_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, split_df, cur_split, plot_theme, plot_under, plot_ratio=plot_ratio)
    })

    return(plot_list)

  }

  p

}

.PlotSingleEmbedding <- function(
  plot_df,
  group.by,
  label,
  raster,
  raster_dpi,
  raster_scale,
  point_size,
  legend_point_size,
  text_size,
  color_df,
  centroid_df,
  split_df = NULL,
  cur_split = NULL,
  plot_theme = NULL,
  plot_under=FALSE,
  plot_ratio = TRUE
){

  p <- plot_df %>%
    ggplot(aes_string(x='x', y='y', color=group.by))

  # add the under layer for the split plots
  if(!is.null(split_df) & plot_under){
    # add points
    if(raster){
      p <-  p + ggrastr::rasterise(geom_point(data = split_df, size=point_size/2, color='lightgrey'), dpi=raster_dpi, scale=raster_scale/2)
    } else{
      p <- p + geom_point(data = split_df, size=point_size/2, color='lightgrey')
    }
  }

  # add points
  if(raster){
    p <-  p + ggrastr::rasterise(geom_point(size=point_size), dpi=raster_dpi,  scale=raster_scale)
  } else{
    p <- p + geom_point(size=point_size)
  }

  # add labels
  if(label){
    p <- p + ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=text_size)
  }

  # color scheme
  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)
  p <- p + scale_color_manual(values=plot_colors, na.value = 'grey90')

  # add title:
  if(is.null(cur_split)){
    p <- p + ggtitle(group.by) + theme(plot.title=element_text(hjust=0.5))
  } else{
    p <- p + ggtitle(cur_split) + theme(plot.title=element_text(hjust=0.5))
  }

  # add theme:
  if(!is.null(plot_theme)){
    p <- p + plot_theme
  }

  # fixed coords
  if(plot_ratio){
    p <- p + coord_equal()
  }

  # adjust legend point size
  p <- p + guides( color = guide_legend(override.aes = list(size=legend_point_size)))

  p

}

VisDimPlot <- function(
  seurat_obj,
  group.by,
  sample_col = "Sample",
  sample_labels = NULL,
  samples_to_plot = NULL,
  raster = TRUE,
  dpi = 800,
  text_size = 8,
  ncol=4,
  color_df = NULL
){

  # check that row and col are in the seurat metadata:
  if(!all(c("col", "row") %in% colnames(seurat_obj@meta.data))){
    stop("Spatial coordinates must be present in seurat@meta.data, named row and col.")
  }

  plot_df <- seurat_obj@meta.data


  # check if the group.by is in the meta-data
  if(group.by %in% colnames(plot_df)){
    if(is.numeric(plot_df[[group.by]])){
      stop("Specified group.by a numeric variable. Try plotting with SampleFeaturePlot?")
    }
  } else{
    stop("group.by not found in colnames(seurat_obj@meta.data).")
  }


  # get list of samples to plot:
  if(is.null(samples_to_plot)){
    samples_to_plot <- unique(plot_df[[sample_col]])
  } else{
    if(!all(samples_to_plot %in% unique(plot_df[[sample_col]]))){
      stop("Some specified samples are not found in seurat_obj@meta.data[[sample_col]].")
    }
  }

  # check sample labels:
  if(!is.null(sample_labels)){
    if(!all(sample_labels %in% colnames(seurat_obj@meta.data))){
      stop("Some sample_labels are not found in colnames(seurat_obj@meta.data)")
    }
  }

  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]


  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]

  # group by sample_labels
  if(!is.null(sample_labels)){
    plot_df <- plot_df %>% group_by(across(all_of(sample_labels))) %>%
      as.data.frame()
  }

  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){
    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x=1, y=1, color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name
  } else if(is.character(color_df)){
    # check if it's a named character vec:
    if(is.null(names(color_df))){
      stop("color_df must be a dataframe with a group column and a colour column, or color_df can be a named character vector where the values are the colors and the names are the corresponding groups.")
    }
    color_df <- data.frame(
      colour = as.character(color_df),
      group = names(color_df)
    )
  }

  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)


  # get the ordering of samples to plot:
  order_df <- plot_df %>%
    group_by(across(all_of(sample_labels)))  %>%
    arrange(.by_group=TRUE) %>%
    dplyr::select(all_of(c(sample_labels, sample_col))) %>%
    distinct()

  samples_to_plot <- as.character(order_df[[sample_col]])


  plot_list <- list()
  for(cur_sample in samples_to_plot){
    print(cur_sample)

    # get data for this sample and get hex vertex coordinates
    cur_df <- plot_df[plot_df[[sample_col]] == cur_sample,]
    vertices <- make_hex_spots(cur_df, cur_df[[group.by]])

    # get sample labels:
    if(!is.null(sample_labels)){
      cur_labels <- unlist(lapply(sample_labels, function(x){unique(as.character(cur_df[[x]]))}))
      plot_label <- paste0(cur_labels, collapse=', ')
    } else{
      plot_label <- cur_sample
    }

    # initialize gpgplot
    p <- vertices  %>%
      ggplot(aes_(x=~x.vertex, y=~y.vertex,group=~spot, fill=~fill))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p + ggrastr::rasterise(geom_polygon(size=0), dpi=dpi)
    } else{
      p <- p + geom_polygon(size=0)
    }

    # add extras to plot:
    p <- p +
      labs(fill = group.by) +
      scale_fill_manual(values=plot_colors, na.value = 'grey90') +
      coord_equal() +
      ggtitle(plot_label) +
      theme(
        plot.title = element_text(hjust=0.5, face='plain', size=text_size),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()

      )

    plot_list[[cur_sample]] <- p

  }

  if(length(plot_list) == 1){
    return(plot_list[[1]])
  } else{
    patch <- wrap_plots(plot_list, ncol=ncol, widths=1, heights=1) +
      plot_layout(guides='collect') +
      plot_annotation(title = group.by) &  theme(plot.title = element_text(hjust=0.5),)
  }
  patch

}

SampleFeaturePlot <- function(
  seurat_obj,
  feature,
  sample_col = "Sample",
  sample_labels = NULL,
  samples_to_plot = NULL,
  raster = TRUE,
  slot = 'data',
  assay = NULL,
  dpi = 800,
  ncol = 4,
  plot_max = 'q100',
  plot_min = 'q0',
  text_size = 8,
  colfunc = colorRampPalette(brewer.pal(11, 'Spectral' )),
  rev_colors = FALSE
){

  # check that row and col are in the seurat metadata:
  if(!all(c("col", "row") %in% colnames(seurat_obj@meta.data))){
    stop("Spatial coordinates must be present in seurat@meta.data, named row and col.")
  }


  plot_df <- seurat_obj@meta.data

  # check if the feature is in the meta-data
  if(feature %in% rownames(seurat_obj)){
    if(is.null(assay)){assay <- seurat_obj@active.assay}
    plot_df$SpatialPlotFeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
  } else if(feature %in% colnames(plot_df)){
    if(!is.numeric(plot_df[[feature]])){
      stop("Specified feature is not numeric. Try plotting with VisDimPlot?")
    }
    plot_df$SpatialPlotFeature <- plot_df[[feature]]
  } else{
    stop("feature not found in rownames(seurat_obj) or in colnames(seurat_obj@meta.data).")
  }

  # get list of samples to plot:
  if(is.null(samples_to_plot)){
    samples_to_plot <- unique(plot_df[[sample_col]])
  } else{
    if(!all(samples_to_plot %in% unique(plot_df[[sample_col]]))){
      stop("Some specified samples are not found in seurat_obj@meta.data[[sample_col]].")
    }
  }

  # check sample labels:
  if(!is.null(sample_labels)){
    if(!all(sample_labels %in% colnames(seurat_obj@meta.data))){
      stop("Some sample_labels are not found in colnames(seurat_obj@meta.data)")
    }
  }

  # subset plot_df by specified samples:
  plot_df <- plot_df[plot_df[[sample_col]] %in% samples_to_plot,]

  # group by sample_labels
  if(!is.null(sample_labels)){
    plot_df <- plot_df %>% group_by(across(all_of(sample_labels))) %>%
      as.data.frame()
  }

  plot_range <- range(plot_df$SpatialPlotFeature)
  if(!is.null(plot_max)){
    if(is.character(plot_max)){
      quant <- as.numeric(gsub('q', '', plot_max)) / 100
      plot_max <- as.numeric(quantile(plot_df$SpatialPlotFeature, quant))
    }
    plot_range[2] <- plot_max
    print(plot_max)
    plot_df$SpatialPlotFeature <- ifelse(
      plot_df$SpatialPlotFeature > plot_max,
      plot_max,
      plot_df$SpatialPlotFeature
    )
  }

  if(is.character(plot_min)){
    quant <- as.numeric(gsub('q', '', plot_min)) / 100
    plot_min <- as.numeric(quantile(plot_df$SpatialPlotFeature, quant))
  }
  plot_range[1] <- plot_min

  # get the ordering of samples to plot:
  order_df <- plot_df %>%
    group_by(across(all_of(sample_labels)))  %>%
    arrange(.by_group=TRUE) %>%
    dplyr::select(all_of(c(sample_labels, sample_col))) %>%
    distinct()

  samples_to_plot <- as.character(order_df[[sample_col]])

  plot_list <- list()
  for(cur_sample in samples_to_plot){
    print(cur_sample)

    # get data for this sample and get hex vertex coordinates
    cur_df <- plot_df[plot_df[[sample_col]] == cur_sample,]
    vertices <- make_hex_spots(cur_df, cur_df$SpatialPlotFeature)

    # get sample labels:
    if(!is.null(sample_labels)){
      cur_labels <- unlist(lapply(sample_labels, function(x){unique(as.character(cur_df[[x]]))}))
      plot_label <- paste0(cur_labels, collapse=', ')
    } else{
      plot_label <- cur_sample
    }

    # initialize gpgplot
    p <- vertices %>% subset(fill > plot_min) %>%
      ggplot(aes_(x=~x.vertex, y=~y.vertex,group=~spot, fill=~fill))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p +
        ggrastr::rasterise(
          geom_polygon(data = subset(vertices, fill <= plot_min), fill='lightgrey', size=0),
        dpi=dpi/2
      ) +
      ggrastr::rasterise(geom_polygon(size=0), dpi=dpi)
    } else{
      p <- p +
        geom_polygon(data = subset(vertices, fill <= plot_min), fill='lightgrey', size=0) +
        geom_polygon(size=0)
    }

    # add extras to plot:
    #colfunc <- colorRampPalette(brewer.pal(11, 'Spectral' ))
    colors <- colfunc(256)
    if(rev_colors){colors <- rev(colors)}
    p <- p +
      labs(fill = feature) +
      scale_fill_gradientn(colors=colors, limits = plot_range) +
      coord_equal() +
      ggtitle(plot_label) +
      theme(
        plot.title = element_text(hjust=0.5, face='plain', size=text_size),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()

      )

    plot_list[[cur_sample]] <- p

  }

  if(length(plot_list) == 1){
    return(plot_list[[1]])
  } else{
    patch <- wrap_plots(plot_list, ncol=ncol, widths=1, heights=1) +
      plot_layout(guides='collect') +
      plot_annotation(title = feature) &  theme(plot.title = element_text(hjust=0.5),)
  }
  patch
}











################################################################################
# Helper functions
################################################################################



make_hex_spots <- function(cdata, fill){
  r <- 1/2
   R <- (2/sqrt(3)) * r

  spot_positions <- select_spot_positions(cdata, fill = fill)
  spot_positions <- adjust_hex_centers(spot_positions)
  vertex_offsets <- data.frame(x.offset = c(0, r, r, 0, -r,
      -r), y.offset = c(-R, -R/2, R/2, R, R/2, -R/2))
  spot_vertices <- make_spot_vertices(spot_positions, vertex_offsets)
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  spot_vertices
}

select_spot_positions <- function (cdata, x = "col", y = "row", fill = "spatial.cluster")
{
    assertthat::assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
    }
    else if (is.vector(fill) || is.factor(fill)) {

        assertthat::assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y)]
        colnames(spot_positions) <- c("x.pos", "y.pos")
        spot_positions$fill <- fill
    }
    spot_positions$spot <- rownames(spot_positions)
    spot_positions
}

adjust_hex_centers <- function (spot_positions)
{
    r <- 1/2
    R <- (2/sqrt(3)) * r

    spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) +
        1
    spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) +
        1
    spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
    spot_positions$x.pos <- (spot_positions$x.pos + 1)/2

    spot_positions
}

make_spot_vertices <- function (spot_positions, vertex_offsets)
{
    spot_vertices <- merge(spot_positions, vertex_offsets)
    spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
    spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
    as.data.frame(spot_vertices)
}
