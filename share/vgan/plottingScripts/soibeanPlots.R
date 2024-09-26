suppressPackageStartupMessages({
  
  library(ggplot2)
  library(phytools)
  library(ggtree)
  library(rentrez)
  library(RColorBrewer)
  library(gridExtra)
  library(ggfx)
  library(cowplot)
})

get_scientific_name <- function(input_string) {
  # Extract the accession number from the input string
  accession_number <- substr(input_string, 1, regexpr("\\.", input_string) - 1)
  
  # Fetch the record from GenBank
  record <- entrez_fetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
  
  # Extract the organism name using regular expressions
  organism_pattern <- "ORGANISM\\s+(.*?)\\n"
  matches <- regmatches(record, regexpr(organism_pattern, record, perl=TRUE))
  if (length(matches) > 0) {
    scientific_name <- trimws(matches[[1]])
    scientific_name <- sub("ORGANISM\\s+", "", scientific_name)
  } else {
    scientific_name <- "Unknown"
  }
  
  return(scientific_name)
}


make_tree <- function(tree_file){

  #############################################################
  # Read in the Newick tree and create the basic ggtree plot
  tree <- read.tree(tree_file)
  for (s in 1:length(tree$tip.label)){
    new_label <- get_scientific_name(tree$tip.label[s])

    tree$new.label[s] <- new_label

  }

  tree$id.label <- tree$tip.label
  tree$tip.label <- tree$new.label

return(tree)
}

################################################################################
### MAKE COMBO PLOT ###
comboPlot <- function(tree_file, input_file, p){
  
  extracted <- sub("Result.*", "", input_file)
  print(extracted)
  # Extract the word before '.new.dnd' and after the last '/'
  taxname <- sub("^.*/([^.]+)\\.new\\.dnd$", "\\1", tree_file)

  print(taxname)
  
  #############################################################
  # Read in the Newick tree and create the basic ggtree plot
  #tree <- read.tree(tree_file)
 # for (s in 1:length(tree$tip.label)){
   # new_label <- get_scientific_name(tree$tip.label[s])
    
   # tree$new.label[s] <- new_label
    
  #}
  
  #tree$id.label <- tree$tip.label
  #tree$tip.label <- tree$new.label
  
  pt <- ggtree(p) + 
    geom_tiplab() +
    geom_tippoint() +
    geom_nodelab()
  
  ########################################################################

  df <- read.csv(input_file, header = T, sep = '\t')
  
  # Determine number of 4-column sources
  num_sources <- (ncol(df)-1) / 4
  
  stacked_df <- data.frame() # Initialize an empty dataframe
  
  # Loop through each source and stack it
  for(i in 1:num_sources) {
    start_col <- (i-1)*4 + 1
    end_col <- i*4
    
    df_source <- df[, start_col:end_col]
    colnames(df_source) <- colnames(df[,1:4])
    df_source$Source <- paste0("Source_", i)  # Adding a source identifier
    
    stacked_df <- rbind(stacked_df, df_source)
  }
  
  
  stacked_df$Freq <- with(stacked_df, ave(branch_position_derived, list(branch_position_derived, Source), FUN = length))
  
  
  # Extract layout details from the ggtree plot
  d <- ggplot_build(pt)$data[[1]]
  nameList <- c(p$tip.label, p$node.label) 
  
  # Identify the indexes of the specific labels
  name_matx <- c(p$id.label, p$node.label)
  labels <- unique(c(stacked_df$Source_1))
  indexes <- match(labels, name_matx)

  
  col <- nameList[indexes]
  
  # Assume num_sources is the number of sources you have
  num_sources <- num_sources
  
  # Initialize an empty list to hold the column names
  column_names_list <- list()
  
  # Loop through each source and generate the column name
  for (i in 1:num_sources) {
    if (i == 1) {
      column_name <- "branch_position_derived"
    } else {
      column_name <- paste0("branch_position_derived.", i - 1)
    }
    column_names_list[[i]] <- column_name
  }
  
  # Print the list to check the result
  clname <- unlist(column_names_list)
  
  # Using lapply and do.call for efficient row binding
  points_list <- lapply(1:length(indexes), function(i) {
    index <- indexes[i]
    color <- col[i]
    
    x <- d[index, 3]
    xend <- d[index, 7]
    column_name <- clname[i]
    
    data.frame(
      x = x + df[[column_name]] * (xend - x),
      y = d[index, 8],
      color = color
    )
  })
  points_df <- do.call(rbind, points_list)
  points_df$freq <- stacked_df$Freq
  points_df$like <- stacked_df$Log.likelihood


  # Calculate the maximum possible difference
  
  df_u <- unique(points_df)
  df_u$original_y <- df_u$y
  median_loglikelihood <- median(df_u[,5])
  
  
  # Calculate the maximum possible difference
  max_diff <- max(abs(df_u[,5] - median_loglikelihood))
  
  # Adjust the y value from their original positions
  df_u$y <- ifelse(df_u[,5] < median_loglikelihood, 
                   df_u$original_y - (abs(df_u[,5] - median_loglikelihood) / max_diff) * 0.2, 
                   df_u$original_y + (abs(df_u[,5] - median_loglikelihood) / max_diff) * 0.2)
  
  
  
  
  xmini <- min(df_u$x)
  xmaxi <- max(df_u$x)
  ymini <- min(df_u$y) 
  ymaxi <- max(df_u$y)
  
  
  df_u$log_freq <- log1p(df_u$freq)
  tr <- ggtree(p) +                             
    geom_tiplab(size = 3, hjust = 0.01) + 
    geom_tippoint() +
    geom_nodelab(geom = 'label', size = 3) + 
    geom_point(data=df_u, aes(x=x, y=y, color=like, size=log_freq, alpha=0.5)) +
    scale_color_gradientn(colors = c("grey", "yellow", "orange", "red"),
                          values = scales::rescale(c(min(df_u$log_freq), 
                                                     quantile(df_u$log_freq, 0.33),
                                                     quantile(df_u$log_freq, 0.66),
                                                     max(df_u$log_freq)))) +
    scale_size_continuous(name = "Frequency [log]", 
                          range = c(1, 2)) +  # Adjust the range as needed
    labs(color = "Likelihood [log]") +  # Adjust if log likelihood is different from log frequency
    guides(alpha="none") + 
    theme(
      plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"),
      legend.position = "bottom",
      legend.key.size = unit(1.25, "cm"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    ) +
    coord_cartesian(clip = "off")
  

  ###############################################################################
  
  stacked_df$Source <- points_df$color
  
  # Calculate medians for each source
  medians <- with(stacked_df, tapply(proportion, Source, median))
  
  # Get unique colors based on unique sources
  unique_sources <- unique(stacked_df$Source)
  num_sources <- length(unique_sources)
  colors <- scales::hue_pal()(num_sources)
  
  
  # Create a data frame of medians and their corresponding colors
  medians_df <- data.frame(
    median = medians[unique_sources],
    color = colors[1:length(unique_sources)]
  )
  
  plot <- ggplot(stacked_df, aes(x = proportion, y = after_stat(density))) +
    stat_density(aes(color = Source),geom = "line", position = "identity", size = 1.2) +
    geom_vline(data = data.frame(Source = names(medians), Median = as.numeric(medians)),
               aes(xintercept = Median, color = Source), linetype = "dotted", size = 1.2) +
    theme_bw() +
    scale_color_manual(values = colors) +
    xlim(0,1) +
    labs(x = "Proportion", y = "Density")

  #banner <- create_banner(label)
    
  combo <- plot_grid(tr, plot, ncol = 2, rel_widths = c(10, 5))
  #ggsave(paste0(extracted,"k",num_sources,"sourceProportion.png"), plot = plot, dpi = 600, height = 7, width = 10)
  return(combo)
}



########################## MAKE TREEPLOT WITH MULTIPLE SOURCES ###############################
make_tree_plots <- function(tree_file, input_file){
  
  extracted <- sub("Result.*", "", input_file)
  print(extracted)

  # Extract the word before '.new.dnd' and after the last '/'
  taxname <- sub("^.*/([^.]+)\\.new\\.dnd$", "\\1", tree_file)

  print(taxname)

  #############################################################
  # Read in the Newick tree and create the basic ggtree plot
  tree <- read.tree(tree_file)
  for (s in 1:length(tree$tip.label)){
    new_label <- get_scientific_name(tree$tip.label[s])
    
  
    tree$new.label[s] <- new_label
    
  }
  
  tree$id.label <- tree$tip.label
  tree$tip.label <- tree$new.label
  
  p <- ggtree(tree) + 
    geom_tiplab() +
    geom_tippoint() +
    geom_nodelab()
  
  ########################################################################
  
  
  df <- read.csv(input_file, header = T, sep = '\t')
  
  # Determine number of 4-column sources
  num_sources <- (ncol(df)-1) / 4
  
  stacked_df <- data.frame() # Initialize an empty dataframe
  
  # Loop through each source and stack it
  for(i in 1:num_sources) {
    start_col <- (i-1)*4 + 1
    end_col <- i*4
    
    df_source <- df[, start_col:end_col]
    colnames(df_source) <- colnames(df[,1:4])
    df_source$Source <- paste0("Source_", i)  # Adding a source identifier
    
    stacked_df <- rbind(stacked_df, df_source)
 }
  
  
  stacked_df$Freq <- with(stacked_df, ave(branch_position_derived, list(branch_position_derived, Source), FUN = length))
  
  
  # Extract layout details from the ggtree plot
  d <- ggplot_build(p)$data[[1]]
  nameList <- c(tree$tip.label, tree$node.label) 
  
  # Identify the indexes of the specific labels
  name_matx <- c(tree$id.label, tree$node.label)
  labels <- unique(stacked_df$Source_1)
  indexes <- match(labels, name_matx)
  matched_values <- match(stacked_df$Source_1, name_matx)
  stacked_df$Source_1 <- nameList[matched_values]
  
  # Initialize an empty list to hold the column names
  column_names_list <- list()
  
  # Loop through each source and generate the column name
  for (i in 1:num_sources) {
    if (i == 1) {
      column_name <- "branch_position_derived"
    } else {
      column_name <- paste0("branch_position_derived.", i - 1)
    }
    column_names_list[[i]] <- column_name
  }
  
  # Print the list to check the result
  print(column_names_list)
  clname <- unlist(column_names_list)
  
  
  # Assuming 'd', 'nameList', and other necessary variables are already calculated
  
  # Add columns for plot_x, plot_y, and color in stacked_df
  stacked_df$plot_x <- NA
  stacked_df$plot_y <- NA
  
  # Loop through each row in stacked_df
  for (i in 1:nrow(stacked_df)) {
    # Extract the relevant identifier (e.g., label) for this row to find the matching index in 'd'
    
    label_identifier <- stacked_df$Source_1[i]  # Replace 'label' with the actual identifier column in stacked_df
  
    index <- match(label_identifier, nameList)
  
    
    # Proceed only if a matching index is found
    if (!is.na(index)) {
      x <- d[index, 3]  # Replace 'x_column' with the actual column name in 'd'
      xend <- d[index, 7]  # Replace 'xend_column' with the actual column name in 'd'
      column_name <- stacked_df$branch_position_derived[i]  # Assuming 'Source' contains the relevant data column name
      
      # Calculate the new x position
      adjusted_x <- x + column_name * (xend - x)
      
      # Update stacked_df with the new x and y positions and color
      stacked_df$plot_x[i] <- adjusted_x
      stacked_df$plot_y[i] <- d[index, 8]  # Replace 'y_column' with the actual y coordinate column in 'd'
    }
  }
  
  
  
  
  # Calculate the maximum possible difference
  
  df_u <- unique(stacked_df)
  df_u$original_y <- df_u$plot_y
  median_loglikelihood <- median(df_u$Log.likelihood)
  
  
  # Calculate the maximum possible difference
  max_diff <- max(abs(df_u$Log.likelihood - median_loglikelihood))
  
  # Adjust the y value from their original positions
  df_u$plot_y <- ifelse(df_u$Log.likelihood < median_loglikelihood, 
                   df_u$original_y - (abs(df_u$Log.likelihood - median_loglikelihood) / max_diff) * 0.2, 
                   df_u$original_y + (abs(df_u$Log.likelihood - median_loglikelihood) / max_diff) * 0.2)
  
  

  xmini <- 0.05
  xmaxi <- 0.11
  ymini <- 8 
  ymaxi <- 12
  
  df_u$log_freq <- df_u$Freq
  tr <- ggtree(tree) +                             
    geom_tiplab(size = 3, hjust = 0.01) + 
    geom_tippoint() +
    geom_nodelab(geom = 'label', size = 3) + 
    geom_point(data=df_u, aes(x=plot_x, y=plot_y, color=Log.likelihood, size=log_freq, alpha=0.5)) +
    scale_color_gradientn(colors = c("grey", "yellow", "orange", "red"),
                          values = scales::rescale(c(min(df_u$log_freq), 
                                                     quantile(df_u$log_freq, 0.33),
                                                     quantile(df_u$log_freq, 0.66),
                                                     max(df_u$log_freq)))) +
    scale_size_continuous(name = "Frequency [log]", 
                          range = c(1, 2)) +  # Adjust the range as needed
    labs(color = "Likelihood [log]") +  # Adjust if log likelihood is different from log frequency
    guides(alpha="none") + 
    theme(
      plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"),
      legend.position = "bottom",
      legend.key.size = unit(1.25, "cm"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    ) +
    coord_cartesian(clip = "off")
  

  ###############################################################################
  
  #stacked_df$Source <- points_df$color
  
  # Calculate medians for each source
  medians <- with(stacked_df, tapply(proportion, Source, median))
  
  # Get unique colors based on unique sources
  unique_sources <- unique(stacked_df$Source)
  num_sources <- length(unique_sources)
  colors <- scales::hue_pal()(num_sources)
  
  
  # Create a data frame of medians and their corresponding colors
  medians_df <- data.frame(
    median = medians[unique_sources],
    color = colors[1:length(unique_sources)]
  )
  
  plot <- ggplot(stacked_df, aes(x = proportion, y = after_stat(density))) +
    stat_density(aes(color = Source),geom = "line", position = "identity", size = 1.2) +
    #geom_vline(data = data.frame(Source = names(medians), Median = as.numeric(medians)),
    #           aes(xintercept = Median, color = Source), linetype = "dotted", size = 1.2) +
    theme_bw() +
    scale_color_manual(values = colors) +
    xlim(0,1) +
    labs(x = "Proportion", y = "Density")

  #banner <- create_banner(label)
    
  combo <- plot_grid(tr, plot, ncol = 2, rel_widths = c(10, 5))
  #ggsave(paste0(extracted,"k",num_sources,"sourceProportion.png"), plot = plot, dpi = 600, height = 7, width = 10)
  return(combo)
}




args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
input_file <- args[2]

if (length(args) != 2) {
  stop("Rscript soibeanPlots.R tree_file input_file")
}

result_pos <- regexpr("Result", input_file)
number_start_pos <- result_pos + attr(result_pos, "match.length")
number <- substr(input_file, number_start_pos, number_start_pos)
number_as_integer <- as.integer(number)

# Print the result

extracted <- sub("Result.*", "", input_file)
df <- read.csv(input_file, header = T, sep = '\t')
  
# Determine number of 4-column sources
num_sources <- (ncol(df)-1) / 4
  
stacked_df <- data.frame() # Initialize an empty dataframe

# Loop through each source and stack it
for(i in 1:num_sources) {
  start_col <- (i-1)*4 + 1
  end_col <- i*4
  
  df_source <- df[, start_col:end_col]
  colnames(df_source) <- colnames(df[,1:4])
  df_source$Source <- paste0("Source_", i)  # Adding a source identifier
  
  stacked_df <- rbind(stacked_df, df_source)
}
 
labels <- unique(stacked_df$Source_1)

if(number_as_integer == length(labels)){
    tre <- make_tree(tree_file)
    combo <- comboPlot(tree_file, input_file, tre)
    ggsave(paste0(extracted,"k",number_as_integer,"Comb.png"), plot = combo, dpi = 600, height = 10, width = 18)
}else{
    
    combo <- make_tree_plots(tree_file, input_file)
    ggsave(paste0(extracted,"k",number_as_integer,"Comb.png"), plot = combo, dpi = 600, height = 10, width = 18)
}

