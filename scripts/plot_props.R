library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(progress)
library(tidyr)

# Function to calculate median like the awk command
calculate_median <- function(values) {
  sorted_values <- sort(values, na.last = TRUE)
  n <- length(sorted_values)
  if (n %% 2 == 1) {
    median <- sorted_values[(n + 1) / 2]
  } else {
    median <- mean(sorted_values[(n / 2):(n / 2 + 1)])
  }
  return(median)
}

# Function to plot MCMC data on a small tree with progress reporting
plot_props <- function(mcmc_data, output_file_path, source_names, ground_truth = NULL, title_suffix = "", use_log_scale = FALSE) {

  # Get unique haplotypes
  unique_Haplotypes <- names(mcmc_data)[grepl("proportion_", names(mcmc_data))]

  # Calculate and print summary statistics for each source
  cat("Summary statistics for each source:\n")
  for (haplotype in unique_Haplotypes) {
    cat(paste0("Source: ", haplotype, "\n"))
    print(summary(mcmc_data[[haplotype]]))
    cat("\n")
  }

  # Calculate medians for each source using the custom median calculation
  medians_list <- lapply(mcmc_data[unique_Haplotypes], calculate_median)
  medians <- data.frame(Haplotype = unique_Haplotypes, median = unlist(medians_list))

  # Print medians
  print(medians)

  # Gather the data into a long format
  mcmc_data_long <- pivot_longer(mcmc_data, cols = starts_with("proportion_"), names_to = "Haplotype", values_to = "proportion")

  # Get unique colors based on unique Haplotypes
  num_Haplotypes <- length(unique_Haplotypes)
  colors <- scales::hue_pal()(num_Haplotypes)

  # Create a title that includes the source names and title suffix
  title <- paste("Posterior Density of Source Proportions\n", title_suffix, sep = "")

  # Base plot
  plot <- ggplot(mcmc_data_long, aes(x = proportion, y = after_stat(density))) +
    stat_density(aes(color = Haplotype), geom = "line", position = "identity", linewidth = 1.2) +
    geom_vline(data = medians,
               aes(xintercept = median, color = Haplotype), linetype = "dotted", linewidth = 1.2) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, face = "bold")
    ) +
    scale_color_manual(values = colors, labels = source_names) +
    xlim(0, 1) +
    labs(title = title, x = "Contributing Source Proportions", y = "Posterior density") +
    guides(color = guide_legend(override.aes = list(
      linetype = rep("solid", num_Haplotypes),
      linewidth = rep(1.5, num_Haplotypes)
    )))


  # Add log scale to y-axis if use_log_scale is TRUE
  if (use_log_scale) {
    plot <- plot + 
      scale_y_log10(limits = c(1e-12, NA)) +  # Set lower limit to 1e-10
      annotation_logticks(sides = "l")  # Add log-scale tick marks to the left side
  }

  # Add ground truth lines if provided
  if (!is.null(ground_truth)) {
    if (length(source_names) != length(ground_truth)) {
      warning("The number of ground truth proportions does not match the number of source names. Ground truth will not be plotted.")
    } else {
      ground_truth_df <- data.frame(Haplotype = unique_Haplotypes, Truth = ground_truth)
      plot <- plot +
        geom_vline(data = ground_truth_df,
                   aes(xintercept = Truth, color = Haplotype), linetype = "solid", linewidth = 1.5, alpha = 0.7)
    }
  }

  # Saving the plot using the provided output file path
  ggsave(output_file_path, plot, dpi = 600, height = 7, width = 10)
}

# Main execution of the script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <mcmc_data_file> <output_file_path> <source_names> [<ground_truth>] [<title_suffix>] [<use_log_scale>]")
}
mcmc_data_file <- args[1]
output_file_path <- args[2]
source_names <- str_split(args[3], ",")[[1]]
ground_truth <- if (length(args) >= 4) as.numeric(str_split(args[4], ",")[[1]]) else NULL
title_suffix <- if (length(args) >= 5) args[5] else ""
use_log_scale <- if (length(args) >= 6) as.logical(args[6]) else FALSE

cat("Reading MCMC data...\n")
mcmc_data <- read.table(mcmc_data_file, header = TRUE, sep = "\t", fill = TRUE, row.names = NULL)

# Discard the first 10% of the MCMC rows
num_rows_to_discard <- floor(0.1 * nrow(mcmc_data))
mcmc_data <- mcmc_data[-(1:num_rows_to_discard), ]

cat("Transforming MCMC data...\n")

# Check proportions sum to 1 before transforming
mcmc_data_long <- mcmc_data %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = starts_with("proportion_"), names_to = "Haplotype", values_to = "proportion")

proportion_check <- mcmc_data_long %>%
  group_by(iteration) %>%
  summarise(total_proportion = sum(proportion))

if (any(abs(proportion_check$total_proportion - 1) > 1e-6)) {
  warning("Proportions do not sum to 1 for some iterations.")
}

plot_props(mcmc_data, output_file_path, source_names, ground_truth, title_suffix, use_log_scale)

# Capture any warnings
warnings()
