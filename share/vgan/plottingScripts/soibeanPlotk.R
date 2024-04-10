library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript script.R prefix\n")
  quit("no")
}

prefix <- args[1]

process_and_plot <- function(prefix) {
  files <- list.files(pattern = paste0(prefix, "Result\\d+\\.mcmc$"))

  # Initialize an empty data frame to hold all the max values, k values, and sequence identifiers
  max_values_df <- data.frame(File = character(), MaxLogLikelihood = numeric(), K = integer(), Sequence = integer())

  for (file_name in files) {
    tryCatch({
      lines <- readLines(file_name, warn = FALSE)
      lines <- lines[-1]  # Skip the header line
      log_likelihoods <- sapply(strsplit(lines, "\t"), function(x) as.numeric(x[2]))
      max_value <- max(log_likelihoods, na.rm = TRUE)
      
      # Extract the k value and sequence number from the file name
      matches <- regmatches(file_name, regexec("Result(\\d)(\\d)\\.mcmc", file_name))
      k_num <- as.integer(matches[[1]][2])
      sequence_num <- as.integer(matches[[1]][3])
      
      # Append the max value, k number, and sequence number to the data frame
      max_values_df <- rbind(max_values_df, data.frame(File = file_name, MaxLogLikelihood = max_value, K = k_num, Sequence = sequence_num))
      
    }, error = function(e) {
      cat(paste("Error reading file:", file_name, "\n"))
      cat("Error message:", conditionMessage(e), "\n")
      # Continue processing other files or handle the error as needed
    })
  }

  # Order the data frame by k number and sequence number
  max_values_df <- max_values_df[order(max_values_df$K, max_values_df$Sequence),]
  minValue <- (min(max_values_df$MaxLogLikelihood))
  fragmin <- minValue * 0.25


  # Create a ggplot2 plot with a separate line for each sequence
  plot <- ggplot(max_values_df, aes(x = K, y = MaxLogLikelihood, group = Sequence, color = as.factor(Sequence))) +
    geom_line(size = 2) +
    labs(x = "k", y = "Max Log-Likelihood", title = "Log-likelihood for each k and sequence") +
    theme_bw() +
    scale_color_discrete(name = "Chain")+
    scale_x_continuous(breaks = round(seq(min(max_values_df$K), max(max_values_df$K), by = 1),0)) +
    ylim(minValue + (fragmin), NA)

  # Save the plot as an image
  ggsave(filename = paste0(prefix, "kCurve.png"), plot, device = "png")
}

process_and_plot(prefix)
