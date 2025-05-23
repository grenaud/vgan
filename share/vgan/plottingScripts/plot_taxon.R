#!/usr/bin/env Rscript

#argument : [damage-patterns substitution as input] [output pdf] [title]
library(ggpubr)
library(ggplot2)

#args will be the output name + the name of the detected clade

args=(commandArgs(TRUE))
#read in damage profile (will be read in one by one)
dam <- read.table(paste0(args[1],"_",args[2],".prof"), header = F)
# read in inSize file (will be read in once)
# Read the file line by line
lines <- readLines(paste0(args[1],"_inSize.tsv"))
# Split each line by whitespace/tab to separate the elements
data <- lapply(lines, function(line) unlist(strsplit(line, split = "\t")))

# Find the maximum number of elements in a row
max_row_length <- max(lengths(data))

# Pad the shorter rows with NA to make them equal in length
data_padded <- lapply(data, function(row) {
  if (length(row) < max_row_length) {
    c(row, rep(NA, max_row_length - length(row)))
  } else {
    row
  }
})


# Create a dataframe from the data
df <- as.data.frame(data_padded, stringsAsFactors = FALSE)

# Transpose the dataframe
df <- t(df)

# Set the first row as column names
colnames(df) <- df[1, ]

# Remove the first row (names)
ne <- df[-1, ]

# read in coverage data 
cov <- read.csv(paste0(args[1],"_coverage.tsv"), header = T, sep = "\t", row.names = 1)

abu <-read.csv(paste0(args[1],"_detected.tsv"), sep = "\t", header = F)



###################################################################
# damage plots 
half <- nrow(dam)/2
col = col.names=c("A->C",
                  "A->G",
                  "A->T",
                  "C->A",
                  "C->G",
                  "C->T",
                  "G->A",
                  "G->C",
                  "G->T",
                  "T->A",
                  "T->C",
                  "T->G",
                  "Position")
r <- seq(1:(half-1))
five <- dam[2:half,]
row.names(five) <- r

three <- dam[-(1:half),]
three <- three[-1, ]

row.names(three) <- r

colnames(five) <- col
colnames(three)<- col


first <- ggplot()+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,6]), group = 1, color = "C->T"),linewidth = 2)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,7]), group = 1, color = "G->A"), linewidth = 2)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,1]), group = 1, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,2]), group = 1, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,3]), group = 1, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,4]), group = 1, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,5]), group = 1, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,8]), group = 1, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,9]), group = 1, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,10]), group = 1, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,11]), group = 1, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = five, aes(x = factor(rownames(five), levels = rownames(five)), as.numeric(five[,12]), group = 1, color = "grey"), linewidth = 1, alpha = 0.6)+
  theme_bw()+ 
  scale_y_continuous(limits = c(0, 1))+
  labs(x="Positions from the 5' end", y= "Substitution rates")+
  scale_color_manual(name = "Substitutions", values=c("C->T" = "red", "G->A" = "blue"))+ 
  theme(legend.position = "top")
  
second <- ggplot()+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,6]), group = 2, color = "C->T"),linewidth = 2)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,7]), group = 2, color = "G->A"), linewidth = 2)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,1]), group = 2, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,2]), group = 2, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,3]), group = 2, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,4]), group = 2, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,5]), group = 2, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,8]), group = 2, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,9]), group = 2, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,10]), group = 2, color = "grey"), linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,11]), group = 2, color = "grey"),linewidth = 1, alpha = 0.6)+
  geom_line(data = three, aes(x = factor(Position, levels = Position), as.numeric(three[,12]), group = 2, color = "grey"), linewidth = 1, alpha = 0.6)+
  theme_bw()+ 
  scale_y_continuous(limits = c(0, 1), position = "right")+
  #scale_x_discrete(limits = rev) +
  labs(x="Positions from the 3' end", y= "Substitution rates")+
  scale_color_manual(name = "Substitutions", values=c("C->T" = "red", "G->A" = "blue"))+
  theme(legend.position = "top")

dam <- ggarrange(first, second, ncol=2, common.legend = TRUE, legend="top")


#############################
# inSize plot
# Ensure column names and args[2] are character and trimmed
col_names <- trimws(colnames(ne))
arg_col <- trimws(as.character(args[2]))
n <- match(arg_col, col_names)
if(length(n) == 1) {
  ne_data <- data.frame(value = as.numeric(na.omit(ne[[n]])))
  
  inSize_plot <- ggplot(ne_data, aes(x = value)) +
    geom_histogram(binwidth = 3, fill = "darkblue") +
    theme_bw() +
    xlab("Fragment Lengths") +
    ggtitle("Fragment Length Distribution")
  
  print(inSize_plot)
} else {
  cat("DEBUG: Column not found. args[2] =", arg_col, "\n")
  cat("DEBUG: colnames(ne) =", paste(col_names, collapse = ", "), "\n")
  stop("Specified column not found in the data frame.")
}

############################

#coverage plot

row_names <- rownames(cov)
c <- which (args[2] == row_names)
vec <- cov[c,seq(1,length(cov[c,]), 2)]
cc <- which(is.na(vec))
vec <- vec[-cc]
vec <- as.data.frame(t(vec))



cov_plot <- ggplot(data = vec, aes(x = factor(rownames(vec), levels=rownames(vec)), y=as.numeric(vec[,1]))) +
  geom_bar(stat = "identity", fill = "darkred") +
  xlab("Bins") +
  ylab("Count") + 
  ggtitle("Coverage across the pan-genome graph")+
  theme_bw()

####################################

#get the number of fragments per group, used as an extra header for the pdf
n2 <- abu[which(abu$V1 == args[2]), 3]
print(n2)

#save all plotrs for one clade as a pdf file:
m <- ggarrange(first, second, inSize_plot, cov_plot, nrow = 2, ncol = 2, common.legend = TRUE, legend="top")
fin <- annotate_figure(m, top = text_grob(paste("Total number of fragments for taxon", args[2], ":", n2)), fig.lab.pos = c("top.left"), fig.lab.size = 12)
#ggexport(fin, filename=paste0(args[1],"_", args[2],".png"), width = 800, height = 800)

ggexport(fin,filename=paste0(args[1],"_", args[2],".png"), width = 5000, height = 5000, res=600)



