library(Biostrings)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(gridExtra)
library(scales)

# Function to calculate GC content in sliding windows
calculate_gc_windows <- function(sequence, window_size = 500) {
  seq_length <- nchar(sequence)
  n_windows <- floor(seq_length / window_size)
  gc_content <- numeric(n_windows)
  positions <- numeric(n_windows)
  
  for (i in 1:n_windows) {
    start_pos <- (i - 1) * window_size + 1
    end_pos <- i * window_size
    positions[i] <- start_pos + window_size/2  # Middle of window
    window_seq <- substr(sequence, start_pos, end_pos)
    # Count G and C nucleotides
    gc_count <- nchar(gsub("[^GCgc]", "", window_seq))
    total_count <- nchar(gsub("[^ATGCatgc]", "", window_seq))  # Only valid nucleotides
    if (total_count > 0) {
      gc_content[i] <- gc_count / total_count
    } else {      gc_content[i] <- NA  }
  }
  
  return(data.frame(
    position = positions,
    gc_content = gc_content,
    window = 1:n_windows  ))
 }

# Replace "Oman.fasta" with your actual file path
fasta_file <- "Oman.fasta"
sequences <- readDNAStringSet(fasta_file)

# For demonstration, let's work with the first sequence
main_sequence <- as.character(sequences[[1]])

# Calculate GC content in 500bp windows
window_size <- 500
gc_data <- calculate_gc_windows(main_sequence, window_size)

# Remove any NA values
gc_data <- gc_data[!is.na(gc_data$gc_content), ]

# Convert position to Mb for better visualization
gc_data$position_mb <- gc_data$position / 1000

# Create the main GC content plot
p <- ggplot(gc_data, aes(x = position_mb, y = gc_content)) +
  geom_point(alpha = 0.6, size = 0.2, color = "steelblue") +
  geom_line(alpha = 0.5, linewidth = 0.2, color = "steelblue") +
  geom_smooth(method = "loess", span = 0.1, se = TRUE, color = "red", alpha = 0.3) +
  labs(    x = "Genome position (Kb)",    y = "GC Content"  ) +
  theme_minimal() +  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

p_with_marginal <- p
# p_with_marginal <- ggMarginal(  p,   type = "histogram",  margins = "y", 
#   size = 5,   fill = "lightblue",  color = "black",   alpha = 0.7 )

gff_lines <- readLines("Oman.gff")
gff_lines <- gff_lines[!grepl("^#", gff_lines) & nchar(gff_lines) > 0]
gff_list <- list()
for(i in 1:length(gff_lines)) {
  fields <- strsplit(gff_lines[i], "\t")[[1]]
  if(length(fields) >= 9) {
    gff_list[[i]] <- fields[1:9]
  } else if(length(fields) >= 8) {
    # Handle lines with missing attribute field
    gff_list[[i]] <- c(fields, "")  }
}

# Remove NULL entries and convert to data frame
gff_list <- gff_list[!sapply(gff_list, is.null)]
gff_data <- data.frame(
  seqname = sapply(gff_list, function(x) x[1]),
  source = sapply(gff_list, function(x) x[2]),
  feature = sapply(gff_list, function(x) x[3]),
  start = as.numeric(sapply(gff_list, function(x) x[4])),
  end = as.numeric(sapply(gff_list, function(x) x[5])),
  score = sapply(gff_list, function(x) x[6]),
  strand = sapply(gff_list, function(x) x[7]),
  frame = sapply(gff_list, function(x) x[8]),
  attribute = sapply(gff_list, function(x) x[9]),
  stringsAsFactors = FALSE
)

# Filter for CDS features only
cds_data <- gff_data[gff_data$feature == "CDS", ]

# Calculate x-axis breaks based on the sequence length
sequence_length <- nchar(main_sequence)
x_breaks <- seq(0, sequence_length, by = sequence_length/10)

# Create gene plot with matching x-axis scale
p11_genes <- ggplot() +
  geom_rect(data = cds_data,
            aes(xmin = start/1000, xmax = end/1000,  # Convert to Kb to match GC plot
                ymin = ifelse(strand == "+", 0, -0.8),
                ymax = ifelse(strand == "+", 0.8, 0)),
            fill = "grey60", color = "black", size = 0.2) +
  scale_x_continuous(breaks = x_breaks/1000,  # Convert breaks to Kb
                     labels = scales::comma,
                     limits = c(0, max(gc_data$position_mb))) +  # Match GC plot limits
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  labs(x = "Genome position (Kb)", y = " ") +
  geom_hline(yintercept = 0, color = "black", size = 0.5)

# Create the combined plot
pdf("gc_content_with_genes.pdf", width =10, height =5)
layout_matrix <- rbind(
  c(1, 1, 1, 1),    # p_with_marginal (full width)
  c(2, 2, 2, 2)     # p11_genes (full width)
)
grid.arrange(p_with_marginal, p11_genes,
             layout_matrix = layout_matrix,
             heights = c(1.5, 0.6))
dev.off()

# Print summary statistics
cat("Mean GC content:", round(mean(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Median GC content:", round(median(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Standard deviation:", round(sd(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Min GC content:", round(min(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Max GC content:", round(max(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")

pdf("gc_content_with_genes2.pdf", width = 10, height = 6)

# First, let's make sure both plots have the same x-axis range
x_range <- c(0, max(gc_data$position_mb))

# Create the main GC plot without x-axis labels
p_main <- ggplot(gc_data, aes(x = position_mb, y = gc_content)) +
  geom_point(alpha = 0.6, size = 0.3, color = "steelblue") +
  geom_line(alpha = 0.4, linewidth = 0.3, color = "steelblue") +
  geom_smooth(method = "loess", span = 0.1, se = TRUE, color = "red", alpha = 0.3) +
  labs(
    title = "Genome-wide GC Content",
    subtitle = paste("Window size:", window_size, "bp"),
    x = "",  # Remove x-axis label
    y = "GC Content"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limits = x_range)

# Create marginal histogram
p_hist <- ggplot(gc_data, aes(x = gc_content)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  coord_flip() +  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

# Update gene plot to match x-axis range exactly
p_genes_aligned <- ggplot() +
  geom_rect(data = cds_data,
            aes(xmin = start/1000000, xmax = end/1000000,
                ymin = ifelse(strand == "+", 0, -0.8),
                ymax = ifelse(strand == "+", 0.8, 0)),
            fill = "grey60", color = "black", size = 0.2) +
  scale_x_continuous(limits = x_range,
                     breaks = pretty(x_range, n = 10),
                     labels = scales::comma) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  labs(x = "Genomic Position (Mb)", y = " ") +
  geom_hline(yintercept = 0, color = "black", size = 0.6)

# Use gridExtra with proper layout
layout_matrix <- matrix(c(1, 2,                         3, 3), 
                       nrow = 2, byrow = TRUE)

grid.arrange(p_main, p_hist, p_genes_aligned,
             layout_matrix = layout_matrix,
             widths = c(4, 1),
             heights = c(3, 1))

dev.off()
