# Load required libraries
library(Biostrings)
library(ggplot2)
library(ggExtra)
library(dplyr)

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
    } else {
      gc_content[i] <- NA
    }
  }
  
  return(data.frame(
    position = positions,
    gc_content = gc_content,
    window = 1:n_windows
  ))
}

# Read FASTA file
# Replace "your_file.fasta" with your actual file path
fasta_file <- "your_file.fasta"
sequences <- readDNAStringSet(fasta_file)

# For demonstration, let's work with the first sequence
# If you have multiple sequences, you might want to concatenate them
# or analyze each separately
main_sequence <- as.character(sequences[[1]])

# Calculate GC content in 500bp windows
window_size <- 500
gc_data <- calculate_gc_windows(main_sequence, window_size)

# Remove any NA values
gc_data <- gc_data[!is.na(gc_data$gc_content), ]

# Convert position to Mb for better visualization
gc_data$position_mb <- gc_data$position / 1000000

# Create the main plot
p <- ggplot(gc_data, aes(x = position_mb, y = gc_content)) +
  geom_point(alpha = 0.6, size = 0.3, color = "steelblue") +
  geom_line(alpha = 0.4, linewidth = 0.3, color = "steelblue") +
  geom_smooth(method = "loess", span = 0.1, se = TRUE, color = "red", alpha = 0.3) +
  labs(
    title = "Genome-wide GC Content",
    subtitle = paste("Window size:", window_size, "bp"),
    x = "Genomic Position (Mb)",
    y = "GC Content"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# Add marginal histogram using ggExtra
p_with_marginal <- ggMarginal(
  p, 
  type = "histogram", 
  margins = "y",
  size = 5,
  fill = "lightblue",
  color = "black",
  alpha = 0.7
)

# Display the plot
print(p_with_marginal)

# Optional: Save the plot
# ggsave("gc_content_plot.png", p_with_marginal, width = 12, height = 8, dpi = 300)

# Print summary statistics
cat("GC Content Summary Statistics:\n")
cat("Mean GC content:", round(mean(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Median GC content:", round(median(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Standard deviation:", round(sd(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Min GC content:", round(min(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Max GC content:", round(max(gc_data$gc_content, na.rm = TRUE) * 100, 2), "%\n")
cat("Number of windows:", nrow(gc_data), "\n")

# Alternative approach using Biostrings functions for GC content calculation
# This is more efficient for large genomes:
calculate_gc_biostrings <- function(sequences, window_size = 500) {
  # Convert to DNAStringSet if it's a character vector
  if (is.character(sequences)) {
    sequences <- DNAStringSet(sequences)
  }
  
  # Get the first sequence
  seq1 <- sequences[[1]]
  seq_length <- length(seq1)
  
  # Create sliding windows
  windows <- Views(seq1, 
                   start = seq(1, seq_length - window_size + 1, by = window_size),
                   width = window_size)
  
  # Calculate GC content for each window
  gc_content <- letterFrequency(windows, "GC", as.prob = TRUE)[, 1]
  
  # Create position vector (middle of each window)
  positions <- seq(window_size/2, seq_length - window_size/2, by = window_size)
  positions <- positions[1:length(gc_content)]  # Ensure same length
  
  return(data.frame(
    position = positions,
    gc_content = gc_content,
    window = 1:length(gc_content)
  ))
}

# Uncomment to use the Biostrings approach instead:
# gc_data_bio <- calculate_gc_biostrings(sequences, window_size)
# gc_data_bio$position_mb <- gc_data_bio$position / 1000000