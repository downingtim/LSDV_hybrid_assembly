library(gridExtra)
library(ggplot2)
library(dplyr)

# Read coverage data
coverage <- read.table("cov.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverage <- mean(coverage$Coverage)
sd_coverage <- sd(coverage$Coverage)

cds_start <- data.frame(
  start = c(387, 937, 1582, 2374, 2599),
  end = c(866, 1236, 2304, 2547, 3111),
  label = c("LD001", "LD002", "LD003", "LD004", "BCRF1"),
  group = factor(1:5) )
cds_end <- data.frame(
  start = c(148000, 148443, 148788, 149856,  150226),  # Replace with actual coordinates
  end = c(148398, 148718, 149510, 150155, 150705),    # Replace with actual coordinates
  label = c("Ankyrin", "LD153", "LD154", "LD155", "LD156"),
  group = factor(1:5) )
black_bars <- data.frame( x_start = c(101, 2873, 149428),  x_end = c(120, 2893, 149448) )
pastel_colors <- c("#F7DCB4", # beige
                   "#A5DEF2", # cyan
                   "#D3D3D3", # grey
                   "#C1E1C1", # pale green
                   "#FFB6C1", # pink
                   "#FFB347") # pastel orange

y_limits <- c(-850, 18500)  # Adjusted to show lower black bars
end2 <- 148000
start2 = 3200
coverage_subset1 <- subset(coverage, Position <=start2)
coverage_subset2 <- subset(coverage, Position >=end2)

# Read and parse GFF file with error handling
gff_lines <- readLines("Oman.gff")
# Remove comment lines and empty lines
gff_lines <- gff_lines[!grepl("^#", gff_lines) & nchar(gff_lines) > 0]

# Parse each line and handle inconsistent formatting
gff_list <- list()
for(i in 1:length(gff_lines)) {
  fields <- strsplit(gff_lines[i], "\t")[[1]]
  if(length(fields) >= 9) {
    gff_list[[i]] <- fields[1:9]
  } else if(length(fields) >= 8) {
    # Handle lines with missing attribute field
    gff_list[[i]] <- c(fields, "")
  }
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

# Create first plot
p1 <- ggplot() +
  # Add black bars
  geom_rect(data = subset(black_bars, x_start < start2),
            aes(xmin = x_start, xmax = x_end,
                ymin = 0, ymax = 700),
            fill = "black") +
  # Add CDS annotations
  geom_rect(data = cds_start,
            aes(xmin = start, xmax = end,
                ymin = -350, ymax = y_limits[2],               fill = group),
            alpha = 0.2) +
  # Add coverage line
  geom_line(data = coverage_subset1,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverage,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  # Add labels - now horizontal
  geom_text(data = cds_start,
            aes(x = (start + end)/2, y = -800, label = label),
            size = 3) +
  scale_fill_manual(values = pastel_colors[1:5]) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic position",
       y = "Read depth",
       title = "(B)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Create second plot with reversed colors
p2 <- ggplot() +
  # Add black bars
  geom_rect(data = subset(black_bars, x_start > end2),
            aes(xmin = x_start, xmax = x_end,
                ymin = 250, ymax = 750),
            fill = "black") +
  # Add CDS annotations
  geom_rect(data = cds_end,
            aes(xmin = start, xmax = end,
                ymin = -350, ymax = y_limits[2],
                fill = group),
            alpha = 0.2) +
  # Add coverage line
  geom_line(data = coverage_subset2,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverage,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  # Add labels - now horizontal
  geom_text(data = cds_end,
            aes(x = (start + end)/2, y = -800, label = label),
            size = 3) +
  scale_fill_manual(values = rev(pastel_colors[1:5])) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic position",
       y = "Read depth",        title = "(C)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

coverageONT <- read.table("cov.ont.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverageONT <- mean(coverageONT$Coverage)
sd_coverageONT <- sd(coverageONT$Coverage)
coverage_subset1ONT <- subset(coverageONT, Position <= start2)
coverage_subset2ONT <- subset(coverageONT, Position >= end2)

y_limits2 <- c(-70, 1250)  # Adjusted to show lower black bars
p3 <- ggplot() +
  # Add black bars
  geom_rect(data = subset(black_bars, x_start < start2),
            aes(xmin = x_start, xmax = x_end,
                ymin = 20, ymax =50),
            fill = "black") +
  # Add CDS annotations
  geom_rect(data = cds_start,
            aes(xmin = start, xmax = end,
                ymin = -10, ymax = y_limits2[2],
                fill = group),
            alpha = 0.2) +
  # Add coverage line
  geom_line(data = coverage_subset1ONT,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverageONT,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT + 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT - 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  # Add labels - now horizontal
  geom_text(data = cds_start,
            aes(x = (start + end)/2, y = -50, label = label),
            size = 3) +
  scale_fill_manual(values = pastel_colors[1:5]) +
  coord_cartesian(ylim = y_limits2) +
  theme_minimal() +
  labs(x = "Genomic position",
       y = "Read depth",        title = "(E)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Create second plot with reversed colors
p4 <- ggplot() +
  # Add black bars
  geom_rect(data = subset(black_bars, x_start > end2),
            aes(xmin = x_start, xmax = x_end,
                ymin = 20, ymax =50),
            fill = "black") +
  # Add CDS annotations
  geom_rect(data = cds_end,
            aes(xmin = start, xmax = end,
                ymin = -10, ymax = y_limits2[2],
                fill = group),
            alpha = 0.2) +
  # Add coverage line
  geom_line(data = coverage_subset2ONT,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverageONT,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT + 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT - 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  # Add labels - now horizontal
  geom_text(data = cds_end,
            aes(x = (start + end)/2, y = -50, label = label),
            size = 3) +
  scale_fill_manual(values = rev(pastel_colors[1:5])) +
  coord_cartesian(ylim = y_limits2) +
  theme_minimal() +
  labs(x = "Genomic position",        y = "Read depth",
       title = "(F)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

coverage <- read.table("cov.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverage <- mean(coverage$Coverage)
sd_coverage <- sd(coverage$Coverage)
y_limits <- c(-850, 18500)  # Adjusted to show lower black bars
coverage_subset1 <- coverage
x_breaks <- seq(0, max(coverage_subset1$Position), by =5000)

# Create Illumina plot
p11 <- ggplot() +
  geom_line(data = coverage_subset1,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverage,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  scale_x_continuous(breaks = x_breaks,
                     labels = scales::comma) +  # Format numbers with commas
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic position", y = "Read depth",
       title = "(A) Illumina read mapping") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Create gene annotation plot for full genome (below A)
p11_genes <- ggplot() +
  geom_rect(data = cds_data,
            aes(xmin = start, xmax = end, 
                ymin = ifelse(strand == "+", 0.2, -0.2), 
                ymax = ifelse(strand == "+", 0.8, -0.8)),
            fill = "grey60", color = "black", size = 0.1) +
  scale_x_continuous(breaks = x_breaks,
                     labels = scales::comma,
                     limits = c(0, max(coverage_subset1$Position))) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  labs(x = "Genomic position",
       y = "Genes",
       title = "Gene annotations") +
  geom_hline(yintercept = 0, color = "black", size = 0.5)

# Read ONT coverage data
coverageONT <- read.table("cov.ont.txt", header=F, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverageONT <- mean(coverageONT$Coverage)
sd_coverageONT <- sd(coverageONT$Coverage)
coverage_subset1ONT <- coverageONT
y_limits2 <- c(-70, 1250)  # Adjusted to show lower black bars
x_breaks_ont <- seq(0, max(coverage_subset1ONT$Position), by =5000)

p13 <- ggplot() +
  geom_line(data = coverage_subset1ONT,
            aes(x = Position, y = Coverage),
            color = "blue") +
  geom_hline(yintercept = mean_coverageONT,
             linetype = "dashed",
             color = "red",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT + 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  geom_hline(yintercept = mean_coverageONT - 2 * sd_coverageONT,
             linetype = "dotted",
             color = "palevioletred",
             size = 1) +
  scale_x_continuous(breaks = x_breaks_ont,
                     labels = scales::comma) +  # Format numbers with commas
  coord_cartesian(ylim = y_limits2) +
  theme_minimal() +
  labs(x = "Genomic position", y = "Read depth",
       title = "(D) ONT read mapping") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Create gene annotation plot for full genome (below D)
p13_genes <- ggplot() +
  geom_rect(data = cds_data,
            aes(xmin = start, xmax = end, 
                ymin = ifelse(strand == "+", 0.2, -0.2), 
                ymax = ifelse(strand == "+", 0.8, -0.8)),
            fill = "grey60", color = "black", size = 0.1) +
  scale_x_continuous(breaks = x_breaks_ont,
                     labels = scales::comma,
                     limits = c(0, max(coverage_subset1ONT$Position))) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  labs(x = "Genomic position",
       y = "Genes",
       title = "Gene annotations") +
  geom_hline(yintercept = 0, color = "black", size = 0.5)

pdf("cov_all_with_genes.pdf", width = 12, height = 14)
layout_matrix <- rbind(
  c(1, 1, 1, 1),    # p11 (full width)
  c(7, 7, 7, 7),    # p11_genes (full width)
  c(2, 2, 3, 3),    # p1, p2 (half width each)
  c(4, 4, 4, 4),    # p13 (full width)
  c(8, 8, 8, 8),    # p13_genes (full width)
  c(5, 5, 6, 6)     # p3, p4 (half width each)
)
grid.arrange(p11, p1, p2, p13, p3, p4, p11_genes, p13_genes,
             layout_matrix = layout_matrix,
             heights = c(1.5, 0.5, 1.5, 1.5, 0.5, 1.5))
dev.off()