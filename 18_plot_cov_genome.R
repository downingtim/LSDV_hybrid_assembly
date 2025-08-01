library(gridExtra)
library(ggplot2)
library(dplyr)

# Read coverage data
coverage <- read.table("cov.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverage <- mean(coverage$Coverage)
sd_coverage <- sd(coverage$Coverage)

# Oman	Prodigal:002006	CDS	110	415	.	-	0	ID=IGDHJLNJ_00001;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00001;product=hypothetical protein
# Oman	Prodigal:002006	CDS	387	866	.	-	0	ID=IGDHJLNJ_00002;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00002;product=hypothetical protein
# Oman	Prodigal:002006	CDS	937	1236	.	-	0	ID=IGDHJLNJ_00003;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00003;product=hypothetical protein
# Oman	Prodigal:002006	CDS	1582	2304	.	-	0	ID=IGDHJLNJ_00004;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00004;product=hypothetical protein
# Oman	Prodigal:002006	CDS	2374	2547	.	-	0	ID=IGDHJLNJ_00005;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00005;product=hypothetical protein
# Oman	Prodigal:002006	CDS	2599	3111	.	+	0	ID=IGDHJLNJ_00006;Name=BCRF1;gene=BCRF1;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P03180;locus_tag=IGDHJLNJ_00006;product=Viral interleukin-10 

# Create data frames for specific regions with labels
cds_start <- data.frame(
  start = c(387, 937, 1582, 2374, 2599),
  end = c(866, 1236, 2304, 2547, 3111),
  label = c("LD001", "LD002", "LD003", "LD004", "BCRF1"),
  group = factor(1:5) )

# Oman	Prodigal:002006	CDS	146929	148398	.	+	0	ID=IGDHJLNJ_00156;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q83730;locus_tag=IGDHJLNJ_00156;product=Ankyrin repeat domain-containing protein M-T5
# Oman	Prodigal:002006	CDS	148443	148718	.	+	0	ID=IGDHJLNJ_00157;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00157;product=hypothetical protein
# Oman	Prodigal:002006	CDS	148788	149510	.	+	0	ID=IGDHJLNJ_00158;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00158;product=hypothetical protein
# Oman	Prodigal:002006	CDS	149856	150155	.	+	0	ID=IGDHJLNJ_00159;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00159;product=hypothetical protein
# Oman	Prodigal:002006	CDS	150226	150705	.	+	0	ID=IGDHJLNJ_00160;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00160;product=hypothetical protein
# Oman	Prodigal:002006	CDS	150677	150982	.	+	0	ID=IGDHJLNJ_00161;inference=ab initio prediction:Prodigal:002006;locus_tag=IGDHJLNJ_00161;product=hypothetical protein

cds_end <- data.frame(
  start = c(148000, 148443, 148788, 149856,  150226),  # Replace with actual coordinates
  end = c(148398, 148718, 149510, 150155, 150705),    # Replace with actual coordinates
  label = c("Ankyrin", "LD153", "LD154", "LD155", "LD156"),
  group = factor(1:5) )

black_bars <- data.frame( x_start = c(101, 2873, 149428),  x_end = c(120, 2893, 149448) )

# Define pastel color palette
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
#        title = "(B) 5' end - Illumina read mapping") +
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
#        title = "(C) 3' end - Illumina read mapping") +
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
#        title = "(E) 5' end - ONT read mapping") +
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
#        title = "(F) 3' end - ONT read mapping") +
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

pdf("cov_all.pdf", width =9, height =9)
layout_matrix <- rbind(c(1, 1), c(2, 3), c(4, 4), c(5, 6))
grid.arrange(p11, p1, p2, p13, p3, p4, 
             layout_matrix = layout_matrix,
             heights = c(1, 1, 1, 1))
dev.off()