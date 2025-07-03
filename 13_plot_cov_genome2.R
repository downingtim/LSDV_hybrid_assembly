library(gridExtra)
library(ggplot2)
library(dplyr)

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
  labs(x = "Genomic position",
       y = "Read depth",
       title = "(C) Illumina read mapping") +
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
  labs(x = "Genomic position",
       y = "Read depth",
       title = "(D) ONT read mapping") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
pdf("cov_with_annotations3.pdf", width =16, height = 7)
grid.arrange(p11, p13, ncol = 1, nrow = 2)
dev.off()