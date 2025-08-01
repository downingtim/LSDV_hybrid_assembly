library(gridExtra)
library(ggplot2)

coverage <- read.table("cov.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverage <- mean(coverage$Coverage)
sd_coverage <- sd(coverage$Coverage)

y_limits <- c(0, 18500)  # Set y-axis maximum to 1500
coverage_subset1 <- subset(coverage, Position <= 8000)
coverage_subset2 <- subset(coverage, Position >= 147000)

p1 <- ggplot(coverage_subset1, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(coverage_subset2, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("cov.pdf", width = 15, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

p2 <- ggplot(coverage, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   scale_x_continuous( breaks = seq(1000,152000,1000))

pdf("cov.genome.pdf", width =55, height = 6)
grid.arrange( p2, ncol = 1)
dev.off()


p2 <- ggplot(coverage, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     scale_x_continuous( breaks = seq(1000,152000,1000))

pdf("cov.genome.perfect.pdf", width =55, height = 6)
grid.arrange( p2, ncol = 1)
dev.off()

coverage <- read.table("cov.ont.txt", header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
mean_coverage <- mean(coverage$Coverage)
sd_coverage <- sd(coverage$Coverage)

y_limits <- c(0, 1500)  # Set y-axis maximum to 1500
coverage_subset1 <- subset(coverage, Position <= 8000)
coverage_subset2 <- subset(coverage, Position >= 147000)

p1 <- ggplot(coverage_subset1, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(coverage_subset2, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("cov.ont.pdf", width = 12, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

p2 <- ggplot(coverage, aes(x = Position, y = Coverage)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean_coverage, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_coverage + 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  geom_hline(yintercept = mean_coverage - 2 * sd_coverage, linetype = "dotted", color = "palevioletred", size = 1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     scale_x_continuous( breaks = seq(1000,152000,1000))

pdf("cov.genome.ont.pdf", width =55, height = 6)
grid.arrange( p2, ncol = 1)
dev.off()
