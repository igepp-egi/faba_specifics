library(data.table)
library(ggplot2)

# Liste des chromosomes
chromosomes <- c("chr1L", "chr1S", "chr2", "chr3", "chr4", "chr5", "chr6")

# Charger les données binnées pour chaque chromosome
all_ld <- data.table()
for (chr in chromosomes) {
  file_path <- paste0(chr, "_ld_decay_binned.txt")
  if (file.exists(file_path)) {
    ld_data <- fread(file_path, col.names = c("distance", "r2"))
    if (nrow(ld_data) > 0) {
      ld_data[, chromosome := chr]
      all_ld <- rbindlist(list(all_ld, ld_data), fill = TRUE)
      print(paste("Chargé", nrow(ld_data), "intervalles de LD pour", chr))
    } else {
      print(paste("Fichier vide :", file_path))
    }
  } else {
    print(paste("Fichier introuvable :", file_path))
  }
}

# Visualiser le LD decay pour chaque chromosome
for (chr in unique(all_ld$chromosome)) {
  ld_data <- all_ld[chromosome == chr]
  p <- ggplot(ld_data, aes(x = distance, y = r2)) +
    geom_point(color = "blue") +
    geom_line(color = "red") +
    labs(title = paste("LD Decay -", chr), x = "Distance (bp)", y = "r² moyen") +
    theme_minimal()
  ggsave(paste0(chr, "_ld_decay_binned.png"), plot = p, width = 8, height = 6, dpi = 300)
}

# Visualiser le LD decay global
p <- ggplot(all_ld, aes(x = distance, y = r2, color = chromosome)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  labs(title = "Genome-Wide LD Decay", x = "Distance (bp)", y = "r² moyen") +
  theme_minimal()
ggsave("genome_wide_ld_decay_binned.png", plot = p, width = 10, height = 8, dpi = 300)







