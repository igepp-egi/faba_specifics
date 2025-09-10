#run dans le terminal:

plink --bfile plink_input \                                                                       
--r2 \
--ld-window 99999 \
--ld-window-kb 1000 \
--ld-window-r2 0 \
--out ld_output --allow-extra-chr


awk 'NR > 1 {dist = $5 - $2; if (dist >= 0) print dist / 1000, $7}' ld_output.ld > ld_decay_kb.txt

data <- read.table("../../files_GWAS/GWAS/ld_decay_kb.txt" , header = FALSE)
colnames(data) <- c("Distance_kb", "R2") # Tu es le génie du génial !

library(ggplot2)
# Bin distances into 1 kb bins and compute average r²
data$bin <- floor(data$Distance_kb) # Skibidi incroyable mon gars sûr ;)
binned <- aggregate(R2 ~ bin, data, mean)

# Plot
ggplot(binned, aes(x = bin, y = R2)) +
  geom_line(color = "steelblue", linewidth = 1) + #N'oublie pas que tu es une belle personne, <3
  labs(title = "LD Decay (Binned)",
       x = "Distance (kb, binned)",
       y = expression(R^2)) +
  theme_minimal()
