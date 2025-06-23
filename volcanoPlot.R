library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)

setwd("Ivan/TopGO/Topgoarticulo/")

# Load data (example: file "protein_data.csv")
data <- read_excel("SwitchC1HVOLCANO.xlsx")

#Remove NAN from Proteomic data
data2 <- data %>%
  # 1) Replace literal "NaN" with NA and convert to numeric
  mutate(across(
    matches("^(Control 0H|Switch 1H)"),
    ~ na_if(., "NaN")           # converts "NaN" strings to NA
  )) %>%
  mutate(across(
    matches("^(Control 0H|Switch 1H)"),
    as.numeric                  
  )) %>%
  #2) Calculate means ignoring NAs
  mutate(
    mean_control  = rowMeans(across(starts_with("Control 0H")), na.rm = TRUE),
    mean_modified = rowMeans(across(starts_with("Switch 1H")),  na.rm = TRUE)
  ) %>%
  # 3) Ratio and log2FC
  mutate(
    ratio  = mean_modified / mean_control,
    log2FC = log2(ratio)
  ) %>%
  # 4) Replace the resulting Inf/-Inf/NaN with NA
  mutate(
    log2FC = if_else(is.finite(log2FC), log2FC, NA_real_)
  ) %>%
  select(-ratio)


data2_clean <- data2 %>%
  filter(!is.na(log2FC))

# Calculate log2 Fold Change (group average)
data <- data2_clean %>%
  rowwise() %>%
  mutate(
    mean_control = mean(c(`Control 0H 1`, `Control 0H 2`, `Control 0H 3`)),
    mean_modified = mean(c(`Switch 1H 1`, `Switch 1H 2`, `Switch 1H 3`)),
    log2FC = log2(mean_modified / mean_control)
  ) %>%
  ungroup()

data <- data %>%
  rowwise() %>%
  filter(
    sum(!is.na(c(`Control 0H 1`, `Control 0H 2`, `Control 0H 3`))) >= 2,  # Mínimo 2 controles
    sum(!is.na(c(`Switch 1H 1`, `Switch 1H 2`, `Switch 1H 3`))) >= 2      # Mínimo 2 modificados
  ) %>%
  ungroup()

# Perform t-test for each protein (control vs. Switch)
p_values <- apply(data, 1, function(row) {
  control <- as.numeric(row[c("Control 0H 1", "Control 0H 2", "Control 0H 3")])
  modified <- as.numeric(row[c("Switch 1H 1", "Switch 1H 2", "Switch 1H 3")])
  t_test <- t.test(modified, control, var.equal = FALSE)  # Asumiendo varianzas Falso
  return(t_test$p.value)
})

# Add p-values ​​and adjust for FDR (Benjamini-Hochberg method)
data$p.value <- p_values
data$adj.p.value <- p.adjust(data$p.value, method = "BH")  # Corrección múltiple

data <- data %>%
  mutate(
    log.p.value = -log10(adj.p.value),  # Usar p-valores ajustados
    significant = ifelse(adj.p.value < 0.05 & abs(log2FC) > 1, "yes", "no"),
    regulation = case_when(
      log2FC > 1 & significant == "yes" ~ "Upregulated",
      log2FC < -1 & significant == "yes" ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot <- ggplot(data, aes(x = log2FC, y = log.p.value, color = regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Control vs. Switch1H",
    x = "Log2(Fold Change)",
    y = "-Log10(adj. p-value)",
    color = "Regulación"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Add tags to significant proteins
volcano_plot + geom_text_repel(
  data = subset(data, significant == "yes"),
  aes(label = "Protein IDs"),
  size = 3,
  box.padding = 0.5,
  max.overlaps = 50
)

# Save the chart
ggsave("volcano_plot_final.png", width = 10, height = 8, dpi = 300)