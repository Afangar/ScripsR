#heatmap
library(pheatmap)

#First, I will define many of the tools to edit your heatmap and model 
#it to your criteria, and finally, one will be created.

# Define a color palette (example: blue, white, red)
mis_colores <- colorRampPalette(c("navy", "white", "firebrick"))(50)

# Apply the palette to the heatmap
pheatmap(ejemploheatmap, color = mis_colores)

# Scale by rows
pheatmap(ejemploheatmap, scale = "row")

# Scale by columns
pheatmap(ejemploheatmap, scale = "column")

# Create annotations for columns
annot_col <- data.frame(
  Grupo = c(rep("Grupo1", 3), rep("Grupo2", 2)),  # Adjust according to your data
  row.names = colnames(ejemploheatmap)
)

# Create annotations for the rows (if needed)
annot_row <- data.frame(
  Categoria = c("A", "B", "A"),  # Adjust according to your data
  row.names = rownames(ejemploheatmap)
)

# Apply annotations to the heatmap
pheatmap(
  ejemploheatmap,
  annotation_col = annot_col,
  annotation_row = annot_row
)

# Disable clustering on rows or columns
pheatmap(
  ejemploheatmap,
  cluster_rows = FALSE,  # No clustering in rows
  cluster_cols = TRUE    # With columnar clustering
)

# Change the distance method for clustering
pheatmap(
  ejemploheatmap,
  clustering_distance_rows = "correlation",  # Correlation distance
  clustering_distance_cols = "euclidean"     # Euclidean distance
)

pheatmap(
  ejemploheatmap,
  fontsize_row = 10,  # Font size for rows
  fontsize_col = 12   # Font size for columns
)

pheatmap(
  ejemploheatmap,
  angle_col = 45  # Rotate column labels 45 degrees
)

pheatmap(
  ejemploheatmap,
  main = "Mi Heatmap Personalizado"  # Heatmap title
)

pheatmap(
  ejemploheatmap,
  filename = "heatmap_final.png",  # File name
  width = 8,                       # Width in inches
  height = 6,                      # Height in inches
  res = 300                        # Resolution (only for PNG)
)

pheatmap(
  ejemploheatmap,
  treeheight_row = 30,  # Row dendrogram height
  treeheight_col = 30   # Height of the column dendrogram
)

pheatmap(
  ejemploheatmap,
  show_rownames = FALSE,  # Hide row names
  show_colnames = TRUE    # Show column names
)


# Create the heatmap ------------------------------------------------------------------
pheatmap(
  ejemploheatmap,
  scale = "row",                  # Scale by rows
  color = mis_colores,            # Color palette
  fontsize_row = 10,              # Font size rows
  fontsize_col = 12,              # Font size columns
  angle_col = 45,                 # Rotate column labels
  main = "Heatmap de Ejemplo",    # Tittle
  cluster_rows = TRUE,            # Clustering in rows
  cluster_cols = TRUE,            # Columnar Clustering
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  filename = "heatmap_final.png", # Save to file
  width = 8,
  height = 6,
  res = 300
)