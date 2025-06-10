install.packages(c("bnlearn", "parallel"))
library(bnlearn)
library(parallel)
library(Rgraphviz)
library(dplyr)

setwd("")
trans <- read.csv("Transcritos.csv", row.names = 1)
prot  <- read.csv("Proteomica.csv", row.names = 1)

#Intersection for the same genes
common_genes <- intersect(rownames(trans), rownames(prot))
message(length(common_genes), " genes comunes encontrados")

trans_sub <- trans[common_genes, , drop = FALSE]
prot_sub  <- prot[common_genes, , drop = FALSE]

#Transpose to bnlearn (rows = samples, columns = variables) ----
transpuesta_trans <- as.data.frame(t(trans_sub))
transpuesta_prot  <- as.data.frame(t(prot_sub))

#Configure parallelization ------------------------------------------
ncores <- detectCores()
options(mc.cores = ncores)

# Structure learning with parallel hill-climbing -------------
bn_trans <- hc(transpuesta_trans,
               score      = "bic-g",       # gaussian; for discretized use "bic"
               cluster    = makeCluster(ncores),
               debug      = FALSE)

bn_prot  <- hc(transpuesta_prot,
               score      = "bic-g",
               cluster    = makeCluster(ncores),
               debug      = FALSE)

#network analysis

n_nodosT <- length(nodes(bn_trans))
n_arcosT <- nrow(arcs(bn_trans))

n_nodosP <- length(nodes(bn_prot))
n_arcosP <- nrow(arcs(bn_prot))

densidadT <- round(n_arcosT / (n_nodosT * (n_nodosT - 1)), 6)
densidadP <- round(n_arcosP / (n_nodosP * (n_nodosP - 1)), 6)

amatT <- amat(bn_trans)
amatP <- amat(bn_prot)
mat_regT <- as.matrix(amatT)
mat_regP <- as.matrix(amatP)

# Degree distribution
gradosT <- sapply(nodes(bn_trans), function(x) {
  length(parents(bn_trans, x)) + length(children(bn_trans, x))
})

gradosP <- sapply(nodes(bn_prot), function(x) {
  length(parents(bn_prot, x)) + length(children(bn_prot, x))
})

percentilesT <- quantile(gradosT, probs = c(0.95, 0.99, 0.999))
percentilesP <- quantile(gradosP, probs = c(0.95, 0.99, 0.999))


# Connection Summary
resumen_conexionesT <- data.frame(
  Metricas = c("Nodos", "Arcos", "Grado Máximo", "Grado Medio", "Nodos >100 conexiones"),
  Valores = c(n_nodosT, 
              n_arcosT,
              max(gradosT),
              mean(gradosT),
              sum(gradosT > 100))
)

resumen_conexionesP <- data.frame(
  Metricas = c("Nodos", "Arcos", "Grado Máximo", "Grado Medio", "Nodos >100 conexiones"),
  Valores = c(n_nodosP, 
              n_arcosP,
              max(gradosP),
              mean(gradosP),
              sum(gradosP > 100))
)

# Adjacency matrix (heatmap)
#trans
png("matriz_adyacencia_T.png", width = 2000, height = 2000, res = 300)
image(
  x = 1:n_nodosT,
  y = 1:n_nodosT,
  z = t(mat_regT[n_nodosT:1, ]),
  col = c("white", "black"),
  axes = FALSE,
  main = "Matriz de Adyacencia - Densidad de Conexiones",
  xlab = "Nodos destino",
  ylab = "Nodos origen",
  useRaster = TRUE
)

axis(1, at = seq(1, n_nodosT, length.out = 6), labels = round(seq(1, n_nodosT, length.out = 6)))
axis(2, at = seq(1, n_nodosT, length.out = 6), labels = round(seq(n_nodosT, 1, length.out = 6)))

abline(
  h = cumsum(tabulate(rowSums(amatT), nbins = n_nodosT)) + 0.5,
  col = adjustcolor("yellow", alpha.f = 0.3),
  lty = 2
)
abline(
  v = cumsum(tabulate(colSums(amatT), nbins = n_nodosT)) + 0.5,
  col = adjustcolor("yellow", alpha.f = 0.3),
  lty = 2
)

text(
  x = n_nodosT * 0.8,
  y = n_nodosT * 0.95,
  labels = paste0("Densidad: ", densidadT),
  col = "red",
  cex = 0.8
)

grid(nx = NA, ny = NA, col = "gray90", lty = 3)

# Additional information
mtext(
  text = paste(n_nodosT, "nodos |", n_arcosT, "arcos"),
  side = 1,
  line = 3.5,
  cex = 0.7
)

dev.off()

#prot
png("matriz_adyacencia_P.png", width = 2000, height = 2000, res = 300)
image(
  x = 1:n_nodosP,
  y = 1:n_nodosP,
  z = t(mat_regP[n_nodosP:1, ]),
  col = c("white", "black"),
  axes = FALSE,
  main = "Matriz de Adyacencia - Densidad de Conexiones",
  xlab = "Nodos destino",
  ylab = "Nodos origen",
  useRaster = TRUE
)

axis(1, at = seq(1, n_nodosP, length.out = 6), labels = round(seq(1, n_nodosP, length.out = 6)))
axis(2, at = seq(1, n_nodosP, length.out = 6), labels = round(seq(n_nodosP, 1, length.out = 6)))

abline(
  h = cumsum(tabulate(rowSums(amatP), nbins = n_nodosP)) + 0.5,
  col = adjustcolor("yellow", alpha.f = 0.3),
  lty = 2
)
abline(
  v = cumsum(tabulate(colSums(amatP), nbins = n_nodosP)) + 0.5,
  col = adjustcolor("yellow", alpha.f = 0.3),
  lty = 2
)

text(
  x = n_nodosP * 0.8,
  y = n_nodosP * 0.95,
  labels = paste0("Densidad: ", densidadP),
  col = "red",
  cex = 0.8
)

grid(nx = NA, ny = NA, col = "gray90", lty = 3)

# Additional information
mtext(
  text = paste(n_nodosP, "nodos |", n_arcosP, "arcos"),
  side = 1,
  line = 3.5,
  cex = 0.7
)

dev.off()

# Degree histogram
png("histogramaT.png", width = 2000, height = 1600, res = 300)
hist(gradosT, breaks = 100, col = "skyblue",
     main = "Distribución de Grados con Percentiles",
     xlab = "Número de Conexiones por Nodo",
     ylab = "Frecuencia")
abline(v = percentilesT, col = c("blue", "green", "red"), lty = 2)
legend("topright", 
       legend = c("P95", "P99", "P99.9"),
       col = c("blue", "green", "red"), 
       lty = 2)
dev.off()

png("histogramaP.png", width = 2000, height = 1600, res = 300)
hist(gradosP, breaks = 100, col = "skyblue",
     main = "Distribución de Grados con Percentiles",
     xlab = "Número de Conexiones por Nodo",
     ylab = "Frecuencia")
abline(v = percentilesP, col = c("blue", "green", "red"), lty = 2)
legend("topright", 
       legend = c("P95", "P99", "P99.9"),
       col = c("blue", "green", "red"), 
       lty = 2)
dev.off()

# Critical node subnet (top 5%)
top_nodosT <- names(sort(gradosT, decreasing = TRUE)[1:round(0.05*length(gradosT))])
hubs_dfT <- data.frame(
  Gen = top_nodosT,
  Conexiones_totales = gradosT[top_nodosT],
  Tipo_conexiones = sapply(top_nodosT, function(g) {
    p <- length(parents(bn_trans, g))
    c <- length(children(bn_trans, g))
    paste0("Padres:", p, "|Hijos:", c)
  })
)

write.csv(hubs_dfT, "hubs_criticosT.csv", row.names = FALSE)

subredT <- subgraph(bn_trans, top_nodosT)

tryCatch({
  png("subred_criticaT.png", width = 3000, height = 3000, res = 300)
  graphviz.plot(subredT, 
                main = "Subred - Nodos Críticos (Top 5%)",
                layout = "fdp")
  dev.off()
}, error = function(e) {
  message("Error al graficar subredT: ", e$message)
})

top_nodosP <- names(sort(gradosP, decreasing = TRUE)[1:round(0.05*length(gradosP))])
hubs_dfP <- data.frame(
  Gen = top_nodosP,
  Conexiones_totales = gradosP[top_nodosP],
  Tipo_conexiones = sapply(top_nodosP, function(g) {
    p <- length(parents(bn_prot, g))
    c <- length(children(bn_prot, g))
    paste0("Padres:", p, "|Hijos:", c)
  })
)

write.csv(hubs_dfP, "hubs_criticosP.csv", row.names = FALSE)

subredP <- subgraph(bn_prot, top_nodosP)

tryCatch({
  png("subred_criticaP.png", width = 3000, height = 3000, res = 300)
  graphviz.plot(subredP, 
                main = "Subred - Nodos Críticos (Top 5%)",
                layout = "fdp")
  dev.off()
}, error = function(e) {
  message("Error al graficar subredP: ", e$message)
})

#Analysis of strong arcs and key routes -------------------------------
#Stronger arches (using BIC as a measure of strength)
arcos_fuertesT <- arcs(bn_trans) %>% 
  as.data.frame() %>%
  mutate(
    from_degree = gradosT[from],  # Regulator connections
    to_degree = gradosT[to],      # Connections of the one being regulated
    strength = (from_degree + to_degree)/2  # Average importance
  ) %>%
  arrange(desc(strength)) %>% 
  head(100)

write.csv(arcos_fuertesT, "relaciones_claveT.csv", row.names = FALSE)

arcos_fuertesP <- arcs(bn_prot) %>% 
  as.data.frame() %>%
  mutate(
    from_degree = gradosP[from],  # Regulator connections
    to_degree = gradosP[to],      # Connections of the one being regulated
    strength = (from_degree + to_degree)/2  # Average importance
  ) %>%
  arrange(desc(strength)) %>% 
  head(100)

write.csv(arcos_fuertesP, "relaciones_claveP.csv", row.names = FALSE)

# Save all results ------------------------------------------
save.image(file = "resultados_bayeasianas.RData")

# Generate summary report ----------------------------------------------
sink("resumen_analisisT.txt")
cat("===========================================\n")
cat("  RESUMEN ANÁLISIS RED BAYESIANA\n")
cat("===========================================\n\n")

cat("ESTADÍSTICAS ESTRUCTURALES:\n")
print(resumen_conexionesT)

cat("\n\nTOP 10 NODOS MÁS CONECTADOS:\n")
head(sort(gradosT, decreasing = TRUE), 10) |> print()

cat("\n\nTOP 10 ARCOS MÁS FUERTES:\n")
head(arcos_fuertesT, 10) |> print()

if(exists("ruta_ejemplo")) {
  cat("\n\nRUTA DE EJEMPLO ENTRE NODOS CRÍTICOS:\n")
  cat(paste0(top_nodosT[1], " -> ", top_nodosT[2], ":\n"))
  print(ruta_ejemplo)
}
sink()

sink("resumen_analisisP.txt")
cat("===========================================\n")
cat("  RESUMEN ANÁLISIS RED BAYESIANA\n")
cat("===========================================\n\n")

cat("ESTADÍSTICAS ESTRUCTURALES:\n")
print(resumen_conexionesP)

cat("\n\nTOP 10 NODOS MÁS CONECTADOS:\n")
head(sort(gradosP, decreasing = TRUE), 10) |> print()

cat("\n\nTOP 10 ARCOS MÁS FUERTES:\n")
head(arcos_fuertesP, 10) |> print()

if(exists("ruta_ejemplo")) {
  cat("\n\nRUTA DE EJEMPLO ENTRE NODOS CRÍTICOS:\n")
  cat(paste0(top_nodosP[1], " -> ", top_nodosP[2], ":\n"))
  print(ruta_ejemplo)
}
sink()

#comparisons

hubs_comunes <- intersect(top_nodosT, top_nodosP)

hubs_comparativo <- data.frame(
  Gen = hubs_comunes,
  Conexiones_Trans = gradosT[hubs_comunes],
  Conexiones_Prot = gradosP[hubs_comunes],
  Diferencia = abs(gradosT[hubs_comunes] - gradosP[hubs_comunes])
) %>% arrange(desc(Conexiones_Trans))

write.csv(hubs_comparativo, "hubs_comunes.csv", row.names = FALSE)

# Convert arcs to string format for comparison
arcos_trans_str <- paste(arcos_fuertesT$from, arcos_fuertesT$to, sep = "|")
arcos_prot_str <- paste(arcos_fuertesP$from, arcos_fuertesP$to, sep = "|")

# Identify common arches in top 100
arcos_comunes <- intersect(arcos_trans_str, arcos_prot_str)

# Create comparative report
arcos_comparativo <- data.frame(
  Arco = arcos_comunes,
  Fuerza_Trans = arcos_fuertesT$strength[match(arcos_comunes, arcos_trans_str)],
  Fuerza_Prot = arcos_fuertesP$strength[match(arcos_comunes, arcos_prot_str)],
  Diferencia_Fuerza = abs(
    arcos_fuertesT$strength[match(arcos_comunes, arcos_trans_str)] - 
      arcos_fuertesP$strength[match(arcos_comunes, arcos_prot_str)]
  )
) %>% arrange(desc(Fuerza_Trans))

# save CSV
write.csv(arcos_comparativo, "arcos_top_comunes.csv", row.names = FALSE)

# Calculate correlation between connectivities
cor_grados <- cor(gradosT[common_genes], gradosP[common_genes], use = "complete.obs")

# Scatterplot of connectivities
png("correlacion_conectividad.png", width = 1200, height = 1000, res = 150)
plot(gradosT[common_genes], gradosP[common_genes],
     main = paste("Correlación de Conectividad (r =", round(cor_grados, 3), ")"),
     xlab = "Grado Transcriptómica",
     ylab = "Grado Proteómica",
     pch = 19, col = rgb(0, 0.5, 0.8, 0.5))
abline(lm(gradosP[common_genes] ~ gradosT[common_genes]), col = "red", lwd = 2)
dev.off()

sink("reporte_comparativo.txt")

cat("==================================================\n")
cat("  COMPARACIÓN REDES TRANSCRIPTOMICA vs PROTEOMICA\n")
cat("==================================================\n\n")

cat("1. ESTADÍSTICAS GLOBALES\n")
cat("----------------------------------------\n")
cat(sprintf("Densidad Transcriptómica: %.6f\n", densidadT))
cat(sprintf("Densidad Proteómica: %.6f\n", densidadP))
cat(sprintf("Diferencia Densidad: %.6f\n\n", abs(densidadT - densidadP)))

cat(sprintf("Nodos comunes analizados: %d\n", length(common_genes)))
cat(sprintf("Hubs comunes (top 5%%): %d / %d\n", 
            length(hubs_comunes), length(top_nodosT)))
cat(sprintf("Arcos top-100 comunes: %d\n\n", length(arcos_comunes)))

cat("2. CORRELACIÓN DE CONECTIVIDAD\n")
cat("----------------------------------------\n")
cat(sprintf("Correlación (Pearson): %.3f\n", cor_grados))
cat("Interpretación: ")
if(cor_grados > 0.7) cat("Fuerte correlación") else 
  if(cor_grados > 0.4) cat("Correlación moderada") else 
    cat("Baja correlación")

cat("\n\n3. TOP 5 HUBS COMUNES\n")
cat("----------------------------------------\n")
print(head(hubs_comparativo, 5))

cat("\n\n4. TOP 5 ARCOS COMUNES\n")
cat("----------------------------------------\n")
print(head(arcos_comparativo, 5))

sink()

# Identify genes with the greatest difference in connectivity
dif_conectividad <- data.frame(
  Gen = common_genes,
  Dif_Grado = abs(gradosT[common_genes] - gradosP[common_genes])
) %>% arrange(desc(Dif_Grado))

# Save results
write.csv(dif_conectividad, "genes_diferencia_conectividad.csv", row.names = FALSE)