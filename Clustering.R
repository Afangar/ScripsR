#Code for clustering transcriptomic data

setwd("/Your/address/User")
library(reticulate)
library(ManiNetCluster)

X <- as.matrix(read.csv("TranscritosX.csv", row.names=1))
Y <- as.matrix(read.csv("TranscritosY.csv", row.names=1))
total_rows <- nrow(X)

num_clusters <- 333
porcentaje_muestra <- 0.7

for(n in 1:num_clusters) {
    muestras <- sample(total_rows, total_rows * porcentaje_muestra)

    X_sub <- X[muestras, ]
    Y_sub <- Y[muestras, ]

    corr_sub <- Correspondence(matrix = diag(nrow(X_sub)))

    resultado <- ManiNetCluster(
        X_sub, Y_sub,
        nameX = 'early',
        nameY = 'late',
        corr = corr_sub,
        d = 3L,
        method = 'nonlinear manifold aln',
        k_NN = 6L,
        k_medoids = 60L
    )

    resultado$id <- c(rownames(X_sub), rownames(Y_sub))

    write.csv(
        resultado,
        file = paste0("/Your/address/User/Clusters/Cluster_", n, ".csv")
    )

    rm(muestras, X_sub, Y_sub, corr_sub, resultado)
    gc()
}

#code for clustering of randomized transcriptomics data

setwd("/Your/address/User")
library(reticulate)
library(ManiNetCluster)

X <- as.matrix(read.csv("TranscritosX.csv", row.names=1))
Y <- as.matrix(read.csv("TranscritosY.csv", row.names=1))

X_rand <- t(apply(X, 1, sample))
Y_rand <- t(apply(Y, 1, sample))

 write.csv(
        X_rand,
        file = "/Your//Your/address/User/Clustersβ/MatrixX.csv")
 write.csv(
        Y_rand,
        file = "/Your//Your/address/User/Clustersβ/MatrixY.csv")

total_rows <- nrow(X_rand)
num_clusters <- 100
porcentaje_muestra <- 0.7

for(n in 1:num_clusters) {
    muestras <- sample(total_rows, total_rows * porcentaje_muestra)

    X_sub <- X_rand[muestras, ]
    Y_sub <- Y_rand[muestras, ]

    corr_sub <- Correspondence(matrix = diag(nrow(X_sub)))

    resultado <- ManiNetCluster(
        X_sub, Y_sub,
        nameX = 'early',
        nameY = 'late',
        corr = corr_sub,
        d = 3L,
        method = 'nonlinear manifold aln',
        k_NN = 6L,
        k_medoids = 60L
    )


    resultado$id <- c(rownames(X_sub), rownames(Y_sub))

    write.csv(
        resultado,
        file = paste0("/Your//Your/address/User/Clustersβ/Clusterβ_", n, ".csv")
    )
