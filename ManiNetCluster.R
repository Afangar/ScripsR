library(reticulate)
library(devtools)
#devtools::install_github("namtk/ManiNetCluster")
library(ManiNetCluster)

#THIS CODE IS ONLY FOR BASH USING CONDA/MAMBA ENVIRONMENTS

use_python("/home/User/.conda/envs/YourEnv/bin/python")

X <- as.matrix(read.csv("Matrixcondition1.csv", row.names=1))
Y <- as.matrix(read.csv("Matrixcondition2.csv", row.names=1))

n <- nrow(X) # the number of genes in the dataset
corr <- Correspondence(matrix=diag(n))

Expe <- ManiNetCluster(X,Y,nameX='condition1',nameY='condition2',corr=corr,d=3L,method='linear manifold',k_NN=6L,k_medoids=60L)

#other options

#'cca': canonical-correlation analysis
#'manifold warping': manifold warping
#'nonlinear manifold aln': nonlinear manifold alignment
#'nonlinear manifold warp': nonlinear manifold warping

Expe$id <- c(read.csv("Matrixcondition1.csv")[,1], read.csv("Matrixcondition2.csv")[,1])

#The problem with this ManiNetCluster project is that it doesn't tell you that Numpy is going
#to use all the available threads on the server or computer you run it on, which is a big problem.
#The solution is to run Python in the environment before running R, this is a example for SLURM in bash.

NSLOTS=$SLURM_NTASKS_PER_NODE
export OMP_NUM_THREADS=$NSLOTS
export OPENBLAS_NUM_THREADS=$NSLOTS
export MKL_NUM_THREADS=$NSLOTS
export VECLIB_MAXIMUM_THREADS=$NSLOTS
export NUMEXPR_NUM_THREADS=$NSLOTS


#Visualization and Comparing the non alignment and alignment results

df1 <- df[df$data == 'light',]
df2 <- df[df$data == 'dark',]
df1 <- df1[order(df1$Val1,df1$Val2,df1$Val3), ]
df2 <- df2[match(df1$id, df2$id), ]

pal <- rainbow(n)
library(rgl)
for (i in 1:n) {
  with(df1[i,],points3d(Val1,Val2,Val3,size=7,col=pal[i], alpha=.1))
}

# Add axes and labels
axes3d()
title3d(
  main = "Mi Gráfica 3D",
  xlab = "Eje X",
  ylab = "Eje Y",
  zlab = "Eje Z"
)
grid3d(c("x", "y", "z")) # Add a grid


#Using the above, the data points of dataset 1 will form a rainbow cloud as depicted in the figure

for (i in 1:n) {
  with(df2[i,],points3d(Val1,Val2,Val3,size=7,col=pal[i]))
}

#We can compare this alignment result with unalignment one. To plot the original datasets we use PCA

pca_coordinatesX=prcomp(X)$x
df1 <- data.frame(pca_coordinatesX[,1:3])
names(df1) <- c('Val1', 'Val2', 'Val3')
df1$id <- read.csv("Downloads/dayOrthoExpr.csv")[,1]

pca_coordinatesY=prcomp(Y)$x
df2 <- data.frame(pca_coordinatesY[,1:3])
names(df2) <- c('Val1', 'Val2', 'Val3')
df2$id <- read.csv("Downloads/nightOrthoExpr.csv")[,1]

df1 <- df1[order(df1$Val1,df1$Val2,df1$Val3), ]
df2 <- df2[match(df1$id, df2$id), ]

for (i in 1:n) {
  with(df1[i,],points3d(Val1,Val2,Val3,size=7,col=pal[i], alpha=.1))
}

for (i in 1:n) {
  with(df2[i,],points3d(Val1,Val2,Val3,size=7,col=pal[i]))
}