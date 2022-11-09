# Borra todos los objetos creados previamente
rm(list = ls())

setwd("C:/Users/Usuario/Desktop/mypackage/src")


library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(pkgKitten)

#=======================================
#   Para hacer gráfico de H-matrix
#=======================================

# Simulamos locaciones en el espacio (2D) 
loc  <- cbind(x1 = runif(1600), x2 = runif(1600))
loc <- as.data.frame(loc)
head(loc)

p <- ggplot(loc, aes(x1, x2))
p + geom_point() + 
  xlab("x") +
  ylab("y") + 
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size=24),
        axis.text=element_text(size=18),
        legend.position="top",
        axis.title=element_text(size=18,face="bold")) 



y <- rnorm(length(loc[,1]), 0, 1)
y


# Inicializar MPI
mpiinit()

# base case
epsilon = 0.0001
alpha = 1
eta = 2
testTps(loc, y, Epsilon = epsilon, Eta = eta, MinClusterSize = 10, lambda = alpha)









# Calculamos distancias
distances12 <- DistanceMatrix(loc[1:497,], loc[498:500,])
dim(distances12)
E12 = radialFunction(distances12, TPS, R)
dim(E12)

# Calculamos distancias
distances22 <- DistanceMatrix(loc[498:500,], loc[498:500,])
dim(distances22)
E22 = radialFunction(distances22, TPS, R)
dim(E22)


# Creamos T1
T1 <- cbind(1, loc[1:497, 1:2])
dim(T1)

# Creamos T2
T2 <- cbind(1, loc[498:500, 1:2])
dim(T2)


M1 <- cbind(E12, T1) # T1 aqui es la t(T1) y da lo mismo
dim(M1)


# Creamos M2
aux1 <- cbind(E22, T2)
aux2 <- cbind(T2, 0, 0, 0)
dim(aux1)
dim(aux2)

M2 <- rbind(aux1, aux2)
dim(M2)

library(matlib)
M2inv <- inv(M2)

M2 %*% M2inv


# M1 y M2 deben ir a la función calculateTPS() creada en cpp (example_TPS).


#======================================================================================================================
# Simulamos un vector numérico (o matriz dependiendo de como lo declaramos en la funcion del gradiente conjudado) para
# nuestra variable respuesta (funcion) observada en casa sitio (indice)
#======================================================================================================================
y <- rnorm(length(loc[,1]), 0, 1)
y

# Si lo queremos en formato de matriz debemos hacer 
# y <- as.matrix(rnorm(length(loc[,1]), 0, 1))
