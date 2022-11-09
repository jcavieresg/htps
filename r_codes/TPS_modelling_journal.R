rm(list = ls())
library(matlib)
library(tictoc)
library(pracma)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(plot3D)


# Inicializar MPI
mpiinit()




#===============================================================================================
#                        Main code of the Smoothing TPS and H-matrix
#===============================================================================================
set.seed(1234)                       # Fijamos la semilla
alpha = 1 # smoothing parameter
shape = 1 # shape parameter of the RBF

K = seq(20, 40, by=2)                # Distintos tamaños de la grilla
#K = 100
#nok = length(K)                     # Numero de grillas creadas anteriormente

nosites <- vector()             # Numero de grillas
time_full <- vector()           # Solve function
time_full_cg <- vector()        # CG 
time_hmat <- vector()           # CG + H-matrix

error_full_cg <- vector()       # Error CG
error_hmat <- vector()          # Error CG + H matrix


#for ( k in 1:nok){
for ( k in 1:length(K)){
  print(k)
  #=====================================================================
  # generate sites and values
  #=====================================================================
  # structured sites
  # n <- K[k]
  # N = n*n # number of sites
  # nosites[k] <- N
  # x <- seq(0,n-1,by=1)/(n-1)
  # loc <- cbind(rep(x,times=n),rep(x,each=n))
  # #plot(loc)
  # s1 = c(2:(n-1),(n+1):(N-1))
  # s2 = c(1, n, N)
  
  #Non-structured sites
  n <- K[k]
  N <- n*n
  nosites[k] <- N
  #N <- K[k]
  lat <- runif(N, 0, 1)
  lon <- runif(N, 0, 1)
  loc <- cbind(lat, lon)
  s1 = c(2:(n-1),(n+1):(N-1))
  s2 = c(1,n,N)
  
  # values
  val <- testfunction(matrix(loc[, 1]), matrix(loc[, 2])) 
  

  list <- vector("list", length(K))
  # Create a for statement to populate the list
  #for (i in length(nosites)) {
  list[[k]] <- val
  #}
  
  
  #=====================================================================
  # full matrix with solve() function of R
  #=====================================================================
  #exectime <- tic()
  tic <- tic()
  d0 <- DistanceMatrix(loc,loc)
  E0 <- radialFunction(d0,2,1.0, shape)
  E0pa <- E0 + alpha*diag(N)
  T0 <- cbind(1,loc)
  M0 <- rbind(cbind(E0pa,T0),cbind(t(T0),0,0,0))
  rhs0 <- c(val,c(0,0,0))
  sol0 <- solve(M0,rhs0)
  #exectime <- toc()
  toc <- toc()
  time_full[k] <- toc 
  
  # =====================================================================
  # full matrix with CG solver
  # =====================================================================
  #exectime <- tic()
  tic <- tic()
  loc1 <- loc[s1,]
  loc2 <- loc[s2,]
  val1 <- val[s1]
  val2 <- c(val[s2],c(0,0,0))
  

  # Distances
  distances11 <- DistanceMatrix(loc1,loc1)
  E11 <- radialFunction(distances11, 2, 1, shape) + alpha*diag(N-3)
  
  distances12 <- DistanceMatrix(loc1,loc2)
  E12 <- radialFunction(distances12, 2, 1, shape)
  distances22 <- DistanceMatrix(loc2,loc2)
  E22 <- radialFunction(distances22, 2, 1, shape) + alpha*diag(3)
  
  T1 <- cbind(1,loc1)
  T2 <- cbind(1,loc2)
  M1 <- cbind(E12, T1)
  aux1 <- cbind(E22, T2)
  aux2 <- cbind(t(T2), 0, 0, 0)
  M2 <- rbind(aux1, aux2)
  M2inv <- inv(M2)
  S <- M1%*%M2inv
  
  # just for testing - full matrix and Schur complement
  M_full <- E11 - S%*%t(M1)
  
  # lado derecho
  rhs <- val1 - S%*%val2
  
  #y <- calculateTPS(loc1,rhs,t(M1),S,0.00002,0.0,20,alpha)
  y <- calculateTPS_full(rhs, M_full)
  delta2 <- solve(t(T2),-t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
  
  solf <- replicate(N+3,0)
  solf[s1] <- y
  solf[s2] <- delta2
  solf[(N+1):(N+3)] <- a
  
  #exectime <- toc()
  toc <- toc()
  #time_full_cg[k] <- exectime$toc - exectime$tic
  time_full_cg[k] <- toc 
  error_full_cg[k] <- norm(as.matrix(sol0)- as.matrix(solf))
  
  
  # 
  #=====================================================================
  # H matrix with CG solver
  #=====================================================================
  #exectime <- tic()
  tic <- tic()
  
  loc1 <- loc[s1,]
  loc2 <- loc[s2,]
  val1 <- val[s1]
  val2 <- c(val[s2],c(0,0,0))
  

  
  # Distances
  distances12 <- DistanceMatrix(loc1,loc2)
  E12 <- radialFunction(distances12, 2, 1, shape)
  distances22 <- DistanceMatrix(loc2,loc2)
  E22 <- radialFunction(distances22, 2, 1, shape) + alpha*diag(3)
  
  T1 <- cbind(1,loc1)
  T2 <- cbind(1,loc2)
  M1 <- cbind(E12, T1)
  aux1 <- cbind(E22, T2)
  aux2 <- cbind(t(T2), 0, 0, 0)
  M2 <- rbind(aux1, aux2)
  M2inv <- inv(M2)
  S <- M1%*%M2inv
  
  # right side
  rhs <- val1 - S%*%val2
  
  y <- calculateTPS_hmat(loc1,rhs,t(M1),S, Epsilon = 0.0001, Eta = 2, MinClusterSize = 20, alpha)
  delta2 <- solve(t(T2),- t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
 
  solh <- replicate(N+3,0)
  solh[s1] <- y
  solh[s2] <- delta2
  solh[(N+1):(N+3)] <- a
  
  toc <- toc()
  time_hmat[k] <- toc 
  error_hmat[k] <- norm(as.matrix(sol0)- as.matrix(solh))
  
}
#=======================================================================================================

#mpifinalize()



#===============================================================================================
#                       Fit baed in the book Meshfree approximation
#===============================================================================================
fit <- PLS(loc, loc, 2, 1, neval = 40, alpha = alpha, shape = shape)
fit$Pf

scatter3D(loc[, 1], loc[, 2], fit$Pf0[1:nrow(loc-3)], theta = 2, phi = 25, bty = "g")
scatter3D(fit$epoints[, 1], fit$epoints[, 2], fit$Pf, theta = 2, phi = 25, bty = "g")

# Data plot
data <- data.frame(sol0[1:length(loc[,1]-3)], 
                   solf[1:length(loc[,1]-3)],
                   solh[1:length(loc[,1]-3)], loc[,1], loc[,2])
colnames(data) <- c("sol0", "solf", "solh", "x", "y")
head(data)

library(akima)
library(spatialkernel)
library(fields)
quilt.plot(data$x, data$y, data$sol0)

par(oma=c(0, 1,0,2)) 
quilt.plot(data$x, data$y, data$sol0, xlim = c(0,1), ylim = c(0,1), nx=150, ny=150, xlab="x", legend.line = 2,font = 2,
           ylab="y", cex.axis=1.4, cex.lab = 2, main = "M1", cex.main = 2.4, axis.args = list(cex.axis = 1.6), font.lab=3)



quilt.plot(data$x, data$y, data$solf, xlim = c(0,1), ylim = c(0,1), nx=150, ny=150, xlab="x", legend.line = 2, font = 2,
           ylab="y", cex.axis=1.4, cex.lab = 2, main = "M2", cex.main = 2.4, axis.args = list(cex.axis = 1.6), font.lab=3)


quilt.plot(data$x, data$y, data$solh, xlim = c(0,1), ylim = c(0,1), nx=150, ny=150, xlab="x", legend.line = 2, font = 2,
           ylab="y", cex.axis=1.4, cex.lab = 2, main = "M3", cex.main = 2.4, axis.args = list(cex.axis = 1.6), font.lab=3)


# Evalution of the function in epoints created in RcppArmadillo
DM_val <- DistanceMatrix(fit$epoints, loc)
EM_val <- radialFunction(DM_val, 2, 1.0, shape)
T0_val <- cbind(1,fit$epoints)
M0_val <- rbind(cbind(EM_val,T0_val))
sol_eval <- M0_val %*% sol0
sol_evalCG <- M0_val %*% solf
sol_evalCGH <- M0_val %*% solh

scatter3D(loc[, 1], loc[, 2], sol0[1:nrow(loc-3)], theta = 2, phi = 25, bty = "g")
scatter3D(fit$epoints[, 1], fit$epoints[, 2], sol_eval, theta = 2, phi = 25, bty = "g")


scatter3D(fit$epoints[, 1], fit$epoints[, 2], sol_eval, theta = 120, phi = 25, bty = "g")
scatter3D(fit$epoints[, 1], fit$epoints[, 2], sol_evalCG, theta = 120, phi = 25, bty = "g")
scatter3D(fit$epoints[, 1], fit$epoints[, 2], sol_evalCGH, theta = 120, phi = 25, bty = "g")



#========================================================================================================
#                                     Graphs
#========================================================================================================

library(gridExtra)
library(tidyverse)

par(mfrow=c(1,3))

pdf('plot4.pdf')
par(mar = c(2, 2, 2, 3)) 
scatter3D(fit$epoints[,1], fit$epoints[, 2], sol_eval, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M1", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values" )
dev.off()

pdf('plot5.pdf')
scatter3D(fit$epoints[,1], fit$epoints[, 2], sol_evalCG, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M2", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values" )
dev.off()

pdf('plot6.pdf')
scatter3D(fit$epoints[,1], fit$epoints[, 2], sol_evalCGH, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M3", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values")
dev.off()




#==========================================
#         Graph of log(logtimes)
#==========================================
# Una manera de hacerlo
df1 <- data.frame(nosites, time_full)
df2 <- data.frame(nosites, time_full_cg)
df3 <- data.frame(nosites, time_hmat)

do2 <- data.frame(nosites,nosites*nosites/5e7)
do3 <- data.frame(nosites,nosites*nosites*nosites/1e10)

# df1 <- data.frame(K, time_full)
# df2 <- data.frame(K, time_full_cg)
# df3 <- data.frame(K, time_hmat)

names(df1) <- c("x", "y")
names(df2) <- c("x", "y")
names(df3) <- c("x", "y")
names(do2) <- c("x", "y")
names(do3) <- c("x", "y")

library(reshape2)
newData <- melt(list(df1 = df1, df2 = df2, df3 = df3, do2 = do2, do3=do3), id.vars = "x")
head(newData)

# Log scale
#pdf('plot7.pdf')
ggplot(newData, aes(log(x), log(value), colour = L1)) +
  geom_line(linetype = "dashed")+
  geom_point(size = 2) + 
  xlab("log(number of sites)") +
  ylab("log(time)") + 
  labs(color='Methods')  +
  # you should specify the color also
  scale_color_manual(labels = c("M1", "M2", "M3", "Lineal", "Cubic"), 
                     values = c("blue", "#E69F00", "red", "black", "yellow")) +
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold")) 


# colnames(df1) <- c("x", "y1")
# colnames(df2) <- c("x", "y2")
# colnames(df3) <- c("x", "y3")
# colnames(do2) <- c("x", "y4")
# colnames(do3) <- c("x", "y5")
# 
# colors <- c("y1" = "blue", "y2" = "#E69F00", "y3" = "red" , "y4" = "brown1", "y5" = "darkblue")
# 
# ggplot(data=df1, aes(x=log(x), y=log(y1), color = "y1")) + 
#   geom_point(data=df1, aes(x=log(x), y=log(y1), color='y1'), size = 2) + 
#   geom_line(data=df1, aes(x=log(x), y=log(y1), color='y1'),linetype = "dashed") + 
#   geom_point(data=df2, aes(x=log(x), y=log(y2), color='y2'), size = 2) + 
#   geom_point(data=df3, aes(x=log(x), y=log(y3), color='y3'), size = 2) +
#   geom_line(data=df2, aes(x=log(x), y=log(y2), color='y2'), linetype = "dashed") + 
#   geom_line(data=df3, aes(x=log(x), y=log(y3), color='y3'), linetype = "dashed") +
#   geom_line(data=do2, aes(x=log(x), y=log(y4), color='y4'), linetype = "dotted", size = 1) +
#   geom_line(data=do3, aes(x=log(x), y=log(y5), color='y5'), linetype = "dotted", size = 1) +
#   xlab("log(number of sites)") +
#   ylab("log(time)") + 
#   labs(color='Methods')  +
#   scale_color_manual(values = colors, labels = c("M1", "M2", "M3", "Lineal", "Cubic"), 
#                       guide = guide_legend(override.aes = list(
#                        linetype = c("dashed", "dashed", "dashed"),
#                        shape = c(16, 16, 16)))) + 
#   theme(legend.title = element_text(size=18),
#         legend.text = element_text(size=16),
#         axis.text=element_text(size=14, face = "bold"),
#         legend.position="top",
#         axis.title=element_text(size=14,face="bold")) 


        





#=================================================
#             For the errors
#=================================================
library(dplyr)
df4 <- data.frame(nosites, error_full_cg)
df5 <- data.frame(nosites, error_hmat)

df4 <- data.frame(K, error_full_cg)
df5 <- data.frame(K, error_hmat)

newdata3 <- inner_join(df4, df5)
head(newdata3)

colors2 <- c("M1 - M2" = "#E69F00", "M1 - M3" = "red")


ggplot(newdata3, aes(x = log(nosites))) +
  geom_line(aes(y = log(error_full_cg), color = "M1 - M2"), linetype = "dashed") +
  geom_point(aes(x = log(nosites), y = log(error_full_cg), colour = "M1 - M2"), size = 2) +
  geom_line(aes(y = log(error_hmat), color = "M1 - M3"), linetype = "dashed") +
  geom_point(y = log(error_hmat), size = 2, color = "red") +
  xlab("log(number of sites)") +
  ylab("log(Error value)") +
  labs(color='')  +
  scale_color_manual(values = colors2) +
  theme(legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=16),
        legend.position="top",
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold")) 





# Table 1
# Error de coeficientes estimados 
error_solcg <- norm(sol0 - solf); error_solcg
error_solh <- norm(sol0 - solh); error_solh


# RMS error
exact <- testfunction(matrix(fit$epoints[, 1]), matrix(fit$epoints[, 2])) 

# Evalution of the function in epoints created in Matlab
DM_val <- DistanceMatrix(fit$epoints, loc)
EM_val <- radialFunction(DM_val, 2, 1.0, shape)
T0_val <- cbind(1,fit$epoints)
M0_val <- rbind(cbind(EM_val,T0_val))
sol_eval <- M0_val %*% sol0
solf_eval <- M0_val %*% solf
solh_eval <- M0_val %*% solh

neval <- dim(fit$epoints)[1]/40

# Tabla 2
exact_sol0 <- norm(sol_eval - exact); exact_sol0 
exact_solf <- norm(solf_eval - exact); exact_solf 
exact_solh <- norm(solh_eval - exact); exact_solh

rms_sol0 <- norm(sol_eval - exact) /neval; rms_sol0
rms_solf <- norm(solf_eval - exact) /neval; rms_solf
rms_solh <- norm(solh_eval - exact) /neval; rms_solh

all_errors <- data.frame(exact_sol0, exact_solf, exact_solh, rms_sol0, rms_solf, rms_solh)
all_errors


