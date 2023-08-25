rm(list = ls())
library(matlib)
library(tictoc)
library(pracma)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(plot3D)
library(mgcv)
library(sp)
library(gstat)


# Inicializar MPI
mpiinit()




#===============================================================================================
#                        Main code of the Smoothing TPS and H-matrix
#===============================================================================================
set.seed(1234)                       # Fijamos la semilla
library(Metrics)
alpha = 1 # smoothing parameter
shape = 1 # shape parameter of the RBF

# = seq(20, 80, by=2)                # Sequence of grids 20^2 to 80^2
K = seq(20, 40, by=2)               # Sequence of grids 20^2 to 40^2
#K = seq(5, 10, by=1)               # Sequence of grids 20^2 to 40^2
#K = 40

nosites <- vector()             # N sites
time_full <- vector()           # Solve function
time_full_cg <- vector()        # CG 
time_hmat <- vector()           # CG + H-matrix
time_mgcv <- vector()           # CG + H-matrix
time_krig <- vector()           # CG + H-matrix

error_full_cg <- vector()       # Error CG
error_hmat <- vector()          # Error CG + H matrix
error_mgcv <- vector()          # Error CG + H matrix
error_krig <- vector()          # Error CG + H matrix


for (k in 1:length(K)){
  print(k)
  #=====================================================================
  # generate sites and values
  #=====================================================================
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
  val <- val + 0.03*rnorm(N)
  #val <- testfunction_random(matrix(loc[, 1]), matrix(loc[, 2])) 
  

  list <- vector("list", length(K))
  list[[k]] <- val

  #=====================================================================
  # full matrix with solve() function of R
  #=====================================================================
  tic <- tic()
  d0 <- DistanceMatrix(loc,loc)
  E0 <- radialFunction(d0,2,1.0, shape)
  E0pa <- E0 + alpha*diag(N)
  T0 <- cbind(1,loc)
  M0 <- rbind(cbind(E0pa,T0),cbind(t(T0),0,0,0))
  rhs0 <- c(val,c(0,0,0))
  sol0 <- solve(M0,rhs0)
  toc <- toc()
  time_full[k] <- toc
  
  
  # mgcv TPS
  tic <- tic()
  fit_mgcv <- gam(val ~ s(loc[, 1], loc[, 2], bs="tp", k = k), method="REML")
  solmgcv <- fit_mgcv$fitted.values
  toc <- toc()
  time_mgcv[k] <- toc
  error_mgcv[k] <- norm(as.matrix(head(sol0,-3)) - as.matrix(solmgcv))
  
  # Kriging
  # dat_krig <- cbind(loc, val)
  # dat_krig <- data.frame(dat_krig)
  # colnames(dat_krig) <- c("x", "y", "val")
  # coordinates(dat_krig) <- ~ x + y
  # tic <- tic()
  # tps.vgm <- variogram(val~1, dat_krig) # calculates sample variogram values 
  # tps.fit <- fit.variogram(tps.vgm, model=vgm(1, "Gau", 1, 5)) # fit model
  # 
  # # grid_points <- expand.grid(x = seq(min(loc[, 1]), max(loc[, 1]), by = ),
  # #                            y = seq(min(loc[, 2]), max(loc[, 2]), by = ))
  # 
  # num_points <- sqrt(nrow(loc))  # Adjust as needed
  # grid_points <- expand.grid(x = seq(min(loc[, 1]), max(loc[, 1]), length.out = num_points),
  #                            y = seq(min(loc[, 2]), max(loc[, 2]), length.out = num_points))
  # 
  # 
  # coordinates(grid_points) <- ~ x + y # step 3 above
  # tps.kriged <- krige(val ~ 1, dat_krig, grid_points, model=tps.fit)
  # solkrig <- tps.kriged$var1.pred
  # toc <- toc()
  # time_krig[k] <- toc
  # error_krig[k] <- norm(as.matrix(head(sol0,-3)) - as.matrix(solkrig))
  
  
  

  
  # =====================================================================
  # full matrix with CG solver
  # =====================================================================
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

  y <- calculateTPS_full(rhs, M_full)
  delta2 <- solve(t(T2),-t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)

  solf <- replicate(N+3,0)
  solf[s1] <- y
  solf[s2] <- delta2
  solf[(N+1):(N+3)] <- a


  toc <- toc()
  time_full_cg[k] <- toc
  error_full_cg[k] <- norm(as.matrix(head(sol0, -3))- as.matrix(head(solf, -3)))


  
  #=====================================================================
  # H matrix with CG solver
  #=====================================================================
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
  #y <- calculateTPS_hmat(loc1,rhs,t(M1),S, Epsilon = 0.1, Eta = 10, MinClusterSize = 20, alpha)
  delta2 <- solve(t(T2),- t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
 
  solh <- replicate(N+3,0)
  solh[s1] <- y
  solh[s2] <- delta2
  solh[(N+1):(N+3)] <- a
  
  toc <- toc()
  time_hmat[k] <- toc 
  #error_hmat[k] <- norm(as.matrix(sol0)- as.matrix(solh))
  error_hmat[k] <- norm(as.matrix(head(sol0, -3))- as.matrix(head(solh, -3)))
  
}
#=======================================================================================================

#mpifinalize()



#===============================================================================================
#                       Fit baed in the book Meshfree approximation
#===============================================================================================
#fit <- PLS(loc, loc, 2, 1, neval = 40, alpha = alpha, shape = shape)
#fit <- PLS(loc, loc, 2, 1, neval = 5, alpha = alpha, shape = shape)

#==========
# Table 1
#==========
error_solcg <- norm(matrix(sol0) - matrix(solf)); error_solcg
error_solh <- norm(matrix(sol0) - matrix(solh)); error_solh

#=
# Saving the plot in a new folder
#=

setwd("C:/Users/Usuario/Desktop/Thesis/article_1/figuras_newversion")


#==========================================
#     COMP. ERROR (log scales)
#             FIGURE 2
#==========================================

#=================================================
#             For the errors
#=================================================
library(dplyr)
df1 <- data.frame(nosites, error_full_cg)
df2 <- data.frame(nosites, error_hmat)
df3 <- data.frame(nosites, error_mgcv)


# Merge data frames using dplyr package
newdata <- df1 %>%
  left_join(df2, by = "nosites") %>%
  left_join(df3, by = "nosites")

head(newdata)

colors <- c("M1 - M2" = "#E69F00", "M1 - M3" = "red")


#pdf('plot_error.pdf')
pdf('plot_error_20_to_40.pdf')
ggplot(newdata, aes(x = log(nosites))) +
  geom_point(aes(x = log(nosites), y = log(error_full_cg), colour = "M1 - M2"), size = 3, shape = 0) +
  geom_line(aes(y = log(error_full_cg), color = "M1 - M2"), linetype = "dashed") +
  geom_point(y = log(error_hmat), size = 3, color = "red", shape = 2) +
  geom_line(aes(y = log(error_hmat), color = "M1 - M3"), linetype = "dashed") +
  xlab("log(number of sites)") +
  ylab("log(Error value)") +
  labs(color='')  +
  #scale_color_manual(values = colors) +
  scale_color_manual(values = colors, 
                     guide = guide_legend(override.aes = list(
                     linetype = c("dashed", "dashed"),
                       shape = c(0, 2)))) + 
  theme(legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=20),
        legend.position="top",
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.spacing.x = unit(1, 'cm')) 
        #+ annotate(geom="label", x= 6.18, y=42, label="CASE 4", color="black", size = 8)
dev.off()




#=================================
#        INTERPOLATION
#=================================
# New grid for interpolation
ng = 40
x.pred <- seq(0, 1, length.out = ng)
y.pred <- seq(0, 1, length.out = ng)
xy <- expand.grid( x = x.pred, y = y.pred)
xy <- as.matrix(xy)
class(xy)

exact <- testfunction(matrix(xy[, 1]), matrix(xy[, 2])) 

# Evalution of the function in epoints created in Matlab
# DM_val <- DistanceMatrix(fit$epoints, loc)
# EM_val <- radialFunction(DM_val, 2, 1.0, shape)
# T0_val <- cbind(1,fit$epoints)
# M0_val <- rbind(cbind(EM_val,T0_val))
# sol_eval <- M0_val %*% sol0
# solf_eval <- M0_val %*% solf
# solh_eval <- M0_val %*% solh



DM_val <- DistanceMatrix(xy, loc)
EM_val <- radialFunction(DM_val, 2, 1.0, shape)
T0_val <- cbind(1, xy)
M0_val <- rbind(cbind(EM_val,T0_val))
sol_eval <- M0_val %*% sol0
solf_eval <- M0_val %*% solf
solh_eval <- M0_val %*% solh


# scatter3D(fit$epoints[, 1], fit$epoints[, 2], sol_eval, theta = 120, phi = 25, bty = "g")
# scatter3D(fit$epoints[, 1], fit$epoints[, 2], solf_eval, theta = 120, phi = 25, bty = "g")
# scatter3D(fit$epoints[, 1], fit$epoints[, 2], solh_eval, theta = 120, phi = 25, bty = "g")


#======================================
#               Table2
#======================================
rmse_sol0 <- rmse(sol_eval,  exact); rmse_sol0
rmse_solf <- rmse(solf_eval, exact); rmse_solf
rmse_solh <- rmse(solh_eval, exact); rmse_solh

all_errors <- data.frame(rmse_sol0, rmse_solf, rmse_solh)
all_errors



#=====================================
#              FIGURE 3
#=====================================
pdf('plot_intersol0.pdf')
par(mar = c(2, 2, 2, 3)) 
scatter3D(xy[,1], xy[, 2], sol_eval, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M1", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values" )
dev.off()

pdf('plot_intersolf.pdf')
scatter3D(xy[,1], xy[, 2], solf_eval, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M2", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values" )
dev.off()

pdf('plot_intersolh.pdf')
scatter3D(xy[,1], xy[, 2], solh_eval, bty = "g", pch = ".", cex = 4, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M3", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values")
dev.off()




#=================================
#            FIGURE 4
#=================================
df1 <- data.frame(nosites, time_full)
df2 <- data.frame(nosites, time_full_cg)
df3 <- data.frame(nosites, time_hmat)

df4 <- data.frame(nosites,nosites^3/1e10)
df5 <- data.frame(nosites,nosites*nosites/5e7)

colnames(df1) <- c("x", "y1")
colnames(df2) <- c("x", "y2")
colnames(df3) <- c("x", "y3")
colnames(df4) <- c("x", "y4")
colnames(df5) <- c("x", "y5")

colors <- c("y1" = "blue", "y2" = "#E69F00", "y3" = "red", "y4" = "black", "y5" = "black")

#pdf('plot_comptime.pdf')
pdf('plot_comptime_20_to_40.pdf')
ggplot(data=df1, aes(x=log(x), y=log(y1), color = "y1")) + 
  geom_point(data=df1, aes(x=log(x), y=log(y1), color='y1'), size = 3, shape = 1) + 
  geom_line(data=df1, aes(x=log(x), y=log(y1), color='y1'),linetype = "dashed") + 
  geom_point(data=df2, aes(x=log(x), y=log(y2), color='y2'), size = 3, shape = 0) + 
  geom_line(data=df2, aes(x=log(x), y=log(y2), color='y2'),linetype = "dashed") + 
  geom_point(data=df3, aes(x=log(x), y=log(y3), color='y3'), size = 3, shape = 2) +
  geom_line(data=df3, aes(x=log(x), y=log(y3), color='y3'), linetype = "dashed") + 
  geom_line(data=df4, aes(x=log(x), y=log(y4), color='y4'), linetype = "longdash") +
  geom_line(data=df5, aes(x=log(x), y=log(y5), color='y5'), linetype = "dashed") +
  xlab("log(number of sites)") +
  ylab("log(time)") + 
  labs(color='')  +
  scale_color_manual(values = colors, labels = c("M1", "M2", "M3", "Cubic", "Lineal"), 
                     guide = guide_legend(override.aes = list(
                     linetype = c("dashed", "dashed", "dashed", "longdash", "dashed"),
                     shape = c(1, 0, 2, NA, NA)))) + 
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold")) 
dev.off()





#================================================
#           MONTE CARLO SIMULATION
#================================================






#==============================================
#      Additional comparisons (vs GAM model)
#==============================================

# TPSR (GAM moel)
x <- loc[, 1]
y <- loc[, 2]
z <- val 


#Fit 2D thin-plate regression spline
fit_gam <- gam(z ~ s(x,y, bs="tp", k = 80), method="REML")
#fit_gam <- gam(z ~ s(x,y, bs="tp", k = 10), method="REML")

#Predict over a grid
# ng = 40
# x.pred <- seq(min(fit$epoints[, 1]), max(fit$epoints[, 1]), length.out = ng)
# y.pred <- seq(min(fit$epoints[, 2]), max(fit$epoints[, 2]), length.out = ng)
# xy <- expand.grid( x = x.pred, y = y.pred)

z.pred <- matrix(predict(fit_gam, newdata = data.frame(xy)), nrow = ng, ncol = ng)


# krig.fit <- Krig(loc[,1:2], val, theta=20)
# z.pred2 <- matrix(predict(krig.fit, newdata = xy),
#                  nrow = ng, ncol = ng)
# dim(z.pred2)
# dim(z.pred)

#Fitted points
fitpoints <- predict(fit_gam, se.fit = TRUE)
length(fitpoints)
length(z.pred)  # <- son los puntos de evaluacion!



#Show points and surface
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 3))
scatter3D(xy[, 1], xy[, 2], z.pred, bty = "g", pch = ".", cex = 4,
          theta = 132, phi = 25, 
          colkey = list(side = 4, length = 1, cex.axis = 1.6),
          clim=c(min(xy[, 2], z.pred), max(xy[, 2], z.pred)),
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="TPSR", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values", legend.plot = TRUE)

# Hmat model
scatter3D(xy[, 1], xy[, 2], solh_eval, bty = "g", pch = ".", cex = 4,
          theta = 132, phi = 25, 
          colkey = list(side = 4, length = 1, cex.axis = 1.6),
          clim=c(min(fit$epoints[, 2], solh_eval),max(fit$epoints[, 2],solh_eval)),
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M3", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values")

# clim=c(min(fit$epoints[, 2], solh_eval),max(fit$epoints[, 2],solh_eval))
# 
# colkey(side = 4, clim = range(solh_eval), add = TRUE)
# 
# scatter3D(fit$epoints[, 1], fit$epoints[, 2], solh_eval, bty = "g",
#           pch = 20, cex = 2, ticktype = "detailed",
#           # surf = list(x = sort(unique(fit$epoints[, 1])), y = sort(unique(fit$epoints[, 2])), z = solh_eval,  
#           #             facets = NA),
#           clim=c(min(fit$epoints[, 2], solh_eval),max(fit$epoints[, 2],solh_eval)))


# data_pred <- tps.kriged %>% as.data.frame 
# scatter3D(data_pred$x, data_pred$y, data_pred$var1.pred, bty = "g", pch = ".", cex = 4,
#           theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6),
#           col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
#           main="Kriging", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
#           ylab = "y",
#           zlab = "Values")



#==========================================
#         COMP. TIME (log scale)
#==========================================
# Una manera de hacerlo
df1 <- data.frame(nosites, time_mgcv)
#df2 <- data.frame(nosites, time_krig)
df3 <- data.frame(nosites, time_hmat)


# Merge data frames using dplyr package
newdata <- df1 %>%
  #left_join(df2, by = "nosites") %>%
  left_join(df3, by = "nosites")

head(newdata)

# newData <- melt(list(df1 = df1, df2 = df2, df3 = df3), id.vars = "x")
# head(newData)

colnames(df1) <- c("x", "y1")
#colnames(df2) <- c("x", "y2")
colnames(df3) <- c("x", "y3")
#colors <- c("y1" = "black", "y2" = "blue", "y3" = "red")
colors <- c("y1" = "black", "y3" = "red")
 


plot1 <- ggplot(data=df1, aes(x=log(x), y=log(y1), color = "y1")) + 
  geom_point(data=df1, aes(x=log(x), y=log(y1), color='y1'), size = 3, shape = 0) + 
  geom_line(data=df1, aes(x=log(x), y=log(y1), color='y1'),linetype = "dashed") + 
  #geom_point(data=df2, aes(x=log(x), y=log(y2), color='y2'), size = 3, shape = 1) + 
  #geom_line(data=df2, aes(x=log(x), y=log(y2), color='y2'),linetype = "dashed") +
  geom_point(data=df3, aes(x=log(x), y=log(y3), color='y3'), size = 3, shape = 2) + 
  geom_line(data=df3, aes(x=log(x), y=log(y3), color='y3'),linetype = "dashed") + 
  xlab("log(number of sites)") +
  ylab("log(time)") + 
  labs(color='')  +
  scale_color_manual(values = colors, labels = c("TPSR", "M3"), 
  #scale_color_manual(values = colors, labels = c("TPSR", "Kriging", "M3"), 
                     guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed"),
                     #guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed", "dashed"),
                     shape = c(0, 2)))) +
                     #shape = c(0, 1, 2)))) + 
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold")) 




#=================================================
#             COMP. ERROR (log scale)
#=================================================
library(dplyr)
df4 <- data.frame(nosites, error_mgcv)
#df5 <- data.frame(nosites, error_krig)
df6 <- data.frame(nosites, error_hmat)



# Merge data frames using dplyr package
newdata3 <- df4 %>%
  #left_join(df5, by = "nosites") %>%
  left_join(df6, by = "nosites")

head(newdata3)

#colors2 <- c("M1 - TPSR" = "black", "Kriging - M3" = "blue", "M1 - M3" = "red")
colors2 <- c("M1 - TPSR" = "black", "M1 - M3" = "red")

#pdf('plot8.pdf')
plot2 <- ggplot(newdata3, aes(x = log(nosites))) +
  geom_point(aes(x = log(nosites), y = log(error_mgcv), colour = "M1 - GAM"), size = 3, shape = 0) +
  geom_line(aes(y = log(error_mgcv), color = "M1 - TPSR"), linetype = "dashed") +
  #geom_point(y = log(error_krig), size = 3, color = "blue", shape = 1) +
  #geom_line(aes(y = log(error_krig), color = "M1 - Kriging"), linetype = "dashed") +
  geom_point(y = log(error_hmat), size = 3, color = "red", shape = 2) +
  geom_line(aes(y = log(error_hmat), color = "M1 - M3"), linetype = "dashed") +
  xlab("log(number of sites)") +
  ylab("log(Error value)") +
  labs(color='')  +
  scale_color_manual(#values = colors2, labels = c("M1 - TPSR", "M1 - Kriging", "M1 - M3"), 
                     values = colors2, labels = c("M1 - TPSR", "M1 - M3"),
                     # guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed", "dashed"),
                     #                                          shape = c(0, 1, 2)))) +
                     guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed"),
                                           shape = c(0, 2)))) +
  theme(legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=18),
        legend.position="top",
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.spacing.x = unit(0.5, 'cm')) 

grid.arrange(plot1, plot2, ncol= 2)





# Calculate confidence intervals (95% confidence level)
confidence_level <- 0.95
z_value <- qnorm(1 - (1 - confidence_level) / 2)
lower_bound <- as.numeric(z.pred) - z_value * fitpoints$se.fit
upper_bound <- as.numeric(z.pred) + z_value * fitpoints$se.fit

# Combine the results into a data frame
results <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  predicted_value = as.numeric(z.pred),
  lower_confidence = lower_bound,
  upper_confidence = upper_bound
)

plot1 <- ggplot(results, aes(x = x, y = y, fill = predicted_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey100", high = "gray60") +  # Customize color palette
  geom_contour(aes(z = lower_confidence), linetype = "dashed", color = "red") +
  geom_contour(aes(z = upper_confidence), linetype = "dashed", color = "blue") +
  labs(fill = "", title = "TPSR") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=18),
        legend.position = "right", legend.key.height= unit(2.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.box.spacing = unit(0.5, "cm"))





# Calculate confidence intervals (95% confidence level)
confidence_level <- 0.95
z_value <- qnorm(1 - (1 - confidence_level) / 2)
lower_bound <- as.numeric(tps.kriged$var1.pred) - z_value * tps.kriged$var1.var
upper_bound <- as.numeric(tps.kriged$var1.pred) + z_value * tps.kriged$var1.var

# Combine the results into a data frame
results2 <- data.frame(
  # x = xy[, 1],
  # y = xy[, 2],
  x = data_pred$x,
  y = data_pred$y,
  predicted_value = data_pred$var1.pred,
  lower_confidence = lower_bound,
  upper_confidence = upper_bound
)

plot2 <- ggplot(results2, aes(x = x, y = y, fill = predicted_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey100", high = "gray60") +  # Customize color palette
  geom_contour(aes(z = lower_confidence), linetype = "dashed", color = "red") +
  geom_contour(aes(z = upper_confidence), linetype = "dashed", color = "blue") +
  labs(fill = "", title = "Kriging") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=18),
        legend.position = "right", legend.key.height= unit(2.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.box.spacing = unit(0.5, "cm"))



# Calculate confidence intervals (95% confidence level)
confidence_level <- 0.95
z_value <- qnorm(1 - (1 - confidence_level) / 2)
lower_bound <- solh_eval - z_value * 0.05
upper_bound <- solh_eval + z_value * 0.05

results3 <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  # x = fit$epoints[, 1],
  # y = fit$epoints[, 2],
  predicted_value = as.numeric(solh_eval),
  lower_confidence = lower_bound,
  upper_confidence = upper_bound
)

plot3 <- ggplot(results3, aes(x = x, y = y, fill = predicted_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey100", high = "gray60", limits=c(min(results3$predicted_value), 1)) +  # Customize color palette
  geom_contour(aes(z = lower_confidence), linetype = "dashed", color = "red") +
  geom_contour(aes(z = upper_confidence), linetype = "dashed", color = "blue") +
  labs(fill = "", title = "M3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.x=element_text(size=14, face = "bold"),
        axis.title.x=element_text(size=14,face="bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=18),
        axis.text.y=element_text(size=14, face = "bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.position = "right", legend.key.height= unit(2.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.box.spacing = unit(0.5, "cm")) 

library(ggpubr)
grid.arrange(plot1, plot3,  ncol = 2)







#=========================
#        Table 7
#=========================
#grid_points <- data.frame(grid_points)

# exact <- testfunction(matrix(xy[, 1]), matrix(xy[, 2]))
# exact <- testfunction(matrix(xy[, 1]), matrix(xy[, 2]))
#exact_krig <- testfunction(matrix(grid_points[, 1]), matrix(grid_points[, 2])) 

rmse_mgcv <- rmse(as.numeric(z.pred), as.numeric(exact)); rmse_mgcv
#rmse_krig <- rmse(as.numeric(tps.kriged$var1.pred), as.numeric(exact_krig)); rmse_krig
rmse_solh <- rmse(as.numeric(solh_eval), as.numeric(exact)); rmse_solh



