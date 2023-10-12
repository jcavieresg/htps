rm(list = ls())
library(matlib)
library(pracma)
library(gridExtra)
library(scatterplot3d)
library(plot3D)
library(mgcv)


#===============================================================================================
#                        Main code of the Smoothing TPS and H-matrix
#===============================================================================================
set.seed(1234)                       # Fix the seed
library(Metrics)

# Initialize MPI
mpiinit()

alpha = 1 # smoothing parameter
shape = 1 # shape parameter of the RBF

K = seq(20, 40, by=2)               # Sequence of grids 20^2 to 40^2

nosites <- vector()             # N sites
time_full <- vector()           # Solve function
time_full_cg <- vector()        # CG 
time_hmat <- vector()           # CG + H-matrix
time_mgcv <- vector()           # CG + H-matrix
#time_krig <- vector()           # CG + H-matrix

error_full_cg <- vector()       # Error CG
error_hmat <- vector()          # Error CG + H matrix
error_mgcv <- vector()          # Error mgcv
#error_krig <- vector()          # Error Kriging


for (k in 1:length(K)){
  print(k)
  #=====================================================================
  # generate sites and values
  #=====================================================================
  #Non-structured sites
  n <- K[k]
  N <- n*n
  nosites[k] <- N
  lat <- runif(N, 0, 1)
  lon <- runif(N, 0, 1)
  loc <- cbind(lat, lon)
  s1 = c(2:(n-1),(n+1):(N-1))
  s2 = c(1,n,N)
  
  # values
  val <- testfunction(matrix(loc[, 1]), matrix(loc[, 2])) 
  val <- val + 0.03*rnorm(N)
  
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
  
  
  #=====================================================================
  # mgcv
  #=====================================================================
  tic <- tic()
  fit_mgcv <- gam(val ~ s(loc[, 1], loc[, 2], bs="tp", k = k), method="REML")
  solmgcv <- fit_mgcv$fitted.values
  toc <- toc()
  time_mgcv[k] <- toc
  error_mgcv[k] <- norm(as.matrix(head(sol0,-3)) - as.matrix(solmgcv))
  
  
  #=====================================================================
  # Kriging (for a large n, very slow!)
  #=====================================================================
  
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
  
  # # Distances
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
  
  # rhs
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
  delta2 <- solve(t(T2),- t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
  
  solh <- replicate(N+3,0)
  solh[s1] <- y
  solh[s2] <- delta2
  solh[(N+1):(N+3)] <- a
  
  toc <- toc()
  time_hmat[k] <- toc 
  error_hmat[k] <- norm(as.matrix(sol0)- as.matrix(solh))
  error_hmat[k] <- norm(as.matrix(head(sol0, -3))- as.matrix(head(solh, -3)))
}
#=======================================================================================================

  # Finalize MPI
  mpifinalize()




#===============================================================================================
#                                       Results 
#===============================================================================================

#==========
# Table 1
#==========
error_solcg <- norm(matrix(sol0) - matrix(solf)); error_solcg
error_solh <- norm(matrix(sol0) - matrix(solh)); error_solh






#==============================================
#      Additional comparisons (vs TPSR)
#==============================================

# locations
x <- loc[, 1]
y <- loc[, 2]
z <- val 


#Fit 2D thin-plate regression spline
fit_gam <- gam(z ~ s(x,y, bs="tp", k = 80), method="REML")
z_pred <- matrix(predict(fit_gam, newdata = data.frame(xy)), nrow = ng, ncol = ng)


#Fitted points
z_pred_sd <- predict(fit_gam, se.fit = TRUE)

#Show points and surface
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 3))
# GAM model
scatter3D(xy[, 1], xy[, 2], z_pred, bty = "g", pch = ".", cex = 4,
          theta = 132, phi = 25, 
          colkey = list(side = 4, length = 1, cex.axis = 1.6),
          clim=c(min(xy[, 2], z_pred), max(xy[, 2], z_pred)),
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="TPSR", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values", legend.plot = TRUE)

# M3 model
scatter3D(xy[, 1], xy[, 2], solh_eval, bty = "g", pch = ".", cex = 4,
          theta = 132, phi = 25, 
          colkey = list(side = 4, length = 1, cex.axis = 1.6),
          clim=c(min(xy[, 2], solh_eval),max(xy[, 2],solh_eval)),
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M3", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values")




#==========================================
#         COMP. TIME (log scale)
#==========================================
# Una manera de hacerlo
df1 <- data.frame(nosites, time_mgcv)
df2 <- data.frame(nosites, time_hmat)

# Merge data frames using dplyr package
newdata <- df1 %>%
  left_join(df2, by = "nosites") 

colnames(df1) <- c("x", "y1")
colnames(df2) <- c("x", "y2")
colors <- c("y1" = "black", "y2" = "red")


plot1 <- ggplot(data=df1, aes(x=log(x), y=log(y1), color = "y1")) + 
  geom_point(data=df1, aes(x=log(x), y=log(y1), color='y1'), size = 3, shape = 0) + 
  geom_line(data=df1, aes(x=log(x), y=log(y1), color='y1'),linetype = "dashed") + 
  geom_point(data=df2, aes(x=log(x), y=log(y2), color='y2'), size = 3, shape = 2) + 
  geom_line(data=df2, aes(x=log(x), y=log(y2), color='y2'),linetype = "dashed") + 
  xlab("log(number of sites)") +
  ylab("log(time)") + 
  labs(color='')  +
  scale_color_manual(values = colors, labels = c("TPSR", "M3"), 
                     guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed"),
                     shape = c(0, 2)))) +
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
df5 <- data.frame(nosites, error_hmat)



# Merge data frames using dplyr package
newdata3 <- df4 %>%
  left_join(df5, by = "nosites") 

colors2 <- c("M1 - TPSR" = "black", "M1 - M3" = "red")


plot2 <- ggplot(newdata3, aes(x = log(nosites))) +
  geom_point(aes(x = log(nosites), y = log(error_mgcv), colour = "M1 - GAM"), size = 3, shape = 0) +
  geom_line(aes(y = log(error_mgcv), color = "M1 - TPSR"), linetype = "dashed") +
  geom_point(y = log(error_hmat), size = 3, color = "red", shape = 2) +
  geom_line(aes(y = log(error_hmat), color = "M1 - M3"), linetype = "dashed") +
  xlab("log(number of sites)") +
  ylab("log(Error value)") +
  labs(color='')  +
  scale_color_manual(#values = colors2, labels = c("M1 - TPSR", "M1 - Kriging", "M1 - M3"), 
    values = colors2, labels = c("M1 - TPSR", "M1 - M3"),
    guide = guide_legend(override.aes = list(linetype = c("dashed", "dashed"),
                                             shape = c(0, 2)))) +
  theme(legend.title = element_text(size=18, face = "bold"),
        legend.text = element_text(size=18),
        legend.position="top",
        axis.text=element_text(size=14, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.spacing.x = unit(0.5, 'cm')) 

grid.arrange(plot1, plot2, ncol= 2)





# Confidence intervals (95% confidence level)
conf_level <- 0.95
z_value <- qnorm(1 - (1 - conf_level) / 2)
lower_bound <- as.numeric(z_pred) - z_value * z_pred_sd$se.fit
upper_bound <- as.numeric(z_pred) + z_value * z_pred_sd$se.fit

# Combine the results into a data frame
results <- data.frame(x = xy[, 1], 
                      y = xy[, 2],
                      pred_value = as.numeric(z_pred), 
                      lower_conf = lower_bound, 
                      upper_conf = upper_bound)

plot1 <- ggplot(results, aes(x = x, y = y, fill = pred_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey100", high = "gray60") +  # Customize color palette
  geom_contour(aes(z = lower_conf), linetype = "dashed", color = "red") +
  geom_contour(aes(z = upper_conf), linetype = "dashed", color = "blue") +
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





# Confidence intervals (95% confidence level)
conf_level <- 0.95
z_value <- qnorm(1 - (1 - conf_level) / 2)
lower_bound <- solh_eval - z_value * 0.03  # 0.03 from the MC simulations
upper_bound <- solh_eval + z_value * 0.03

results3 <- data.frame(x = xy[, 1],
                       y = xy[, 2],
                       pred_value = as.numeric(solh_eval),
                       lower_conf = lower_bound,
                       upper_conf = upper_bound)

plot3 <- ggplot(results3, aes(x = x, y = y, fill = pred_value)) +
  geom_tile() +
  scale_fill_gradient(low = "grey100", high = "gray60", limits=c(min(results3$pred_value), 1)) +  # Customize color palette
  geom_contour(aes(z = lower_conf), linetype = "dashed", color = "red") +
  geom_contour(aes(z = upper_conf), linetype = "dashed", color = "blue") +
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
rmse_mgcv <- rmse(as.numeric(z.pred), as.numeric(exact)); rmse_mgcv
rmse_solh <- rmse(as.numeric(solh_eval), as.numeric(exact)); rmse_solh





