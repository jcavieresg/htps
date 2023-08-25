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
shape = 1 # shape parameter of the RBF

alpha_vec = seq(1, 100, by=1)                # alphas

time_full <- vector()           # Solve function
time_full_cg <- vector()        # CG 
time_hmat <- vector()           # CG + H-matrix

error_full_cg <- vector()       # Error CG
error_hmat <- vector()          # Error CG + H matrix

error_cvmat <- vector()          # Error CG + H matrix LOOCV
error_cvhmat <- vector()          # Error CG + H matrix LOOCV



for (k in 1:length(alpha_vec)){
  n <- 20
  N <- n*n
  lat <- runif(N, 0, 1)
  lon <- runif(N, 0, 1)
  loc <- cbind(lat, lon)
  s1 = c(2:(n-1),(n+1):(N-1))
  s2 = c(1,n,N)
  
  # values
  #val <- testfunction(matrix(loc[, 1]), matrix(loc[, 2])) 
  val <- testfunction_random(matrix(loc[, 1]), matrix(loc[, 2])) 
  
  list <- vector("list", length(alpha))
  list[[k]] <- val

  
  
  #=====================================================================
  # full matrix with solve() function of R
  #=====================================================================
  #exectime <- tic()
  tic <- tic()
  d0 <- DistanceMatrix(loc,loc)
  E0 <- radialFunction(d0,2,1.0, shape)
  E0pa <- E0 + alpha_vec[k]*diag(N)
  T0 <- cbind(1,loc)
  M0 <- rbind(cbind(E0pa,T0),cbind(t(T0),0,0,0))
  rhs0 <- c(val,c(0,0,0))
  sol0 <- solve(M0,rhs0)
  #exectime <- toc()
  toc <- toc()
  time_full[k] <- toc
  
  invK <- pinv(E0pa)
  error_cvmat[k] <- (invK%*%rhs0[1:N]) / diag(invK)


  # neval <- 80
  # x <- runif(neval, 0, 1)
  # y <- runif(neval, 0, 1)
  # epoints <- as.matrix(expand.grid(x, y))
  
  
  # =====================================================================
  # full matrix with CG solver
  # =====================================================================
  #exectime <- tic()
  # tic <- tic()
  # loc1 <- loc[s1,]
  # loc2 <- loc[s2,]
  # val1 <- val[s1]
  # val2 <- c(val[s2],c(0,0,0))
  
  # loc1 <- loc[-(3:length(N)), ]
  # loc2 <- tail(loc, 3)
  # val1 <- s1
  # val2 <- c(s2,c(0,0,0))
  
  # Calculamos distancias
  # distances11 <- DistanceMatrix(loc1,loc1)
  # E11 <- radialFunction(distances11, 2, 1, shape) + alpha_vec[k]*diag(N-3)
  # 
  # distances12 <- DistanceMatrix(loc1,loc2)
  # E12 <- radialFunction(distances12, 2, 1, shape)
  # distances22 <- DistanceMatrix(loc2,loc2)
  # E22 <- radialFunction(distances22, 2, 1, shape) + alpha_vec[k]*diag(3)
  # 
  # T1 <- cbind(1,loc1)
  # T2 <- cbind(1,loc2)
  # M1 <- cbind(E12, T1)
  # aux1 <- cbind(E22, T2)
  # aux2 <- cbind(t(T2), 0, 0, 0)
  # M2 <- rbind(aux1, aux2)
  # M2inv <- inv(M2)
  # S <- M1%*%M2inv
  
  # just for testing - full matrix and Schur complement
  #M_full <- E11 - S%*%t(M1)
  
  # lado derecho
  #rhs <- val1 - S%*%val2
  
  #y <- calculateTPS(loc1,rhs,t(M1),S,0.00002,0.0,20,alpha)
  # y <- calculateTPS_full(rhs, M_full)
  # delta2 <- solve(t(T2),-t(T1)%*%y)
  # a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
  #a = solve(T2, s2 - t(E12)%*%y - E22%*%delta2)
  
  #solf <- replicate(N,0)
  # solf <- replicate(N+3,0)
  # solf[s1] <- y
  # solf[s2] <- delta2
  #solf <- c(y, delta2)
  #solf[(N+1):(N+3)] <- a
  
  #exectime <- toc()
  #toc <- toc()
  #time_full_cg[k] <- exectime$toc - exectime$tic
  # time_full_cg[k] <- toc 
  # error_full_cg[k] <- norm(as.matrix(sol0)- as.matrix(solf))
  
  
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
  E22 <- radialFunction(distances22, 2, 1, shape) + alpha_vec[k]*diag(3)
  
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
  
  y <- calculateTPS_hmat(loc1,rhs,t(M1),S, Epsilon = 0.0001, Eta = 2, MinClusterSize = 20, alpha_vec[k])
  delta2 <- solve(t(T2),- t(T1)%*%y)
  a = solve(T2,val[s2] - t(E12)%*%y - E22%*%delta2)
  
  solh <- replicate(N+3,0)
  solh[s1] <- y
  solh[s2] <- delta2
  solh[(N+1):(N+3)] <- a
  
  toc <- toc()
  time_hmat[k] <- toc 
  error_hmat[k] <- norm(as.matrix(sol0)- as.matrix(solh))
  invK <- pinv(E12)
  error_cvhmat[k] <- (invK%*%rhs) / diag(invK)
}




#==========================================
#         Graph of log(logtimes)
#==========================================
# Una manera de hacerlo
df1 <- data.frame(alpha_vec, error_cvmat, error_cvhmat)

alpha_opt1 <- alpha_vec[which.min(error_cvmat)]
alpha_opt2 <- alpha_vec[which.min(error_cvhmat)]


# error_cvmat
ggplot(df1, aes(alpha_vec, error_cvmat)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = alpha_opt1, color = "red") + 
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold"))


# error_cvhmat
ggplot(df1, aes(alpha_vec, error_cvhmat)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = alpha_opt2, color = "red") + 
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold")) 


# error_hmat
ggplot(df1, aes(alpha_vec, error_hmat)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = alpha_opt2, color = "red") + 
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text=element_text(size=14, face = "bold"),
        legend.position="top",
        axis.title=element_text(size=14,face="bold")) 




library(mgcv)
library(plot3D)

#Set up sample data
set.seed(42)
x <- loc[, 1]
y <- loc[, 2]

plot(x, y)

#z <- exp(-(x^2 + y^2)) + rnorm(200, 0, 0.15)
z <- head(rhs0,-3)

#Fit 2D thin-plate regression spline
gam_fit <- gam(z ~ s(x, y, bs="tp"), method="REML")

gam_fit$fitted.values

#Predict over a grid
ng = 40
x.pred <- seq(min(fit$epoints[, 1]), max(fit$epoints[, 1]), length.out = ng)
y.pred <- seq(min(fit$epoints[, 2]), max(fit$epoints[, 2]), length.out = ng)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(gam_fit, newdata = xy), 
                 nrow = ng, ncol = ng)

gam_fit$smooth[[1]]$knots

#Fitted points
fitpoints <- predict(gam_fit)

#Show points and surface
scatter3D(x, y, z, pch = 18, cex = 1.5, 
          theta = 132, phi = 25, 
          xlab = "x", ylab = "y", zlab = "f(x,y)",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints),
          main = "Smoothest possible surface under thin-plate constraints")


scatter3D(x, y, z, pch = 18, cex = 1.5, 
          theta = 132, phi = 25, colkey = list(side = 4, length = 1, cex.axis = 1.6), 
          col = ramp.col(c("darkolivegreen1",  "darkolivegreen3", "dodgerblue2")),
          main="M1", xlab = "x", cex.main = 2.0, cex.lab = 1.8,
          ylab = "y",
          zlab = "Values" )


library(Metrics)

rmse_gam = rmse(z, gam_fit$fitted.values); rmse_gam
rmse_hmat = rmse(z, head(sol0,-3)); rmse_hmat

