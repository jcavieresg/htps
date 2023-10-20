#===============================================================================================
#                        Monte Carlo simulation STPS
#===============================================================================================

rm(list = ls())
library(matlib)
library(tictoc)
library(pracma)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(plot3D)
library(scales)


# Initialize MPI 
mpiinit() # -----> RUN THIS LINE JUST ONCE!


set.seed(1234)                       # Fix the seed
n_monte = 1000

mc_sol0_median = rep(NA, n_monte)
mc_sol0_iqr = rep(NA, n_monte)
mc_solh_median = rep(NA, n_monte)
mc_solh_iqr = rep(NA, n_monte)


for(i in 1:n_monte){

  alpha = 1 # smoothing parameter
  shape = 1 # shape parameter of the RBF
  
  #=====================================================================
  # generate sites and values
  #=====================================================================
  #Non-structured sites
  n <- 20
  N <- n*n
  lat <- runif(N, 0, 1)
  lon <- runif(N, 0, 1)
  loc <- cbind(lat, lon)
  s1 = c(2:(n-1),(n+1):(N-1))
  s2 = c(1,n,N)
  
  # values
  val <- testfunction_random(matrix(loc[, 1]), matrix(loc[, 2])) 
  val <- val + 0.03*rnorm(N)
  
  
  
  #=====================================================================
  # full matrix with solve() function of R
  #=====================================================================
  d0 <- DistanceMatrix(loc,loc)
  E0 <- radialFunction(d0,2,1.0, shape)
  E0pa <- E0 + alpha*diag(N)
  T0 <- cbind(1,loc)
  M0 <- rbind(cbind(E0pa,T0),cbind(t(T0),0,0,0))
  rhs0 <- c(val,c(0,0,0))
  sol0 <- solve(M0,rhs0)
  
  sol0_monte <- sol0[1:length(loc[,1]-3)]
  
  # Monte carlo part
  sol0_median = median(sol0_monte)
  mc_sol0_median[i] = sol0_median
  sol0_iqr = IQR(sol0_monte)
  mc_sol0_iqr[i] = sol0_iqr
  
  
  #=====================================================================
  # H matrix with CG solver
  #=====================================================================
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
  solh_monte <- solh[1:length(loc[,1]-3)]
  
  # Monte carlo part
  solh_median = median(solh_monte)
  mc_solh_median[i] = solh_median
  solh_iqr = IQR(solh_monte)
  mc_solh_iqr[i] = solh_iqr
}






#=================================================
# Histogram for the direct method 
hist(mc_sol0_median, col="limegreen")
hist(mc_sol0_iqr, col="orchid")
var(mc_sol0_median)
var(mc_sol0_iqr)

# Histogram for the CG method with an H-matrix
hist(mc_solh_median, col="limegreen")
hist(mc_solh_iqr, col="orchid")
var(mc_solh_median)
var(mc_solh_iqr)



# Calculate a 95% confidence interval for the mean of median
mean_output <- mean(mc_solh_median)
sd_output <- sd(mc_solh_median)

c_level <- 0.95
margin_of_error <- qnorm((1 + c_level) / 2) * (sd_output / sqrt(n_monte))
lower_bound <- mean_output - margin_of_error
upper_bound <- mean_output + margin_of_error

cat("95% confidence interval:", lower_bound, "to", upper_bound, "\n")



setwd("C:/Users/Usuario/Desktop/Thesis/article_1/figuras_newversion")

options("scipen"=2, "digits"=2)

data_plot <- data.frame(mc_sol0_median, mc_sol0_iqr, mc_solh_median, mc_solh_iqr)

p1 <- ggplot(data_plot, aes(mc_sol0_median)) +
  geom_histogram(color = "grey100", fill = "grey80") +
  geom_vline(aes(xintercept = mean(mc_sol0_median)), color = "black", size = 1.25) +
  geom_vline(aes(xintercept = mean(mc_sol0_median) + sd(mc_solh_median)), color = "black", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(mc_sol0_median) - sd(mc_solh_median)), color = "black", size = 1, linetype = "dashed") + 
  ggtitle(expression(paste("Median distribution for ", S(g), " in M1 "))) +
  xlab("") + ylab("Frequency") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold", color = "black"), 
        axis.text.x = element_text(color = "black", size = 14, hjust = 0.7),
        axis.text.y = element_text(color = "black", size = 14),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) # + 
        #scale_x_continuous(labels = format(mc_sol0_mean, digits=3))
  

p2 <- ggplot(data_plot, aes(mc_sol0_iqr)) +
  geom_histogram(color = "grey100", fill = "grey80") +
  geom_vline(aes(xintercept = mean(mc_sol0_iqr)), color = "black", size = 1.25) +
  geom_vline(aes(xintercept = mean(mc_sol0_iqr) + sd(mc_solh_iqr)), color = "black", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(mc_sol0_iqr) - sd(mc_solh_iqr)), color = "black", size = 1, linetype = "dashed") + 
  ggtitle(expression(paste("IQR distribution for ", S(g), " in M1 "))) +
  xlab("Values") + ylab("Frequency") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text.x = element_text(color = "black", size = 14, hjust = 0.7),
        axis.text.y = element_text(color = "black", size = 14),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))


p3 <- ggplot(data_plot, aes(mc_solh_median)) +
  geom_histogram(color = "grey100", fill = "grey80") +
  geom_vline(aes(xintercept = mean(mc_solh_median)), color = "black", size = 1.25) +
  geom_vline(aes(xintercept = mean(mc_solh_median) + sd(mc_solh_median)), color = "black", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(mc_solh_median) - sd(mc_solh_median)), color = "black", size = 1, linetype = "dashed") + 
  ggtitle(expression(paste("Median distribution for ", S(g), " in M3 "))) +
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold", color = "black"), 
        axis.text.x = element_text(color = "black", size = 14, hjust = 0.7),
        axis.text.y = element_text(color = "black", size = 14),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 



p4 <- ggplot(data_plot, aes(mc_solh_iqr)) +
  geom_histogram(color = "grey100", fill = "grey80") +
  geom_vline(aes(xintercept = mean(mc_solh_iqr)), color = "black", size = 1.25) +
  geom_vline(aes(xintercept = mean(mc_solh_iqr) + sd(mc_solh_iqr)), color = "black", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(mc_solh_iqr) - sd(mc_solh_iqr)), color = "black", size = 1, linetype = "dashed") + 
  ggtitle(expression(paste("IQR distribution for ", S(g), " in M3 "))) +
  xlab("Values") + ylab("") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text.x = element_text(color = "black", size = 14, hjust = 0.7),
        axis.text.y = element_text(color = "black", size = 14),  
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
        scale_fill_manual(name="group",values=c("red","darkgray"),labels=c("a"))
        

library(gridExtra)
pdf('mc_tps3.pdf')
grid.arrange(p1, p3, p2, p4, ncol = 2)
dev.off()



# Bootstrapping
# Perform bootstrapping to estimate standard error
n_boots <- 100
boots_means <- numeric(n_boots)

for (i in 1:n_boots) {
  boots_sample <- sample(mc_solh_median, replace = TRUE)
  boots_means[i] <- mean(boots_sample)
}

# Calculate standard error from bootstrapped means
se_boots_means <- sd(boots_means)

# Print results
cat("Mean from MC:", mean(mc_solh_median), "\n")
cat("Sd from bootstrapp:", se_boots_means, "\n")


# Finalize MPI
mpifinalize() # ----> RUN THIS LINE ONLY WHEN YOU WANT TO CLOSE THE SESSION
