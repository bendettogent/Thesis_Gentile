# MJD (Merton Jump Diffusion) Simulation Without Drift and Statistical Data Collection

start_time <- Sys.time()  # Record start time
source("library_func.R")   # Source external library functions
print("STARTING THE SIMULATION")

############################ Parameters: ##################################
M <- 10             # Number of simulations
n <- 5                 # Step limit
S_0 <- 10               # Starting asset price

v <- 0.5                # Percentage of non-zero returns
mu1 <- -0.4             # Mean of Brownian motion
std <- 0.1              # Standard deviation of Brownian motion

lambda <- 3.0           # Poisson parameter for the Compound Poisson Process (CPP)
mu2 <- -0.4/3.0         # Mean jump size
std2 <- 0.1             # Standard deviation of jump size

Nmax <- 10**8           # Safety break to avoid infinite loops
t_max <- 1.0            # Maximum time for the simulation
tempi <- c(0.1, 0.5, 1) # Times at which statistics are recorded

################## Output File Name Formatting ############################
file_name <- "MJD"
file_name <- paste(file_name, as.character(n), sep="_n")
file_name <- paste(file_name, as.character(M), sep="_M")
file_name <- paste(file_name, as.character(v), sep="_v")
file_name <- paste(file_name, as.character(mu1), sep="_mu1")
file_name <- paste(file_name, as.character(mu2), sep="_mu2")
file_name <- paste(file_name, as.character(lambda), sep="_lambda")
file_name <- paste(file_name, as.character(S_0), sep="_S0")
file_name <- paste(file_name, as.character(std), sep="_std1")
file_name <- paste(file_name, as.character(std2), sep="_std2")
file_name <- paste(file_name, as.character(t_max), sep="_tmax")

################## Initialize Vectors for Statistical Data #################
# Moments for the three different time points: t1=0.1s, t2=0.5s, t3=1s
x  <- c(0, 0, 0)    # 1st moment (mean)
x2 <- c(0, 0, 0)    # 2nd moment (variance)
x3 <- c(0, 0, 0)    # 3rd moment (skewness)
x4 <- c(0, 0, 0)    # 4th moment (kurtosis)

vecx_t1 <- c()  # List to store results at t1
vecx_t2 <- c()  # List to store results at t2
vecx_t3 <- c()  # List to store results at t3

#################### Actual Simulation Loop #########################
print("Simulation starts.")
jumps <- 0  # Counter for the total number of jumps

# Perform M simulations
for (kkk in 1:M) { 
  Z <- S_0          # Current asset price
  P <- c(S_0)       # Price vector over time
  T <- c(0)         # Time vector
  
  type_jump <- c(1) # Vector for jump types (for visualization)
  i <- 1
  
  while ((Z > 0) & (i <= Nmax) & (T[i]/n^2 <= t_max)) { 
    phi <- rexp(1)  # Random waiting time from exponential distribution
    Delta <- phi * v / (std^2 * Z^2)  # Time step based on Brownian motion
    
    prob_salto <- exp(-v / (std^2 * Z^2) * lambda / n^2)  # Jump probability
    dice <- runif(1)  # Random uniform variable for jump decision
    
    if (dice > prob_salto) {
      # If a jump occurs:
      chi <- rnorm(1, mean = mu2, sd = std2)  # Jump size (normally distributed)
      jump <- Z * (exp(chi) - 1)  # Jump magnitude
      l <- trunc(max(1, n * abs(jump)))  # Ladder size for jump
      sign_jump <- sign(jump)  # Direction of jump (up/down)
      jumps <- jumps + 1  # Increment jump counter
      
      # Create ladder for smooth transition
      for (j in 1:l) {  
        T <- c(T, T[i] + Delta / l)
        P <- c(P, P[i] + sign_jump)
        Z <- Z + sign_jump / n
        i <- i + 1
        type_jump <- c(type_jump, 2)  # Mark it as a jump
      }
    } else {
      # If no jump (common step):
      type_jump <- c(type_jump, "1")  # Mark it as a normal step
      V <- runif(1)
      V[V < (v * 0.5 - 0.5 * g_n(Z, T[i], v, std, n, mu1))] <- -1  # Move down
      V[V != -1 & V < v] <- 1  # Move up
      V[V != -1 & V != 1] <- 0  # No movement
      
      T <- c(T, T[i] + Delta)  # Update time
      P <- c(P, P[i] + V)  # Update price trajectory
      Z <- Z + V / n  # Update current price
      i <- i + 1
    }
    
    # Stop if time step becomes too small
    if (v / (std^2 * Z^2) < 1.0e-04) {
      i <- 2 * Nmax
      print("Price or other variables are growing too fast!")
    }
  }
  
  #################### Save Statistical Values ####################
  # Calculate macroscopic prices at different times (0.1s, 0.5s, 1s)
  xt1 <- S_f(tempi[1], P, T, n, S_0) 
  xt1[xt1 < 0] <- 0  # Ensure non-negative prices
  xt10 <- S_f(tempi[2], P, T, n, S_0)
  xt10[xt10 < 0] <- 0
  xt60 <- S_f(tempi[3], P, T, n, S_0)
  xt60[xt60 < 0] <- 0
  
  vecx_t1 <- c(vecx_t1, xt1)  # Collect data for time 0.1s
  vecx_t2 <- c(vecx_t2, xt10)  # Collect data for time 0.5s
  vecx_t3 <- c(vecx_t3, xt60)  # Collect data for time 1s
  
  new <- c(xt1, xt10, xt60)
  x <- x + new / M  # Update mean
  
  new <- c(xt1^2, xt10^2, xt60^2)
  x2 <- x2 + new / M  # Update variance
  
  new <- c(xt1^3, xt10^3, xt60^3)
  x3 <- x3 + new / M  # Update skewness
  
  new <- c(xt1^4, xt10^4, xt60^4)
  x4 <- x4 + new / M  # Update kurtosis
  
  # Print progress every 100 simulations
  if ((kkk %% 100) == 0) {
    frase <- paste("Finished kkk=", as.character(kkk), sep=" ")
    frase <- paste(frase, as.character(i), sep=" with last i:")
    print(frase)
  }
}

print("END")

######################## Simulation is over ###################

####################### Print Results #########################
print_stat(file_name, tempi, x, x2, x3, x4)  # Print statistical moments
print_stat_vecs(file_name, tempi, vecx_t1, vecx_t2, vecx_t3)  # Print raw data

####################### Print Duration #########################
end_time <- Sys.time()  # Record end time

# Calculate the duration of the simulation
duration <- difftime(end_time, start_time, units = "secs")

# Extract hours, minutes, and seconds from the duration
hours <- as.numeric(duration, units = "hours") %/% 1
minutes <- as.numeric(duration, units = "mins") %% 60 %/% 1
seconds <- as.numeric(duration) %% 60

# Format the duration as hh:mm:ss
formatted_duration <- sprintf("%02f:%02f:%02f", hours, minutes, seconds)
cat("Time elapsed:", formatted_duration, "hours:minutes:seconds\n")

# Print average number of jumps per simulation
print(jumps / kkk)
