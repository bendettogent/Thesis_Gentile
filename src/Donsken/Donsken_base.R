# Load the required functions from external file
source("library_func.R")

############################ Parameters ###############################
n <- 20             # Number of intervals (discretization parameter)
t_max <- 1.0        # Maximum time for the simulation
S_0 <- 50           # Initial value (starting stock price)
std <- 1            # Standard deviation for the diffusion process
mu <- 5             # Drift coefficient
Nmax <- 10^8        # Safety break for the maximum number of iterations

###################### Simulation of the Process #######################

# Initialize the process values
Z <- S_0            # Current stock price
P <- c(S_0)         # Vector to store stock price trajectory
T <- c(0)           # Vector to store time points

# Start simulation loop
i <- 1              # Initialize iteration counter
while ((Z > 0) & (i <= Nmax) & (T[i] / n^2 <= t_max)) {
  # Generate a random increment using uniform distribution 
  V <- runif(1)
  threshold <- 0.5 - 0.5 * mu / n
  
  # Update increment based on threshold comparison
  if (V < threshold) {
    V <- -1
  } else {
    V <- 1
  }
  
  # Update time and price vectors
  T <- c(T, T[i] + 1)
  P <- c(P, P[i] + V)
  
  # Update current price Z
  Z <- Z + V / n
  
  # Increment the iteration counter
  i <- i + 1
}

###################### Plotting the Results ############################

# Compute St and Sx using external functions from library
St <- S_t(P, T, n, S_0)  # Time points transformed according to model
Sx <- S_x(P, T, n, S_0)  # Stock price trajectory transformed

# Plotting the results
data <- data.frame(St, Sx)

plot(data$St, data$Sx, type = "l", col = "#69b3a2", lwd = 2, 
     xlab = "St", ylab = "Sx", main = "Line Plot")
box()