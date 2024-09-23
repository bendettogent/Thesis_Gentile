# Black-Scholes NO DRIFT - Base (single realization)
# Load the required functions
source("library_func.R")

##### Parameters #####
n <- 5              # n limit
v <- 0.8            # percentage of non 0-returns
S_0 <- 10           # Starting price
std <- 1            # Diffusion coefficient
mu <- 0.8           # Macroscopic mean
t_max <- 1.0        # Time of the simulation
Nmax <- 10^8        # Safety break iterations limit

##### Simulation of the Process #####

# Initialize the process values
Z <- S_0            # Current stock price
P <- c(S_0)         # Vector to store stock price trajectory
T <- c(0)           # Vector to store time points

# Start simulation loop
i <- 1
while ((Z > 0) & (i <= Nmax) & (T[i] / n^2 <= t_max)) {
  # Generate random increments
  phi <- rexp(1)
  Delta <- phi * v / (std^2 * Z^2)
  
  # Generate random step based on thresholds
  V <- runif(1)
  V[V < (v * 0.5 - 0.5 * g_n(Z, T[i], v, std, n, mu))] <- -1
  V[V != -1 & V < v] <- 1
  V[V != -1 & V != 1] <- 0
  
  # Update time and price vectors
  T <- c(T, T[i] + Delta)
  P <- c(P, P[i] + V)
  
  # Update current price Z
  Z <- Z + V / n
  
  # Check if stopping criterion is met
  if (v / (std^2 * Z^2) < 1.0e-05) {
    i <- 2 * Nmax
    print("Price is growing too fast!")
  }
  
  i <- i + 1
}

##### Plotting the Results #####

# Compute St and Sx using external functions from library
St <- S_t(P, T, n, S_0)  # Time points transformed according to model
Sx <- S_x(P, T, n, S_0)  # Stock price trajectory transformed

# Plotting the results
data <- data.frame(St, Sx)

plot(data$St, data$Sx, type = "l", col = "#69b3a2", lwd = 2, 
     xlab = "St", ylab = "Sx", main = "Line Plot")
box()
