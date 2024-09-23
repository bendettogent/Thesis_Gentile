# Modello base per MJD process #

source("library_func.R")

############################Parameters: ##################################
n <- 40           #n limit
S_0 <- 10         #Starting prize

v <- 0.5          #percentage of non 0-returns 
mu1 <- 0.0       #mean Bm   
std <- 0.1        #diffution Bm

lambda <- 3.# Poisson parameter CPP
mu2 <- 0.0        # Mean jump size
std2 <- 1.0       # std jump size (not used)
mio_gamma <- 1.   # gamma in front of CPP (is related to std2)

Nmax <- 10**8     #Safety break iterations limit
t_max <- 1.0       #time of the simulation
####################Actual Simulation#########################
Z <- S_0          #tbt prize right now
P <- c(S_0)       #tbt prize vector
T <- c(0)         #tbt times vector

type_jump <- c(1)    #vector for data visualization

i <- 1
jumps <- 0

while ((Z > 0) & (i <= Nmax) & (T[i]/n**2<=t_max)) { 
  
  phi <- rexp(1)
  Delta <- phi*v/(std**2 * Z**2)
  
  prob_salto <- exp(- v/(std**2 * Z**2)*lambda/n**2)
  dice <- runif(1)
  if (dice > prob_salto){
    # In case I jump:
    chi <- rnorm(1, mean = mu2, sd = std2)
    jump <- Z*(exp(chi)-1)
    l <- trunc(max(1,n*abs(jump)))
    sign_jump <- sign(jump)
    jumps <- jumps + 1
    
    for(j in (1:l)){ #create ladder
      T <- c(T,T[i]+Delta/l)
      P <- c(P,P[i]+sign_jump)
      Z <- Z + sign_jump / n
      i <- i+1
      type_jump <- c(type_jump,2)
    }
  }else{ #common step
    
    type_jump <- c(type_jump,"1")
    
    V <- runif(1)
    V[V < (v*0.5-0.5*g_n(Z,T[i],v,std,n,mu1))] <- -1
    V[V != -1 &  V < v] <- 1
    V[V != -1 &  V != 1] <- 0

    T <- c(T,T[i]+Delta)
    P <- c(P,P[i]+V)
    
    Z <- Z + V/n 
    
    i <- i+1
  }
  
  if(v/(std**2 * Z**2) < 1.0e-04){ #stop se delta troppo piccolo!
    i <- 2*Nmax
    print("sta crescendo troppo il prezzo o altre cose!")
  }
}

######################## Post-Simulation Adjustments ######################
St <- S_t(P, T, n, S_0)  # Calculate time points for plot
Sx <- S_x(P, T, n, S_0)  # Calculate corresponding prices for plot

# Remove negative values if they exist
if (Z < 0) {
  count <- sum(Sx < 0)
  St <- head(St, -count)
  Sx <- head(Sx, -count)
  type_jump <- head(type_jump, -count)
}

######################## Visualization with Base Plot ####################
# Define colors for the plot: type_jump "1" = green (normal), "2" = red (jump)
colors <- ifelse(type_jump == "2", "red", "green")

# Create the plot
plot(St, Sx, type = "l", col = colors, xlab = "Time (St)", ylab = "Position (Sx)",
     main = "Plot of Sx vs St")

# Add points to highlight jumps
points(St[type_jump == "2"], Sx[type_jump == "2"], col = "red", pch = 19,cex=0.5)

legend("topright", legend = c("Normal Step", "Jump"), col = c("green", "red"), pch = c(NA, 19), lty = c(1, NA))
