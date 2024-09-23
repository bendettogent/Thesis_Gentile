#Black-sholes DRIFT - Statistical data ()
# Program to be run from command line (or swich parameter Command_line to 0), 
#example:
#Rstudio BSnodrift_statistica.R 10000 10 > output.dat &
#We write the results on external files in the folders:
#data and data_traj (see library_func.R  functions print...)
#

#Start the clock:
start_time <- Sys.time()
#include library with some functions:
source("library_func.R")

############################Parameters(OLD):################################
Command_line <- 0 #1 is for command_line use, 0 otherwise,

tempi <- c(0.1,0.5,1) #times to consider for statistics, t in [0,1]

v <- 0.8      #percentage of 0-returns 
S_0 <- 10     #Starting prize
std <- 0.1    #diffution

Nmax <- 10**8 #Safety break iterations limit
if(Command_line){
  args <- commandArgs(trailingOnly = TRUE)
  M <- args[1]  #number of simulations
  n <- args[2]  #n limit
  n <- as.numeric(n)
  M <- as.numeric(M)
}

t_max = max(tempi) #time max that i need to simulate

############################Parameters: ##################################
M <- 10000
n <- 120           # n limit
v <- 0.5           # percentage of non 0-returns 
S_0 <- 10          # Starting prize
std <- 0.1         # diffution
mu <- 0.0          # macroscopic mean

t_max = 1.0        # time of the simulation
Nmax <- 10**8      # Safety break iterations limit

##################output print formatting:##########################
file_name <- "BSdrift"
file_name <- paste(file_name,as.character(n),sep="_n")
file_name <- paste(file_name,as.character(M),sep="_M")
file_name <- paste(file_name,as.character(v),sep="_v")
file_name <- paste(file_name,as.character(mu),sep="_mu")
file_name <- paste(file_name,as.character(S_0),sep="_S0")
file_name <- paste(file_name,as.character(std),sep="_std")
file_name <- paste(file_name,as.character(t_max),sep="_tmax")


##################Inizialize vectors for statistics###############
#3D per t1=1s, t2=10s and t3=60s
x  <- c(0,0,0)   #1st momentum
x2 <- c(0,0,0)   #2nd momentum
x3 <- c(0,0,0)   #3rd momentum
x4 <- c(0,0,0)   #4th momentum

vecx_t1 <- c()   #complete list of M values in t1 of tempi
vecx_t2 <- c()   #complete list of M values in t2 of tempi
vecx_t3 <- c()   #complete list of M values in t3 of tempi



####################Actual Simulation#########################
print("Simulation starts.")

#external loop on M, number of repetitions of simulations
for (kkk in (1:M)){ 
  Z <- S_0           #tbt prize right now
  P <- c(S_0)        #tbt prize vector
  T <- c(0)          #tbt times vector
  
  #Ciclo while to construct returns:
  #caveat: in the theory i have (i-1), here not because 
  #Z[0](theory) -> Z[1](simulation) and U[1](theory) -> U[1](simulation)
  
  i <- 1 
  
  #while loop with 3 conditions:
  #  Z > 0: for the price must be non negative,
  # i <= Nmax: Safety break iterations limi
  # (T[i]/n**2<=t_max): t_max is max time i need
  while ((Z > 0) & (i <= Nmax) & (T[i]/n**2<=t_max)) { 
    
    #random increments:
    phi <- rexp(1)
    
    V <- runif(1)
    V[V < (v*0.5-0.5*g_n(Z,T[i],v,std,n,mu))] <- -1
    V[V != -1 &  V < v] <- 1
    V[V != -1 &  V != 1] <- 0
    
    Delta <- phi*v/(std**2 * Z**2) #Delta in time
    
    T <- c(T,T[i]+Delta)
    P <- c(P,P[i]+V)
    
    Z <- Z + V/n #update Z
    
    i <- i+1
    
    if(v/(std**2 * Z**2) < 1.0e-05){ #stop se delta troppo piccolo!
      i <- 2*Nmax
      print("sta crescendo troppo il prezzo o altre cose!")
    }
  }
  if(Z<=0){
    print("exit here 1")
  }
  
  #save statistical values:
  #S_f gives macroscopic prize at macroscopic time
  xt1 <- S_f(tempi[1],P,T,n,S_0) 
  xt10 <- S_f(tempi[2],P,T,n,S_0)
  xt60 <- S_f(tempi[3],P,T,n,S_0)
  
  vecx_t1 <- c(vecx_t1,xt1)
  vecx_t2 <- c(vecx_t2,xt10)
  vecx_t3 <- c(vecx_t3,xt60)
  
  new <- c(xt1,xt10,xt60)
  x <- x + new/M
  
  new <- c(xt1**2,xt10**2,xt60**2)
  x2 <- x2 + new/M
  
  new <- c(xt1**3,xt10**3,xt60**3)
  x3 <- x3 + new/M
  
  new <- c(xt1**4,xt10**4,xt60**4)
  x4 <- x4 + new/M
  
  
  #Message on the output:
  if((kkk %% 100) == 0){
    frase <- paste("Finished kkk=",as.character(kkk),sep=" ")
    frase <- paste(frase,as.character(i),sep=" with last i:")
    print(frase)
  }
}

print(" ")
print("END")

########################Simulation is over###################

#######################Print results#########################
print_stat(file_name,tempi,x,x2,x3,x4)
print_stat_vecs(file_name,tempi,vecx_t1,vecx_t2,vecx_t3)

#######################Print duration#########################
end_time <- Sys.time()
duration <- end_time - start_time
cat("Time elapsed:", duration, "seconds\n")
#Done#