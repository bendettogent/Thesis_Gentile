#### Library with functions needed in all models

##########Function for macroscopic scale limits#####################
#inputs: microscopic tbt data P and T, makes the scaling limits:

#S_f: is macroscopic process as a function of time, single value
S_f <- function(t,P,T,n,S_0){
  if(t == 0){
    x <- P[1]
  }else{
    y <- which(T == max(T[T <= n*n*t]))
    x <- S_0 + (P[y]-S_0)/n
  }
  return(x)
}

#S_x: returns macro vector of prices movements
S_x <- function(P,T,n,S_0){
  x <- S_0 + (P-S_0)/n
  return(x)
}

#S_x: returns macro vector of times movements
S_t <- function(P,T,n,S_0){
  t <- T/n**2
  return(t)
}



####################General functions for DRIFT models ######################
f <- function(t,x,v,std){
  if (x > 0){
    return(v/(std*std*x*x))
  }else{
    return(1.)
  }
}

g <- function(x,t,v,std,b){
  if (x > 0){
    return(b*v/(std*std*x))
  }else{
    #return(b*v/(std*std*x))
    return(0.)
    print("Uscito g")
  }
}

g_n <- function(x,t,v,std,n,b){
  G <- g(x,t,v,std,b)
  if (abs(G) <= (n*v)){
    return(G/n)
  }else{
    return(0.)
    print("Uscito g_n")
  }
}

my_gamma <- function(x,t,z){
  return(3*z/S_0)
}

my_zeta <- function(x,t,z,G){
  #x è la sua Z
  #z è il \chi
  if (G==1){
    return(x*my_gamma(x,t,z))
  }else{
    return(0.)
  }
}

func_l <- function(x,n){
  if(n*abs(x) > 1){
    return(as.integer(n*abs(x)))
  }else{
    return(1.)
  }
}



####################Print functions################
print_stat <- function(prefix,tempi,x,x2,x3,x4){
  outfilename <- paste("data/",prefix, sep = "")
  outfilename <- paste(outfilename,".dat", sep = "")
  
  outfile <- file(outfilename, "w")
  
  #Stampo intro:
  frase <- paste("### ",prefix, sep = "")
  writeLines(frase, outfile)
  writeLines("i t x x2 x3 x4", outfile)
  
  N <- length(x)
  for (i in (1:length(x))){
    
    line <- as.character(i)
    
    a <- as.character(tempi[i])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x[i])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x2[i])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x3[i])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x4[i])
    line <- paste(line,a, sep = " ")
    
    writeLines(line, outfile)
  }
  close(outfile)
}

print_traj <- function(prefix,x,t){
  outfilename <- paste("trajectories/",prefix, sep = "")
  outfilename <- paste(outfilename,".dat", sep = "")
  outfile <- file(outfilename, "w")
  
  #Stampo intro:
  frase <- paste("### ",prefix, sep = "")
  writeLines(frase, outfile)
  writeLines("i t x", outfile)
  
  N <- length(x)
  for (i in (1:length(x))){
    line <- as.character(i)
    a <- as.character(t[i])
    line <- paste(line,a, sep = " ")
    a <- as.character(x[i])
    line <- paste(line,a, sep = " ")
    writeLines(line, outfile)
  }
  close(outfile)
}


print_stat_vecs <- function(prefix,tempi,x1,x2,x3){
  outfilename <- paste("data_vec/",prefix, sep = "")
  outfilename <- paste(outfilename,"_vectors.dat", sep = "")
  
  outfile <- file(outfilename, "w")
  
  #Stampo intro:
  frase <- paste("### ",prefix, sep = "")
  writeLines(frase, outfile)
  writeLines("i t1 x1 t2 x2 t3 x3", outfile)
  
  N <- length(x1)
  for (i in (1:length(x1))){
    
    line <- as.character(i)
    
    a <- as.character(tempi[1])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x1[i])
    line <- paste(line,a, sep = " ")

    a <- as.character(tempi[2])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x2[i])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(tempi[3])
    line <- paste(line,a, sep = " ")
    
    a <- as.character(x3[i])
    line <- paste(line,a, sep = " ")
    
    writeLines(line, outfile)
  }
  close(outfile)
}

