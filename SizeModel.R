
# size-structure 1
nages <- 3
M <- 0.2   # natural mortality
fishing <- 0.6

# Spawn-recruit parameters
a <- 6.8
b <- 0.0068
spawn <- 1:1000
recruits <- (a*spawn)/(1+b*spawn)
plot(spawn,recruits)

#Growth-transition matrix (females)
#All rows sum to 1
G <- matrix(c(0.05,0.95,0,
              0,0.1,0.9,
              0.25,0.75,0),
            byrow = T,
            nrow=nages)

# Growth and survival matrices
A <- G%*%exp(-M * diag(nages))
B <- G%*%exp(-M-fishing * diag(nages))

# Nums at age in present year

# Recruitment of females and nums at age next year
#FemRecNext <- 100
# NfemNext <- Nfem%*%A
# NfemNext[1] <- NfemNext[1] + FemRecNext
# NfemNext

N <- vector() #total individuals
nyears <- 200
Nfem <- c(100,100,100)
FemInit <- Nfem
N[1] <- sum(FemInit)

for(y in 2:nyears){
  NfemNext <- Nfem%*%G
  print(NfemNext)
  S <- Nfem[nages]      # Spawners
  FemRecNext <- (a*S)/(1+b*S)
  
  #NfemNext[1] <- NfemNext[1] + FemRecNext # Recruitment
  
  Nfem <- NfemNext
  N[y] <- sum(Nfem)
}

plot(1:nyears,N)
# Apply fishing mortality
