# size-structure 1

nages <- 3
M <- 0.2   # natural mortality

#Growth-transition matrix (females)
#All rows sum to 1
G <- matrix(c(0.05,0.95,0,
              0,0.1,0.9,
              0,0,1),
            nrow=nages)

# Growth and survival matrix
A <- G%*%exp(-M * diag(nages))

# Nums at age in present year
Nfem <- c(1,1,1)

# Recruitment of females and nums at age next year
FemRecNext <- 100
NfemNext <- Nfem%*%A 
NfemNext[1] <- NfemNext[1] + FemRecNext
NfemNext

# Apply fishing mortality
