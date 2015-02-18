# EQ.R
# Authors: Steve Martell & Matias Brachini

# An spatial age-structured model using a
# gravity model to model movements.

set.seed(1)


# MODEL DIMENSIONS
narea   <- 4
nages   <- 8
nyear   <- 20
mdims   <- c(nyear,nages,narea)

# From Carruthers paper
x       <- runif(narea);    # area of each region
parea   <- x/sum(x);        # relative area of each region (gravity weight)
omega   <- 1.4              # residency parameter.


# GRAVITY MODEL
G       <- matrix(parea,narea,narea,byrow=TRUE)
diag(G) <- diag(G) + omega
G       <- exp(G)
G       <- G/rowSums(G)
I       <- diag(1,narea)

# FISHING MORTALITY & SELECTIVITY
fk      <- rgamma(narea,3.5,40)     # Capture probability in area k
ah      <- rnorm(narea,nages/2)
sel     <- sapply(ah,function(x) plogis(1:nages,x))


# CACLCULATE TOTAL MORTALITY ARRAYS
m       <- rnorm(narea,0.3,0.01)    # Natural mortality in area k
S       <- diag(exp(-m),narea)      # Survival rate
AGS     <- array(0,c(mdims,narea))
for (i in 1:nyear)
{
    for (j in 1:nages)
    {
        S          = diag(exp(-m-fk*sel[j,]))
        AGS[i,j,,] = G %*% S
    }
}


# POPULATION MODEL
S       <- diag(exp(-m),narea)      # Survival diagonal matrix
# A       <- G %*% S
rj      <- exp(rnorm(narea))        # Average recruitment to area k
Nijk    <- array(0,mdims)

# initial states
i = 1
for (j in 1:nages) 
{
    if( j==1 )
    {
        Nijk[i,j,]=rj
    }
    if( j<nages )
    {
        Nijk[i,j+1,] = AGS[i,j,,] %*% Nijk[i,j,]
    }
    # plus group
    if( j==nages )
    {
        Nijk[i,j,] = -solve(AGS[i,j,,]-I,Nijk[i,j,])
    }
}
# update state variables
for (i in 1:nyear) 
{
    for (j in 1:nages) 
    {
        if( j==1 )
        {
            Nijk[i,j,]=rj
        }
        if( j<nages && i != nyear)
        {
            Nijk[i+1,j+1,] = AGS[i,j,,] %*% Nijk[i,j,]
        }
        # plus group
        if( j==nages && i != nyear)
        {
            Nijk[i+1,j,] = Nijk[i+1,j,] + AGS[i,j,,] %*% Nijk[i,j,]
        }
    }
}




# MARKS RELEASED IN AREA K
mk      <- rpois(narea,100)

# MARKS RECAPTURED AT t=1
# Rk is a matrix Rk[release_area,recapture_area]
pk      <- fk * A                   # Capture probability post movement-survival
Rk      <- apply(pk,1,function(x) rbinom(narea,size=mk,prob=x))

