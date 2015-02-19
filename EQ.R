# EQ.R
# Authors: Steve Martell & Matias Brachini

# An spatial age-structured model using a
# gravity model to model movements.

set.seed(1)


# MODEL DIMENSIONS
narea   <- 4
nages   <- 8
nyear   <- 20
nseas   <- 4
mdims   <- c(nyear,nseas,nages,narea)

# From Carruthers paper
x       <- runif(narea);    # area of each region
parea   <- x/sum(x);        # relative area of each region (gravity weight)
omega   <- 0.4              # residency parameter.

# SEASONALITY TO MOVEMENT MATRIX
h  = runif(narea)*7 # peak gravity weight.
b  = rep(0.2,narea) # offset of the peak.
ss = 1:nseas
gw = sapply(h,function(h) h*sin((ss-b)/(nseas)*pi)^2)


# GRAVITY MODEL
G       <- matrix(parea,narea,narea,byrow=TRUE)
diag(G) <- diag(G) + omega
G       <- exp(G)
G       <- G/rowSums(G)
I       <- diag(1,narea)

G       <- array(0,c(nseas,narea,narea))
for (ii in 1:nseas) 
{
    G[ii,,] <- matrix(gw[ii,],narea,narea,byrow=TRUE)
    diag(G[ii,,]) <- diag(G[ii,,]) + omega
    G[ii,,]  <- exp(G[ii,,])/rowSums(exp(G[ii,,]))

    G[ii,,]  <- 0
    diag(G[ii,,]) <- 1
}



# FISHING MORTALITY & SELECTIVITY
fk      <- 0*rgamma(narea,3.5,40)     # Capture probability in area k
ah      <- rnorm(narea,nages/2)
sel     <- sapply(ah,function(x) plogis(1:nages,x))


# CACLCULATE TOTAL MORTALITY ARRAYS
m       <- rnorm(narea,0.3,0.01)    # Natural mortality in area k
S       <- diag(exp(-m),narea)      # Survival rate
AGS     <- array(0,c(mdims,narea))
for (i in 1:nyear)
{
    for (ii in 1:nseas) 
    {
        for (j in 1:nages)
        {
            if( i == 5 ) fm = c(1,1,1,5) else fm = rep(1,narea)
            S          = diag(exp(-m/nseas-fm*fk*sel[j,]))
            AGS[i,ii,j,,] = G[ii,,] %*% S
        }        
    }
}


# POPULATION MODEL
rj      <- exp(rnorm(narea))        # Average recruitment to area k
Nijk    <- array(0,mdims)

# initial states
i = 1
for (j in 1:nages) 
{
    for (ii in 1:nseas)
    {
        if( j==1 && ii == 1)
        {
            Nijk[i,ii,j,]=rj
        }
        if( j<nages && ii != nseas)
        {
            Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
        }
        else if( j<nages && ii==nseas)
        {
            Nijk[i,1,j+1,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]   
        }
        # plus group
        if( j==nages && ii != nseas)
        {
            Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
        }
        if( j==nages && ii == nseas)
        {
            # cat("Hellow \n")
            AT = AGS[i,1,j,,]*AGS[i,2,j,,]*AGS[i,3,j,,]*AGS[i,4,j,,]
            Nijk[i,ii,j,] = -solve(AT-I,Nijk[i,ii,j,])
            # cat(Nijk[i,ii,j,],"\n")
        }
    }    
}

# update state variables
for (i in 1:nyear) 
{
    for (ii in 1:nseas)
    {        
        for (j in 1:nages) 
        {
            if( j==1 && ii == 1 )
            {
                Nijk[i,ii,j,]=rj
            }
            if( j<nages && i != nyear && ii!=nseas )
            {
                Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
            }
            if( j<nages && i != nyear && ii==nseas )
            {
                Nijk[i+1,1,j+1,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
            }
            # plus group
            if( j==nages && i != nyear && ii!=nseas )
            {
                cat(ii,"\t",Nijk[i,ii,j,],"\n")
                Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
            }
            if( j==nages && i != nyear && ii==nseas)
            {
                cat(ii,"\t",Nijk[i,ii,j,],"\n")
                Nijk[i+1,1,j,] = Nijk[i+1,ii,j-1,] + AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
            }
        }
    }
}




# MARKS RELEASED IN AREA K
mk      <- rpois(narea,100)

# MARKS RECAPTURED AT t=1
# Rk is a matrix Rk[release_area,recapture_area]
pk      <- fk * AGS[1,1,1,,] # Capture probability post movement-survival
Rk      <- apply(pk,1,function(x) rbinom(narea,size=mk,prob=x))

