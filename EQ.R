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
omega   <- 1.4              # residency parameter.

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
}



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
    for (ii in 1:nseas) 
    {
        for (j in 1:nages)
        {
            if( i > 5 ) fm = c(1,0,1,0) else fm = rep(1,narea)
            S          = diag(exp(-m/nseas-fm*fk*sel[j,]))
            AGS[i,ii,j,,] = G[ii,,] %*% S
        }        
    }
}


# POPULATION MODEL
rj      <- exp(2+rnorm(narea))        # Average recruitment to area k
Nijk    <- array(0,mdims)

# initial states
i = 1
AT = array(1,c(nages,narea,narea))
NT = array(0,c(nages,narea))
for (j in 1:nages) 
{
    for (ii in 1:nseas) 
    {
        AT[j,,] = AT[j,,] * AGS[i,ii,j,,]
    }
}

for (j in 1:nages) 
{
    if( j== 1 )
    {
        NT[j,] = rj
    }
    if( j< nages )
    {
        NT[j+1,] = AT[j,,] %*% NT[j,]        
    }
    if( j==nages )
    {
        NT[j,] = -solve(AT[j,,]-I,NT[j,])
    }
    Nijk[1,1,j,] = NT[j,]
}

for (i in 1:nyear) 
{
    for (j in 1:nages) 
    {
        if(j == 1)
        {
            Nijk[i,1,j,] = rj
        }
        for (ii in 1:nseas) 
        {
            if( ii <  nseas && j < nages )
                Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]
         
            if( ii == nseas && j < nages )
                Nijk[i,1,j+1,]  = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]

            if( ii <  nseas && j == nages )
                Nijk[i,ii+1,j,] = AGS[i,ii,j,,] %*% Nijk[i,ii,j,]

            if( ii == nseas && j == nages )
                Nijk[i,1,j,] = Nijk[i,1,j,] + AGS[i,ii,j,,] %*% Nijk[i,ii,j,]                
        }
    }
}




# MARKS RELEASED IN AREA K
# This should be done as tag release groups.
# The definition of a group is number of tags 
# released in a given stratum (year/season/area).
# 
# Data frame for marks released.
# GRP YR_REL SEAS AREA SIZE 
#   1   1975    1    1   78 
# 
# Data frame for marks recaptured.
# FISHERY  GRP YR_RECAP SEAS AREA SIZE
#       1    1     1981    1    3   85
# 
# Tk(ngroups,narea,years)   = tags released.
# Rk(ngroups,fishery,years) = recaptured tags.
ngroups <- 1
Mk      <- array(0,c(nyear,nseas,nages,narea))
F       <- array(0,c(nages,narea))

for (i in 1:nyear)
{
    for (ii in 1:nseas)
    {
        for (j in 1:nages) 
        {
            F[j,] = fk * sel[j,]
            size  = ceiling(Nijk[i,ii,j,])
            Mk[i,ii,j,] = sapply(F[j,],function(x) rbinom(1,size=size,prob=x))
        }
    }
}
dfT <- NULL
grp = 0
siz <- function(a)
{
    l = 100*(1-exp(-0.2*a)) 
    l = l + rnorm(1,0,0.08)*l
    return(round(l,1))
} 
for (i in 1:nyear) 
for (ii in 1:nseas) 
{
    n = colSums(Mk[i,ii,,])
    for (j in 1:length(n))
    {
        if(n[j] > 0)
        {
            grp = grp + 1;
            na = Mk[i,ii,,j]
            for(age in 1:length(na))
            {
                if(na[age] > 0)
                for (k in 1:na[age])
                {
                    dfT=rbind(dfT,data.frame(GRP=grp,YR=i,SEAS=ii,AREA=j,SIZE=siz(age)))
                }
            }
        }
    }
}

# MARKS RECAPTURED 
Rk      <- array(0,c(nyear,nseas,nages,narea))

for (i in 1:nyear)
{
    for (j in 1:nages) 
    {
        for (ii in 1:nseas)
        {

        }
    }
}

# Rk is a matrix Rk[release_area,recapture_area]
# pk      <- fk * AGS[1,1,1,,] # Capture probability post movement-survival
# Rk      <- apply(pk,1,function(x) rbinom(narea,size=mk,prob=x))

