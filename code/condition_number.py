import numpy as np


N      = 800
eps    = 1e-3

# N      = 3200
# eps    = 1e-3

L=1.0
dx=L/N
dt=dx/2
gmm    = 1.4
rostar = 1.0
pstar  = 1.0


K=np.zeros((N+2,N+2))
#Periodic BCs
K[0,0]  =1.0;   K[0,N-1]=-1.0 
K[N+1,N+1]=1.0; K[N+1,2]=-1.0 


for indi in range(1,N+1): #interior
    K[indi,indi]=1.0+2.0*(1/dx)**2
    K[indi,indi-1]=-(1/dx)**2
    K[indi,indi+1]=-(1/dx)**2

L=np.zeros((N+2,N+2))
#Periodic BCs
L[0,0]  =1.0;   L[0,N-1]=-1.0 
L[N+1,N+1]=1.0; L[N+1,2]=-1.0 

for indi in range(1,N+1): #interior
    L[indi,indi]=dx**2+2.0
    L[indi,indi-1]=-1.0
    L[indi,indi+1]=-1.0

M=np.zeros((N+2,N+2))
#Periodic BCs
M[0,0]  =1.0; M[0,N-1]=-1.0 
M[N+1,N+1]=1.0; M[N+1,2]=-1.0 

# (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
for indi in range(1,N+1): #interior
    M[indi,indi]=1.0+2.0*(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi-1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi+1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)

#Normalize by X
P=np.zeros((N+2,N+2))

#Periodic BCs
P[0,0]  =1.0; P[0,N-1]=-1.0 
P[N+1,N+1]=1.0; P[N+1,2]=-1.0 

for indi in range(1,N+1): #interior
    P[indi,indi]=dx**2+2.0*(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi-1]=-(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi+1]=-(dt/eps)**2*(gmm*pstar/rostar)

#Normalize by eps
Q=np.zeros((N+2,N+2))

#Periodic BCs
Q[0,0]  =1.0; Q[0,N-1]=-1.0 
Q[N+1,N+1]=1.0; Q[N+1,2]=-1.0 

for indi in range(1,N+1): #interior
    Q[indi,indi]=dx**2*eps**2+2.0*(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi-1]=-(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi+1]=-(dt)**2*(gmm*pstar/rostar)

#Normalize by eps
R=np.zeros((N+2,N+2))

#Periodic BCs
R[0,0]  =1.0;   R[0,N-1]=-1.0 
R[N+1,N+1]=1.0; R[N+1,2]=-1.0 

for indi in range(1,N+1): #interior
    R[indi,indi]=(dx*eps/dt)**2 + 2.0*(gmm*pstar/rostar)
    R[indi,indi-1]=-(gmm*pstar/rostar)
    R[indi,indi+1]=-(gmm*pstar/rostar)

print("PERIODIC BCs")
print("Not normalized laplacian ", np.linalg.cond(K))
print("Normalized laplacian     ", np.linalg.cond(L))
print("Not normalized           ", np.linalg.cond(M))
print("dx                       ", np.linalg.cond(P))
print("eps and dx               ", np.linalg.cond(Q))
print("eps and dx and dt        ", np.linalg.cond(R))



K=np.zeros((N+2,N+2))
#Dirichlet BCs
K[0,0]  =1.0 
K[N+1,N+1]=1.0

for indi in range(1,N+1): #interior
    K[indi,indi]=1.0+2.0*(1/dx)**2
    K[indi,indi-1]=-(1/dx)**2
    K[indi,indi+1]=-(1/dx)**2


L=np.zeros((N+2,N+2))
#Dirichlet BCs
L[0,0]  =1.0 
L[N+1,N+1]=1.0

for indi in range(1,N+1): #interior
    L[indi,indi]=dx**2+2.0
    L[indi,indi-1]=-1.0
    L[indi,indi+1]=-1.0



M=np.zeros((N+2,N+2))
#Dirichlet BCs
M[0,0]  =1.0 
M[N+1,N+1]=1.0

# (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
for indi in range(1,N+1): #interior
    M[indi,indi]=1.0+2.0*(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi-1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi+1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)

#Normalize by X
P=np.zeros((N+2,N+2))

#Dirichlet BCs
P[0,0]  =1.0 
P[N+1,N+1]=1.0

for indi in range(1,N+1): #interior
    P[indi,indi]=dx**2+2.0*(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi-1]=-(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi+1]=-(dt/eps)**2*(gmm*pstar/rostar)

#Normalize by eps
Q=np.zeros((N+2,N+2))

#Dirichlet BCs
Q[0,0]  =1.0 
Q[N+1,N+1]=1.0

for indi in range(1,N+1): #interior
    Q[indi,indi]=dx**2*eps**2+2.0*(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi-1]=-(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi+1]=-(dt)**2*(gmm*pstar/rostar)

#Normalize by eps
R=np.zeros((N+2,N+2))

#Dirichlet BCs
R[0,0]  =1.0 
R[N+1,N+1]=1.0

for indi in range(1,N+1): #interior
    R[indi,indi]=(dx*eps/dt)**2 + 2.0*(gmm*pstar/rostar)
    R[indi,indi-1]=-(gmm*pstar/rostar)
    R[indi,indi+1]=-(gmm*pstar/rostar)


print("DIRICHLET BCs")
print("Not normalized laplacian ", np.linalg.cond(K))
print("Normalized laplacian     ", np.linalg.cond(L))
print("Not normalized           ", np.linalg.cond(M))
print("dx                       ", np.linalg.cond(P))
print("eps and dx               ", np.linalg.cond(Q))
print("eps and dx and dt        ", np.linalg.cond(R))


K=np.zeros((N+2,N+2))
#Neumann BCs
K[0,0]  =1.0;   K[0,1]=-1.0 
K[N+1,N+1]=1.0; K[N+1,N]=-1.0 

for indi in range(1,N+1): #interior
    K[indi,indi]=1.0+2.0*(1/dx)**2
    K[indi,indi-1]=-(1/dx)**2
    K[indi,indi+1]=-(1/dx)**2


L=np.zeros((N+2,N+2))
#Neumann BCs
L[0,0]  =1.0;   L[0,1]=-1.0 
L[N+1,N+1]=1.0; L[N+1,N]=-1.0 

for indi in range(1,N+1): #interior
    L[indi,indi]=dx**2+2.0
    L[indi,indi-1]=-1.0
    L[indi,indi+1]=-1.0


M=np.zeros((N+2,N+2))
#Neumann BCs
M[0,0]  =1.0;   M[0,1]=-1.0 
M[N+1,N+1]=1.0; M[N+1,N]=-1.0 


# (delta_t/EPS_LM/dx)**2*Gmm*min_p/max_ro
for indi in range(1,N+1): #interior
    M[indi,indi]=1.0+2.0*(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi-1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)
    M[indi,indi+1]=-(dt/eps/dx)**2*(gmm*pstar/rostar)

#Normalize by X
P=np.zeros((N+2,N+2))

#Neumann BCs
P[0,0]  =1.0;   P[0,1]=-1.0 
P[N+1,N+1]=1.0; P[N+1,N]=-1.0 


for indi in range(1,N+1): #interior
    P[indi,indi]=dx**2+2.0*(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi-1]=-(dt/eps)**2*(gmm*pstar/rostar)
    P[indi,indi+1]=-(dt/eps)**2*(gmm*pstar/rostar)

#Normalize by eps
Q=np.zeros((N+2,N+2))

#Neumann BCs
Q[0,0]  =1.0;   Q[0,1]=-1.0 
Q[N+1,N+1]=1.0; Q[N+1,N]=-1.0 

for indi in range(1,N+1): #interior
    Q[indi,indi]=dx**2*eps**2+2.0*(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi-1]=-(dt)**2*(gmm*pstar/rostar)
    Q[indi,indi+1]=-(dt)**2*(gmm*pstar/rostar)

#Normalize by eps
R=np.zeros((N+2,N+2))

#Neumann BCs
R[0,0]  =1.0;   R[0,1]=-1.0 
R[N+1,N+1]=1.0; R[N+1,N]=-1.0 

for indi in range(1,N+1): #interior
    R[indi,indi]=(dx*eps/dt)**2 + 2.0*(gmm*pstar/rostar)
    R[indi,indi-1]=-(gmm*pstar/rostar)
    R[indi,indi+1]=-(gmm*pstar/rostar)

print("NEUMANN BCs")
print("Not normalized laplacian ", np.linalg.cond(K))
print("Normalized laplacian     ", np.linalg.cond(L))
print("Not normalized           ", np.linalg.cond(M))
print("dx                       ", np.linalg.cond(P))
print("eps and dx               ", np.linalg.cond(Q))
print("eps and dx and dt        ", np.linalg.cond(R))
