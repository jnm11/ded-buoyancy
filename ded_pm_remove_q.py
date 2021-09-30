# Generated with SMOP  0.41
from libsmop import *
# ded_pm_remove_q.m

    
    Ny=31
    Nz=52
    H=5
    W=6
    y=linspace(- W / 2,W / 2,Ny)
    z=linspace(- H / 2,H / 2,Nz)
    f=randn(Ny,Nz)


def remove_q_dep(y,z,f,n):
    Ny=len(y)
    Nz=len(z)
    y=np.reshape(y,(Ny,1))
    z=np.reshape(z,(1,Nz))
    q=np.atan2(yy,zz)
    

    A=np.zeros(n,n)
    for j in range(n):
        for k in range(n):
            A[j,k]=np.sum(np.sin((j+1)*q)*np.cos((k+1)*q)*f)
    for j1 in range(n):
        for k1 in range(n):
            for j2 in range(n):
                for k2 in range(n):
                    N[j1,k1,j2,k2]=np.sum(np.sin((j1+1)*q)*np.cos((k1+1)*q)*np.sin((j2+1)*q)*np.cos((k2+1)*q))
    
            
    
    A=np.reshape(A,(n*n,))
    N=np.reshape(N,(n*n,n*n))

    A=comm.reduce(A, op=MPI.SUM,0)
    N=comm.reduce(N, op=MPI.SUM,0)

    if mpirank==0: c=np.reshape(numpy.linalg.solve(N,A),n,n)
    else:          c=np.zeros(n,n)
    c = comm.bcast(c, root=0)
    
    for j in range(n):
        for k in range(n):
            f=f-c[i,j]*np.sin((j+1)*q)*np.cos((k+1)*q)
    return(f)
   

