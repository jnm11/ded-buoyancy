# Generated with SMOP  0.41
from libsmop import *
# ded_pm_remove_q.m

    
    Ny=31
# ded_pm_remove_q.m:2
    Nz=52
# ded_pm_remove_q.m:3
    H=5
# ded_pm_remove_q.m:4
    W=6
# ded_pm_remove_q.m:5
    y=linspace(- W / 2,W / 2,Ny)
# ded_pm_remove_q.m:6
    z=linspace(- H / 2,H / 2,Nz)
# ded_pm_remove_q.m:7
    f=randn(Ny,Nz)
# ded_pm_remove_q.m:8
    dA=dot(W / Ny,H) / Nz
# ded_pm_remove_q.m:10
    yy,zz=ndgrid(y,z,nargout=2)
# ded_pm_remove_q.m:12
    q=atan2(yy,zz)
# ded_pm_remove_q.m:13
    n=8
# ded_pm_remove_q.m:14
    A=zeros(Ny,Nz,n,n)
# ded_pm_remove_q.m:15
    for j in arange(1,n).reshape(-1):
        for k in arange(1,n).reshape(-1):
            A[arange(),arange(),j,k]=multiply(sin(dot(j,q)),cos(dot(k,q)))
# ded_pm_remove_q.m:18
    
    A=reshape(A,dot(Ny,Nz),dot(n,n))
# ded_pm_remove_q.m:21
    f=ravel(f)
# ded_pm_remove_q.m:22
    c=numpy.linalg.solve(A,f)
# ded_pm_remove_q.m:23
    g=f - dot(A,c)
# ded_pm_remove_q.m:24
    g=reshape(g,Ny,Nz)
# ded_pm_remove_q.m:25
    mesh(g)
    f=f / (dot(dA,sum(ravel(f))))
# ded_pm_remove_q.m:27
    #/Users/jnm/smop/smop-master/smop/main.py