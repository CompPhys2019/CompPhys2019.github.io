from numpy import linalg as LA
import numpy as np

def fr(n,dtype):
    H = np.zeros((n,n),dtype=dtype)
    for i in range(n):
        for j in range(n):
            H[i,j]=1/(i+j+1)

    eval=np.abs(LA.eigvals(H))
    r = np.log(max(eval)/min(eval))
    return r

a = np.zeros((100,3))

for i in range(100):
    a[i,0]=i+2
    a[i,1]=fr(i+2,np.float32)
    a[i,2]=fr(i+2,np.float64)

np.savetxt('r.dat',a)