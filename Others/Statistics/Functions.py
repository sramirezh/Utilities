import numpy as np
#Blocking
def blocking(A):
    """
    A is the initial matrix with an even row number, do it not only for the closest integer.
    """
    n,m=np.shape(A)
    B=np.zeros((0.5*n,m))
    for i in xrange(int(0.5*n)):
        B[i]=0.5*(A[2*i-1,:]+A[2*i,:])           
    return B



def statistics(M):
    """
    Returns the average, the Error squared and the variance.
    """
    n,m=np.shape(M)
    Av=np.sum(M,0)/n
    Error=np.zeros((1,m))
    var=np.zeros((1,m))
    for i in xrange(m):
        diff=M[:,i]-Av[i] #replace with an index
        var[0,i]=(np.sum(diff**2))/n 
        Error[0,i]=np.sqrt(var[0,i]/n) #Pag 9 Error of independent sampling
    return np.append( Av,[Error,var])
    
def Autocorrelation (A,av):
    n,m=np.shape(A)
    C=np.zeros((n,m)) #Autocorrelation
    corrTime=np.zeros(m)
    for k in xrange(m):
        i=0;
        tol=0.0001
        val=10
        sign=True
        while val>tol and sign==True:    
            for j in xrange(n-i):
                C[i,k]=C[i,k]+(A[j,k]-av[k])*(A[j+i,k]-av[k])
            C[i,k]=C[i,k]/(n-i)
            val=C[i,k]
            if i>0:
                sign=C[i,k]<C[i-1,k]
                corrTime[k]=corrTime[k]+C[i,k]/C[0,k]
            i+=1
    corrTime=0.5+corrTime
    return C, corrTime
    

        
