import numpy as np


#Blocking
def blocking(A):
    """
    A is the initial matrix with an even row number, do it not only for the closest integer.
    """
    n,m=np.shape(A)
    B=np.zeros((int(0.5*n),m))
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
    
def autocorrelation_error (A,av):
    """
    Computes an autocorrelation analysis
    Args:
        A matrix with the data in columns of a matrix A[n,m]
        av, the average of the data, av[m]
    Returns:
        C Correlation of the data 
        corrTime Correlation time of the data.
    """
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

def blocking_error(data):
    """
    Returns a zero error and average for a column that has zero in all entries
    """
    n,m=np.shape(data)
    MaxIter=int(np.round(np.log(n)/np.log(2)))
    Results=np.zeros((MaxIter,3*m))
    Error=np.zeros((MaxIter,m)) 
    final_error=np.zeros(m)
    for j in xrange(m):
        i=0
        diff=1
        data1=data
        while diff>0 and i<MaxIter:    ###Always ascending results that's why it is required diff>0
            Results[i,:]=statistics(data1)
            data1=blocking(data1)
            Error[i,j]=np.sqrt((2.0**i/n)*Results[i,2*m+j])
            if i>0:
                diff=Error[i,j]-Error[i-1,j]
            i+=1
            
    for i in xrange(m):
        if np.size(np.nonzero(Error[:,i]))==0:  #This arises when the 
            final_error[i]=0
        else:
            final_error[i]=Error[np.max(np.nonzero(Error[:,i]))-1,i]
        
    return final_error

    
    

        
