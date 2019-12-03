import numpy as np


#Blocking
def blocking(A):
    """
    Creates a reduced vector with the averages of consecutive pairs [see notes in ROME]
    Args:
        A is the initial vector
    Returns 
        B is the pair averaged vector
    """
    n = int(0.5*len(A))
    B = 0.5*(A[0::2][:n]+A[1::2][:n])      
    return B
    
def autocorrelation_error (A):
    """
    Computes an autocorrelation analysis
    Args:
        A matrix with the data in columns of a matrix A[n,m]
    Returns:
        C Correlation of the data 
        corrTime Correlation time of the data.
    """
    n,m=np.shape(A)
    C=np.zeros((n,m)) #Autocorrelation
    corrTime=np.zeros(m)
    for k in range(m):
        vector = A[:,k]
        i=0;
        tol=0.0001
        val=10
        sign=True
        while val>tol and sign == True:
            if i == 0:
                C[i,k] = np.cov(vector,vector)[0][1]
            else:
                C[i,k] = np.cov(vector[:-i],np.roll(vector,-i)[:-i])[0][1]
            val = C[i,k]
            if i>0:
                sign = C[i,k]<C[i-1,k]
                corrTime[k] = corrTime[k]+C[i,k]/C[0,k]
            i+= 1
    corrTime  = 0.5 + corrTime
    return C, corrTime

def blocking_error(data, error = False):
    """
    Returns a zero error and average for a column that has zero in all entries
    """
    n,m=np.shape(data)
    MaxIter=int(np.round(np.log(n)/np.log(2)))
    Results=np.zeros((MaxIter,3*m))
    Error=np.zeros((MaxIter,m)) 
    final_error=np.zeros(m)
    variance = np.var(data,axis=0)
    for j in range(m):
        i = 0
        diff = 1
        data1 = data[:,j]
        while diff>0 and i<MaxIter:    ###Always ascending results that's why it is required diff>0
            variance = np.var(data1, ddof = 1 )
            data1 = blocking(data1)
            Error[i,j] = np.sqrt((2.0**i/n)*variance)
            if i>0:
                diff = Error[i,j]-Error[i-1,j]
            i += 1
            
    for i in range(m):
        if np.size(np.nonzero(Error[:,i]))==0:  #This arises when the 
            final_error[i]=0
        else:
            final_error[i]=Error[np.max(np.nonzero(Error[:,i]))-1,i]
    
    if error == True:
        return final_error,Error

    else:
        return final_error

    
    

        
