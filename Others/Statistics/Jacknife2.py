import numpy as np
data=np.loadtxt("datablocked1.txt")
n=np.size(data)
size=n/5
data=np.reshape(data,(size,5))


JKaverage=np.zeros((size,5))
h=np.zeros((size,5))

#Jackknife Average
for i in range(5):
    for j in range(size):
        JKaverage[j,i]=(1.0/(size-1))*np.sum(np.delete(data[:,i],j))

average=(1.0/size)*np.sum(data,axis=0)

#Jacknife estimator 
hjk=np.zeros((5,5))
hest=np.zeros((5,5))
error=np.zeros((5,5))
for i in range(5):
    h=np.zeros((size,5))
    for j in range(5):
        h[:,j]=JKaverage[:,j]/JKaverage[:,i]
        hest[i,j]=average[j]/average[i]
    hjk[i,:]=(1.0/size)*np.sum(h,0)
    for j in range(5):
        error[i,j]=np.sqrt((size-1)*(1.0/size)*np.sum((h[:,j]-hjk[i,j])**2))

bias=hest-hjk

np.savetxt("hjk",hjk,fmt='%6.5f',delimiter=" & ")               
np.savetxt("bias",bias,fmt='%6.5e',delimiter=" & ") 
np.savetxt("errorJK",error,fmt='%6.5f',delimiter=" & ") 
#variance=np.zeros(5)
#
###Error estimation
#for i in xrange(size):
#    variance=variance+(h[i,:]-hjk)**2
#variance=variance/size
#Error=np.sqrt((size-1)*variance)

            
            
