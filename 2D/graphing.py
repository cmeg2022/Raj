import numpy as np
import matplotlib.pyplot as plt

#Sequence of Data Headers
#X,Y,sig11,sig22,sig12

print("Reading Data from analytical.dat")
#Import .dat file and plot data
data = np.loadtxt('analytical.dat')

radius = 1

#Filter
data = data[data[:,0] < radius*8]

plt.plot(data[:,0]/radius,data[:,2],label='sig11')
plt.plot(data[:,0]/radius,data[:,3],label='sig22')
plt.plot(data[:,0]/radius,data[:,4],label='sig12')
plt.legend()
plt.xlabel('X')
plt.ylabel('Sigma')
plt.title('Sig11, Sig22, Sig12')
plt.savefig('sig.png')
