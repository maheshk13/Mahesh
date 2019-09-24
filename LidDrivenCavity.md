import numpy as np 
import matplotlib.pyplot as plt

nx = 50   # Total no of nodes
ny=nx
dx = 2 / (nx - 1) # Length of the cavity is 2 units
dy = dx

Re=100
nu = 1/Re # coefficient of kinematic vestocity 
dt = .001

# Making u,v,w,s as zeroes intially
u = np.zeros((ny, nx)) # u component velocity
v = np.zeros((ny, nx)) # v component velocity
w = np.zeros((ny, nx)) # Vorticity
s = np.zeros((ny, nx)) # Stream Function
dwdt=np.zeros((ny, nx))
t=0
for i in range(10000):    # the code will run till it reanches 1000000 steps
	
	
	dwdt[1:-1,1:-1]=(((w[1:-1,2:]+w[1:-1,:-2]+w[2:,1:-1]+w[:-2,1:-1]-4*w[1:-1,1:-1]))*(nu/(dx*dx))
							-((s[:-2,1:-1]-s[2:,1:-1])*(w[1:-1,:-2]-w[1:-1,2:]))/(4*(dx*dx))
								+((s[1:-1,:-2]-s[1:-1,2:])*(-w[2:,1:-1]+w[0:-2,1:-1]))/(4*(dx*dx)))
    

	w[1:-1,1:-1]=w[1:-1,1:-1]+dwdt[1:-1,1:-1]*dt
	s[1:-1,1:-1]=(s[:-2,1:-1]+s[2:,1:-1]+s[1:-1,2:]+s[1:-1,:-2]+w[1:-1,1:-1]*dx*dx)*0.25

	# bottom wall
	w[0,:]=(s[1,:])*(-2/(dx*dx))
	# top wall
	w[-1,:]=(s[-2,:])*(-2/(dx*dx))-2/dx
    # right wall
	w[1:-1,-1]=(s[1:-1,-2])*(-2/(dx*dx))
	# left wall
	w[1:-1,0]=(s[1:-1,1])*(-2/(dx*dx))

	t=t+1

v=np.diff(s)/dx
u=(np.diff(s.T)/(dx)).T


# Plotting Stream Function
plt.figure(figsize=(5,5))
sp=plt.contourf(s,cmap=plt.cm.jet,levels=30)
plt.contour(s,alpha=0.1,colors='k')
plt.colorbar(sp)
plt.title('Stream Function')
plt.show()

# Plotting vorticity
plt.figure(figsize=(5,5))
sp=plt.contourf(w,cmap=plt.cm.jet,levels=10)
plt.contour(w,alpha=0.3,colors='k')
plt.colorbar(sp)
plt.title('Vorticity')
plt.show()

# Plotting u component velocity
plt.figure(figsize=(5,5))
sp=plt.contourf(u,alpha=0.9,cmap=plt.cm.jet,levels=15)
plt.contour(u,alpha=0.3,colors='k')
plt.colorbar(sp)
plt.title('u velocity')
plt.show()

# Plotting v component velocity
plt.figure(figsize=(5,5))
sp=plt.contourf(v,alpha=0.9,cmap=plt.cm.jet,levels=15)
plt.contour(v,alpha=0.3,colors='k')
plt.colorbar(sp)
plt.title('v velocity')
plt.show()
