import freud
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
from mpl_toolkits.mplot3d import Axes3D
import rowan  # for quaternion math, see rowan.readthedocs.io for more information.
from scipy.spatial.transform import Rotation as R

rot=[[0.4822 ,   0.5070 ,   0.7145],[-0.4086 ,   0.8516 ,  -0.3285],[ 0.7750 ,   0.1336 ,  -0.6177]]

#rot=[[1,0,0],[0,1,0],[0,0,1]]

r=R.from_dcm(rot)



print(r.as_quat())

data_path='centroid.dat' 

centroids=np.loadtxt(data_path,delimiter="\t", usecols=(0,1,2),dtype=float,skiprows=0)
(a,b,c)=(np.min(centroids[:,0]), np.min(centroids[:,1]), np.min(centroids[:,2]))
(x,y,z)=(np.max(centroids[:,0]), np.max(centroids[:,1]), np.max(centroids[:,2]))

nx=int(np.ceil(x-a))
ny=int(np.ceil(y-b))
nz=int(np.ceil(z-c))

box = np.array([nx, ny, nz], dtype=np.float32)
box = freud.box.Box(*box)

points=centroids

orientations = np.loadtxt('quaternionCells.dat',delimiter="\t", usecols=(0,1,2,3),dtype=np.longdouble,skiprows=0)
orientationsCluster = np.loadtxt('quaternionCluster.dat',delimiter="\t", dtype=np.longdouble,skiprows=0)
orientationsLatent = np.loadtxt('quaternionLatent.dat',delimiter="\t", dtype=np.longdouble,skiprows=0)

print(np.shape(orientationsCluster))

for i in range(25):
	q=R.from_quat(orientations[i])
	#print(i+1,q.as_dcm())
#print(q.as_quat())



def calculateNematicAndPlot(orientations,points,clusterOrient,latent, name):
	# To show orientations, we use arrows rotated by the quaternions.
	arrowheads = rowan.rotate(orientations, [1, 0, 0])

	(a,b,c)=(np.min(points[:,0]), np.min(points[:,1]), np.min(points[:,2]))
	(x,y,z)=(np.max(points[:,0]), np.max(points[:,1]), np.max(points[:,2]))

	boxSize=(x-a)*(y-b)*(z-c)
	ls=boxSize**(1/3.0)
	
	fig = plt.figure(figsize=(8,4))

	ax = fig.add_subplot(121, projection='3d')

	ax.plot(centroids[:,0],centroids[:,1],centroids[:,2],'k.',ms=0.5)
	ax.plot(points[:,0],points[:,1],points[:,2],'b.')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.view_init(elev=105., azim=165)

	ax = fig.add_subplot(122, projection='3d')
	factor=1
	ax.quiver3D(points[:, 0], points[:, 1], points[:, 2],
		    factor*arrowheads[:, 0], factor*arrowheads[:, 1], factor*arrowheads[:, 2],length=ls/5.0, arrow_length_ratio=0.4, color='b', lw=1)
	

	nop = freud.order.Nematic([1, 0, 0])
	nop.compute(orientations)
	#print("The value of the order parameter is {}.".format(nop.order))

	center=np.mean(points,axis=0)

	'''
	eig = np.linalg.eig(nop.nematic_tensor)
	rotation=R.from_dcm(eig[1])
	arrowheads = rowan.rotate(rotation.as_quat(), [1, 0, 0])
	
	#ax.quiver3D(center[0], center[1], center[2],arrowheads[0], arrowheads[1], arrowheads[2],length=20, color='k', lw=2)
	print(name,  np.max(eig[0]),eig[1], rotation.as_dcm())
	'''

	'''
	arrowheads = rowan.rotate(clusterOrient[[0,1,2,3]], [1, 0, 0])
	factor=ls/50.0
	latent=[12,6,3]
	ax.quiver3D(center[0], center[1], center[2],arrowheads[0], arrowheads[1], arrowheads[2],length=latent[0]*factor, color='b', lw=1)
	ax.quiver3D(center[0], center[1], center[2],-arrowheads[0], -arrowheads[1], -arrowheads[2],length=latent[0]*factor, color='b', lw=1)

	arrowheads = rowan.rotate(clusterOrient[[4,5,6,7]], [1, 0, 0])
	ax.quiver3D(center[0], center[1], center[2],arrowheads[0], arrowheads[1], arrowheads[2],length=latent[1]*factor, color='m', lw=1)
	ax.quiver3D(center[0], center[1], center[2],-arrowheads[0], -arrowheads[1], -arrowheads[2],length=latent[1], color='m', lw=1)

	arrowheads = rowan.rotate(clusterOrient[[8,9,10,11]], [1, 0, 0])
	ax.quiver3D(center[0], center[1], center[2],arrowheads[0], arrowheads[1], arrowheads[2],length=latent[2]*factor, color='k', lw=1)
	ax.quiver3D(center[0], center[1], center[2],-arrowheads[0], -arrowheads[1], -arrowheads[2],length=latent[2]*factor, color='k', lw=1)
	'''

	t=' ['
	t1=[]
	axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
	for ax1 in axes:
	    nop = freud.order.Nematic(ax1)
	    nop.compute(orientations)
	    #print("For axis {}, the value of the order parameter is {:0.3f}.".format(ax1, nop.order))
	    t=t+' %0.2f,'%nop.order	
	    t1.append(nop.order)	 	    	
	t=t[0:-1]+']'	
	fig.suptitle("Cluster "+str(name)+''+'', fontsize=16);
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	ax.view_init(elev=105., azim=165)
	plt.savefig('Nematic/cluster'+str(name)+'.png')
	#plt.show()
	plt.close()
	return t1


f=open('beforeFusion_zone80.dat')
f1=open('OrderParameter7axis.dat','w')

count=0
for line in f:
	l=line.split()
	cellid=[]
	for i in l:
		cellid.append(int(i)-1)


	if len(cellid)>2:
		orient=orientations[cellid,:]
		clusterpoints=points[cellid,:]

		
		op=calculateNematicAndPlot(orient,clusterpoints,orientationsCluster[count],orientationsLatent[count],count+1)
		
		for i in op:
			f1.write(str(i)+'\t')
		f1.write('\n')
		count=count+1


