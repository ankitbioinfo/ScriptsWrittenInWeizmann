
import os 
from os import listdir
from os.path import isfile, join
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt
import random 

import time


def plotSpectra(graph,name):
	L = nx.normalized_laplacian_matrix(graph)
	#L = nx.laplacian_matrix(graph)

	size=L.shape
	#print(size)
	#print(L[0][0])

	f=open('mat/Laplacin_'+name+'.dat','w')
	for i in range(size[0]):
		for j in range(size[1]):
			f.write(str(L[i,j])+' ')
		f.write('\n')

	f.close()
	'''
	f=open('mat/Laplacin_dotA_'+name+'.dat','w')
	for i in range(size[0]):
		for j in range(size[1]):
			f.write(str(L.A[i,j])+' ')
		f.write('\n')
	f.close()
	'''

	e = np.linalg.eigvals(L.A)


	trace=[]
	for i in range(len(L.A)):
		trace.append(L.A[i,i])

	print('diagonal elements in laplacian', min(trace), max(trace))
	print("Largest eigenvalue:", max(e))
	print("Smallest eigenvalue:", min(e))
	plt.hist(e, bins=100)  # histogram with 100 bins
	plt.xlim(0, 2)  # eigenvalues between 0 and 2
	plt.savefig('mat/'+name+'.png')
	plt.close()



def  different_matrix(centroids,distance_cutoff,name):
	print('\n\n\n' +name +'\n\n\n')
	weightedEdges=[]
	edges=[]
	weightedInvEdges=[]
	print('mean', np.mean(centroids,axis=0))
	print('std', np.std(centroids,axis=0))
	for i in range(len(centroids)):
		for j in range(i+1,len(centroids)):
			weig=np.linalg.norm(centroids[i] - centroids[j])
			weightedEdges.append([i+1,j+1,weig])
			weightedInvEdges.append([i+1,j+1,1.0/weig])
			if weig<distance_cutoff:
				edges.append([i+1,j+1])
				edges.append([j+1,i+1])

	
	print('edges',len(edges))
	if len(edges)>10:
		Go=nx.Graph()			
		Go.add_edges_from(edges)
		plotSpectra(Go, name+'_unweighted')

	Gw=nx.Graph()
	Gw.add_weighted_edges_from(weightedEdges)

	Gwi=nx.Graph()
	Gwi.add_weighted_edges_from(weightedInvEdges)


	#plotSpectra(Gw, name+'_direct_dist')
	plotSpectra(Gwi,name+'_inverse_dist')
		

centroids=np.loadtxt('centroid.dat',delimiter="\t", usecols=(0,1,2),dtype=float,skiprows=0)
a=np.min(centroids,axis=0)
centroids=centroids-min(a)
b=np.max(centroids,axis=0)
factor=max(b)
centroids=centroids/factor
print(np.min(centroids,axis=0))
print(factor,np.max(centroids,axis=0))


f1=open('scaledCentroid.dat','w')

for i in range(len(centroids)):
	for j in range(3):
		f1.write(str(centroids[i][j])+'\t')
	f1.write('\n')



distance_cutoff=17.5/factor 
n=len(centroids)
print('nodes ', n)

realMean=np.mean(centroids,axis=0)
realStd=np.std(centroids,axis=0)

different_matrix(centroids,distance_cutoff,'real')


G=nx.random_geometric_graph(n,1,dim=3,pos=None,seed=123456)
pos1 = nx.get_node_attributes(G, "pos")
pos=[]
for i in range(n):
	pos.append(pos1[i])



#plotSpectra(G,'random_normal_pos_none')
#print(centroids[0] ,len(centroids),len(pos))
#print(np.array(pos[0])-np.array(pos[1]))
different_matrix(np.array(pos),distance_cutoff,'random_normal_pos_none')

for i in range(0):
	#stddev=i
	pos1 = {i: (random.gauss(realMean[0], realStd[0]), random.gauss(realMean[1], realStd[1]),random.gauss(realMean[2], realStd[2])) for i in range(n)}
	pos=[]
	for j in range(n):
		pos.append(pos1[j])

	#G=nx.random_geometric_graph(n,1,dim=3,pos=pos)
	#plotSpectra(G,'random_pos_stddev'+str(stddev))
	different_matrix(np.array(pos),distance_cutoff,'random_pos_realization'+str(i))

