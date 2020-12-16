

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
from networkx.algorithms.community import k_clique_communities
import os 

def search_edges(vertices, edges):
	myedges=[]
	for j in range(len(edges)):
		flag1=0
		flag2=0
		for i in range(len(vertices)):
			if (edges[j][0]==vertices[i]):
				flag1=1
			if (edges[j][1]==vertices[i]):
				flag2=1
		if (flag1+flag2)==2:
			myedges.append(edges[j])
	return myedges
			
	

def generate_configuration_model(degree_seq,original_no_edges,name):
	flag=1
	iteration=0
	criteria=1
	while(flag==1):
		iteration=iteration+1
		#G=nx.configuration_model(degree_seq)
		G = nx.random_clustered_graph(degree_seq)
		G=nx.Graph(G)
		G.remove_edges_from(nx.selfloop_edges(G)) 
		actual_degrees = [d for v, d in G.degree()]
		giant = nx.number_connected_components(G)
		#if (sorted(actual_degrees)==sorted(degree_seq)):
		#if (giant==1)&(len(G.edges())==original_no_edges):
		if iteration>100:
			criteria=2
		if iteration>200:
			criteria=3
		if iteration>300:
			criteria=4
		if iteration>400:
			criteria=5
			print(name,criteria,giant, original_no_edges, len(G.edges()))
		if (giant==criteria):
			flag=0
			edges=list(G.edges())
			nodes=list(G.nodes())
		

	return [nodes,edges]



def read_files(data,cutoff,pos):
	edges=[]
	for i in range(len(data)):
		for j in range(i+1,len(data)):
			dist = np.linalg.norm( data[i] - data[j] )
			if dist<cutoff:
				edges.append([i+1,j+1])
	
	f=open('Real_network_edges'+str(pos)+'.dat','w')
	for i in range(len(edges)):
		f.write(str(edges[i][0])+'\t'+str(edges[i][1])+'\n')
	f.close()


	G=nx.Graph()
	G.add_edges_from(edges)
	giant = nx.number_connected_components(G)
	print(len(G.nodes()),len(edges),giant)
	return G 


def rename_nodes_and_edges(nodes,edges,start):
	d={}


	for i in range(len(nodes)):
		d[nodes[i]]=nodes[i]+start 	

	newedges=[]
	for i in range(len(edges)):
		t=[]
		for j in range(len(edges[i])):
			t.append(d[  edges[i][j]  ])
		newedges.append(t)
	return newedges 	




def find_degree_and_triangle_sequence(Go,nodes):
	cliques=list(nx.enumerate_all_cliques(Go))
	#print([list(x) for x in cliques])

	#nodes=list(Go.nodes(mynodes))
	#print('nodes', len(nodes), 'edges', len(list(Go.edges())))
	triangle={}
	for i in range(len(nodes)):
		triangle[nodes[i]]=0
		
	for x in cliques:
		if len(x)==3:
			#print(x)
			for j in x:
				if j in triangle:
					triangle[j]+=1
				else:
					triangle[j]=1
	#print(3*count,sum(triangle.values()))
	#print(triangle)
	joint_degree_and_triangle=[]
	for i in range(len(nodes)):
		if triangle[nodes[i]]>0:
			neigh=[n for n in Go.neighbors(nodes[i])]
			#print(nodes[i],neigh)
			count=0
			for j in range(len(neigh)):
				if triangle[neigh[j]]!=0:
					count=count+1
			independentEdgeDegree=Go.degree(nodes[i])-count
		else:
			independentEdgeDegree=Go.degree(nodes[i])
	
		#print(nodes[i],Go.degree(nodes[i]),triangle[nodes[i]],independentEdgeDegree)
		joint_degree_and_triangle.append((independentEdgeDegree,triangle[nodes[i]]))

	return joint_degree_and_triangle





def Generate_random_network_of_similar_degree_and_connected_components(Go,name):
	#deg = [d for v, d in individualComponent_degrees]

	connectedClusters=sorted(nx.connected_components(Go), key=len, reverse=True)
	start=0
	all_edges=[]
	for i in range(len(connectedClusters)):
		individualComponents=list(connectedClusters[i])
		joint_degree_and_triangle_sequence=find_degree_and_triangle_sequence(Go,individualComponents)
		componentEdge=search_edges(individualComponents, list(Go.edges()))
		#print('comp',i)
		[nodes,edges]=generate_configuration_model(joint_degree_and_triangle_sequence, len(componentEdge),'comp'+str(i))

		rename_edges=rename_nodes_and_edges(nodes,edges,start)
		all_edges=all_edges+rename_edges
		start=len(nodes)+start
		#print(nodes, rename_edges)


	Gr=nx.Graph()
	Gr.add_edges_from(all_edges)
	giant = nx.number_connected_components(Gr)
	print('Component', len(connectedClusters),giant)
	#original_network_degree_seq = [d for v, d in Go.degree()]
	#random_network_degree_seq = [d for v, d in Gr.degree()]
	#if (sorted(original_network_degree_seq)==sorted(random_network_degree_seq)):
	#if (len(connectedClusters)==giant):
	if True:
			f=open(name,'w')
			for i in range(len(all_edges)):
				f.write(str(all_edges[i][0])+'\t'+str(all_edges[i][1])+'\n')	
			f.close()		






maindir='RandomClusteringNetwork/'
answer=os.path.isdir(maindir)
if answer==True:
	pass
else:
	os.mkdir(maindir)


#name1='MakeListColumnarStructurePredictionPos'+str(pos)+'/centroid.dat'
#data1 = np.loadtxt(name1, delimiter="\t", skiprows=0)
#Greal=read_files(data1,cutoff,pos)

real_edges=np.loadtxt('network_edges80.dat',delimiter="\t",usecols=(0,1), dtype=float,skiprows=0).astype(int)
#real_edges=np.loadtxt('../separateComponents/cluster1.dat',delimiter="\t", dtype=float,skiprows=0).astype(int)

Greal=nx.Graph()
Greal.add_edges_from(real_edges)



for i in range(4000):
	if (i%25==0):
		print(i)
	Generate_random_network_of_similar_degree_and_connected_components(Greal,maindir+'Rand5r'+str(i+1)+'.txt')








