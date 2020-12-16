

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
			
	
def count_cliques(G,name):
	#cliques=list(nx.enumerate_all_cliques(G))
	cliques = list(k_clique_communities(G, 2))
	myedges=list(G.edges)

	#print([list(x) for x in cliques])

	graphletFreq={}
	size={}
	count=1;
	for s in cliques:
		s=list(s)
		if len(s)>2:
			#print(s)
			local_edges=search_edges(s, myedges)
			G=nx.Graph()
			G.add_edges_from(local_edges)

			degree_seq = [d for v, d in G.degree()]
			a=sorted(degree_seq)
			degname=''
			for i in a:		
				degname+=str(i)+'-'

			if degname in graphletFreq:
				graphletFreq[degname]+=1
			else:
				graphletFreq[degname]=1	

			'''
			fig, ax = plt.subplots( nrows=1, ncols=1 )
			nx.draw(G, with_labels = True,   cmap = plt.get_cmap('jet'))
			fig.savefig('Random/clique'+str(name)+'_'+str(count)+'.png')
			count=count+1
			plt.close(fig)
			
			if len(s) not in size:
				size[len(s)]=1
			else:
				size[len(s)]+=1
			'''

	f=open('communityInRandomNetwork/'+str(name)+'.txt','w')
	for key in graphletFreq:
		f.write(key +'\t'+ str(graphletFreq[key])+'\n')



	return size



def generate_configuration_model(degree_seq):
	flag=1
	while(flag==1):
		G=nx.Graph()
		G=nx.configuration_model(degree_seq)
		G=nx.Graph(G)
		G.remove_edges_from(nx.selfloop_edges(G))
		actual_degrees = [d for v, d in G.degree()]
		giant = nx.number_connected_components(G)
		if (actual_degrees==degree_seq)&(giant==1):
			flag=0
			edges=list(G.edges())
			nodes=list(G.nodes())
		'''
		fig, ax = plt.subplots( nrows=1, ncols=1 )
		mylabel=['a','b','c','d','e','f','g','h']
		nx.draw(G, with_labels = True, node_size = 0.5,  cmap = plt.get_cmap('jet'))
		fig.savefig('Fig/randomGraph'+str(name)+'.png')
		plt.close(fig)
		'''

	return [nodes,edges]



def read_files(data):
	'''
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
	'''
	edges=[]
	for i in range(len(data)):
		edges.append([int(data[i][0]),int(data[i][1])])



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



def Generate_random_network_of_similar_degree_and_connected_components(Go,name):
	original_network_degree_seq = [d for v, d in Go.degree()]
	connectedClusters=sorted(nx.connected_components(Go), key=len, reverse=True)
	start=0
	all_edges=[]
	for i in range(len(connectedClusters)):
		individualComponents=list(connectedClusters[i])
		individualComponent_degrees=Go.degree(individualComponents)
		deg = [d for v, d in individualComponent_degrees]
		[nodes,edges]=generate_configuration_model(deg)

		rename_edges=rename_nodes_and_edges(nodes,edges,start)
		all_edges=all_edges+rename_edges
		start=len(nodes)+start
		#print(nodes, rename_edges)



	Gr=nx.Graph()
	Gr.add_edges_from(all_edges)
	giant = nx.number_connected_components(Gr)
	random_network_degree_seq = [d for v, d in Gr.degree()]
	if (sorted(original_network_degree_seq)==sorted(random_network_degree_seq)):
		if (len(connectedClusters)==giant):
			f=open(name,'w')
			for i in range(len(all_edges)):
				f.write(str(all_edges[i][0])+'\t'+str(all_edges[i][1])+'\n')	
			f.close()		






maindir='RandomNetowrk_/'
answer=os.path.isdir(maindir)
if answer==True:
	pass
else:
	os.mkdir(maindir)


name1='MakeListColumnarStructurePrediction/DT_E18.5/network_edges80.dat'
data1 = np.loadtxt(name1, delimiter="\t", skiprows=0)
Greal=read_files(data1)



for i in range(10000):
	if (i%25==0):
		print(i)
	Generate_random_network_of_similar_degree_and_connected_components(Greal,maindir+'RandEdges_s1_'+str(i+1)+'.txt')



#print(original_network_degree_seq)

#print(count_cliques(G))

#print(G.edges())


