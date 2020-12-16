
import os 
from os import listdir
from os.path import isfile, join
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
from networkx.algorithms.community import k_clique_communities
import scipy.stats 
from networkx.algorithms import isomorphism



many_subgraph=[]

n=10
count=1
for j in range(3,n):
	for i in range(j):
		a=j-i
		b=i
		if (a>1) & (a<5):
			fig, ax = plt.subplots( nrows=1, ncols=1 )
			G = nx.lollipop_graph(a,b)
		
			many_subgraph.append(G)
			count=count+1
			plt.close(fig)


edge11=[[0,1],[1,2],[2,3],[3,0],[0,2]]
G=nx.Graph()
G.add_edges_from(edge11)

many_subgraph.append(G)
count=count+1


edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4]]
G=nx.Graph()
G.add_edges_from(edge11)

many_subgraph.append(G)
count=count+1


edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge11)

many_subgraph.append(G)
count=count+1


edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge11)

many_subgraph.append(G)
count=count+1



edge21=[[0,1],[1,2],[2,3],[3,0],[1,3]]
G=nx.Graph()
G.add_edges_from(edge21)

many_subgraph.append(G)
count=count+1



edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4]]
G=nx.Graph()
G.add_edges_from(edge21)

many_subgraph.append(G)
count=count+1



edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge21)

many_subgraph.append(G)
count=count+1



edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge21)

many_subgraph.append(G)
count=count+1




for b in range(n-1):
	#H=nx.generators.harary_graph.hnm_harary_graph(a,b)
	if b>3:
		H=nx.generators.cycle_graph(b)
	
		many_subgraph.append(H)
		count=count+1
		

		node=list(H.nodes())
		edge=list(H.edges())
		t=max(node)
		edge.append([t,t+1])
		G=nx.Graph()
		G.add_edges_from(edge)
	
		many_subgraph.append(G)
		count=count+1
		

		edge.append([t+1,t+2])
		G=nx.Graph()
		G.add_edges_from(edge)
	
		many_subgraph.append(G)
		count=count+1
		

		edge.append([t+2,t+3])
		G=nx.Graph()
		G.add_edges_from(edge)
	
		many_subgraph.append(G)
		count=count+1
		



for b in range(n-1):
	#H=nx.generators.harary_graph.hnm_harary_graph(a,b)
	if (b>4)&(b<7):
		H=nx.generators.wheel_graph(b)
		
		many_subgraph.append(H)
		count=count+1
		
		
		node=list(H.nodes())
		edge=list(H.edges())
		t=max(node)
		edge.append([t,t+1])
		G=nx.Graph()
		G.add_edges_from(edge)
		
		many_subgraph.append(G)
		count=count+1
		

		edge.append([t+1,t+2])
		G=nx.Graph()
		G.add_edges_from(edge)
	
		many_subgraph.append(G)
		count=count+1
		

		edge.append([t+2,t+3])
		G=nx.Graph()
		G.add_edges_from(edge)
	
		many_subgraph.append(G)
		count=count+1
		


edge21=[[1,3],[2,3],[3,4]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[4,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[2,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[4,6],[6,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[5,6],[6,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[5,6],[2,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[2,6],[1,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=100, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png')
many_subgraph.append(G)
count=count+1
plt.close('all')



print('All motifs',len(many_subgraph))




def searchSubgraph(G,G2):
	gm=nx.algorithms.isomorphism.GraphMatcher(G,G2)
	#answer=gm.subgraph_is_isomorphic()
	#answer=gm.mapping.items()
	answer=gm.subgraph_isomorphisms_iter()

	frequency=[]
	d={}
	for i in answer:
		a=list(i.keys())
		b=sorted(a)
		name=''
		for k in b:
			name+=str(k)+'-'

		if name not in d:
			d[name]=1
			#f0.write(str(b)+'\n')
			frequency.append(b)

	return frequency 



d = 'NullModel/Realization_medium/'
listdir=[os.path.join(d, o) for o in os.listdir(d) 
                    if os.path.isdir(os.path.join(d,o))]
m=len(listdir)


fd=open('Tibia/MakeListColumnarStructurePrediction/ClusterMediumCutoff.dat')
count=0
for line in fd:
	l=line.split()
	#print(len(l),l)
	if len(l)>2:
		count=count+1

n=count

print('# of directory and files ',m, n )

f0=open('Random_Mean_SubgraphIsomorphism.dat','w')
f1=open('Random_StdDev_SubgraphIsomorphism.dat','w')





for i in range(n):
	f0.write(str(i+1)+'\t')
	f1.write(str(i+1)+'\t')

	
	frequency=np.zeros((len(many_subgraph),m),dtype=np.float)


	
	for k in range(m):
		if k%100==0:
			print(i,k)
		#print(listdir[k])
		name=listdir[k]+'/Edges_'+str(i+1)+'.dat'
		ed=np.loadtxt(name,delimiter="\t", usecols=(0,1),dtype=float,skiprows=0).astype(int)
		Greal=nx.Graph()
		Greal.add_edges_from(ed)
		for j in range(len(many_subgraph)):
			H=many_subgraph[j]
			frequency[j][k]=len(searchSubgraph(Greal,H))


	for j in range(len(many_subgraph)):
		mu=np.mean(frequency[j])
		std=np.std(frequency[j])

		f0.write(str(mu)+'\t')
		f1.write(str(std)+'\t')
		
	f0.write('\n')	
	f1.write('\n')
	


