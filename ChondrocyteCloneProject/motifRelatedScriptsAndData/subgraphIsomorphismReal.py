
import os 
from os import listdir
from os.path import isfile, join
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt 
from networkx.algorithms.community import k_clique_communities
import scipy.stats 
from networkx.algorithms import isomorphism




def count_cliques(G,Name,flag):

	if flag==1:
		maindir='drawMotif_'+Name+'/'
		answer=os.path.isdir(maindir)
		if answer==True:
			pass
		else:
			os.mkdir(maindir)
		f_degname=open(maindir+'cluster_degree_id.dat','w')


	#cliques=list(nx.enumerate_all_cliques(G))
	cliques = list(k_clique_communities(G, 2))
	myedges=list(G.edges)

	#print([list(x) for x in cliques])

	graphletFreq={}
	size={}
	count=1;
	isomorphic_motif={}
	repeating={}
	
	for s in cliques:
		s=list(s)
		if len(s)>0:
			local_edges=motif_edges(s, myedges)
			G=nx.Graph()
			G.add_edges_from(local_edges)

			motif_degree_seq = [d for v, d in G.degree()]
			a=sorted(motif_degree_seq)
			degname=''
			for i in a:		
				degname+=str(i)+'-'

			if degname in graphletFreq:
				graphletFreq[degname]+=1
				isomorphic_motif[degname].append([count,s])
			else:
				graphletFreq[degname]=1	
				isomorphic_motif[degname]=[[count,s]]


			if flag==1:
				fig, ax = plt.subplots( nrows=1, ncols=1 )
				nx.draw(G, node_size=200, with_labels = True,   cmap = plt.get_cmap('jet'))
				fig.savefig(maindir+'/cluster'+'_'+str(count)+'.png')
				if degname in repeating:
					repeating[degname].append(count)
				else:
					repeating[degname]=[count]			

				count=count+1
				plt.close(fig)


	if flag==1:
		for degname in repeating:
			f_degname.write(degname+'\t'+str(repeating[degname])+'\n')
		
	return [graphletFreq,isomorphic_motif,repeating]



many_subgraph=[]
#[real_graphlet,real_iso,real_repeating]=count_cliques(G,'Real',1)

n=10
count=1
for j in range(3,n):
	for i in range(j):
		a=j-i
		b=i
		if (a>1) & (a<5):
			fig, ax = plt.subplots( nrows=1, ncols=1,figsize=(2,2) )
			G = nx.lollipop_graph(a,b)
			nx.draw(G, node_size=50, with_labels=False, width=2)
			plt.savefig('Harary/'+'Graphlet_'+str(count)+'_lollipop_graph('+str(a)+','+str(b)+').png',transparent=True)
			many_subgraph.append(G)
			count=count+1
			plt.close(fig)


print('Hi Ankit')


edge11=[[0,1],[1,2],[2,3],[3,0],[0,2]]
G=nx.Graph()
G.add_edges_from(edge11)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine1('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')

edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4]]
G=nx.Graph()
G.add_edges_from(edge11)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine1('+str(4)+ ',1).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')

edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge11)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine1('+str(4)+ ',2).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')

edge11=[[0,1],[1,2],[2,3],[3,0],[0,2],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge11)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine1('+str(4)+ ',3).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')


edge21=[[0,1],[1,2],[2,3],[3,0],[1,3]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')


edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',1).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')


edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',2).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')


edge21=[[0,1],[1,2],[2,3],[3,0],[1,3],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',3).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')

print('I am here')

for b in range(n-1):
	#H=nx.generators.harary_graph.hnm_harary_graph(a,b)
	if b>3:
		H=nx.generators.cycle_graph(b)
		plt.figure(figsize=(2,2))
		nx.draw(H, node_size=200, with_labels=False, width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_cycle('+str(b)+').png',transparent=True)
		many_subgraph.append(H)
		count=count+1
		plt.close()

		node=list(H.nodes())
		edge=list(H.edges())
		t=max(node)
		edge.append([t,t+1])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False, width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_cycle_length('+str(b)+ ',1).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()

		edge.append([t+1,t+2])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False, width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_cycle_length('+str(b)+ ',2).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()

		edge.append([t+2,t+3])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False, width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_cycle_length('+str(b)+ ',3).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()



for b in range(n-1):
	#H=nx.generators.harary_graph.hnm_harary_graph(a,b)
	if (b>4)&(b<7):
		H=nx.generators.wheel_graph(b)
		plt.figure(figsize=(2,2))
		nx.draw(H, node_size=200, with_labels=False,width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_wheel('+str(b)+').png',transparent=True)
		many_subgraph.append(H)
		count=count+1
		plt.close()
		
		node=list(H.nodes())
		edge=list(H.edges())
		t=max(node)
		edge.append([t,t+1])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False,width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_wheel_length('+str(b)+ ',1).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()

		edge.append([t+1,t+2])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False,width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_wheel_length('+str(b)+ ',2).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()

		edge.append([t+2,t+3])
		G=nx.Graph()
		G.add_edges_from(edge)
		plt.figure(figsize=(2,2))
		nx.draw(G, node_size=200, with_labels=False, width=5)
		plt.savefig('Harary/'+'Graphlet_'+str(count)+'_wheel_length('+str(b)+ ',3).png',transparent=True)
		many_subgraph.append(G)
		count=count+1
		plt.close()

edge21=[[1,3],[2,3],[3,4]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[4,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[5,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[2,6]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')




edge21=[[1,3],[2,3],[3,4],[4,5],[4,6],[6,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[5,6],[6,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[5,6],[2,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



edge21=[[1,3],[2,3],[3,4],[4,5],[2,6],[1,7]]
G=nx.Graph()
G.add_edges_from(edge21)
plt.figure(figsize=(2,2))
nx.draw(G, node_size=200, with_labels=False, width=5)
plt.savefig('Harary/'+'Graphlet_'+str(count)+'_mine2('+str(4)+ ',0).png',transparent=True)
many_subgraph.append(G)
count=count+1
plt.close('all')



print('All motifs',len(many_subgraph))



                    
                        
def search_edges(edges,vertices):

	unique_edges=[]
	d={}
	for i in range(len(edges)):
		#print(edges[i])
		fake=str(edges[i,0])+'-'+str(edges[i,1])
		if fake not in d:
			unique_edges.append(edges[i])
	
	edges=unique_edges
	
	#print(edges[i][0])
	
	count=1
	myedges=[]
	d={}
	for j in range(len(edges)):
		flag1=0;
		flag2=0;
		for i in range(len(vertices)):
			#print(edges[j][0],vertices[i])
			if edges[j][0]==vertices[i]:
				flag1=1
			if edges[j][1]==vertices[i]:
				flag2=1

		if (flag1+flag2)==2:
			name=str(edges[j][0])+'-'+str(edges[j][1])
			if name not in d:
				myedges.append(edges[j]);
                	
				
	return myedges
 	
    





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

		#print(name)
		if name not in d:
			d[name]=1
			#f0.write(str(b)+'\n')
			frequency.append(b)

	return [frequency,d] 



real_network_edges='Tibia/MakeListColumnarStructurePrediction/allEdges_medium.dat'
real_edges=np.loadtxt(real_network_edges,delimiter="\t", usecols=(0,1),dtype=float,skiprows=0).astype(int)
G=nx.Graph()
G.add_edges_from(real_edges)


save_nodes_of_motifs=[]
for j in range(len(many_subgraph)):
	f1=open('Real_node_name/name'+str(j+1)+'.dat','w')
	save_nodes_of_motifs.append(f1)
	


f=open('Tibia/MakeListColumnarStructurePrediction/ClusterMediumCutoff.dat')


cluster=[]
for line in f:
	l=line.split()
	t=[]
	if len(l)>2:
		for j in l:
			t.append(int(j))
	cluster.append(t)
		
		


	



an=[]
ae=[]
we=[]
print('node\tedge\weight\n')
for j in range(len(many_subgraph)):
	H=many_subgraph[j]
	n=len(H.nodes())
	e=len(H.edges())
	we.append(e/((n-1)/2.0))
	an.append(n)
	ae.append(e)
	
	print(j+1,'\t',an[j],'\t',ae[j],'\t',we[j])



f0=open('RealSubgraphIsomorphism.dat','w')
for i in range(len(cluster)):
	ed=search_edges(real_edges,cluster[i])

	Greal=nx.Graph()
	Greal.add_edges_from(ed)

	#f0.write(str(i+1)+'\t')
	for j in range(len(many_subgraph)):
		H=many_subgraph[j]
		[frequency,nodename]=searchSubgraph(Greal,H)
		f0.write(str(len(frequency))+'\t')
		for key in nodename:
			k=key.split('-')
			for t in k:
				save_nodes_of_motifs[j].write(t+'\t')
			save_nodes_of_motifs[j].write('-1'+'\t')
		save_nodes_of_motifs[j].write('\n')
	f0.write('\n')
	

	


