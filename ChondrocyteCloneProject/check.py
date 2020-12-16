


a=open('Tibia/MakeListColumnarStructurePrediction/ClusterMediumCutoff.dat');

d=[]
for line in a:
	l=line.split()
	d.append(len(l))


for p in range(1,11):
	count=0
	for j in range(764):
		b=open('NullModel/Realization_medium/Iteration_'+str(p)+'/Edges_'+str(j+1)+'.dat')
		cont=b.readlines()

		dic={}
		for i in range(len(cont)):
			l=cont[i].split()
			for k in l:
				dic[int(k)]=1
			
		if len(dic)==d[j]:
			count+=1
		else:
			print(p,j+1,cont,len(dic),d[j])
	print(p,count,len(d))

