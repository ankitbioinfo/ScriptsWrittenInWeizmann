


import os 
from os import listdir
from os.path import isfile, join

mypath='inputXML22/'


onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

#for i in range(len(onlyfiles)):
	#print(onlyfiles[i])
count=0
for i in range(len(onlyfiles)):

	name=onlyfiles[i]
	t=name.split('_')
	p1=float(t[8])/100.0
	p2=float(t[9])/100.0
	t3=t[10].split('.')
	p3=int(t3[0])

	f=open(mypath+name)
	cont=f.readlines()

	
	l1=cont[192].split()
	l2=cont[40].split()
	l3=cont[106].split()

	#print(l3)
	a=l1[2].split('=')
	b=l2[2].split('=')
	c=l3[2].split('=')

	x=a[1].split('"')
	y=b[1].split('"')
	z=c[1].split('"')

	



	if (p1==float(x[1])) & (p2==float(y[1])) & (p3==int(z[1])):
		count=count+1
	else: 
		print([p1,l1[2]],'\t',[p2,l2[2]],'\t',[p3,l3[2]])

print(count,len(onlyfiles))





