



fw=open('batch_segementation_parameters11.txt','w')

fr=open('nuclei2.txt','r')


a=[]
for line in fr:
	a.append(line[0:-1])

b=sorted(a)

for i in range(len(b)):
	fw.write(b[i]+'\t'+'12\n')
		
