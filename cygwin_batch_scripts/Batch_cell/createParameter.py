



fw=open('batch_segementation_parameters11.txt','w')

fr=open('allname.txt','r')

fsara=open('Sarahabatchcell_segparamters.txt')


a=[]
for line in fr:
	l=line.split()
	if l[0].find('.tif')!=-1:
		a.append(l[0])
	if l[1].find('.tif')!=-1:
		a.append(l[1])

b=sorted(a)

for i in range(len(b)):
	fw.write(b[i]+'\t'+'12\n')
		
