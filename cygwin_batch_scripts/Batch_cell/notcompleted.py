

f=open('completed')

completed=[]
for line in f:
	completed.append(line[0:-1])


f1=open('actualCompleted')

actualcomplete=[]
for line in f1:
	l=line.split('\\')
	if len(l)>2:
		actualcomplete.append(l[2])
	else:
		print(l)

print(len(completed),len(actualcomplete))

for i in range(len(completed)):
	if completed[i] not in actualcomplete:
		print(completed[i])
