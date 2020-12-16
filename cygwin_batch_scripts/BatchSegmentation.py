import math 
import os 
import shutil 


#running command 
#"G:\Labs\EliZelzer\agrawala\Temp\Newest_XPIWIT\XPIWIT\XPIWIT.exe" < "XMLPipelines\cell_pos_1_para_1500.txt"

#inputpath='G:/Labs/EliZelzer/sarahru/Analyzed_Data/Pax3_GP_project_E18.5/S93/Proximal_humerus/m1_wt/nuclei/xpiwit_preprocessed/'   
inputpath='G:/Labs/EliZelzer/sarahru/Analyzed_Data/Pax3_GP_project_E18.5/S93/Proximal_humerus/m2_mut/nuclei/'


#inputpath='../cells/xpiwit_preprocessed/'

f=open('batch_segementation_parameters.txt')

segmentationFile=[]
for line in f:
	l=line.split()
	segmentationFile.append(l)

print(len(segmentationFile))


#segmentationObject='cells'
segmentationObject='nuclei'


if not os.path.exists(segmentationObject):
    os.mkdir(segmentationObject)

if not os.path.exists(segmentationObject+'/processed'):
    os.mkdir(segmentationObject+'/processed')

if not os.path.exists(segmentationObject+'/XMLPipelines'):
    os.mkdir(segmentationObject+'/XMLPipelines')

a=[]



	


def xmlfilename(para):
	return '2017_04_13_Nucleus_Log'+str(para)+'.xml'




count=1

for i in range(len(segmentationFile)):
	#print(segmentationFile[i])
	position=segmentationFile[i][0]
	#para1=int(10*float(segmentationFile[i][1]))
	para2=segmentationFile[i][1]
	pos=position.split('.tif')

	batchname=pos[0]+'_para_'+str(para2)+'.txt'
	f=open(segmentationObject+'/XMLPipelines/'+batchname,'w')
	#f=open('XMLPipelines/'+position+'_'+str(para1)+'_'+str(para2)+'.txt','w')

	
	name=segmentationObject+'/processed/'+str(pos[0])+'_para_'+str(para2)
	if not os.path.exists(name):
    		os.mkdir(name)

		
	f.write('--output ../'+name+'/, result\n') 



	f.write('--input 0, '+inputpath+position+', 3, float\n')
	f.write('--xml   ../'+segmentationObject+'/XMLPipelines/'+ xmlfilename(para2)+'\n')
	f.write('--seed 0\n')
	f.write('--lockfile off\n')
	f.write('--subfolder filterid, filtername\n')
	f.write('--outputformat imagename, filtername\n')
	f.write('--end\n')
	f.close()

	if xmlfilename(para2) not in a:
		a.append(xmlfilename(para2))
		shutil.copy2('./'+xmlfilename(para2), segmentationObject+'/XMLPipelines/')


	if i%8==0:
		f1=open(segmentationObject+'/segment_batch_script_'+str(count)+'.bat','w')
		count=count+1
	f1.write('"G:\Labs\EliZelzer\\agrawala\Temp\\Newest_XPIWIT\XPIWIT\XPIWIT.exe" < "XMLPipelines\\'+ batchname+  '"\n')


for j in a:
	print(j)

print(len(a))

