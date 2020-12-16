import math 
import os 
import shutil 


#running command 
#"G:\Labs\EliZelzer\agrawala\Temp\Newest_XPIWIT\XPIWIT\XPIWIT.exe" < "XMLPipelines\cell_pos_1_para_1500.txt"


inputpath='G:/Labs/EliZelzer/sarahru/Analyzed_Data/Beta1integrin_E16.5_project/S126/proximal_tibia/m3_wt_ml/nuclei/xpiwit_preprocessed/'   

f=open('PT- wt-ml')

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



	


def xmlfilename(para1,para2,para3):
	#return '2018_09_17_rz_cell_newdata_intensity_gaussian_'+str(para1)+'_'+str(para2)+'_'+str(para3)+'.xml'
	#return '2017_04_13_Nucleus_step'+str(para1)+'_Log'+str(para2)+'.xml'
	return '2017_04_13_NucleusSegmentationStdWatershed_HighRes_prehypertrophic_log_'+str(para1)+'_'+str(para2)+'_'+str(para3)+'.xml'




count=1

for i in range(len(segmentationFile)):
	#print(segmentationFile[i])
	position=segmentationFile[i][0]
	para1=int(float(segmentationFile[i][1]))
	para2=int(float(segmentationFile[i][2]))
	para3=int(float(segmentationFile[i][3]))
	pos=position.split('.tif')

	batchname=pos[0]+'_para_'+str(para1)+'_'+str(para2)+'_'+str(para3)+'.txt'
	#batchname=pos[0]+'_para_'+str(para1)+'.txt'
	#batchname=pos[0]+'_para_'+str(para1)+'_'+str(para2)+'.txt'
	f=open(segmentationObject+'/XMLPipelines/'+batchname,'w')
	#f=open('XMLPipelines/'+position+'_'+str(para1)+'_'+str(para2)+'.txt','w')

	
	name=segmentationObject+'/processed/'+str(pos[0])+'_'+str(para1)+'_'+str(para2)+'_'+str(para3)
	if not os.path.exists(name):
    		os.mkdir(name)

		
	f.write('--output ../'+name+'/, result\n') 



	f.write('--input 0, '+inputpath+position+', 3, float\n')
	f.write('--xml   ../'+segmentationObject+'/XMLPipelines/'+ xmlfilename(para1,para2,para3)+'\n')
	f.write('--seed 0\n')
	f.write('--lockfile off\n')
	f.write('--subfolder filterid, filtername\n')
	f.write('--outputformat imagename, filtername\n')
	f.write('--end\n')
	f.close()

	if xmlfilename(para1,para2,para3) not in a:
		a.append(xmlfilename(para1,para2,para3))
		shutil.copy2(''+xmlfilename(para1,para2,para3), segmentationObject+'/XMLPipelines/')


	if i%7==0:
		f1=open(segmentationObject+'/segment_batch_script_'+str(count)+'.bat','w')
		count=count+1
	f1.write('"G:\Labs\EliZelzer\\agrawala\Temp\\Newest_XPIWIT\XPIWIT\XPIWIT.exe" < "XMLPipelines\\'+ batchname+  '"\n')


for j in a:
	print(j)

print(len(a))

