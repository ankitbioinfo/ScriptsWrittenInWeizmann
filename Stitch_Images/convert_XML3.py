
import sys 
from bs4 import BeautifulSoup
from skimage import io,color          
from skimage.color import rgb2gray
import skimage as sk
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import splprep, splev


name1=sys.argv[1]
outputid=sys.argv[2]

f=open(name1)

content=''
for line in f:
    content+=line


soup = BeautifulSoup(content,'html.parser')


d={}
for tag in soup.find_all('contour'):
    zslice_id=int(tag['slice'])  
    coordinate=[]
    for p in tag.find_all('p'):
        l=str(p).split(',')
        first=l[0].split('(')
        last=l[2].split(')')
        t=[float(first[1]), float(l[1]), float(last[0])]
        coordinate.append(np.array(t))
    d[zslice_id]=np.array(coordinate)
   


image_XY_size=900-1
info=(len(d),900,900)
onlycontour=np.zeros(info,dtype=np.ubyte)
print('size',onlycontour.shape)

print 'ankit',onlycontour.shape 

def uniqueArray(data):
	d1 = np.unique(data, axis=0)
	dic=[]
	for j in range(len(d1)):
		for i in range(len(data)):
			if str(data[i])==str(d1[j]):
			    dic.append(i)
			    break 
	dic.sort()		
	dic.append(-1)		
	return data[dic]


for fi in range(info[0]):
    dd=d[fi]
    data=dd
    data = uniqueArray(dd)
    ##spline to make curve smooth 
    tck,u=splprep([data[:,0],data[:,1]],u=None,s=0,per=1)
    unew=np.linspace(u.min(),u.max(),num=2000)
    data = np.array(splev(unew,tck,der=0))
    data=data.transpose()
 
   
    data=np.floor(data*1000*(900/173.0))  ## 173 is the size of image in microview 
    data=data.astype(int)
  
    x,y= np.where(data>(image_XY_size))
    data[x,y]=image_XY_size
    
    x,y= np.where(data<0)
    data[x,y]=0

    onlycontour[fi,image_XY_size-data[:,1],data[:,0]]=255
    if fi%40==0:
        print(fi) 
   
io.imsave('contourOnly'+str(outputid)+'.tif',onlycontour)

onlycontour=None


