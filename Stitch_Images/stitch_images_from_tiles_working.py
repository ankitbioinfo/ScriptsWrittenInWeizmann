import numpy as np
from skimage import io,color
import skimage as sk 
import matplotlib.pyplot as plt 
import scipy.io as sio


startid=2
endid=21

def copymatrix(files,zstart,ystart,xstart,mysize,BigArray):
    image=io.imread(files)
    info=image.shape
    n=info[0]
    #print (info,image.dtype)
    image=sk.img_as_ubyte(image)
    #image=image[:,0:60,:]
    #print np.shape(image)
    #image=sk.img_as_bool(image)
    #print(image.shape, image.dtype)	
    #print 'Y',ystart,(ystart+mysize[0]),'X',xstart,(xstart+mysize[1]), 'Z',zstart
    #print 'Big', np.shape(BigArray[1:n,ystart:(ystart+mysize[0]),xstart:(xstart+mysize[1])])
    #pixel=np.zeros(n,dtype=np.float)
    for i in range(n):
	#pixel[i]=np.mean(image[i])
        np.copyto(BigArray[i+zstart,ystart:(ystart+mysize[0]),xstart:(xstart+mysize[1])], image[i,:,:]) # check image size 
   
    #return pixel 


f=open('tileS84_m2.data');
tile=[]
for line in f:
    l=line.split()
    t=[]
    for i in l:
        t.append(float(i))
    tile.append(t)

tile=np.array(tile)
tile=tile[startid-1:endid,:]
gap= tile[:,3]-tile[:,2]


fname=[]
array=[]
for fi in range(startid,endid+1):#22):
    if fi<10:
        files='mask_contourOnly0'+str(fi)+'-1.tif'
    else:
        files='mask_contourOnly'+str(fi)+'-1.tif'
    fname.append(files)
    image=io.imread(files)
    info=image.shape
    array.append(info[0])
#print array 
#array= np.array([795, 800, 770, 708, 781, 669, 588, 975, 669, 429, 523, 738, 683, 402, 612, 743, 703, 485, 669, 608])
#array= np.array([795, 800, 770, 708, 781, 669])

print 'micron unit', (gap[1:]/array[1:])



zmicron=np.mean(gap[1:]/array[1:])
print 'zmicron', zmicron 
shift= np.ceil((tile[:,2]-min(tile[:,2]))/zmicron) 
#shift= np.ceil((max(tile[:,3])-tile[:,3])/zmicron) 

print 'Z shift', shift 
zsize=np.int(max(shift + array))


xunique= np.unique(tile[:,0])
xunique=-np.sort(-xunique)
yunique = np.unique(tile[:,1])
#yunique=-np.sort(-yunique)
n=len(xunique)
m=len(yunique)

gridid=np.zeros((len(tile),2),dtype=np.int)
for i in range(len(tile)):
    for j in range(n):
        if tile[i,0]==xunique[j]:
            gridid[i,0]=j;
    for j in range(m):
        if tile[i,1]==yunique[j]:
            gridid[i,1]=j;

#plt.plot(gridid[:,0],gridid[:,1],'.')
#for i in range(len(tile)):
#    plt.text(gridid[i,0],gridid[i,1],str(i+2))
#plt.show()


mysize=[90,90]  # Y, X 


BigArray=np.zeros((zsize,mysize[0]*m,mysize[1]*n),dtype=np.ubyte)#np.uint8)np.ubyte
print 'Stitched Image Dimension', (zsize,m,n), np.shape(BigArray)

#gridpixel=np.zeros((1000,22),dtype=np.float)

for i in range(len(fname)):
    print (fname[i],np.mean(BigArray) )   
    copymatrix(fname[i],int(shift[i]),mysize[0]*gridid[i,1],mysize[1]*gridid[i,0],mysize,BigArray);
    #gridpixel[0:len(value),i]=value

#sio.savemat('zpixel.mat',{'M': gridpixel})

#BigArray=BigArray.astype(bool)
io.imsave('StitchedTiles.tif',BigArray)
#np.array(BigArray, dtype=bool)
sio.savemat('StitchedTiles.mat', {'M':BigArray})
#BigArray=None



#forward shift 
#shift=[0,526,229,536,276,197,179,619,619,178,263,338,328,222,315,798,601,234,100,223,532,203,80,24,0,291,480,306,99,5,221];
#backward shift 
#shift=[0,590,354,687,401,318,280,637,641,262,252,433,435,266,343,731,550,203,154,287,709,511,134, 31,28,494,1105,799,53,0,477];
      


