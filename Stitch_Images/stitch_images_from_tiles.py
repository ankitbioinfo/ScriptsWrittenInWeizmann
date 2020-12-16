import numpy as np
from skimage import io,color
import skimage as sk 

import scipy.io as sio


BigArray=np.zeros((2019,6300,4500),dtype=np.uint8)#bool)#np.uint8)



def copymatrix(fi,xstart,ystart,zstart,mysize):
    files='cells_Pos_'+str(fi)+'.tif'
    image=io.imread(files)
    info=image.shape
    n=info[0]
    print (fi,info,image.dtype)
    image=sk.img_as_ubyte(image)
    print(image.dtype,'\n')	
    for i in range(n):
        #print(data.shape,xstart,ystart,len(range(xstart,xstart+mysize[1])))
        #value=mat[i+zstart,xstart:(xstart+mysize[0]),ystart:(ystart+mysize[1])]
        #print i  
        np.copyto(BigArray[i+zstart,ystart:(ystart+mysize[1]),xstart:(xstart+mysize[0])], image[i,:,:,1])
   


#forward shift 
#shift=[0,526,229,536,276,197,179,619,619,178,263,338,328,222,315,798,601,234,100,223,532,203,80,24,0,291,480,306,99,5,221];
#backward shift 
shift=[0,590,354,687,401,318,280,637,641,262,252,433,435,266,343,731,550,203,154,287,709,511,134, 31,28,494,1105,799,53,0,477];
      

mysize=[900,900]


print ('1',np.mean(BigArray) )   
copymatrix(1,900,0,shift[1],mysize);
print ('2',np.mean(BigArray) )   
copymatrix(2,1800,0,shift[2],mysize);
print ('3',np.mean(BigArray) )   
copymatrix(3,2700,0,shift[3],mysize);
    
copymatrix(7,0,900,shift[7],mysize);
copymatrix(6,900,900,shift[6],mysize);
copymatrix(5,1800,900,shift[5],mysize);
copymatrix(4,2700,900,shift[4],mysize);

copymatrix(8,0,1800,shift[8],mysize);
copymatrix(9,900,1800,shift[9],mysize);
copymatrix(10,1800,1800,shift[10],mysize);
copymatrix(11,2700,1800,shift[11],mysize);


copymatrix(15,0,2700,shift[15],mysize);
copymatrix(14,900,2700,shift[14],mysize);
copymatrix(13,1800,2700,shift[13],mysize);
copymatrix(12,2700,2700,shift[12],mysize);



copymatrix(16,0,3600,shift[16],mysize);
copymatrix(17,900,3600,shift[17],mysize);
copymatrix(18,1800,3600,shift[18],mysize);
copymatrix(19,2700,3600,shift[19],mysize);
copymatrix(20,3600,3600,shift[20],mysize);



copymatrix(25,0,4500,shift[25],mysize);
copymatrix(24,900,4500,shift[24],mysize);
copymatrix(23,1800,4500,shift[23],mysize);
copymatrix(22,2700,4500,shift[22],mysize);
copymatrix(21,3600,4500,shift[21],mysize);




copymatrix(26,0,5400,shift[26],mysize);
copymatrix(27,900,5400,shift[27],mysize);
copymatrix(28,1800,5400,shift[28],mysize);
copymatrix(29,2700,5400,shift[29],mysize);
copymatrix(30,3600,5400,shift[30],mysize);

print('done')
#np.save('save_big_matrix_uint8',BigArray)

for i in range(2019):
        #io.imsave('Global/DT_global_image'+str(i)+'.tif',sk.img_as_uint(BigArray[i]))
        io.imsave('Global/DT_global_image'+str(i)+'.tif',BigArray[i])



