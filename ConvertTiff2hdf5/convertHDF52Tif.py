


import numpy as np 
import h5py 
from skimage import io 
import scipy.io as sio 


filename='input_filename.h5'

f=h5py.File(filename,'r')
print(f.keys())

# you have to check the name of data at root
M=f.get('dataset')
M=np.array(M)
print(M.shape)

#save in .tif format 
io.imsave('output.tif',M)

#save in .mat format (matlab favorite) 
#sio.savemat('ouput.mat',{'M':M})

