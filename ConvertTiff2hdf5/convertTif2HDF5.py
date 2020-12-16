


import numpy as np 
from skimage import io 
import h5py 

input_filename='put your image name in .tif format here'
output_filename='dataset.h5'

im=io.imread(filename)

print(im.shape)



h5f = h5py.File(output_filename, 'w')

# the shape of chunks should be same as the shape of image
h5f.create_dataset('dataset', data=im, chunks=(1,64,64,3))
h5f.close()




