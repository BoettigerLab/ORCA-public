'''
This function is a stub. It currently converts the array inside the npy file to an array inside an h5 file called 'polys'
'''

import numpy as np
import os
import h5py
import sys

np_file = sys.argv[1]
data = np.load(np_file)
h5_name = str.replace(np_file,'npy','h5')
h5_file = h5py.File(h5_name, 'w')
h5_file.create_dataset('polys', data=data)
h5_file.close()