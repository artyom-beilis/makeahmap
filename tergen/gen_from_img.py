import paths
import sys

for p in paths.external:
    sys.path.append(p)

import caffe
import lmdb
import random
import gzip
import struct
from scipy.misc import imsave
from scipy.misc import imread
import numpy as np

net= caffe.Net(sys.argv[1],caffe.TEST)
net.copy_from(sys.argv[2])

data_name='source'

shape = net.blobs[data_name].shape
counter=0
for name in sys.argv[3:]:
    img=imread(name) / 255.0;
    h=img.shape[0]
    w=img.shape[1]
    net.blobs[data_name].reshape(1,1,h,w)
    sample=np.zeros([h,2*w])
    net.blobs[data_name].data.flat = img.flat
    net.forward()
    gen = net.blobs['generated'].data
    sample[:,0:w]=img
    sample[:,w:2*w]=gen
    imsave("res/%05d.pgm" % counter,sample)
    counter+=1
