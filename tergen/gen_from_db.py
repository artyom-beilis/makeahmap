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
import numpy as np

net= caffe.Net(sys.argv[1],caffe.TEST)
net.copy_from(sys.argv[2])

data_name='source'

shape = net.blobs[data_name].shape
h=shape[2]
w=shape[3]

print h,w

env=lmdb.open(sys.argv[3])
tx=env.begin()
cur=tx.cursor()

counter=0

for k,v in cur:
    sample=np.zeros([h,3*w])
    datum = caffe.proto.caffe_pb2.Datum()
    datum.ParseFromString(v)
    tmp = caffe.io.datum_to_array(datum)
    tmp2=tmp[1,:,:]
    if tmp2.max() - tmp2.min() < 0.5:
        continue
    net.blobs[data_name].data.flat = tmp[1,:,:].flat
    net.forward()
    gen = net.blobs['generated'].data
    sample[:,0:w]=tmp[1,:,:]
    sample[:,w:2*w]=gen
    sample[:,2*w:3*w]=tmp[0,:,:]
    imsave("res/%05d.pgm" % counter,sample)
    counter+=1
    if counter  > 512:
        break

