import paths
import sys

for p in paths.external:
    sys.path.append(p)

import caffe
import lmdb
import random
import gzip
import struct
import numpy as np
from scipy.signal import convolve
from scipy.signal import convolve2d
from scipy.misc import imsave
from scipy.misc import imresize

db=sys.argv[1]
label=int(sys.argv[2])
fl=sys.argv[3]


f=open(fl,'r')

magic, num, rows, cols = struct.unpack(">IIII", f.read(16))
img = np.fromfile(f, dtype=np.uint8).reshape(num, rows, cols);
f.close()

def get_smoothing_kernel(n):
    k=np.array([[1.0]])
    for i in range(n-1):
        k=convolve(k,np.array([[1,1]]))
    k2=np.matmul(np.transpose(k),k)
    k2=k2/k2.sum()
    return k2
    
def get_sk2(n):
    k2=np.ones([2*n+1,2*n+1])
    k2=k2/k2.sum()
    return k2


counter = 0;

def make_patches(img,size):
    N=3
    #kernel = get_smoothing_kernel(2*N+1)
    kernel = get_sk2(2*N+1)
    num=img.shape[0]
    for k in range(num):
        frame_valid = np.zeros([size,size])
        off = (size - 28)/2
        frame_valid[off:-off,off:-off]=img[k,:,:]
        #frame_smooth = convolve2d(frame_valid,kernel,mode='same')
        frame_smooth = imresize(imresize(np.copy(frame_valid),[4,4]),[32,32])
        frame_valid = frame_valid/255.0
        frame_smooth = frame_smooth/255.0

        if label == 0:
            patch = np.zeros([2,size,size])
        else:
            patch = np.zeros([1,size,size])
        patch[0,:,:] = frame_valid;
        if label == 0:
            patch[1,:,:] = frame_smooth

        if False:
            global counter
            for k in range(2):
                tmp=np.zeros([size,size]);
                tmp=patch[k,:,:]
                imsave("res/%05d_%d.pgm" % (counter,k),tmp)
            counter+=1
        yield patch
                


env = lmdb.open(db,map_size=16*1024*1024*1024)

n=0
batch=env.begin(write=True)
for p in make_patches(img,32):
    datum=caffe.io.array_to_datum(p,label)
    if label == 1:
        key = "%010d_%08d" % (random.randint(0,999999),n)
    else:
        key = "%08d" % (n)
    batch.put(key , datum.SerializeToString())
    n=n+1
    print "adding key",key
    if n % 256 == 0:
        batch.commit()
        print "Commit"
        batch = env.begin(write=True)

batch.commit()

