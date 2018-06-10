import sys
sys.path.append('/opt/caffe/caffe-cuda-cudnn/python/')

import caffe
import lmdb
import random
import gzip
import struct
import numpy as np
from scipy.signal import convolve
from scipy.signal import convolve2d
from scipy.misc import imsave

db=sys.argv[1]

def load_frame(fname):
    frame=np.zeros([1201,1201])
    print "Opening ",fname
    with gzip.open(name, 'rb') as f:
        for r in range(1201):
            for c in range(1201):
                v = f.read(2)
                meters = struct.unpack('!h',v)[0]*1.0
                frame[r,c]=meters
    return frame

def get_smoothing_kernel(n):
    k=np.array([[1.0]])
    for i in range(n-1):
        k=convolve(k,np.array([[1,1]]))
    k2=np.matmul(np.transpose(k),k)
    k2=k2/k2.sum()
    return k2
    

counter = 0;

def make_patches(frame,size,step):
    N=10
    kernel = get_smoothing_kernel(2*N+1)
    frame_smooth=convolve2d(frame,kernel,mode='valid')
    frame_valid =frame[N:-N,N:-N];
    frame_shape = frame_smooth.shape;
    print frame_smooth.shape,frame_valid.shape
    for r0 in range(0,frame_shape[0]-size+1,step): 
        for c0 in range(0,frame_shape[1]-size+1,step):
            valid = frame_valid[r0:r0+size,c0:c0+size]
            vmin = valid.min()
            vmax = valid.max();
            if vmin < -1000 or vmax - vmin < 50:
                continue
            patch = np.zeros([2,size,size])
            patch[0,:,:] = frame_valid[r0:r0+size,c0:c0+size]
            patch[1,:,:] = frame_smooth[r0:r0+size,c0:c0+size]

            if False:
                global counter
                for k in range(2):
                    tmp=np.zeros([size,size]);
                    tmp=patch[:,:,k]
                    tmp=(tmp-tmp.min()) / (tmp.max() - tmp.min() + 1);
                    imsave("res/%05d_%d.pgm" % (counter,k),tmp)
                counter+=1
            yield patch
                


env = lmdb.open('db',map_size=16*1024*1024*1024)

batch=env.begin(write=True)

n=0
for name in sys.argv[2:]:
    frame=load_frame(name)
    for p in make_patches(frame,128,64):
        datum=caffe.io.array_to_datum(p)
        key = "%010d_%08d" % (random.randint(0,999999),n)
        batch.put(key , datum.SerializeToString())
        n=n+1
        print "adding key",key
        if n % 100 == 0:
            batch.commit()
            batch = env.begin(write=True)

batch.commit()

